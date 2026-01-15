# ============================================================
# calc_iso_diff.py
#
# Calculate the difference between estimated and measured
# otolith oxygen isotope values (iso_diff).
#
# iso_diff is defined as:
#   iso_diff = iso_est - iso_measured
#
# This script is provided to support reproducibility of the
# isotopic reconstruction analyses presented in the manuscript by 
# Y Gou et al. 
# Revealing thermal habitat use of chum salmon during primary ocean migration.
# submitted to Marine Ecology Progress Series
# ============================================================

import os
import numpy as np
import pandas as pd
import datetime as dt
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature

# Code related to process data of Japan Coastal Ocean Predictability Experiment 2M (JCOPE2M) was not shown here.

# Date processing for Time-averaged calculation of d18Oseawater by monthly scale: center_date ± 14 day
def build_date_window(center_date_str, date_range=14, date_interval=1):
    center = dt.datetime.strptime(center_date_str, "%Y%m%d")
    offsets = list(range(-date_range, date_range + 1, date_interval))
    dates = [(center + dt.timedelta(days=dd)).strftime("%Y%m%d") for dd in offsets]
    l_center = offsets.index(0)
    return dates, l_center


# Read JCOPE2 in 3-dimension coordi
def read_3d_from_dat(filepath, li=866, lj=620, byteorder=">"): #li, lj: grid frame of data 
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"file not found: {filepath}")

    filesize = os.path.getsize(filepath)
    rec_size = li * (lj + 1) # float for one single record
    n_float = filesize // 4
    n_rec_full = n_float // rec_size  # total record = vertical layer number

    dtype = np.dtype(byteorder + "f4")

    field_3d = np.zeros((n_rec_full, lj, li), dtype=np.float32)
    with open(filepath, "rb") as f:
        for k in range(n_rec_full):
            data = np.fromfile(f, dtype=dtype, count=rec_size)
            if data.size != rec_size:
                raise ValueError(
                    f"{rec_size}, {data.size}"
                )
            arr2d = data.reshape((lj + 1, li))
            field_3d[k] = arr2d[:lj, :] 

    return field_3d


#  extract water layer ≤30m based on the fact that juvenile chum salmon preliminarily utilize waters shallower than 30m (Morita et al. 2011)
#  if depth .le. depth_lim, average the upper layer, ELSE return NaN
def shallow_average(field_3d, depth_3d, depth_limit=30.0):
    nz = min(field_3d.shape[0], depth_3d.shape[0])
    F = field_3d[:nz]
    D = depth_3d[:nz]

    if F.shape != D.shape:
        raise ValueError(f"{F.shape} ≠ {D.shape}")

    mask = D <= depth_limit  # True: shallow water
    count = mask.sum(axis=0)
    sum_field = np.where(mask, F, 0.0).sum(axis=0)

    field_shallow = np.full_like(sum_field, np.nan, dtype=np.float32)
    valid = count > 0
    field_shallow[valid] = (sum_field[valid] / count[valid]).astype(np.float32)
    return field_shallow


def compute_iso_sw_single_time_simple(sal2d):
    return (0.44 * sal2d - 15.13).astype(np.float32) # d18Osw construction based on the linear relationship with SSS based on LeGrande & Schmidt, 2006.


def compute_iso_sw_avg_simple(sal_window):
    nt = sal_window.shape[0]
    iso_sum = np.zeros_like(sal_window[0], dtype=np.float64)

    for it in range(nt):
        iso_l = compute_iso_sw_single_time_simple(sal_window[it])
        iso_sum += iso_l

    iso_avg = iso_sum / float(nt)
    return iso_avg.astype(np.float32)


# Domain #

def compute_iso_maps_series_depthlimited(
    hdddir,
    layer_path,
    center_dates,
    iso_values,
    li=866,
    lj=620,
    date_range=14,
    date_interval=1,
    byteorder=">",
    depth_limit=30.0,
    use_bottom_mask=True,
    out_prefix=None,
    make_plots=True,
    region=None,      # (lat_min, lat_max, lon_min, lon_max)
    temp_range=None,  # in Celsius degree ºC
    shading="auto",
):
    if len(center_dates) != len(iso_values):
        raise ValueError("different lenght for center_dates & iso_values")

    # lat lon
    lon1d, lat1d = load_lon_lat(hdddir, li=li, lj=lj)
    lon2d, lat2d = np.meshgrid(lon1d, lat1d)

    # region bounds
    if region is not None:
        lat_min, lat_max, lon_min, lon_max = region
        mask_region = (lat2d >= lat_min) & (lat2d <= lat_max) & \
                      (lon2d >= lon_min) & (lon2d <= lon_max)
        region_bounds = (lat_min, lat_max, lon_min, lon_max)
    else:
        mask_region = None
        region_bounds = None

#  bottom mask
    if use_bottom_mask:
        bottom = load_bottom(hdddir, li=li, lj=lj)
        mask_bottom = bottom > 1.0
    else:
        mask_bottom = None

#  depth
    depth_3d = load_depth_3d(layer_path, li=li, lj=lj)  # (nz, lj, li)
    nz_depth = depth_3d.shape[0]
    if out_prefix is not None:
        out_prefix = os.path.abspath(out_prefix)
        head, tail = os.path.split(out_prefix)

        if out_prefix.endswith(os.sep) or os.path.isdir(out_prefix) or (tail and "." not in tail):
            out_dir = out_prefix.rstrip(os.sep)
            base = "iso_diff"
        else:
            out_dir = head
            base = os.path.splitext(tail)[0]

        os.makedirs(out_dir, exist_ok=True)
    else:
        out_dir = None
        base = None

    all_results = []
    all_rows = []
    n_case = len(center_dates)

    # mapping bg
    def plot_field(field, vmin, vmax, cmap, title, suffix):
        if out_dir is None or base is None:
            return
    
        png_path = os.path.join(out_dir, f"{base}_case{idx:02d}_{date_str}{suffix}")

        proj = ccrs.PlateCarree()
        fig = plt.figure(figsize=(7.5, 4.8))
        ax = plt.axes(projection=proj)

        ax.add_feature( 
            cfeature.LAND.with_scale("10m"),
            facecolor="#666666",   # land colour
            edgecolor="none",
            zorder=1,
        )
        ax.coastlines(resolution="10m", linewidth=0.4, zorder=2)

        ocean_cmap = plt.get_cmap(cmap).copy()
        ocean_cmap.set_bad(color="none") 
    
        im = ax.pcolormesh(
            lon2d, lat2d, field,
            transform=proj,
            shading=shading,
            vmin=vmin, vmax=vmax,
            cmap=ocean_cmap,
            zorder=3,
        )

        if region_bounds is not None:
            lat_min_b, lat_max_b, lon_min_b, lon_max_b = region_bounds
        else:
            lat_min_b, lat_max_b = lat2d.min(), lat2d.max()
            lon_min_b, lon_max_b = lon2d.min(), lon2d.max()
    
        ax.set_extent([lon_min_b, lon_max_b, lat_min_b, lat_max_b], crs=proj)
    
        xticks = np.arange(135, 181, 10)
        yticks = np.arange(35, 56, 5) 
        ax.set_xticks(xticks, crs=proj)
        ax.set_yticks(yticks, crs=proj)
        ax.set_xticklabels([f"{x}°E" for x in xticks])
        ax.set_yticklabels([f"{y}°N" for y in yticks])
    

        cb = fig.colorbar(im, ax=ax, orientation="vertical", pad=0.02)
        cb.set_label(title)
    
        ax.set_xlabel("Longitude")
        ax.set_ylabel("Latitude")
        ax.set_title(title)

        fig.subplots_adjust(left=0.08, right=0.98, bottom=0.12, top=0.95)  
        fig.savefig(png_path, dpi=300)
        plt.close(fig)

    ## main cycle ##
    for idx, (date_str, iso_val) in enumerate(zip(center_dates, iso_values), start=1):
        print(f"\nProcessing {idx}/{n_case} : date={date_str}")

        
        date_list, l_center = build_date_window(date_str, date_range, date_interval)
        nt = len(date_list)

        sal_shallow_window = np.zeros((nt, lj, li), dtype=np.float32)
        temp_shallow_window = np.zeros((nt, lj, li), dtype=np.float32)

        
        for it, d in enumerate(date_list):
            s_path = os.path.join(hdddir, "SSS", f"S_{d}.dat")
            t_path = os.path.join(hdddir, "SST", f"T_{d}.dat")

            sal_3d = read_3d_from_dat(s_path, li=li, lj=lj, byteorder=byteorder)
            temp_3d = read_3d_from_dat(t_path, li=li, lj=lj, byteorder=byteorder)

            nz_s = sal_3d.shape[0]
            if nz_s != temp_3d.shape[0]:
                raise ValueError(f"{d}: sal_3d nz={nz_s} ≠ temp_3d nz={temp_3d.shape[0]}")
            nz_use = min(nz_s, nz_depth)
            sal_use = sal_3d[:nz_use]
            temp_use = temp_3d[:nz_use]
            depth_use = depth_3d[:nz_use]

            sal_shallow_window[it] = shallow_average(sal_use, depth_use, depth_limit=depth_limit)
            temp_shallow_window[it] = shallow_average(temp_use, depth_use, depth_limit=depth_limit)

       
        iso_sw_avg = compute_iso_sw_avg_simple(sal_shallow_window)


        temp_center = temp_shallow_window[l_center]


        if temp_range is not None:
            tmin, tmax = temp_range
            mask_temp = (temp_center >= tmin) & (temp_center <= tmax)
        else:
            mask_temp = None

        # calculation for iso_est, based on temperature dependency of d18Ooto fractionation (Gou et al. 2022)
        iso_est = (-0.186 * temp_center + 3.219 + iso_sw_avg).astype(np.float32)
        iso_diff = (iso_val - iso_est).astype(np.float32)

        if mask_bottom is not None:
            for arr in (iso_sw_avg, iso_est, iso_diff):
                arr[~mask_bottom] = np.nan


        if mask_region is not None:
            for arr in (iso_sw_avg, iso_est, iso_diff):
                arr[~mask_region] = np.nan


        if mask_temp is not None:
            for arr in (iso_sw_avg, iso_est, iso_diff):
                arr[~mask_temp] = np.nan

        # packaging results
        case_result = {
            "case_id": idx,
            "date": date_str,
            "lon2d": lon2d,
            "lat2d": lat2d,
            "iso_sw_avg": iso_sw_avg,
            "iso_est": iso_est,
            "iso_diff": iso_diff,
        }
        all_results.append(case_result)


        if out_dir is not None and base is not None:
            df_case = pd.DataFrame({
                "case_id": np.full(lon2d.size, idx, dtype=int),
                "date": np.full(lon2d.size, date_str),
                "lon": lon2d.ravel(),
                "lat": lat2d.ravel(),
                "temp": temp_center.ravel(),        
                "iso_sw_avg": iso_sw_avg.ravel(),
                "iso_est": iso_est.ravel(),
                "iso_diff": iso_diff.ravel(),
            })

            df_case = df_case.dropna(subset=["iso_diff"])
            all_rows.append(df_case)


        if make_plots and out_dir is not None and base is not None:
            #  constraining degree; the value calculated based on the combined uncertainty of d18Ooto-WT, d18Osw-Salinity and the precision of stable isotope analysis
            plot_field(np.abs(iso_diff), 0, 0.31, "hot",
                       f"|iso_diff| (‰), id_{idx}, {date_str}", "_iso_diff.png")
            plot_field(iso_est, -3, 3, "seismic",
                       f"iso_est (‰), id_{idx}, {date_str}", "_iso_est.png")
            plot_field(iso_sw_avg, -1.5, 1.5, "seismic",
                       f"iso_sw_avg (‰), id_{idx}, {date_str}", "_iso_sw_avg.png")


    if out_dir is not None and base is not None and len(all_rows) > 0:
        df_all = pd.concat(all_rows, ignore_index=True)
        csv_path = os.path.join(out_dir, f"{base}_iso_maps_allcases_depth<={int(depth_limit)}m.csv")
        df_all.to_csv(csv_path, index=False)
        print(f"\nCiallo～(∠・ω< )⌒☆ All results saved: {csv_path}")

    return all_results

hdddir = ""
layer_path = ""

# Input for one individual
dates = []
isos  = []


results = compute_iso_maps_series_depthlimited(
    hdddir=hdddir,
    layer_path=layer_path,
    center_dates=dates,
    iso_values=isos,
    li=866,
    lj=620,
    date_range=14,
    date_interval=1,
    byteorder=">",
    depth_limit=30.0,
    use_bottom_mask=True,
    out_prefix="",  
    make_plots=True,
    region=(35, 55, 135, 180),   # Study area
    temp_range=(3.0, 16.0),      # Temperature threshold with 1ºC extension beyond the known ocean habitat temperature range of 4–15ºC (Beschta et al. 1987, Nagata et al. 2007, Kitada et al. 2023)
)




