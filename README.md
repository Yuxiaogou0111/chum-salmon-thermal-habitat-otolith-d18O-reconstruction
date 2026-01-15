# chum-salmon-thermal-habitat-otolith-d18O-reconstruction
Supplementary datasets for otolith δ18O-based reconstruction

This directory contains one script for calculation and corresponding nine CSV files for otolith-based individual-level reconstructions.

Each CSV file (~40 MB) represents spatial–temporal solutions constrained by otolith δ18O measurements and seawater isoscapes.

All variables are provided at the otolith-layer resolution and are intended to support reproducibility of spatial and thermal reconstruction analyses presented in the manuscript.

This repository hosts supplementary datasets associated with the manuscript:


[Yuxiao Gou et al.] (2026). Revealing thermal habitat use of chum salmon during primary ocean migration. submitted to Marine Ecology Progress Series.

## Description for Dataset##
The repository contains nine CSV files (~40 MB each) representing reconstructed otolith oxygen isotope results for individual cases.

Each file includes the following columns:

- `case_id` – Otolith layer identifier
- `date` – Reconstructed time period
- `lon` – Longitude (°E)
- `lat` – Latitude (°N)
- `temp` – Sea surface temperature (°C)
- `iso_sw_avg` – Seawater oxygen isotope (δ18O, ‰ VSMOW)
- `iso_est` – Estimated otolith oxygen isotope (δ18O, ‰ VPDB)
- `iso_diff` – Difference between estimated and measured otolith δ18O (‰)
