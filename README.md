# Guide to Python code for SOTE model

The code contained in this repository corresponds to the SOTE model, as described in Kramer & Mau (2020). Click [here](https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2020WR027456) for the corresponding article in Water Resources Research.

1. In order to run the code, execute SOTE_script.py

2. Input parameters are as follows:

- Irrigation water salinity (*Cirr*, mmol_c/L)
- Irrigation water sodicity fraction (*Eirr*, nondimensional)
- Initial soil solution salinity (*C_init*, mmol_c/L)
- Initial soil sodicity fraction (*E_init*, nondimensional)
- Initial relative soil water content (*s_init*, nondimensional)
- Soil depth (*depth*, mm)
- Rain water salinity (*Crain*, mmol_c/L)
- Rain water sodicity fraction (*Erain*, nondimensional)
- Probably of rain (*rain_prob*, 1/day)
- Mean rainfall depth (*mean_height*, mm)
- Rainy season length (*days*)
- Minimum evapotranspiration rate (*ET_w*, mm)
- Time step (*dt*, day)
- Simulation length (*years*, integer)
- Rainfall season length (*days*, days)
- Ratio of Irrigation to ET: (*ET_ratio*, integer)
- Number of simulations (*runs*, integer)


Full citation to article:
