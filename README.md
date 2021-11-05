# BioLectorPy
This package is used to plot results from cultivations in a BioLector. This is an ongoing development where functionality is expanded as needed. The input files are the xlsx-spreadsheets that can be exported from the BioLection software and a simple csv-file (see example) that describes the content of each well in the wellplate. 
Current features include:
 - Support of both BioLector and BioLector Pro
 - Plotting of any parameter (pH, DO, biomass) for:
   - Individual wells
   - Each strain / species across medium conditions
   - Multiple species across one or more medium conditions
 - Well-plate plots summarizing results, e.g. in terms of max-growth rate, fold-increase etc.
 - Estimation of growth rates based on Generalized additive models

# Installation
pip install BioLectorPy

# Contribution
Feel free to contribute by rasing issues or PRs, or reach out to me directly at snorre.sulheim@sintef.no

# Credit
This package is developed by me, Snorre Sulheim, at SINTEF Industry, Department of Biotechnology and Nanomedicine, Norway. If you use this package in your work please cite this repository. 

# License
All files available in this repository is subject to the license given in the license file. 
