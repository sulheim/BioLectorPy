#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
name: plot_boiolector_data
date: 26.02.21
author: Snorre Sulheim

A script for testing the genome-scale metabolic models (GEMs) used in the SEP-outcompete project
"""
# from BioLector import biolector
import sys
sys.path.append("../src/BioLector")
import biolector

from pathlib import Path
folder = Path("data")

################### Bioloector experiment week 5, 2021 ######################
"""
Aim: 
- Test 6 strains on four different carbon sources as a strain QC
- Build experience with the Biolector
- Test parameters wrt medium carbon source concentration, stirring, 
  pre-inocculum, pH-control

Performed by:
- Gunn Broli and Snorre Sulheim
"""

# Experiment #1
# Stirring: 800 rpm
# Default concentration 4g/L
# S. coelicolor missing

if 1:
    fn = folder / "2007-0249_OUT-MF-U05-21_20210201_112613_N1_12-50-10 exp1.xlsx"
    wellmap_fn = folder / "wellmap_U05_exp1.csv"
    fig_folder = Path("results/biolector_U05_exp1")
    B = biolector.PlotResults(fn, wellmap_fn, fig_folder, "Biolector U05-1")
    # B.plot_all_strains_growth(fig_format = "svg", gain_id = "Biomass - 30")
    # B.plot_wellplate(fig_format = "svg", gain_id = "Biomass - 30")
    # # B.plot_individual_well_data_all_strains(fig_format = "svg", gain_id = "Biomass - 30")


################### Vaccibody experiment week xx, 2021 ######################
if 1:
    fn = folder / "VAC-MF-U12-21_Oppsett1 ny BL.xlsx"
    wellmap_fn = folder / "wellmap_U05_exp1.csv"
    fig_folder = Path("results/biolector_U12")
    B = biolector.PlotResults(fn, wellmap_fn, fig_folder, "Biolector U05-1", pro_version = True)
    B.print_parameters()
    B.plot_all_strains_growth(fig_format = "svg", gain_id = "Biomass1 Gain=5")
    B.plot_wellplate(fig_format = "svg", gain_id = "Biomass1 Gain=5")
    B.plot_individual_well_data_all_strains(fig_format = "png", gain_id = "Biomass1 Gain=5", DO_range = (0, 100))
    