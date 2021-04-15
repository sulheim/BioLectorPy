#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
name: biolector
date: 01.03.21
author: Snorre Sulheim

Holds the Biolector class use to plot cultivation data from biolector experiments.
"""
import sys
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from pathlib import Path
import copy
import warnings

plt.rcParams.update({'font.size':16})

well_N_to_fig_shape = {
    1: (1, 1),
    2: (1, 2),
    3: (1, 3),
    4: (2, 2),
    5: (2, 3),
    6: (2, 3),
    7: (2, 4),
    8: (2, 4)
}

class PlotResults(object):
    def __init__(self, data_fn, wellmap_fn, fig_folder, experiment_name = None, pro_version = False):
        self._set_pro_version(pro_version)
        self.parse_data(data_fn)
        self.parse_channel_info(data_fn)
        self.parse_wellmap(wellmap_fn)
        self._set_fig_folder(fig_folder)
        self.experiment_name = experiment_name
        self._set_wellplate_shape()

    def _set_pro_version(self, pro_version):
        self.pro_version = pro_version
        if pro_version:
            self.sheet_name = "Raw Data"
            self.well_column_name = "Well"
            self.skiprows = 23
            self.drop_columns = ["Description", "Content"]
            self.keep_columns = ["Well", "Channel"]
            self.pH_key = 'Cali.pH(HP8) Gain=7'
            self.DO_key = 'Cali.DO(PSt3) Gain=7'
        else:
            self.sheet_name = "raw_data"
            self.well_column_name = "WELL No."
            self.skiprows = 22
            self.drop_columns = ["DESCRIPTION", "CONTENT"]
            self.keep_columns = ["WELL No.", "CHANNEL"]
            self.pH_key = 'Cal.pH [-]:FS=1'
            self.DO_key = 'Cal.pO2 [-]:FS=2'

    def _set_fig_folder(self, path):
        self.fig_folder = Path(path)
        self.fig_folder.mkdir(parents=True, exist_ok=True)

    def _set_wellplate_shape(self):
        self.plate_rows = ['A','B','C','D','E','F']
        self.plate_cols = range(1, 9)

    def parse_data(self, data_fn):
        df = pd.read_excel(data_fn, sheet_name=self.sheet_name, skiprows=self.skiprows)
        df.columns = list(df.columns[:4])+list(df.iloc[1, 4:])
        df = df.loc[df[self.well_column_name].notna(), :]
        df.drop(columns = self.drop_columns, inplace = True)
        self.raw_df = df
        self.df = pd.melt(df, id_vars=self.keep_columns, var_name="Time [h]")

    def parse_channel_info(self, data_fn):
        # Parsing channel info from seperate tble for the "standard" BioLector version
        if not self.pro_version:
            df_channel = pd.read_excel(data_fn, sheet_name=self.sheet_name, skiprows=12, usecols = "A, B, G", nrows=5, index_col = 0)
            df_channel["FILTERNAME"] = df_channel["FILTERNAME"].str.strip()
            self.channel_map = df_channel[['FILTERNAME', 'GAIN']].astype(str).agg(' - '.join, axis=1).to_dict()
            parameters = []
            for i, x in enumerate(list(self.df["CHANNEL"])):
                try:
                    parameters.append(self.channel_map[x])
                except KeyError:
                    parameters.append(x)
            self.df["Parameter"] = parameters
        else:
            self.df["Parameter"] = self.df[self.keep_columns[1]].str[4:]
        # print(self.df.head())

    def parse_wellmap(self, wellmap_fn):
        dtype_dict = {"Well": str, "Strain": str, "Medium":str, "Replicate": int}
        self.df_wellmap = pd.read_csv(wellmap_fn, header = 0, dtype = dtype_dict, sep = None)
        self.df_wellmap.Medium.fillna('', inplace = True)
        self.df = self.df.merge(self.df_wellmap, how='inner', left_on=self.well_column_name, right_on='Well')

    def plot_strain_growth(self, strain_name, gain_id = "Biomass - 30", fig_format = "svg", show = False, ci = 95, confidence_range = True):
        # First check the parameter value (gain id)
        if not self.check_parameter(gain_id):
            print("Could not plot strain growth - wrong biomass ID")
            return False

        if confidence_range:
            units = None
            estimator = np.mean
        else:
            units="Replicate"
            estimator = None

        idx = (self.df["Parameter"] == gain_id) & (self.df["Strain"]==strain_name)
        df_i = self.df.loc[idx, :]
        
        fig, ax = plt.subplots(1, figsize = (12, 8))
        sns.lineplot(data = df_i, x = "Time [h]", y = "value", hue = "Medium", ax = ax, ci = ci_def, units = units, estimator = None)
        ax.set(ylabel = gain_id)
        ax.set_title(strain_name)
        if confidence_range:
            fig.text(0.1, 0.01, "Shaded regions display {0}% CI".format(ci))
        fn = self.fig_folder / "{0}_{1}.{2}".format(strain_name, gain_id, fig_format)
        fig.savefig(fn)
        if show:
            plt.show()
        plt.close()
        return True

    def plot_all_strains_growth(self, gain_id = "Biomass - 50", fig_format = "svg", confidence_range = True):
        # First check the parameter value (gain id)
        if not self.check_parameter(gain_id):
            print("Could not plot strain growth")
            return False

        strains = self.df["Strain"].unique()
        for strain in strains:
            if not isinstance(strain, str):
                continue
            print(strain)
            self.plot_strain_growth(strain, gain_id, fig_format, confidence_range = confidence_range)
        return True
    def print_parameters(self):
        print("Possible parameters: {0}".format(self.df["Parameter"].unique()))

    def check_parameter(self, value):
        acceptable_parameters = list(self.df["Parameter"].unique())
        if not value in acceptable_parameters:
            print("The provided parameter {0} is not in the list of acceptable parameters: {1}".format(value, 
                   "; ".join(acceptable_parameters)))
            return False
        else:
            return True

    def plot_wellplate(self, gain_id = "Biomass - 50", fig_format = "svg", show = False, 
                        rolling_mean_window = 5, fold_change_ref = ("min", "max")):
        """
        Provide a wellplate plot showing the fold change.
        - rolling_mean_window: Indicate the window size of a moving average used to smooth the data to reduce noise
        - fold_change_ref: tuple indicating wheter the fold change ration should be computed between max or min values or start and stop values, or a combination of these.
        """
        # First check the parameter value (gain id)
        if not self.check_parameter(gain_id):
            print("Could not plot wellplate")
            return False

        df_i = self.df.loc[self.df["Parameter"]==gain_id, :]
        # Default settings
        
        fig, ax = plt.subplots(1, figsize = (12, 8))
        vmax = 0
        cmap = plt.cm.YlOrBr

        x_arr = []
        y_arr = []
        c_arr = []
        for i, a in enumerate(self.plate_rows):
            for j in self.plate_cols:
                well_id = "{0}0{1}".format(a, j)
                strain = list(df_i.loc[df_i["Well"]==well_id, "Strain"])[0]
                if not isinstance(strain, str):
                    text = "Empty"
                    ax.scatter(j, 6-i, s = 3.2e3, c = "w", vmin = 1, edgecolor = "k", linewidth = 1)
                else:
                    medium = list(df_i.loc[df_i["Well"]==well_id, "Medium"])[0]
                    fold_change = _get_fold_change(df_i.loc[df_i["Well"]==well_id, "value"], 
                                                   rolling_mean_window, fold_change_ref)
                    x_arr.append(j)
                    y_arr.append(6-i)
                    c_arr.append(fold_change)
                    text = "{0}\n{1}".format(strain, medium)
                ax.annotate(text, (j, 6-i), fontsize = 8, ha = "center", va = "center")#, fontweight='bold')
        vmax = np.max(c_arr)
        sc = ax.scatter(x_arr, y_arr, c = c_arr, s = 3.2e3, cmap = cmap, vmin = 1, vmax = vmax)

        ax.set_xlim(0.5, 8.5)
        ax.set_yticks(range(1,7))
        ax.set_yticklabels(self.plate_rows[::-1])
        ax.set_ylim(0.5, 6.5)
        plt.colorbar(sc, label = "Fold change")
        sns.despine()
        fig.suptitle("{0}: {1} gain".format(self.experiment_name, gain_id))

        # Save fig
        fn = self.fig_folder / "wellplate_{0}.{1}".format(gain_id, fig_format)
        fig.savefig(fn)
        if show:
            plt.show()
        plt.close()
        return True

    def plot_individual_well_data(self, strain_name, gain_id = "Biomass - 30", fig_format = "svg", show = False, pH_range = (6,8), DO_range = (50, 100)):
        # First check the parameter value (gain id)
        if not self.check_parameter(gain_id):
            print("Could not individual well data")
            return False

        idx = self.df["Strain"]==strain_name
        df_i = self.df.loc[idx, :]

        fig_shape = _get_figure_shape(df_i)
        # fig_shape = (3,2)
        well_ids = df_i["Well"].unique()
        
        if fig_shape[0] == 1:
            figsize = (16, 6)
        else:
            figsize = (16, 8)        

        fig, axes = plt.subplots(*fig_shape, figsize = figsize, sharex = True, sharey = True)
        axes = axes.flatten()
        k = 0
        for i in range(fig_shape[0]):
            for j in range(fig_shape[1]):
                ax = axes[k]
                if k < len(axes):
                    well_id = well_ids[k]
                    well_idx = df_i["Well"]==well_id
                    df_well = df_i.loc[well_idx, :]

                    print(well_id)

                    # Plot growth
                    biomass_idx = df_well["Parameter"] == gain_id
                    df_biomass = df_well.loc[biomass_idx, :]
                    l1, = ax.plot(df_biomass["Time [h]"], df_biomass["value"], 'k', label = "Biomass", zorder = 10)
                    

                    # Plot pH
                    ax2 = ax.twinx()
                    pH_idx = df_well["Parameter"] == self.pH_key
                    df_pH = df_well.loc[pH_idx, :]
                    l2, = ax2.plot(df_pH["Time [h]"], df_pH["value"], label = "pH", c = "b")

                    # Plot DO
                    ax3 = ax.twinx()
                    ax3.spines["right"].set_position(("axes", 1.3))
                    _make_patch_spines_invisible(ax3)

                    DO_idx = df_well["Parameter"] == self.DO_key
                    df_DO = df_well.loc[DO_idx, :]
                    l3, = ax3.plot(df_DO["Time [h]"], df_DO["value"], label = "pO2", c = "r")

                    ax2.set_ylim(pH_range)
                    ax3.set_ylim(DO_range)

                    if i == fig_shape[0]-1:
                        ax.set_xlabel("Time [h]")
                        ax2.set_xlabel("Time [h]")
                        ax3.set_xlabel("Time [h]")

                    if j == 0:
                        ax.set_ylabel(gain_id+ " gain")

                    if j == fig_shape[1] -1:
                        ax3.spines["right"].set_visible(True)
                        ax2.set_ylabel("pH")
                        ax2.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
                        ax3.set_ylabel("Dissolved oxygen")

                    else:
                        ax2.get_yaxis().set_visible(False)
                        ax3.get_yaxis().set_visible(False)

                    medium = list(df_well["Medium"])[0]

                    ax.set_title("{0}: {1} - {2}".format(well_id, strain_name, medium), fontsize = 12)


                    ax.yaxis.label.set_color(l1.get_color())
                    ax2.yaxis.label.set_color(l2.get_color())
                    ax3.yaxis.label.set_color(l3.get_color())

                    tkw = dict(size=4, width=1.5)
                    ax.tick_params(axis='y', colors=l1.get_color(), **tkw)
                    ax2.tick_params(axis='y', colors=l2.get_color(), **tkw)
                    ax3.tick_params(axis='y', colors=l3.get_color(), **tkw)
                    ax.tick_params(axis='x', **tkw)

                    if k == 0:
                        lines = [l1, l2, l3]
                        fig.legend(lines, [l.get_label() for l in lines], loc= "upper left",
                                     ncol=3)
            
                    ax.set_zorder(10)
                    ax2.set_zorder(2)
                    ax3.set_zorder(1)
                    ax.patch.set_visible(False)
                
                else:
                    ax.set_visible(False)    
                k += 1
                
        fig.tight_layout(pad = 0.3)

        plt.subplots_adjust(right = 0.8, top = 0.85, left = 0.06, bottom = 0.1)
        # sns.lineplot(data = df_i, x = "Time [h]", y = "value", hue = "Medium", ax = ax)
        # ax.set(ylabel = gain_id)
        
        fn = self.fig_folder / "all_data_{0}.{1}".format(strain_name, fig_format)
        fig.savefig(fn)
        if show:
            plt.show()
        plt.close()

    def plot_individual_well_data_all_strains(self, gain_id = "Biomass - 30", fig_format = "svg", pH_range = (6,8), DO_range = (50, 100)):
        print("Plot data for each well")

        strains = self.df["Strain"].unique()
        for strain in strains:
            if not isinstance(strain, str):
                continue
            self.plot_individual_well_data(strain, gain_id = gain_id, fig_format = fig_format,
                 pH_range = pH_range, DO_range = DO_range)

def _make_patch_spines_invisible(ax):
    """
    From https://matplotlib.org/stable/gallery/ticks_and_spines/multiple_yaxis_with_spines.html
    """
    ax.set_frame_on(True)
    ax.patch.set_visible(False)
    for sp in ax.spines.values():
        sp.set_visible(False)

def _get_figure_shape(df):
    well_ids = df["Well"].unique()
    n = len(well_ids)
    return well_N_to_fig_shape[n]
# def _get_figure_shape(df):
#     well_ids = df["Well"].unique()
#     n_well_letters = len(set([x[0] for x in well_ids]))
#     n_well_numbers = len(set([x[-1] for x in well_ids]))
#     return n_well_letters, n_well_numbers

def _get_fold_change(series, rolling_mean_window, fold_change_ref):
    smoothed_values = series.rolling(rolling_mean_window, center=True).mean()

    if fold_change_ref[0].lower() == "min":
        min_value = smoothed_values.min()
    else:
        min_value = smoothed_values[0]

    if fold_change_ref[1].lower() == "max":
        max_value = smoothed_values.max()
    else:
        max_value = smoothed_values[-1]
    return max_value/min_value
        

if __name__ == '__main__':
    
    folder = "../../data/biolector/"
    fn = folder + "/2007-0249_OUT-MF-U05-21_20210201_112613_N1_12-50-10 exp1.xlsx"
    wellmap_fn = folder + "wellmap_U05_exp1.csv"
    fig_folder = "../../results/biolector_U05_exp1"
    B = PlotResults(fn, wellmap_fn, fig_folder, "Biolector U05-1")
    # B.plot_all_strains_growth(fig_format = "svg", gain_id = "Biomass - 30")
    #B.plot_wellplate(fig_format = "svg", gain_id = "Biomass - 30")
    B.plot_individual_well_data(strain_name = "E. coli", gain_id = "Biomass - 30")

