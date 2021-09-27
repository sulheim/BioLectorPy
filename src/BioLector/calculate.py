#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
name: calculate
created: 11.05.21
author: Snorre Sulheim

Holds methods used to estimate statistics from biolector cultivation data
"""
from pathlib import Path
from pygam import LinearGAM, LogisticGAM, s, l, f
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt

def calculate_instant_growth_rates(df, well_id, param = "Biomass - 30", smooth_size = 5, method = "GAM", 
                                   save_fn = None, lag_time = 5):
    """
    Returns a list consisting of growth rate in each time point for a given well.
    :param parameter: Name of the parameter (string), i.e. "Biomass - 30"
    :param well_no: Name of the well (string), i.e. "A01"
    :param B: Object from 'PlotResults' function in biolector.py
    """
    df_i = df.loc[(df.Parameter == param) & (df.Well == well_id), :]
    
    t_arr = np.array(df_i["Time [h]"])
    x_arr = np.array(df_i.value)

    # Smooth the data to reduce noise
    if method == "GAM":
        x_smooth, t_smooth = GAM(t_arr, x_arr, lag_time = lag_time)
    else:
        x_smooth, t_smooth = moving_average(t_arr, x_arr, smooth_size, lag_time = lag_time)
        # x_arr = x_arr[t1:t2]

    if save_fn:
        fig, ax = plt.subplots(1, figsize = (14, 8))
        ax.plot(t_arr, x_arr, label = "Data")
        ax.plot(t_smooth, x_smooth, label = "Fitted data")
        plt.savefig(save_fn)
        plt.close()

    # Calculate growth rate in each time point
    log_xarr = np.log(x_smooth)
    log_dx = log_xarr[1:] - log_xarr[:-1]
    dt = t_smooth[1:] - t_smooth[:-1]
    return log_dx/dt, t_smooth[:-1]

def calculate_max_growth_rate(df, well_id, param = "Biomass - 30", smooth_size = 5, method = "GAM", lag_time = 5):
    """
    Calculate the maximum growth rate based on instant growth rates. A lag-time is included to avoid the weird increase in growth rate often observed during the first few hours in data from the biolector
    """
    instant_growth_rate, t_arr = calculate_instant_growth_rates(df, well_id, param, smooth_size, method, 
                                lag_time = lag_time)

    # discard lag time
    max_gr = np.max(instant_growth_rate[t_arr > lag_time])
    return max_gr

        # for m in mediums:
        #     df_m = df_i.loc[df_i["Medium"] == m]
        #     print(df_m)
        #     print(df_m.columns)
        #     # for w 


def moving_average(t, x, w, lag_time = 5):
    """
    Function to smooth the data to reduce noise.
    :param x: List of data points before noise reduction.
    :param w: Length of sliding window.
    This function will be taking the convolution of the sequence x and a sequence of ones of length w.
    The chosen mode 'valid' means that the convolution product is only given for points where the sequences overlap completely.
    """
    xx = x[t>lag_time]
    tt = t[t>lag_time]
    x_smooth = np.convolve(xx, np.ones(w), 'valid') / w

    t1 = int((w-1)//2)
    t2 = -int(np.ceil((w-1)/2))
    t_smooth = tt[t1:t2]
    return x_smooth, t_smooth

def GAM(t, x, lag_time = 3, plot = False):
    xx = x[t>lag_time]
    tt = t[t>lag_time]
    tt_gs = t[t>lag_time]
    # fit = LinearGAM(n_splines = 10).fit(t,x)
    tt_gs.shape = (tt.shape[0],1)
    fit = LinearGAM(s(0) + l(0)).gridsearch(tt_gs,xx,lam=np.logspace(-3, 3, 11))
    if plot:
        plt.plot(t,x)
        plt.plot(tt,fit.predict(tt), label = "gridsearch")
        plt.legend()
        plt.show()

    return fit.predict(tt), tt

if __name__ == '__main__':
    folder = Path("C:/Users/snorres/git/Outcompete/data/biolector/week15")

    fn = folder / "2007-0249_OUT-MF-U15-21_20210413_101357_N1_10-23-45.xlsx"
    wellmap_fn = folder / "wellmap_U15_exp1.csv"
    fig_folder = "../../results/biolector_U05_exp1"
    B = PlotResults(fn, wellmap_fn, fig_folder, "Biolector U05-1", delimiter = ",")
    print(B.df)