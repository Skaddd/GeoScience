# -*- coding:utf-8 -*-

__projet__ = "ProjetTelecomPython2020"
__nom_fichier__ = "InterEventTime"
__author__ = "Corentin Gouache"
__date__ = "janvier 2020"

import pandas as pd
import math
from scipy import stats
import numpy as np
import matplotlib.pyplot as plt


def GraphInterEventTime(time):

    time.sort()
    Norm_IET = np.diff(time) * np.size(time) / (time.max() - time.min())

    beta = Norm_IET.var() / Norm_IET.mean()
    delta = 4.4 * (1 - 4 * (Norm_IET.mean() / beta - 0.5) ** 2)
    Ms_pourcent = Norm_IET.mean() / beta * 100 + delta

    Tmp = Norm_IET[Norm_IET> 0]
    sample_Norm_IE = np.linspace(Tmp.min(), Tmp.max(), 100)
    density_Norm_IET = stats.gaussian_kde(Norm_IET)

    gamma = Norm_IET.mean() / beta
    C = 1 / (math.gamma(gamma) * beta ** gamma)
    proba = C * sample_Norm_IE ** (gamma - 1) * np.exp(-sample_Norm_IE / beta)

    x_log = np.logspace(np.log10(Tmp.min()), np.log10(Tmp.max()), 100)

    fig, ax = plt.subplots()
    ax.plot(sample_Norm_IE, proba, label='Data')
    ax.plot(x_log, 1 / (2 * x_log), label='OU')
    ax.plot(x_log, np.exp(-1 * x_log) / 2, label='Poisson')
    ax.set_xscale('log')
    ax.set_xlabel('Normalized IET')
    ax.set_yscale('log')
    ax.set_ylabel('Probability')
    ax.legend()
    plt.show()


    if Ms_pourcent <= 100:
        return Ms_pourcent
    elif (1 / beta) * 100 <= 100:
        return 1 / beta


if __name__ == '__main__':

    df = pd.read_csv('Data/ReNaSS_Essaim_Barcelonnette_2003-2015.txt', sep="\s+")
    time = np.array(df['sec'])
    Ms_pourcent = GraphInterEventTime(time)
    print(Ms_pourcent)