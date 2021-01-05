# -*- coding:utf-8 -*-

__projet__ = "ProjetTelecomPython2020"
__nom_fichier__ = "GutenbergRichter"
__author__ = "Corentin Gouache"
__date__ = "janvier 2020"

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


def AB(df, cutoff_mag):

    '''
    Function which creates the graph of the Gutenberg-Richter law and returns the parameters a and b of this law.
    :param df: initial data
    :param timevalue: list of two elements with the minimal and maximal time value (selected with the slider).
    :param magnitudevalue: list of two elements with the minimal and maximal magnitude value (selected with the slider).
    :param option: determines the kind of data you want to return
    :return: 
    option = 1: return the coefficient of the linear regression between the lists of magnitude and the liste of 'log10(N(m>M))'
    option = 2: return the two lists mag (magnitude) and logM (log10(N(m>M)))
    '''

    dmag = 0.1
    mag = np.arange(int(round(cutoff_mag / dmag)), int(round(df['mag'].max() / dmag)) + 1, 1).tolist()
    mag = [i * dmag for i in mag]
    n = len(mag)
    logCumulNumber = []

    for i in range(n):
        filtered_df = df[(df['mag'] >= (mag[i] - dmag / 2))]
        logCumulNumber.append(np.log10(len(filtered_df['mag'])))

    (b, a) = np.polyfit(mag, logCumulNumber, 1)

    y = [a + b * x for x in mag]

    fig, ax = plt.subplots()
    ax.scatter(mag, logCumulNumber, s=10, color='tab:blue', label='Data')
    ax.plot(mag, y, color='tab:red', label='Model')
    ax.set_xlabel('Magnitude')
    ax.set_ylabel('log10(Cumulative #)')
    ax.set_title('y = ' + str(round(a, 2)) + ' ' + str(round(b, 2)) + 'x')
    ax.legend()
    plt.show()

    return a, b

if __name__ == '__main__':

    df = pd.read_csv('Data/ReNaSS_1980-2011_full.txt', sep="\s+")
    a, b = AB(df, 2.2)
    print(a, b)






