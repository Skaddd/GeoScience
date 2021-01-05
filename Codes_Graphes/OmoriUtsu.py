import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


def GraphOmoriUtsu(df, gap, calcul_OU, ms_max_mag, foreshock):

    """
    Fonction qui crée un graphe Omori-Utsu et calcule les K- & p-values
    :param df: données (table pandas)
    :param gap: taille des bars d'histograme (en s)
    :param calcul_OU: calcule de la loi OU (1) ou non (0)
    :param ms_max_mag: 1 = le jour du main shock est considéré comme celui ayant la plus grande magnitude max
                       0 = le jour du main shock est considéré comme celui ayant le plus grand nombre de séisme
    :param foreshock: 1 = analyse des précurseurs
                      0 = analyse des répliques
    :return: objet graphe qui trace le nombre de séisme par jour (ligne) et la magnitude max par jour (point)
             Le titre est constitué de la localisation (jour) du main shock ainsi que des K- et p-values calculées
    """

    nb_gap = int((df['sec'].max() - df['sec'].min()) // gap)
    nt = []
    magmax = []
    title = 'Earthquake number and maximal magnitude per day'

    for i in range(nb_gap + 1):
        dfbis = df[((i * gap + df['sec'].min()) <= df['sec']) & (((i + 1) * gap + df['sec'].min()) > df['sec'])]
        nt.append(len(dfbis))
        magmax.append(dfbis['mag'].max())

    fig, ax1 = plt.subplots()
    x = [i for i in range(nb_gap + 1)]
    # Ordonnée de gauche : # séismes / jour
    color = 'tab:blue'
    ax1.set_xlabel('time (day)')
    ax1.set_ylabel('Earthquake # per day', color=color)
    y = nt
    ax1.plot(x, y, color=color)
    ax1.tick_params(axis='y', labelcolor=color)
    # Ordonnée de droite : magnitude max / jour
    ax2 = ax1.twinx()
    color = 'tab:red'
    ax2.set_ylabel('Max magnitude per day', color=color)
    ax2.scatter(x, magmax, s=10, color=color)
    ax2.tick_params(axis='y', labelcolor=color)

    if calcul_OU == 1:
        if foreshock == 1:
            a, b, idx_ms = RegressionOU_foreshock(1, nt, ms_max_mag, magmax)
            study = 'foreshocks'
        else:
            a, b, idx_ms = RegressionOU(1, nt, ms_max_mag, magmax)
            study = 'aftershocks'

        title = 'Main shock at ' + str(idx_ms) + ' days / p-value (' + study + ') = ' + str(round(-a, 2)) + \
                ' / K-value = ' + str(round(10 ** b))

    ax1.set_title(title)
    fig.tight_layout()
    plt.show()

    if calcul_OU == 1:
        return a, b, idx_ms

    
def RegressionOU(c, nt, ms_max_mag, magmax):

    """
    Fonction qui calcule les coefficients de la régression de la loi OU à partir des répliques
    :param c: c-value
    :param nt: liste des nombres de séisme par jour
    :param ms_max_mag: 1 = le jour du main shock est considéré comme celui ayant la plus grande magnitude max
                       0 = le jour du main shock est considéré comme celui ayant le plus grand nombre de séisme
    :param magmax: liste des magnitudes max par jour
    :return: les 2 coefficients de la régression et la localisation (jour) du main shock
    """

    # Définition de main shock ?
    # 1) la magnitude la plus élevée
    if ms_max_mag == 1:
        idx_ms = magmax.index(max(magmax))
    # 2) Le nombre de séisme le plus élevé
    else:
        idx_ms = nt.index(max(nt))

    t = [(idx_ms + i) for i in range(len(nt) - idx_ms)]
    x = []
    y = []
    for i in range(len(nt) - idx_ms):
        if nt[i + idx_ms] != 0:
            y.append(np.log10(nt[i + idx_ms]))
            x.append(np.log10(c + t[i]))
    [a, b] = np.polyfit(x, y, 1)

    return a, b, idx_ms


def RegressionOU_foreshock(c, nt, ms_max_mag, magmax):
    """
    Fonction qui calcule les coefficients de la régression de la loi OU à partir des précurseurs
    :param c: c-value
    :param nt: liste des nombres de séisme par jour
    :param ms_max_mag: 1 = le jour du main shock est considéré comme celui ayant la plus grande magnitude max
                       0 = le jour du main shock est considéré comme celui ayant le plus grand nombre de séisme
    :param magmax: liste des magnitudes max par jour
    :return: les 2 coefficients de la régression et la localisation (jour) du main shock
    """

    # Définition de main shock ?
    # 1) la magnitude la plus élevée
    if ms_max_mag == 1:
        idx_ms = magmax.index(max(magmax))
    # 2) Le nombre de séisme le plus élevé
    else:
        idx_ms = nt.index(max(nt))

    x = []
    y = []
    for i in range(idx_ms):
        if nt[i] != 0:
            y.append(np.log10(nt[i]))
            x.append(np.log10(c + i))
    [a, b] = np.polyfit(x, y, 1)

    return a, b, idx_ms

                
if __name__ == '__main__':

    df = pd.read_csv('Data/CDSA_EssaimSaintes_2004-2005.txt', sep="\s+")

    # 1 jour = 86 400 sec
    gap = 3600 * 24
    # Calcul de p et K
    GraphOmoriUtsu(df, gap, 0, 0, 0)
