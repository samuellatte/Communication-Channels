"""
Script Dessin

Latte Samuel
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import patches
import matplotlib as mpl
import scipy.stats as st
from scipy.stats import norm


# Il faut s'assurer que ces paramètres soient communs entre le main et le dessin
TX = np.array([40,300])  # Position des émetteurs
resol_x = 1 # résolution en x
resol_y = 1 # résolution en y
Puissance_dBm = np.zeros([300, 80])
Vector_Puissances_dBm = np.zeros([290])
SNR_dB = np.zeros([300, 80])
Vector_SNR_dB = np.zeros([290])
Delays = np.zeros([300, 80])
Vector_Delays = np.zeros([290])
RicedB = np.zeros([300, 80])
Vector_RicedB = np.zeros([290])
distance = np.arange(10.5,300,1)
Tau = np.zeros([6])
Alpha = np.zeros([6])
Tap = np.zeros([6])
TDL = np.zeros([6])
Params = np.zeros([4])
PathLoss = np.zeros([290])
CellRange = np.zeros([290])
Fading = np.zeros([290])
ParamsOH = np.zeros([4])
PathLossOH = np.zeros([290])
CellRangeOH = np.zeros([290])
FadingOH = np.zeros([290])
Gamma = np.zeros([290])
WholeCellCoverage = np.zeros([290])

# On commence par définir les murs et les couleurs :
murs = {'mur1': [np.array([0,0]), np.array([80,0])],
        'mur2': [np.array([20,10]), np.array([20,300])],
        'mur3': [np.array([60,10]), np.array([60,70])],
        'mur4': [np.array([60,80]), np.array([60,180])],
        'mur5': [np.array([60,190]), np.array([60,270])],
        'mur6': [np.array([60,280]), np.array([60,300])],
        'mur7': [np.array([0,10]), np.array([20,10])],
        'mur8': [np.array([60,10]), np.array([80,10])],
        'mur9': [np.array([60,80]), np.array([80,80])],
        'mur10': [np.array([60,70]), np.array([80,70])],
        'mur11': [np.array([60,190]), np.array([80,190])],
        'mur12': [np.array([60,180]), np.array([80,180])],
        'mur13': [np.array([60,280]), np.array([80,280])],
        'mur14': [np.array([60,270]), np.array([80,270])]}

colors_rgb = [
    [0,         0,         0],
    [0,         0,    0.5156],
    [0,         0,    0.5312],
    [0,         0,    0.5469],
    [0,         0,    0.5625],
    [0,         0,    0.5781],
    [0,         0,    0.5938],
    [0,         0,    0.6094],
    [0,         0,    0.6250],
    [0,         0,    0.6406],
    [0,         0,    0.6562],
    [0,         0,    0.6719],
    [0,         0,    0.6875],
    [0,         0,    0.7031],
    [0,         0,    0.7188],
    [0,         0,    0.7344],
    [0,         0,    0.7500],
    [0,         0,    0.7656],
    [0,         0,    0.7812],
    [0,         0,    0.7969],
    [0,         0,    0.8125],
    [0,         0,    0.8281],
    [0,         0,    0.8438],
    [0,         0,    0.8594],
    [0,         0,    0.8750],
    [0,         0,    0.8906],
    [0,         0,    0.9062],
    [0,         0,    0.9219],
    [0,         0,    0.9375],
    [0,         0,    0.9531],
    [0,         0,    0.9688],
    [0,         0,    0.9844],
    [0,         0,    1.0000],
    [0,    0.0156,    1.0000],
    [0,    0.0312,    1.0000],
    [0,    0.0469,    1.0000],
    [0,    0.0625,    1.0000],
    [0,    0.0781,    1.0000],
    [0,    0.0938,    1.0000],
    [0,    0.1094,    1.0000],
    [0,    0.1250,    1.0000],
    [0,    0.1406,    1.0000],
    [0,    0.1562,    1.0000],
    [0,    0.1719,    1.0000],
    [0,    0.1875,    1.0000],
    [0,    0.2031,    1.0000],
    [0,    0.2188,    1.0000],
    [0,    0.2344,    1.0000],
    [0,    0.2500,    1.0000],
    [0,    0.2656,    1.0000],
    [0,    0.2812,    1.0000],
    [0,    0.2969,    1.0000],
    [0,    0.3125,    1.0000],
    [0,    0.3281,    1.0000],
    [0,    0.3438,    1.0000],
    [0,    0.3594,    1.0000],
    [0,    0.3750,    1.0000],
    [0,    0.3906,    1.0000],
    [0,    0.4062,    1.0000],
    [0,    0.4219,    1.0000],
    [0,    0.4375,    1.0000],
    [0,    0.4531,    1.0000],
    [0,    0.4688,    1.0000],
    [0,    0.4844,    1.0000],
    [0,    0.5000,    1.0000],
    [0,    0.5156,    1.0000],
    [0,    0.5312,    1.0000],
    [0,    0.5469,    1.0000],
    [0,    0.5625,    1.0000],
    [0,    0.5781,    1.0000],
    [0,    0.5938,    1.0000],
    [0,    0.6094,    1.0000],
    [0,    0.6250,    1.0000],
    [0,    0.6406,    1.0000],
    [0,    0.6562,    1.0000],
    [0,    0.6719,    1.0000],
    [0,    0.6875,    1.0000],
    [0,    0.7031,    1.0000],
    [0,    0.7188,    1.0000],
    [0,    0.7344,    1.0000],
    [0,    0.7500,    1.0000],
    [0,    0.7656,    1.0000],
    [0,    0.7812,    1.0000],
    [0,    0.7969,    1.0000],
    [0,    0.8125,    1.0000],
    [0,    0.8281,    1.0000],
    [0,    0.8438,    1.0000],
    [0,    0.8594,    1.0000],
    [0,    0.8750,    1.0000],
    [0,    0.8906,    1.0000],
    [0,    0.9062,    1.0000],
    [0,    0.9219,    1.0000],
    [0,    0.9375,    1.0000],
    [0,    0.9531,    1.0000],
    [0,    0.9688,    1.0000],
    [0,    0.9844,    1.0000],
    [0,    1.0000,    1.0000],
    [0.0156,    1.0000,    0.9844],
    [0.0312,    1.0000,    0.9688],
    [0.0469,    1.0000,    0.9531],
    [0.0625,    1.0000,    0.9375],
    [0.0781,    1.0000,    0.9219],
    [0.0938,    1.0000,    0.9062],
    [0.1094,    1.0000,    0.8906],
    [0.1250,    1.0000,    0.8750],
    [0.1406,    1.0000,    0.8594],
    [0.1562,    1.0000,    0.8438],
    [0.1719,    1.0000,    0.8281],
    [0.1875,    1.0000,    0.8125],
    [0.2031,    1.0000,    0.7969],
    [0.2188,    1.0000,    0.7812],
    [0.2344,    1.0000,    0.7656],
    [0.2500,    1.0000,    0.7500],
    [0.2656,    1.0000,    0.7344],
    [0.2812,    1.0000,    0.7188],
    [0.2969,    1.0000,    0.7031],
    [0.3125,    1.0000,    0.6875],
    [0.3281,    1.0000,    0.6719],
    [0.3438,    1.0000,    0.6562],
    [0.3594,    1.0000,    0.6406],
    [0.3750,    1.0000,    0.6250],
    [0.3906,    1.0000,    0.6094],
    [0.4062,    1.0000,    0.5938],
    [0.4219,    1.0000,    0.5781],
    [0.4375,    1.0000,    0.5625],
    [0.4531,    1.0000,    0.5469],
    [0.4688,    1.0000,    0.5312],
    [0.4844,    1.0000,    0.5156],
    [0.5000,    1.0000,    0.5000],
    [0.5156,    1.0000,    0.4844],
    [0.5312,    1.0000,    0.4688],
    [0.5469,    1.0000,    0.4531],
    [0.5625,    1.0000,    0.4375],
    [0.5781,    1.0000,    0.4219],
    [0.5938,    1.0000,    0.4062],
    [0.6094,    1.0000,    0.3906],
    [0.6250,    1.0000,    0.3750],
    [0.6406,    1.0000,    0.3594],
    [0.6562,    1.0000,    0.3438],
    [0.6719,    1.0000,    0.3281],
    [0.6875,    1.0000,    0.3125],
    [0.7031,    1.0000,    0.2969],
    [0.7188,    1.0000,    0.2812],
    [0.7344,    1.0000,    0.2656],
    [0.7500,    1.0000,    0.2500],
    [0.7656,    1.0000,    0.2344],
    [0.7812,    1.0000,    0.2188],
    [0.7969,    1.0000,    0.2031],
    [0.8125,    1.0000,    0.1875],
    [0.8281,    1.0000,    0.1719],
    [0.8438,    1.0000,    0.1562],
    [0.8594,    1.0000,    0.1406],
    [0.8750,    1.0000,    0.1250],
    [0.8906,    1.0000,    0.1094],
    [0.9062,    1.0000,    0.0938],
    [0.9219,    1.0000,    0.0781],
    [0.9375,    1.0000,    0.0625],
    [0.9531,    1.0000,    0.0469],
    [0.9688,    1.0000,    0.0312],
    [0.9844,    1.0000,    0.0156],
    [1.0000,    1.0000,         0],
    [1.0000,    0.9844,         0],
    [1.0000,    0.9688,         0],
    [1.0000,    0.9531,         0],
    [1.0000,    0.9375,         0],
    [1.0000,    0.9219,         0],
    [1.0000,    0.9062,         0],
    [1.0000,    0.8906,         0],
    [1.0000,    0.8750,         0],
    [1.0000,    0.8594,         0],
    [1.0000,    0.8438,         0],
    [1.0000,    0.8281,         0],
    [1.0000,    0.8125,         0],
    [1.0000,    0.7969,         0],
    [1.0000,    0.7812,         0],
    [1.0000,    0.7656,         0],
    [1.0000,    0.7500,         0],
    [1.0000,    0.7344,         0],
    [1.0000,    0.7188,         0],
    [1.0000,    0.7031,         0],
    [1.0000,    0.6875,         0],
    [1.0000,    0.6719,         0],
    [1.0000,    0.6562,         0],
    [1.0000,    0.6406,         0],
    [1.0000,    0.6250,         0],
    [1.0000,    0.6094,         0],
    [1.0000,    0.5938,         0],
    [1.0000,    0.5781,         0],
    [1.0000,    0.5625,         0],
    [1.0000,    0.5469,         0],
    [1.0000,    0.5312,         0],
    [1.0000,    0.5156,         0],
    [1.0000,    0.5000,         0],
    [1.0000,    0.4844,         0],
    [1.0000,    0.4688,         0],
    [1.0000,    0.4531,         0],
    [1.0000,    0.4375,         0],
    [1.0000,    0.4219,         0],
    [1.0000,    0.4062,         0],
    [1.0000,    0.3906,         0],
    [1.0000,    0.3750,         0],
    [1.0000,    0.3594,         0],
    [1.0000,    0.3438,         0],
    [1.0000,    0.3281,         0],
    [1.0000,    0.3125,         0],
    [1.0000,    0.2969,         0],
    [1.0000,    0.2812,         0],
    [1.0000,    0.2656,         0],
    [1.0000,    0.2500,         0],
    [1.0000,    0.2344,         0],
    [1.0000,    0.2188,         0],
    [1.0000,    0.2031,         0],
    [1.0000,    0.1875,         0],
    [1.0000,    0.1719,         0],
    [1.0000,    0.1562,         0],
    [1.0000,    0.1406,         0],
    [1.0000,    0.1250,         0],
    [1.0000,    0.1094,         0],
    [1.0000,    0.0938,         0],
    [1.0000,    0.0781,         0],
    [1.0000,    0.0625,         0],
    [1.0000,    0.0469,         0],
    [1.0000,    0.0312,         0],
    [1.0000,    0.0156,         0],
    [1.0000,         0,         0],
    [0.9844,         0,         0],
    [0.9688,         0,         0],
    [0.9531,         0,         0],
    [0.9375,         0,         0],
    [0.9219,         0,         0],
    [0.9062,         0,         0],
    [0.8906,         0,         0],
    [0.8750,         0,         0],
    [0.8594,         0,         0],
    [0.8438,         0,         0],
    [0.8281,         0,         0],
    [0.8125,         0,         0],
    [0.7969,         0,         0],
    [0.7812,         0,         0],
    [0.7656,         0,         0],
    [0.7500,         0,         0],
    [0.7344,         0,         0],
    [0.7188,         0,         0],
    [0.7031,         0,         0],
    [0.6875,         0,         0],
    [0.6719,         0,         0],
    [0.6562,         0,         0],
    [0.6406,         0,         0],
    [0.6250,         0,         0],
    [0.6094,         0,         0],
    [0.5938,         0,         0],
    [0.5781,         0,         0],
    [0.5625,         0,         0],
    [0.5469,         0,         0],
    [0.5312,         0,         0],
    [0.5156,         0,         0],
    [0.5000,         0,         0],
]

# On définit ensuite les variables globales

def dessine1D(): # Dessine l'avion, les émetteurs, les carrés de couleur, la colorbar
    fig1, ax1 = plt.subplots()
    fig2, ax2 = plt.subplots()
    fig3, ax3 = plt.subplots()
    fig4, ax4 = plt.subplots()
    fig5, ax5 = plt.subplots()
    fig6, ax6 = plt.subplots()
    fig7, ax7 = plt.subplots()
    fig8, ax8 = plt.subplots()
    # fig9, ax9 = plt.subplots()
    # fig10, ax10 = plt.subplots()
    # fig11, ax11 = plt.subplots()
    # fig12, ax12 = plt.subplots()
    # fig13, ax13 = plt.subplots()
    fig14, ax14 = plt.subplots()


    n = Params[0]
    p0 = Params[1]

    ax1.set_title("Received power")
    ax1.set_xlabel("Log(d) [m]")
    ax1.set_ylabel("$P_{RX}$ [dBm]")
    ax1.grid()
    ax1.semilogx(distance, Vector_Puissances_dBm, color = 'blue')
    ax1.semilogx(distance, PathLoss, color = 'green', label = 'Path loss, n = %.2f' % (n))
    ax1.legend(loc="upper right")

    ax2.set_title("SNR at UE, 200MHz BW")
    ax2.set_xlabel("Log(d) [m]")
    ax2.set_ylabel("SNR [dB]")
    ax2.grid()
    ax2.semilogx(distance, Vector_SNR_dB, color='blue')

    ax3.set_title("Delay spread")
    ax3.set_xlabel("d [m]")
    ax3.set_ylabel("Delay [s]")
    ax3.grid()
    ax3.plot(distance, Vector_Delays, color='blue')

    ax4.set_title("Rice factor")
    ax4.set_xlabel("d [m]")
    ax4.set_ylabel("Rice factor [dB]")
    ax4.grid()
    ax4.plot(distance, Vector_RicedB, color='blue')

    ax5.set_title("Cell range as a function of connection probability at cell edge")
    ax5.set_xlabel("Log(d) [m]")
    ax5.set_ylabel("Probability")
    ax5.grid()
    ax5.semilogx(distance, CellRange, color='blue')
    ax5.semilogx(distance, np.ones(290)*0.9, color='red', linestyle='dashed')

    sigmaL = Params[2]
    meanL = Params[3]

    ax6.set_title("Statistical fading")
    ax6.set_xlabel("Log(d) [m]")
    ax6.set_ylabel("Fading [dBm]")
    ax6.grid()
    ax6.semilogx(distance, Fading, color='blue', label='$\mu_L$ = %i, $\sigma_L$ = %.2f' %(meanL, sigmaL))
    ax6.legend(loc="upper left")

    ax7.hist(Fading, density=True, bins=100, label="Fading")
    mn, mx = ax7.set_xlim()
    ax7.set_xlim(mn, mx)
    kde_xs = np.linspace(mn, mx, 300)
    kde = st.gaussian_kde(Fading)
    ax7.plot(kde_xs, kde.pdf(kde_xs), label="PDF")
    ax7.legend(loc="upper left")
    ax7.set_ylabel("Probability")
    ax7.set_xlabel("Data")
    ax7.set_title("Fading histogram")

    x = np.linspace(meanL - 3 * sigmaL, meanL + 3 * sigmaL, 100)

    ax8.set_title("Statistical fading distribution")
    ax8.set_xlabel("Power [dBm]")
    ax8.set_ylabel("Probability")
    ax8.plot(x, norm.pdf(x, meanL, sigmaL), 'red', label='Theoretical, $\mu_L$ = %i, $\sigma_L$ = %.2f' %(meanL, sigmaL))
    ax8.plot(kde_xs, kde.pdf(kde_xs), label="Histogram")
    ax8.legend(loc="upper left")
    ax8.grid()

    # nOH = ParamsOH[0]
    # p0OH = ParamsOH[1]
    #
    # ax9.set_title("Puissance en dBm")
    # ax9.set_xlabel("Distance [m]")
    # ax9.set_ylabel("Puissance [dBm]")
    # ax9.grid()
    # ax9.semilogx(distance, Vector_Puissances_dBm, color='blue', label='n = {}'.format(nOH))
    # ax9.semilogx(distance, PathLossOH, color='red', label='p0 = {}'.format(p0OH))
    # ax9.legend(loc="upper right")
    #
    # ax10.set_title("Cell range en fonction de la probabilité de connexion - OH")
    # ax10.set_xlabel("Distance [m]")
    # ax10.set_ylabel("Probabilité")
    # ax10.grid()
    # ax10.semilogx(distance, CellRangeOH, color='blue')
    #
    # ax11.set_title("Fading - OH")
    # ax11.set_xlabel("Distance [m]")
    # ax11.set_ylabel("Fading [dBm]")
    # ax11.grid()
    # ax11.semilogx(distance, FadingOH, color='blue')
    #
    # ax12.hist(FadingOH, density=True, bins=100, label="Data")
    # mn, mx = ax12.set_xlim()
    # ax12.set_xlim(mn, mx)
    # kde_xs = np.linspace(mn, mx, 300)
    # kde = st.gaussian_kde(FadingOH)
    # ax12.plot(kde_xs, kde.pdf(kde_xs), label="PDF")
    # ax12.legend(loc="upper left")
    # ax12.set_ylabel("Probability")
    # ax12.set_xlabel("Data")
    # ax12.set_title("Histogram - OH")
    #
    # sigmaLOH = ParamsOH[2]
    # meanLOH = ParamsOH[3]
    # x = np.linspace(meanLOH - 3 * sigmaLOH, meanLOH + 3 * sigmaLOH, 100)
    #
    # ax13.set_title("Statistical fading distribution - OH")
    # ax13.set_xlabel("Power [dBm]")
    # ax13.set_ylabel("Probability")
    # ax13.plot(x, norm.pdf(x, meanLOH, sigmaLOH), 'red', label='Theoretical')
    # ax13.plot(kde_xs, kde.pdf(kde_xs), label="Histogram")
    # ax13.legend(loc="upper left")
    # ax13.grid()

    reversedGamma = np.flip(Gamma)
    reversedWholeCellCoverage = np.flip(WholeCellCoverage)

    ax14.set_title("Whole cell coverage probability")
    ax14.set_xlabel("$\gamma$ [dB]")
    ax14.set_ylabel("$F_u$")
    ax14.set_xlim(0, 20)
    ax14.grid()
    ax14.plot(reversedGamma, reversedWholeCellCoverage, color='blue', label = '$\sigma_L$ = %.2f' %(sigmaL))
    ax14.plot(Gamma, np.ones(290) * 0.9, color='red', linestyle='dashed')
    ax14.legend(loc="best")

    plt.show()

def dessine2D(): # Dessine l'avion, les émetteurs, les carrés de couleur, la colorbar
    fig1, ax1 = plt.subplots()
    fig2, ax2 = plt.subplots()
    fig3, ax3 = plt.subplots()
    fig4, ax4 = plt.subplots()
    max_dBm = np.amax(Puissance_dBm)
    min_dBm = np.amin(Puissance_dBm)
    max_SNR = np.amax(SNR_dB)
    min_SNR = np.amin(SNR_dB)
    max_Delay = np.nanmax(Delays)
    min_Delay = np.nanmin(Delays)
    max_Rice = np.nanmax(RicedB)
    min_Rice = np.nanmin(RicedB)

    ligne, col = Puissance_dBm.shape

    # On trace les murs :
    for cle in murs:
        ax1.plot([murs[cle][0][0],murs[cle][1][0]],[murs[cle][0][1],murs[cle][1][1]],'grey')
        ax2.plot([murs[cle][0][0],murs[cle][1][0]],[murs[cle][0][1],murs[cle][1][1]],'grey')
        ax3.plot([murs[cle][0][0],murs[cle][1][0]],[murs[cle][0][1],murs[cle][1][1]],'grey')
        ax4.plot([murs[cle][0][0],murs[cle][1][0]],[murs[cle][0][1],murs[cle][1][1]],'grey')

    # On trace l'émetteurs :
    ax1.plot(TX[0],TX[1], marker="*", color='b', label='Emetteur')
    ax2.plot(TX[0],TX[1], marker="*", color='b', label='Emetteur')
    ax3.plot(TX[0],TX[1], marker="*", color='b', label='Emetteur')
    ax4.plot(TX[0],TX[1], marker="*", color='b', label='Emetteur')

    # On itère sur chaque élément d'une ligne
    for i in np.arange(0, col, 1):
        PuissancedBm_col = Puissance_dBm[:, i]
        SNRdB_col = SNR_dB[:, i]
        Delay_col = Delays[:, i]
        Rice_col = RicedB[:, i]

        for j in np.arange(0, ligne, 1):
            ind_color_Prx = int((len(colors_rgb) * (PuissancedBm_col[j] - min_dBm)) / (max_dBm - min_dBm)) # 0 = noir, 1 = min/bleu, 257 = max/rouge
            ind_color_SNR = int((len(colors_rgb) * (SNRdB_col[j] - min_SNR)) / (max_SNR - min_SNR))
            if np.isnan(Delay_col[j]):
                ind_color_Delay = 0
            else:
                ind_color_Delay = int( ((len(colors_rgb) * (Delay_col[j] - min_Delay)) / (max_Delay - min_Delay))+1 ) # 0 = noir, 1 = min/bleu, 257 = max/rouge
            if np.isnan(Rice_col[j]):
                ind_color_Rice = 0
            else:
                ind_color_Rice = int( ((len(colors_rgb) * (Rice_col[j] - min_Rice)) / (max_Rice - min_Rice))+1 ) # 0 = noir, 1 = min/bleu, 257 = max/ro

            if ind_color_Prx >= 257:
                ind_color_Prx = 256
            if ind_color_SNR >= 257:
                ind_color_SNR = 256
            if ind_color_Delay >= 257:
                ind_color_Delay = 256
            if ind_color_Rice >= 257:
                ind_color_Rice = 256

            ax1.add_patch(patches.Rectangle((i, 299 - j), resol_x, resol_y, edgecolor=None, facecolor=colors_rgb[ind_color_Prx], fill=True))
            ax2.add_patch(patches.Rectangle((i, 299 - j), resol_x, resol_y, edgecolor=None, facecolor=colors_rgb[ind_color_SNR],fill=True))
            ax3.add_patch(patches.Rectangle((i, 299 - j), resol_x, resol_y, edgecolor=None, facecolor=colors_rgb[ind_color_Delay],fill=True))
            ax4.add_patch(patches.Rectangle((i, 299 - j), resol_x, resol_y, edgecolor=None, facecolor=colors_rgb[ind_color_Rice],fill=True))

    ax1.set_title("Received power [dBm]")
    ax2.set_title("SNR [dB]")
    ax3.set_title("Delay spread [s]")
    ax4.set_title("Rice factor [dB]")
    ax1.axis('equal')
    ax2.axis('equal')
    ax3.axis('equal')
    ax4.axis('equal')
    ax1.set_facecolor('black')
    ax2.set_facecolor('black')
    ax3.set_facecolor('black')
    ax4.set_facecolor('black')
    ax1.get_xaxis().set_visible(False)
    ax1.get_yaxis().set_visible(False)
    ax2.get_xaxis().set_visible(False)
    ax2.get_yaxis().set_visible(False)
    ax3.get_xaxis().set_visible(False)
    ax3.get_yaxis().set_visible(False)
    ax4.get_xaxis().set_visible(False)
    ax4.get_yaxis().set_visible(False)

    cmap = plt.get_cmap('jet')
    norm1 = mpl.colors.Normalize(vmin=min_dBm, vmax=max_dBm)
    norm2 = mpl.colors.Normalize(vmin=min_SNR, vmax=max_SNR)
    norm3 = mpl.colors.Normalize(vmin=min_Delay, vmax=max_Delay)
    norm4 = mpl.colors.Normalize(vmin=min_Rice, vmax=max_Rice)
    sm1 = plt.cm.ScalarMappable(cmap=cmap, norm=norm1)
    sm1.set_array([])
    sm2 = plt.cm.ScalarMappable(cmap=cmap, norm=norm2)
    sm2.set_array([])
    sm3 = plt.cm.ScalarMappable(cmap=cmap, norm=norm3)
    sm3.set_array([])
    sm4 = plt.cm.ScalarMappable(cmap=cmap, norm=norm4)
    sm4.set_array([])
    plt.colorbar(sm1, ax=ax1, ticks=np.linspace(min_dBm, max_dBm, 7))
    plt.colorbar(sm2, ax=ax2, ticks=np.linspace(min_SNR, max_SNR, 7))
    plt.colorbar(sm3, ax=ax3, ticks=np.linspace(min_Delay, max_Delay, 7))
    plt.colorbar(sm4, ax=ax4, ticks=np.linspace(max_Rice, min_Rice, 7))
    plt.show()

def dessineCIR():
    fig1, ax1 = plt.subplots()
    fig2, ax2 = plt.subplots()

    # RX = np.array([ [20,5] , [60,5] , [60,75] , [60,185] , [60,275] ])  # Position du récepteur

    ax1.stem(Tau, Alpha, basefmt=" ")
    ax1.set_title("Physical channel model, [60,275]")
    ax1.set_xlabel(r"$\tau$ [s]")
    ax1.set_ylabel(r'$\vert h(\tau,t) \vert$')
    ax1.set_xlim([0,2*10**-6])
    ax1.grid()

    ax2.stem(Tap, TDL, basefmt=" ")
    ax2.set_title("Tapped delay line, [60,275], wideband assumption")
    ax2.set_xlabel("tap")
    ax2.set_ylabel(r'$\vert h_{TDL}(\tau,t) \vert$')
    ax2.grid()

    plt.show()

def main(): # En fonction du nombre de rebond(s) précisé(s), ouvre un fichier texte et consigne les valeurs dans les matrices

    plot = 1

    if plot == 1:
        k = 0
        with open('/Users/slatte/Documents/Comm Channels/ProjectRayTracing/ProjectRayTracing/FileVectorPrxdBm.txt',
                  'r') as FileVectorPrxdBm, \
                open('/Users/slatte/Documents/Comm Channels/ProjectRayTracing/ProjectRayTracing/FileVectorSNR.txt',
                     'r') as FileVectorSNR, \
                open('/Users/slatte/Documents/Comm Channels/ProjectRayTracing/ProjectRayTracing/FileVectorDelaySpread.txt',
                     'r') as FileVectorDelaySpread, \
                open('/Users/slatte/Documents/Comm Channels/ProjectRayTracing/ProjectRayTracing/FileVectorRiceFactor.txt',
                     'r') as FileVectorRiceFactor, \
                open('/Users/slatte/Documents/Comm Channels/ProjectRayTracing/ProjectRayTracing/FilePathLoss.txt',
                    'r') as FilePathLoss, \
                open('/Users/slatte/Documents/Comm Channels/ProjectRayTracing/ProjectRayTracing/FileCellRange.txt',
                     'r') as FileCellRange, \
                open('/Users/slatte/Documents/Comm Channels/ProjectRayTracing/ProjectRayTracing/FileFading.txt',
                     'r') as FileFading, \
                open('/Users/slatte/Documents/Comm Channels/ProjectRayTracing/ProjectRayTracing/FilePathLossOH.txt',
                     'r') as FilePathLossOH, \
                open('/Users/slatte/Documents/Comm Channels/ProjectRayTracing/ProjectRayTracing/FileCellRangeOH.txt',
                     'r') as FileCellRangeOH, \
                open('/Users/slatte/Documents/Comm Channels/ProjectRayTracing/ProjectRayTracing/FileFadingOH.txt',
                     'r') as FileFadingOH, \
                open('/Users/slatte/Documents/Comm Channels/ProjectRayTracing/ProjectRayTracing/FileGamma.txt',
                     'r') as FileGamma, \
                open('/Users/slatte/Documents/Comm Channels/ProjectRayTracing/ProjectRayTracing/FileWholeCellCoverage.txt',
                     'r') as FileWholeCellCoverage:


            for Prx, SNR, Delay, Rice, Path, Cell, Fad, PathOH, CellOH, FadOH, Gam, Cover in zip(FileVectorPrxdBm, FileVectorSNR, FileVectorDelaySpread, FileVectorRiceFactor, FilePathLoss, FileCellRange, FileFading, FilePathLossOH, FileCellRangeOH, FileFadingOH, FileGamma, FileWholeCellCoverage):

                Vector_Puissances_dBm[k] = Prx

                Vector_SNR_dB[k] = SNR

                Vector_Delays[k] = Delay

                Vector_RicedB[k] = Rice

                PathLoss[k] = Path

                CellRange[k] = Cell

                Fading[k] = Fad

                PathLossOH[k] = PathOH

                CellRangeOH[k] = CellOH

                FadingOH[k] = FadOH

                Gamma[k] = Gam

                WholeCellCoverage[k] = Cover

                k += 1

        k = 0
        with open(
                '/Users/slatte/Documents/Comm Channels/ProjectRayTracing/ProjectRayTracing/FileParameters.txt','r') as FileParameters, \
                open(
                    '/Users/slatte/Documents/Comm Channels/ProjectRayTracing/ProjectRayTracing/FileParametersOH.txt',
                    'r') as FileParametersOH:

            for Param, ParamOH in zip(FileParameters, FileParametersOH):

                # print(Param)
                Params[k] = Param

                ParamsOH[k] = ParamOH

                k += 1

        dessine1D()

    if plot == 2:
        k = 0
        j = 0
        with open('/Users/slatte/Documents/Comm Channels/ProjectRayTracing/ProjectRayTracing/FilePrxdBm.txt', 'r') as FilePrxdBm,\
            open('/Users/slatte/Documents/Comm Channels/ProjectRayTracing/ProjectRayTracing/FileSNRdB.txt', 'r') as FileSNRdB,\
            open('/Users/slatte/Documents/Comm Channels/ProjectRayTracing/ProjectRayTracing/FileDelaySpread.txt', 'r') as FileDelaySpread,\
            open('/Users/slatte/Documents/Comm Channels/ProjectRayTracing/ProjectRayTracing/FileRiceFactor.txt', 'r') as FileRiceFactor:

            for Prx,SNR,Delay,Rice in zip(FilePrxdBm,FileSNRdB,FileDelaySpread,FileRiceFactor):
                if k < int(80 / resol_x):

                    if float(Prx) > -45:
                        Prx = -45
                    if float(Prx) < -90:
                        Prx = -90
                    Puissance_dBm[299-j, k] = Prx

                    if float(SNR) > 15:
                        SNR = 15
                    if float(SNR) < -10:
                        SNR = -10
                    SNR_dB[299-j, k] = SNR

                    if float(Delay) > 10:
                        Delay = np.nan
                    Delays[299-j, k] = Delay

                    if float(Rice) < -100:
                        Rice = np.nan
                    RicedB[299-j, k] = Rice

                    j += 1
                    if j == 300:
                        j = 0
                        k += 1
            dessine2D()

    if plot == 3:
        k = 0
        with open('/Users/slatte/Documents/Comm Channels/ProjectRayTracing/ProjectRayTracing/FileTau.txt','r') as FileTau, \
            open('/Users/slatte/Documents/Comm Channels/ProjectRayTracing/ProjectRayTracing/FileAlpha.txt','r') as FileAlpha:

            for tau, alpha in zip(FileTau, FileAlpha):

                    Tau[k] = tau

                    Alpha[k] = alpha

                    k += 1

        k = 0
        with open('/Users/slatte/Documents/Comm Channels/ProjectRayTracing/ProjectRayTracing/FileTap.txt','r') as FileTap, \
                open('/Users/slatte/Documents/Comm Channels/ProjectRayTracing/ProjectRayTracing/FileTDL.txt','r') as FileTDL:

            for tap,tdl in zip(FileTap, FileTDL):

                    # print(tap)
                    Tap[k] = tap

                    TDL[k] = tdl

                    k += 1

        dessineCIR()

main()

