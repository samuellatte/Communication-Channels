import numpy as np
import matplotlib.pyplot as plt
from scipy import constants

# On commence par définir les murs :

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

coins = {'coin1': np.array([20,10]),
        'coin2': np.array([60,10]),
        'coin3': np.array([60,80]),
        'coin4': np.array([60,190]),
        'coin5': np.array([60,280])}

# On définit ensuite les variables

fig, ax = plt.subplots()

def calc_param():
    for cle in murs:
        u = np.array(murs[cle][1]-np.array(murs[cle][0]))
        u = u / np.linalg.norm(u)
        n = np.array([u[1], -u[0]])  # Permet de calculer le vecteur normal de chaque mur
        murs[cle].append(u)  # Position 2 dans la liste
        murs[cle].append(n)  # Position 3 dans la liste


def verif_transmission(TX,RX):
    d = RX - TX  # On définit le vecteur d
    for cle in murs:  # On vérifie quel mur est traversé par le trajet direct
        if murs[cle][2][0] * d[1] != murs[cle][2][1] * d[0]:
            t = (d[1] * (TX[0] - murs[cle][0][0]) - d[0] * (TX[1] - murs[cle][0][1])) / (murs[cle][2][0] * d[1] - murs[cle][2][1] * d[0])
            s = TX - murs[cle][0]
            if 0 < t < np.linalg.norm(murs[cle][1] - murs[cle][0]):  # Vérifie s'il y a une transmission
                if np.sign(np.dot(s, murs[cle][3])) != np.sign(np.dot(RX - murs[cle][0], murs[cle][3])) and np.sign(np.dot(RX - murs[cle][0], murs[cle][3])) != 0 and np.sign(np.dot(s, murs[cle][3])) != 0:
                    return True
    return False


def calc_image(TX,RX,nbr_rebond):
    if nbr_rebond >= 0 : # On calcule seulement le trajet direct
        Transmission = verif_transmission(TX,RX)
        ClosestDist = 10 ** 13
        if not Transmission :
            ax.plot([TX[0], RX[0]], [TX[1], RX[1]], color="b")
        else:
            for c in coins:
                normDiff = np.linalg.norm(RX - coins[c])
                if normDiff < ClosestDist:
                    ClosestDist = normDiff
                    ind = c
            ax.plot([TX[0], coins[ind][0]], [TX[1], coins[ind][1]], color="b")
            ax.plot([coins[ind][0], RX[0]], [coins[ind][1], RX[1]], color="b")

        if 1 <= nbr_rebond : # On calcule pour 1 rebond
            for cle in murs:
                # print('murs 1ère reflection : ',murs[cle])
                s1 = TX - murs[cle][0]
                image1 = TX - 2 * np.dot(s1, murs[cle][3]) * murs[cle][3]
                d1 = RX - image1
                if murs[cle][2][0] * d1[1] != murs[cle][2][1] * d1[0]: # Condition pour éviter que t1 soit infini
                    t1 = (d1[1] * (image1[0] - murs[cle][0][0]) - d1[0] * (image1[1] - murs[cle][0][1])) / (murs[cle][2][0] * d1[1] - murs[cle][2][1] * d1[0])
                    if 0 < t1 < np.linalg.norm(murs[cle][1] - murs[cle][0]) and np.sign(np.dot(s1, murs[cle][3])) == np.sign(np.dot(RX - murs[cle][0], murs[cle][3])): # Si t1 est dans le mur et que la condition avec s et n est vérifiée
                        P1 = murs[cle][0] + t1 * murs[cle][2] # On a donc un point de rélflexion
                        # print('P1 = ',P1)
                        TransmissionTXP1 = verif_transmission(TX, P1)
                        TransmissionP1RX = verif_transmission(P1, RX)
                        if not TransmissionTXP1 and not TransmissionP1RX:
                            ax.plot([TX[0], P1[0]], [TX[1], P1[1]], 'r', lw=1)
                            ax.plot([P1[0], RX[0]], [P1[1], RX[1]], 'r', lw=1)

                    if 2 <= nbr_rebond:
                        for cle2 in murs:  # On itère sur chaque mur
                            if cle2 != cle:  # Condition pour ne pas retomber sur une image2 qui soit TX
                                s2 = image1 - murs[cle2][0]
                                image2 = image1 - 2 * np.dot(s2, murs[cle2][3]) * murs[cle2][3]
                                d2 = RX - image2
                                if murs[cle2][2][0] * d2[1] != murs[cle2][2][1] * d2[0]:  # Condition pour éviter que t2 soit infini
                                    t2 = (d2[1] * (image2[0] - murs[cle2][0][0]) - d2[0] * (image2[1] - murs[cle2][0][1])) / (murs[cle2][2][0] * d2[1] - murs[cle2][2][1] * d2[0])
                                    if 0 < t2 < np.linalg.norm(murs[cle2][1] - murs[cle2][0]) and np.sign(np.dot(s2, murs[cle2][3])) == np.sign(np.dot(RX - murs[cle2][0], murs[cle2][3])):  # Si t2 est dans le mur et que la condition avec s et n est vérifiée
                                        P2 = murs[cle2][0] + t2 * murs[cle2][2]  # On a donc un point de rélflexion
                                        # print('murs 2ème reflection : ', murs[cle2])
                                        # print('image2 = ', image2)
                                        # print('P2 = ', P2)
                                        d1 = P2 - image1
                                        t1 = (d1[1] * (image1[0] - murs[cle][0][0]) - d1[0] * (image1[1] - murs[cle][0][1])) / (murs[cle][2][0] * d1[1] - murs[cle][2][1] * d1[0])
                                        if 0 < t1 < np.linalg.norm(murs[cle][1] - murs[cle][0]) and (np.sign(np.dot(s1, murs[cle][3])) * np.sign(np.dot(P2 - murs[cle][0], murs[cle][3])) >= 0):  # Si t1 est dans le mur et que la condition avec s et n est vérifiée:
                                            # print('t1 : ', t1)
                                            P1 = murs[cle][0] + t1 * murs[cle][2]
                                            # print('P1 = ', P1)
                                            TransmissionTXP1 = verif_transmission(TX, P1)
                                            TransmissionP1P2 = verif_transmission(P1, P2)
                                            TransmissionP2RX = verif_transmission(P2, RX)
                                            if not TransmissionTXP1 and not TransmissionP1P2 and not TransmissionP2RX:
                                                ax.plot([TX[0], P1[0]], [TX[1], P1[1]], 'y', lw=1)
                                                ax.plot([P1[0], P2[0]], [P1[1], P2[1]], 'y', lw=1)
                                                ax.plot([P2[0], RX[0]], [P2[1], RX[1]], 'y', lw=1)


def dessine_murs(TX,RX):
    for cle in murs:
        ax.plot([murs[cle][0][0],murs[cle][1][0]],[murs[cle][0][1],murs[cle][1][1]], 'grey', lw=1 )
    ax.plot(TX[0], TX[1], marker="*", color='r', label='TX')
    ax.plot(RX[0], RX[1], marker="*", color='b', label='RX')
    ax.legend(loc="best")


def main():
    # Partie définition des paramètres :
    TX = np.array([40,300])  # Position de l'émetteur
    RX = np.array([[60,275]])
    # RX = np.array([ [20,5] , [60,5] , [60,75] , [60,185] , [60,275] ])  # Position du récepteur
    nbr_rebond = 2
    calc_param()

    # Partie dessin :
    for rx in RX:
        dessine_murs(TX, rx)
        calc_image(TX,rx,nbr_rebond)
    ax.axis('equal')
    ax.set_facecolor('black')
    ax.set_title("Rays arriving at [60,275]")
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)
    plt.show()

main()