
import numpy as np
import matplotlib.pyplot as plt
from corpus import Corpus
from utils import *


def onebody():

    # Constantes
    h = 0.01 # Step size for runge kutta
    G = .01 # Constante de gravidade  CONDIÇÂO DE EQUILIBRIO ACHADA: .01
    Ma = 2  # Massa do corpo a (estático) CONDIÇÂO DE EQUILIBRIO ACHADA:2.
    Mb = 5  # Massa do corpo b (solto) CONDIÇÂO DE EQUILIBRIO ACHADA:5.
    d0 = 2  # Distancia entre os corpos CONDIÇÂO DE EQUILIBRIO ACHADA:2
    v0 = .01  # velocidade inicial do corpo a no eixo y CONDIÇÂO DE EQUILIBRIO ACHADA: .01
    debbuging = False  # True se quiser ver os prints
    timespace = 200  # Aumentar o número de amostras de tempo

    # Corpus
    A = Corpus(np.array([0., 0.]), np.array([0., 0.]), Ma)  # Corpo imóvel
    B = Corpus(np.array([d0, 0.]), np.array([0., v0]), Mb)  # Corpo móvel

    # Arrays de posições
    bArray = [B.pos.copy()]
    aArray = [A.pos.copy()]

    # Arrays de velocidade
    speedsB = [B.abs_vel]
    speedsA = [A.abs_vel]

    for i in range(int(timespace / h)):

        vdist = B.pos - A.pos  # Vetor distancia
        dist = np.linalg.norm(vdist)  # Modulo
        udist = vdist / dist  # Vetor unitário

        dist, cntp_vel = rungekutta4(gravitational_acel, [dist, B.cntp_vel], i, h,
                                     (G, Ma))  # Distancia no futuro aproximada por Runge kutta

        vvel_cntp = udist * cntp_vel  # Usa o vetor unitário da distancia * o valor absoluto velocidade pra calcular o vetor velocidade cntp

        B.velocity += vvel_cntp

        B.pos += B.velocity  # Soma ao vetor distancia a variação em função da velocidade inicial

        B.abs_vel = np.linalg.norm(B.velocity)
        speedsB.append(B.abs_vel)
        speedsA.append(A.abs_vel)

        bArray.append(B.pos.copy())
        aArray.append(A.pos.copy())
        if debbuging:
            print('Tempo:', i)
            print('Vetor velocidade centripeta:', vvel_cntp)
            print('Vetor velocidade absoluta:', B.velocity)
            print('Posição:', B.pos)
    plot_space(bArray, aArray)
    plot_animate(bArray, aArray, timespace)

    timeframe = np.arange(0, timespace + h, h)
    plot_speeds(speedsB, speedsA, timeframe, False)

onebody()