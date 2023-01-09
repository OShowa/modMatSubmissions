
import numpy as np
import matplotlib.pyplot as plt
from corpus import Corpus
from utils import *


def onebody(): # um corpo solto, um corpo estático (massivo)

    # Constantes
    h = 0.01  # Step size for runge kutta
    G = .01  # Constante de gravidade  CONDIÇÂO DE EQUILIBRIO ACHADA: .01
    Ma = 2  # Massa do corpo a (estático) CONDIÇÂO DE EQUILIBRIO ACHADA:2.
    Mb = 5  # Massa do corpo b (solto) CONDIÇÂO DE EQUILIBRIO ACHADA:5.
    d0 = 2  # Distancia entre os corpos CONDIÇÂO DE EQUILIBRIO ACHADA:2
    v0 = 0.01  # velocidade inicial do corpo a no eixo y CONDIÇÂO DE EQUILIBRIO ACHADA: .01
    debbuging = False  # True se quiser ver os prints
    timespace = 20  # Aumentar o número de amostras de tempo

    print(v0)
    # Corpus
    A = Corpus(np.array([0., 0.]), np.array([0., 0.]), Ma)  # Corpo imóvel
    B = Corpus(np.array([d0, 0.]), np.array([0., v0]), Mb)  # Corpo móvel

    # Arrays de posições
    bArray = [B.pos.copy()]
    aArray = [A.pos.copy()]

    # Arrays de velocidade
    speedsB = [B.abs_vel]
    speedsA = [A.abs_vel]

    maxsp = B.abs_vel
    minsp = B.abs_vel
    for i in range(int(timespace / h)):

        vdist = B.pos - A.pos  # Vetor distancia
        dist = np.linalg.norm(vdist)  # Modulo
        udist = vdist / dist  # Vetor unitário

        dist, B.cntp_vel = rungekutta4(gravitational_acel, [dist, 0], i, h,
                                     (G, Ma))  # Distancia no futuro aproximada por Runge kutta

        vvel_cntp = udist * B.cntp_vel  # Usa o vetor unitário da distancia * o valor absoluto velocidade pra calcular o vetor velocidade cntp

        B.velocity += vvel_cntp

        B.pos += B.velocity # Soma ao vetor distancia a variação em função da velocidade inicial


        B.abs_vel = np.linalg.norm(B.velocity)
        if B.abs_vel > maxsp:
            maxsp = B.abs_vel
        elif B.abs_vel < minsp:
            minsp = B.abs_vel
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

    print("max speed: ", maxsp)
    print("min speed: ", minsp)
    timeframe = range(0, len(speedsA))
    plot_speeds(speedsB, speedsA, timeframe, False)

def twobody(): # dois corpos soltos

    h = 0.01  # Step size for runge kutta
    G = .01  # Constante de gravidade  CONDIÇÂO DE EQUILIBRIO ACHADA: .01
    Ma = 2  # Massa do corpo a CONDIÇÂO DE EQUILIBRIO ACHADA:2.
    Mb = 2  # Massa do corpo b CONDIÇÂO DE EQUILIBRIO ACHADA:5.
    d0 = 2 # Distancia entre os corpos CONDIÇÂO DE EQUILIBRIO ACHADA:2
    v0a = -.005  # velocidade inicial do corpo a no eixo y CONDIÇÂO DE EQUILIBRIO ACHADA: .01
    v0b = .005  # velocidade inicial do corpo b no eixo y
    debbuging = False  # True se quiser ver os prints
    timespace = 8  # Aumentar o número de amostras de tempo

    # Corpos
    A = Corpus(np.array([0., 0.]), np.array([0., v0a]), Ma)
    B = Corpus(np.array([d0, 0.]), np.array([0., v0b]), Mb)

    # Arrays de posições
    bArray = [B.pos.copy()]
    aArray = [A.pos.copy()]

    # Arrays de velocidade
    speedsB = [B.abs_vel]
    speedsA = [A.abs_vel]

    for i in range(int(timespace / h)):

        vdist = B.pos - A.pos
        dist = np.linalg.norm(vdist)
        B.udist = vdist/dist

        A.udist = -B.udist

        # calculando movimento de B no tick
        new_dist_b, B.cntp_vel = rungekutta4(gravitational_acel, [dist, 0], i, h,
                                     (G, Ma))

        B.vvel_cntp = B.udist * B.cntp_vel

        B.velocity += B.vvel_cntp

        B.pos += B.velocity

        # calculando movimento de A no tick
        new_dist_a, A.cntp_vel = rungekutta4(gravitational_acel, [dist, 0], i, h,
                                             (G, Mb))

        A.vvel_cntp = A.udist * A.cntp_vel

        A.velocity += A.vvel_cntp

        A.pos += A.velocity

        # velocidades absolutas
        B.abs_vel = np.linalg.norm(B.velocity)
        A.abs_vel = np.linalg.norm(A.velocity)
        speedsB.append(B.abs_vel)
        speedsA.append(A.abs_vel)

        bArray.append(B.pos.copy())
        aArray.append(A.pos.copy())

    plot_space(bArray, aArray)
    plot_animate(bArray, aArray, timespace)

    timeframe = range(0, len(speedsA))
    plot_speeds(speedsB, speedsA, timeframe, True)




twobody()
