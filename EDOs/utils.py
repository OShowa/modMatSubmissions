import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation


# s'' = -MG/s^2
# Transformar essa equalçao num sistema de equações de primeira ordem

# v = s'
# v' = a = -MG/s^2
# Fazendo a substituição esse é o novo sistema
# em coordenadadas fica
# y = [s(t),v(t)] entao y' = [v(t),-MG/s(t)^2]
def gravitational_acel(d: float, t: float, M: float, G: float):
    return np.array([d[1], -(M * G) / (d[0] ** 2)])  # retorna um array com [velocidade,aceleração] = [s'(t),s''(t)]


def rungekutta4(f, y0, t, h, args=()):
    # n = len(t)
    # y = np.zeros((n, len(y0)))
    y = y0
    # print('y0:',y0)
    # print(h)
    k1 = f(y, t, *args)
    # print('k1:[velocidade:',k1[0],'aceleracao:',k1[1],']')
    k2 = f(y + k1 * h / 2., t + h / 2., *args)
    # print('k2:[velocidade:',k2[0],'aceleracao:',k2[1],']')
    k3 = f(y + k2 * h / 2., t + h / 2., *args)
    # print('k3:[velocidade:',k3[0],'aceleracao:',k3[1],']')
    k4 = f(y + k3 * h, t + h, *args)
    # print('k4:[velocidade:',k4[0],'aceleracao:',k4[1],']')

    y = y + (h / 6.) * (k1 + 2 * k2 + 2 * k3 + k4)
    # print('y:',y)
    return y


def plot_space(listA, listB):

    xA = [listA[i][0] for i in range(0, len(listA))]
    yA = [listA[i][1] for i in range(0, len(listA))]
    xB = [listB[i][0] for i in range(0, len(listB))]
    yB = [listB[i][1] for i in range(0, len(listB))]

    plt.plot(xA, yA, c='blue')
    plt.plot(xB, yB, c='red')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.grid()


def plot_animate(listA, listB, timespace):
    fig, ax = plt.subplots()

    ims = []

    if timespace > 30:
        nframes = 1800
        if nframes > len(listA):
            nframes = len(listA)
    else:
        nframes = int(timespace * 60)

    step = len(listA) // nframes

    for i in range(0, len(listA), step):
        xA, yA = listA[i][0], listA[i][1]
        xB, yB = listB[i][0], listB[i][0]

        impack = []
        impack.append(ax.scatter(xA, yA, c='blue'))
        impack.append(ax.scatter(xB, yB, c='red'))

        ims.append(impack)

    ani = animation.ArtistAnimation(fig, ims, interval=20, blit=True)
    plt.show()


def plot_speeds(listA, listB, timeframe, twobodies=True):

    print(listA[:20])

    plt.plot(timeframe, listA, c='blue')
    if twobodies:
        plt.plot(timeframe, listB, c='red')
    plt.xlabel('Tempo')
    plt.ylabel('Velocidade')
    plt.grid()
    plt.show()



# t = np.linspace(0,100000,1000)
# rungekutta4(gravitational_acel,np.array([20.0,0.0]),t,(10,10))