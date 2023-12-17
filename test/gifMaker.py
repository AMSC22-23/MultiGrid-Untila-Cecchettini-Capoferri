import matplotlib.pyplot as plt
import numpy as np
from os import listdir
import math
from matplotlib.animation import FuncAnimation


def vectorFromFile(fileName):
    with open(fileName,'r') as file:
        content = file.readlines()
    n = int(content[0])
    f = []
    for i in range(n):
        f.append(float(content[i+1]))
    return np.array(f)


def formatSol(n,u):
    out = []
    temp = []

    for i in range(n):
        for j in range(n):
            temp.append(u[i*n+j])
        out.append(temp)
        temp = []

    return np.matrix(out)




def createGif(n,dir):

    fig, ax = plt.subplots()

    def update(frame):
        fileName = dir + str(frame) + '.mtx'
        vecSol = vectorFromFile(fileName)
        N = int(math.sqrt(len(vecSol)))
        mat = formatSol(N,vecSol)

        ax.clear()
        ax.imshow(mat, cmap='viridis')
        ax.set_title(f'Frame {frame}')

    ani = FuncAnimation(fig, update, frames = range(n), interval=500, repeat=True)
    #ani.save('animation.gif', writer='ffmpeg',bitrate=100,fps=1)

    plt.show()
        
def createGif3d(n,dir):
    width = 10.
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')

    def update(frame):
        fileName = dir + str(frame) + '.mtx'
        vecSol = vectorFromFile(fileName)
        N = int(math.sqrt(len(vecSol)))
        mat = formatSol(N,vecSol)

        X = np.arange(0,N)
        Y = np.arange(0,N)

        X,Y = np.meshgrid(X,Y)

        ax.clear()
        
        ax.plot_surface(X,Y,mat, cmap='viridis')
        ax.set_zlim([0,3.5])

    ani = FuncAnimation(fig, update, frames = range(n), interval=500, repeat=True)
    #ani.save('animation.gif', writer='ffmpeg',bitrate=100,fps=2)

    plt.show()


if __name__ == "__main__":

    folder = './output/'
    lista = list(listdir(folder))
    n = len(lista)
    createGif3d(n,folder)
