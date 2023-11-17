import matplotlib.pyplot as plt
import matplotlib.animation as animation
from math import *

def prod(A,b):
    c = []
    for i in range(len(b)):
        c.append(sum([A[i][j] * b[j] for j in range(len(b))]))
    return(c)

def jacobiIteration(x_k,A,b):
    x_k1 = []
    for i in range(len(x_k)):
        C = sum([A[i][j] * x_k[j] for j in range(i)]) + sum([A[i][j] * x_k[j] for j in range(i+1,len(b))])
        x_k1.append((b[i] - C)/A[i][i])
    return(x_k1)

def test(A, b, xi, iteration, niter, xe, titolo):
    asse_x = [float(i) for i in range(0, niter + 1)]
    asse_y = [sqrt(sum([(xi[j] - xe[j]) ** 2 for j in range(len(xe))]))]

    fig, ax = plt.subplots()
    ax.set_xlim([0, len(xe)])
    ax.set_ylim([-1,3])

    xx = [i for i in range(len(xe))]
    fram, = ax.plot(xx, xi)

    def update(frame):
        nonlocal xi
        xi = iteration(xi, A, b)
        norm = sqrt(sum([(xi[j] - xe[j]) ** 2 for j in range(len(xe))]))
        asse_y.append(norm)
        fram.set_ydata(xi)
        return fram,

    ani = animation.FuncAnimation(fig, update, repeat=True, frames=1000,interval=100)

    ani.save('animation.gif', writer='ffmpeg',bitrate=100,fps=150)


def triDiag1(n):
    A = [[0 for j in range(n)] for i in range(n)]
    A[0][0] = 2
    for i in range(1,n):
        A[i][i] = 2
        A[i-1][i] = 1
        A[i][i-1] = 1
    
    return A

n = 50
A = triDiag1(n)
xe = [1. for i in range(n)]
b = prod(A,xe)
x0 = [0. for i in range(n)]
test(A,b,x0,jacobiIteration,4000,xe,"Jacobi")
