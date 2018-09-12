import plotly.offline as py
import plotly.graph_objs as go
from math import sqrt


def normalize(u, v):
    norm = sqrt(u[0]*u[0] + v[0]*v[0])
    u[0] /= norm
    v[0] /= norm


def squared(a):
    return a*a


h = 0.1  # size of the intervals of t

distance = 8  # total t

lineX = [-0.5]  # initial x position

lineY = [2.2]  # initial y position

p = [0.07]  # initial x velocity

q = [-0.49]  # initial y velocity

lineZ = [lineX[0]*lineY[0]]

for i in range(0, int((distance/h)-1)):
    p.append(p[i] - (2*h*lineY[i]*p[i]*q[i]) /
             (1+squared(lineX[i])+squared(lineY[i])))
    q.append(q[i] - (2*h*lineX[i]*p[i]*q[i]) /
             (1+squared(lineX[i])+squared(lineY[i])))
    lineX.append(lineX[i] + h*p[i])
    lineY.append(lineY[i] + h*q[i])
    lineZ.append(lineX[i+1]*lineY[i+1])

trace = go.Scatter3d(
    x=lineX, y=lineY, z=lineZ, line=dict(width=5)
)

surfaceX = [i/100 for i in range(-500, 500)]

surfaceY = [i/100 for i in range(-500, 500)]

surfaceZ = [['nan' for j in surfaceY] for i in surfaceX]

i = 0

for u in surfaceX:
    j = 0
    for v in surfaceY:
        if (squared(u) + squared(v) <= 9):
            surfaceZ[i][j] = u*v
        j += 1
    i += 1

data = [
    go.Surface(x=surfaceX, y=surfaceY, z=surfaceZ, opacity=0.8),
    trace
]

py.plot(data)
