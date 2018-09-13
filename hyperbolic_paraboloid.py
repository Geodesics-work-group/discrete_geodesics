import plotly.offline as py
import plotly.graph_objs as go
from numpy import *

# Geodesic curve

h = 0.1  # size of the intervals of t

distance = 8  # total t

lineX = [-0.5]  # initial x position

lineY = [2.2]  # initial y position

p = [0.07]  # initial x velocity

q = [-0.49]  # initial y velocity

lineZ = [lineX[0]*lineY[0]]

for i in range(0, int((distance/h)-1)):
    p.append(p[i] - (2*h*lineY[i]*p[i]*q[i]) /
             (1+square(lineX[i])+square(lineY[i])))
    q.append(q[i] - (2*h*lineX[i]*p[i]*q[i]) /
             (1+square(lineX[i])+square(lineY[i])))
    lineX.append(lineX[i] + h*p[i])
    lineY.append(lineY[i] + h*q[i])
    lineZ.append(lineX[i+1]*lineY[i+1])

trace = go.Scatter3d(
    x=lineX, y=lineY, z=lineZ, line=dict(width=10), marker=dict(size=6)
)

# Paraboloid

surfaceX = surfaceY = linspace(-3, 3, 300)

surfaceZ = [['nan' for j in surfaceY] for i in surfaceX]

i = 0

for u in surfaceX:
    j = 0
    for v in surfaceY:
        if (square(u) + square(v) <= 9):
            surfaceZ[i][j] = u*v
        j += 1
    i += 1

data = [
    go.Surface(x=surfaceX, y=surfaceY, z=surfaceZ, opacity=0.8),
    trace
]

py.plot(data, filename='paraboloid.html')
