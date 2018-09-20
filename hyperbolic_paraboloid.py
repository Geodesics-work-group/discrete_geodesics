import plotly.offline as py
import plotly.graph_objs as go
from numpy import *


# Geodesic curve


h = 0.5  # size of the intervals of t

distance = 4  # total t

u = [-0.5]  # initial x position

v = [2.2]  # initial y position

p = [0.07]  # initial x velocity

q = [-0.49]  # initial y velocity

lineZ = [u[0]*v[0]]

for i in range(0, int((distance/h))):
    p.append(p[i] - (2*h*v[i]*p[i]*q[i]) /
             (1+square(u[i])+square(v[i])))
    q.append(q[i] - (2*h*u[i]*p[i]*q[i]) /
             (1+square(u[i])+square(v[i])))
    u.append(u[i] + h*p[i])
    v.append(v[i] + h*q[i])
    lineZ.append(u[i+1]*v[i+1])

trace = go.Scatter3d(
    x=u, y=v, z=lineZ, line=dict(width=10), marker=dict(size=6)
)

# Paraboloid

surfaceX = surfaceY = linspace(-3, 3, 300)

surfaceZ = [['nan' for j in surfaceY] for i in surfaceX]

i = 0

for x in surfaceX:
    j = 0
    for y in surfaceY:
        if (square(x) + square(y) <= 9):
            surfaceZ[i][j] = x*y
        j += 1
    i += 1

data = [
    go.Surface(x=surfaceX, y=surfaceY, z=surfaceZ, opacity=0.8),
    trace
]

py.plot(data, filename='paraboloid.html')
