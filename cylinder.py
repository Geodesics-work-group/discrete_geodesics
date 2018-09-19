from numpy import *
import plotly.offline as py
import plotly.graph_objs as go

# Geodesic curve


def normalize(u, v):
    norm = sqrt(square(u[0]) + square(v[0]))
    u[0] = u[0] / norm
    v[0] = v[0] / norm


h = 0.1  # size of the intervals of t

distance = 5  # total t

u = [0]  # initial u

v = [0]  # initial v

p = [12]  # initial x velocity

q = [8]  # initial y velocity

normalize(p, q)

lineX = [sin(u[0])]

lineY = [cos(u[0])]


for i in range(0, int((distance/h))):
    p.append(p[i])
    q.append(q[i])
    u.append(u[i] + h*p[i])
    v.append(v[i] + h*q[i])
    lineX.append(sin(u[i+1]))
    lineY.append(cos(u[i+1]))

trace = go.Scatter3d(
    x=lineX, y=lineY, z=v, line=dict(width=10), marker=dict(size=6)
)


# Sphere

u = linspace(0, 2*pi, 200)
v = linspace(0, 3, 100)
meshV, meshU = meshgrid(v, u)

data = [
    go.Surface(
        x=sin(meshU),
        y=cos(meshU),
        z=meshV,
        opacity=0.8
    ),
    trace
]

fig = go.Figure(data=data)
py.plot(fig, filename='cylinder.html')
