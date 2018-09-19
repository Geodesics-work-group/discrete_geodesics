from numpy import *
import plotly.offline as py
import plotly.graph_objs as go

# Geodesic curve


def normalize(u, v):
    norm = sqrt(square(u[0]) + square(v[0]))
    u[0] = u[0] / norm
    v[0] = v[0] / norm


h = 0.1  # size of the intervals of t

distance = 2*pi  # total t

u = [0]  # initial u

v = [0]  # initial v

p = [12]  # initial x velocity

q = [8]  # initial y velocity

normalize(p, q)

lineX = [cos(v[0])*cos(u[0])]

lineY = [sin(v[0])*cos(u[0])]

lineZ = [sin(u[0])]

for i in range(0, int((distance/h))):
    p.append(p[i] - h*(cos(u[i])*sin(u[i])*square(q[i])))
    q.append(q[i] + h*(2*tan(u[i])*p[i]*q[i]))
    u.append(u[i] + h*p[i])
    v.append(v[i] + h*q[i])
    lineX.append(cos(v[i+1])*cos(u[i+1]))
    lineY.append(sin(v[i+1])*cos(u[i+1]))
    lineZ.append(sin(u[i+1]))

trace = go.Scatter3d(
    x=lineX, y=lineY, z=lineZ, line=dict(width=10), marker=dict(size=6)
)


# Sphere

u = linspace(-pi/2, pi/2, 100)
v = linspace(0, 2*pi, 200)

vGrid, uGrid = meshgrid(u, v)
surfaceX = outer(cos(v), cos(u))
surfaceY = outer(sin(v), cos(u))
surfaceZ = outer(ones(200), sin(u))

data = go.Data([
    go.Surface(
        x=surfaceX,
        y=surfaceY,
        z=surfaceZ,
        opacity=0.8
    ),
    trace
])

fig = go.Figure(data=data)
py.plot(fig, filename='sphere.html')
