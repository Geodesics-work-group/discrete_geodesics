import plotly.offline as py
import plotly.graph_objs as go

x = [i/100 for i in range(-500, 500)]

y = [i/100 for i in range(-500, 500)]

z = [[0 for j in y] for i in x]

i = 0

for u in x:
    j = 0
    for v in y:
        if (u*u + v*v <= 9):
            z[i][j] = u*v
        else:
            z[i][j] = 'nan'
        j += 1
    i += 1

data = [
    go.Surface(x=x, y=y, z=z)
]

py.plot(data)
