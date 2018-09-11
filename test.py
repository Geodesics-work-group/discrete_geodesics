import plotly.offline as py
import plotly.graph_objs as go

x = [i/10 for i in range(-20, 20)]

y = [i/10 for i in range(-20, 20)]

z = [[i*j for j in y] for i in x]

data = [
    go.Surface(x=x, y=y, z=z)
]

py.plot(data)
