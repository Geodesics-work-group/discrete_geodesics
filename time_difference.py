import plotly.offline as py
import plotly.graph_objs as go
from numpy import *
import time


def vectorSum(a, b):
    return (a[0]+b[0], a[1]+b[1], a[2]+b[2], a[3]+b[3])


def scalarMult(s, v):
    return (s*v[0], s*v[1], s*v[2], s*v[3])


def f(u, v, p, q):
    res = (p, q, -(2*v*p*q)/(1+square(u)+square(v)),
           -(2*u*p*q)/(1+square(u)+square(v)))
    return res

def findPosition(i, t):
    j = 0
    while ((j < 1024) and (t[j+1] <= i)): j+=1
    return j



# Geodesic curve


n = 50 # amount of iterations in a fixed point problem

distance = 8  # total t

u0 = -0.5  # initial x position

v0 = 2.2  # initial y position

p0 = 0.07  # initial x velocity

q0 = -0.49  # initial y velocity


lineX = linspace(10, 1000, 100)


FEulerY = []
HeunY = []
MiddleY = []
TrapeziumY = []
BEulerY = []
RKY = []


for x in lineX:

    h = distance/x

    # Forward Euler

    u = [u0]
    v = [v0]
    p = [p0]
    q = [q0]

    start_time = time.perf_counter()

    for k in range(0, int(x-1)):
        derivative = f(u[k], v[k], p[k], q[k])
        u.append(u[k]+h*derivative[0])
        v.append(v[k]+h*derivative[1])
        p.append(p[k]+h*derivative[2])
        q.append(q[k]+h*derivative[3])
    
    elapsed_time = time.perf_counter() - start_time

    FEulerY.append(elapsed_time)

    # Heun
    del u[:]
    del v[:]
    del p[:]
    del q[:]

    u = [u0]
    v = [v0]
    p = [p0]
    q = [q0]
    
    start_time = time.perf_counter()

    for k in range(0, int(x-1)):
        first = f(u[k], v[k], p[k], q[k])
        approx = vectorSum((u[k], v[k], p[k], q[k]), scalarMult(h, first))
        second = f(approx[0], approx[1], approx[2], approx[3])
        u.append(u[k]+(h/2)*(first[0] + second[0]))
        v.append(v[k]+(h/2)*(first[1] + second[1]))
        p.append(p[k]+(h/2)*(first[2] + second[2]))
        q.append(q[k]+(h/2)*(first[3] + second[3]))


    elapsed_time = time.perf_counter() - start_time
        
    HeunY.append(elapsed_time)

    # Runge-Kutta
    del u[:]
    del v[:]
    del p[:]
    del q[:]

    u = [u0]
    v = [v0]
    p = [p0]
    q = [q0]

    start_time = time.perf_counter()

    for k in range(0, int(x-1)):
        k1 = scalarMult(h, f(u[k], v[k], p[k], q[k]))
        k2 = scalarMult(h, f(u[k]+(1/2)*k1[0], v[k]+(1/2)*k1[1], p[k]+(1/2)*k1[2], q[k]+(1/2)*k1[3]))
        k3 = scalarMult(h, f(u[k]+(1/2)*k2[0], v[k]+(1/2)*k2[1], p[k]+(1/2)*k2[2], q[k]+(1/2)*k2[3]))
        k4 = scalarMult(h, f(u[k]+k3[0], v[k]+k3[1], p[k]+k3[2], q[k]+k3[3]))
    
        u.append(u[k]+(1/6)*(k1[0]+2*k2[0]+2*k3[0]+k4[0]))
        v.append(v[k]+(1/6)*(k1[1]+2*k2[1]+2*k3[1]+k4[1]))
        p.append(p[k]+(1/6)*(k1[2]+2*k2[2]+2*k3[2]+k4[2]))
        q.append(q[k]+(1/6)*(k1[3]+2*k2[3]+2*k3[3]+k4[3]))


    elapsed_time = time.perf_counter() - start_time
        
    RKY.append(elapsed_time)

    # Middle point
    del u[:]
    del v[:]
    del p[:]
    del q[:]

    u = [u0]
    v = [v0]
    p = [p0]
    q = [q0]

    start_time = time.perf_counter()

    derivative = f(u[0], v[0], p[0], q[0])
    u.append(u[0]+h*derivative[0])
    v.append(v[0]+h*derivative[1])
    p.append(p[0]+h*derivative[2])
    q.append(q[0]+h*derivative[3])

    for k in range(1, int(x-1)):
        derivative = f(u[k], v[k], p[k], q[k])
        u.append(u[k-1]+2*h*derivative[0])
        v.append(v[k-1]+2*h*derivative[1])
        p.append(p[k-1]+2*h*derivative[2])
        q.append(q[k-1]+2*h*derivative[3])


    elapsed_time = time.perf_counter() - start_time
        
    MiddleY.append(elapsed_time)


    # Trapezium

    del u[:]
    del v[:]
    del p[:]
    del q[:]

    u = [u0]
    v = [v0]
    p = [p0]
    q = [q0]
    
    start_time = time.perf_counter()

    for k in range(0, int(x-1)):
        first = f(u[k], v[k], p[k], q[k])
        approx = vectorSum((u[k], v[k], p[k], q[k]), scalarMult(h, first))
        newU = approx[0]
        newV = approx[1]
        newP = approx[2]
        newQ = approx[3]
        for i in range(0, n):
            newF = f(newU, newV, newP, newQ)
            newU = u[k]+(h/2)*(first[0] + newF[0])
            newV = v[k]+(h/2)*(first[1] + newF[1])
            newP = p[k]+(h/2)*(first[2] + newF[2])
            newQ = q[k]+(h/2)*(first[3] + newF[3])
        
        second = f(newU, newV, newP, newQ) 
        u.append(u[k]+(h/2)*(first[0] + second[0]))
        v.append(v[k]+(h/2)*(first[1] + second[1]))
        p.append(p[k]+(h/2)*(first[2] + second[2]))
        q.append(q[k]+(h/2)*(first[3] + second[3]))    
        
    elapsed_time = time.perf_counter() - start_time
        
    TrapeziumY.append(elapsed_time)

    
    #Backwards Euler

    del u[:]
    del v[:]
    del p[:]
    del q[:]

    u = [u0]
    v = [v0]
    p = [p0]
    q = [q0]

    start_time = time.perf_counter()

    for k in range(0, int(x-1)):
        first = f(u[k], v[k], p[k], q[k])
        approx = vectorSum((u[k], v[k], p[k], q[k]), scalarMult(h, first))
        newU = approx[0]
        newV = approx[1]
        newP = approx[2]
        newQ = approx[3]
        for i in range(0, n):
            newF = f(newU, newV, newP, newQ)
            newU = u[k]+h*(newF[0])
            newV = v[k]+h*(newF[1])
            newP = p[k]+h*(newF[2])
            newQ = q[k]+h*(newF[3])
        
        second = f(newU, newV, newP, newQ) 
        u.append(u[k]+h*second[0])
        v.append(v[k]+h*second[1])
        p.append(p[k]+h*second[2])
        q.append(q[k]+h*second[3])


    elapsed_time = time.perf_counter() - start_time
        
    BEulerY.append(elapsed_time)

#Traces

ForwardEuler = go.Scatter(
    x = lineX,
    y = FEulerY,
    mode = 'lines',
    name = 'Euler hacia adelante'
)

BackEuler = go.Scatter(
    x = lineX,
    y = BEulerY,
    mode = 'lines',
    name = 'Euler hacia atrás'
)

Trapezium = go.Scatter(
    x = lineX,
    y = TrapeziumY,
    mode = 'lines',
    name = 'Trapecio'
)

Heun = go.Scatter(
    x = lineX,
    y = HeunY,
    mode = 'lines',
    name = 'Heun'
)

MiddlePoint = go.Scatter(
    x = lineX,
    y = MiddleY,
    mode = 'lines',
    name = 'Punto medio'
)

RungeK = go.Scatter(
    x = lineX,
    y = RKY,
    mode = 'lines',
    name = 'Runge-Kutta'
)

data = [
    ForwardEuler,
    BackEuler,
    Trapezium,
    Heun,
    MiddlePoint,
    RungeK
]

layout = go.Layout(
    xaxis= dict(
        title='Cantidad de pasos'
    ),
    yaxis = dict(
        title = 'Tiempo (s)',
        type='log',
        autorange = True
    )
)

fig = go.Figure(data=data, layout=layout)

py.plot(fig, filename='difference.html')
