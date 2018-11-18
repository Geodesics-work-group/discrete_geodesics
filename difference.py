import plotly.offline as py
import plotly.graph_objs as go
from numpy import *


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

# Adams-Bashforth-Moulton


ABMU = [u0]
ABMV = [v0]
p = [p0]
q = [q0]
ABMT = [0]

ABMH = distance/1024

for k in range(0, 3):
    k1 = scalarMult(ABMH, f(ABMU[k], ABMV[k], p[k], q[k]))
    k2 = scalarMult(ABMH, f(ABMU[k]+(1/2)*k1[0], ABMV[k]+(1/2)*k1[1], p[k]+(1/2)*k1[2], q[k]+(1/2)*k1[3]))
    k3 = scalarMult(ABMH, f(ABMU[k]+(1/2)*k2[0], ABMV[k]+(1/2)*k2[1], p[k]+(1/2)*k2[2], q[k]+(1/2)*k2[3]))
    k4 = scalarMult(ABMH, f(ABMU[k]+k3[0], ABMV[k]+k3[1], p[k]+k3[2], q[k]+k3[3]))
    
    ABMU.append(ABMU[k]+(1/6)*(k1[0]+2*k2[0]+2*k3[0]+k4[0]))
    ABMV.append(ABMV[k]+(1/6)*(k1[1]+2*k2[1]+2*k3[1]+k4[1]))
    p.append(p[k]+(1/6)*(k1[2]+2*k2[2]+2*k3[2]+k4[2]))
    q.append(q[k]+(1/6)*(k1[3]+2*k2[3]+2*k3[3]+k4[3]))
    ABMT.append(ABMT[k]+ABMH)

for k in range(3, 1024):
    f0 = f(ABMU[k], ABMV[k], p[k], q[k])
    f1 = f(ABMU[k-1], ABMV[k-1], p[k-1], q[k-1])
    f2 = f(ABMU[k-2], ABMV[k-2], p[k-2], q[k-2])
    f3 = f(ABMU[k-3], ABMV[k-3], p[k-3], q[k-3])
    prediction = vectorSum((ABMU[k], ABMV[k], p[k], q[k]), scalarMult(ABMH/24, vectorSum(vectorSum(vectorSum(scalarMult(55, f0), scalarMult(-59, f1)), scalarMult(37, f2)), scalarMult(-9, f3))))
    fplus = f(prediction[0], prediction[1], prediction[2], prediction[3])
    ABMU.append(ABMU[k]+(ABMH/24)*(9*fplus[0]+19*f0[0] - 5*f1[0] + f2[0]))
    ABMV.append(ABMV[k]+(ABMH/24)*(9*fplus[1]+19*f0[1] - 5*f1[1] + f2[1]))
    p.append(p[k]+(ABMH/24)*(9*fplus[2]+19*f0[2] - 5*f1[2] + f2[2]))
    q.append(q[k]+(ABMH/24)*(9*fplus[3]+19*f0[3] - 5*f1[3] + f2[3]))
    ABMT.append(ABMT[k]+ABMH)

for x in lineX:

    h = distance/x

    # Forward Euler
    
    del p[:]
    del q[:]

    u = [u0]
    v = [v0]
    p = [p0]
    q = [q0]
    diff = [0]
    t = 0

    for k in range(0, int(x-1)):
        derivative = f(u[k], v[k], p[k], q[k])
        u.append(u[k]+h*derivative[0])
        v.append(v[k]+h*derivative[1])
        p.append(p[k]+h*derivative[2])
        q.append(q[k]+h*derivative[3])
        t = t+h
        pos = findPosition(t, ABMT)
        trueU = ABMU[pos]+(t-ABMT[pos])*(ABMU[pos+1]-ABMU[pos])/ABMH
        trueV = ABMV[pos]+(t-ABMT[pos])*(ABMV[pos+1]-ABMV[pos])/ABMH
        dist = square(u[k+1]-trueU)+square(v[k+1]-trueV)
        diff.append(dist)
        
    total = 0
    for d in diff:
        total += d
    FEulerY.append(total/x)

    # Heun
    del u[:]
    del v[:]
    del p[:]
    del q[:]
    del diff[:]

    u = [u0]
    v = [v0]
    p = [p0]
    q = [q0]
    diff = [0]
    t = 0

    for k in range(0, int(x-1)):
        first = f(u[k], v[k], p[k], q[k])
        approx = vectorSum((u[k], v[k], p[k], q[k]), scalarMult(h, first))
        second = f(approx[0], approx[1], approx[2], approx[3])
        u.append(u[k]+(h/2)*(first[0] + second[0]))
        v.append(v[k]+(h/2)*(first[1] + second[1]))
        p.append(p[k]+(h/2)*(first[2] + second[2]))
        q.append(q[k]+(h/2)*(first[3] + second[3]))
        t = t+h
        pos = findPosition(t, ABMT)
        trueU = ABMU[pos]+(t-ABMT[pos])*(ABMU[pos+1]-ABMU[pos])/ABMH
        trueV = ABMV[pos]+(t-ABMT[pos])*(ABMV[pos+1]-ABMV[pos])/ABMH
        dist = square(u[k+1]-trueU)+square(v[k+1]-trueV)
        diff.append(dist)


    total = 0
    for d in diff:
        total += d
        
    HeunY.append(total/x)

    # Runge-Kutta
    del u[:]
    del v[:]
    del p[:]
    del q[:]
    del diff[:]

    u = [u0]
    v = [v0]
    p = [p0]
    q = [q0]
    diff = [0]
    t = 0

    for k in range(0, int(x-1)):
        k1 = scalarMult(h, f(u[k], v[k], p[k], q[k]))
        k2 = scalarMult(h, f(u[k]+(1/2)*k1[0], v[k]+(1/2)*k1[1], p[k]+(1/2)*k1[2], q[k]+(1/2)*k1[3]))
        k3 = scalarMult(h, f(u[k]+(1/2)*k2[0], v[k]+(1/2)*k2[1], p[k]+(1/2)*k2[2], q[k]+(1/2)*k2[3]))
        k4 = scalarMult(h, f(u[k]+k3[0], v[k]+k3[1], p[k]+k3[2], q[k]+k3[3]))
    
        u.append(u[k]+(1/6)*(k1[0]+2*k2[0]+2*k3[0]+k4[0]))
        v.append(v[k]+(1/6)*(k1[1]+2*k2[1]+2*k3[1]+k4[1]))
        p.append(p[k]+(1/6)*(k1[2]+2*k2[2]+2*k3[2]+k4[2]))
        q.append(q[k]+(1/6)*(k1[3]+2*k2[3]+2*k3[3]+k4[3]))

        t = t+h
        pos = findPosition(t, ABMT)
        trueU = ABMU[pos]+(t-ABMT[pos])*(ABMU[pos+1]-ABMU[pos])/ABMH
        trueV = ABMV[pos]+(t-ABMT[pos])*(ABMV[pos+1]-ABMV[pos])/ABMH
        dist = square(u[k+1]-trueU)+square(v[k+1]-trueV)
        diff.append(dist)


    total = 0
    for d in diff:
        total += d
        
    RKY.append(total/x)

    # Middle point
    del u[:]
    del v[:]
    del p[:]
    del q[:]
    del diff[:]

    u = [u0]
    v = [v0]
    p = [p0]
    q = [q0]
    diff = [0]
    t = 0

    derivative = f(u[0], v[0], p[0], q[0])
    u.append(u[0]+h*derivative[0])
    v.append(v[0]+h*derivative[1])
    p.append(p[0]+h*derivative[2])
    q.append(q[0]+h*derivative[3])
    t = t+h
    pos = findPosition(t, ABMT)
    trueU = ABMU[pos]+(t-ABMT[pos])*(ABMU[pos+1]-ABMU[pos])/ABMH
    trueV = ABMV[pos]+(t-ABMT[pos])*(ABMV[pos+1]-ABMV[pos])/ABMH
    dist = square(u[1]-trueU)+square(v[1]-trueV)
    diff.append(dist)

    for k in range(1, int(x-1)):
        derivative = f(u[k], v[k], p[k], q[k])
        u.append(u[k-1]+2*h*derivative[0])
        v.append(v[k-1]+2*h*derivative[1])
        p.append(p[k-1]+2*h*derivative[2])
        q.append(q[k-1]+2*h*derivative[3])
        t = t+h
        pos = findPosition(t, ABMT)
        trueU = ABMU[pos]+(t-ABMT[pos])*(ABMU[pos+1]-ABMU[pos])/ABMH
        trueV = ABMV[pos]+(t-ABMT[pos])*(ABMV[pos+1]-ABMV[pos])/ABMH
        dist = square(u[k+1]-trueU)+square(v[k+1]-trueV)
        diff.append(dist)


    total = 0
    for d in diff:
        total += d
        
    MiddleY.append(total/x)


    # Trapezium

    del u[:]
    del v[:]
    del p[:]
    del q[:]
    del diff[:]

    u = [u0]
    v = [v0]
    p = [p0]
    q = [q0]
    diff = [0]
    t = 0

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
        t = t+h
        pos = findPosition(t, ABMT)
        trueU = ABMU[pos]+(t-ABMT[pos])*(ABMU[pos+1]-ABMU[pos])/ABMH
        trueV = ABMV[pos]+(t-ABMT[pos])*(ABMV[pos+1]-ABMV[pos])/ABMH
        dist = square(u[k+1]-trueU)+square(v[k+1]-trueV)
        diff.append(dist)     
        
    total = 0
    for d in diff:
        total += d
        
    TrapeziumY.append(total/x)

    
    #Backwards Euler

    del u[:]
    del v[:]
    del p[:]
    del q[:]
    del diff[:]

    u = [u0]
    v = [v0]
    p = [p0]
    q = [q0]
    diff = [0]
    t = 0

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
        t = t+h
        pos = findPosition(t, ABMT)
        trueU = ABMU[pos]+(t-ABMT[pos])*(ABMU[pos+1]-ABMU[pos])/ABMH
        trueV = ABMV[pos]+(t-ABMT[pos])*(ABMV[pos+1]-ABMV[pos])/ABMH
        dist = square(u[k+1]-trueU)+square(v[k+1]-trueV)
        diff.append(dist)


    total = 0
    for d in diff:
        total += d
        
    BEulerY.append(total/x)

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
        title = 'Error',
        type='log',
        autorange = True
    )
)

fig = go.Figure(data=data, layout=layout)

py.plot(fig, filename='difference.html')
