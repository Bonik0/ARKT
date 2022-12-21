import numpy as np
import matplotlib.pyplot as plt

def at(time):
    return np.pi / 2 - 0.000129 * time

def exp(h):
    return 2.71828 ** h

def func(t, Typegraf):
    step = [[115, 2 * 5849411, 2 * 33798, 2 * 226233], [262, 2339760, 5443, 116573], [467, 453714, 2653, 29188], [937, 131222, 2631, 16258]]
    m0 = 615070
    a = (step[0][3] - step[0][2])/step[0][0]
    F = step[0][1]
    g = 9.8
    s = 20 
    c = 0.0045
    B = 0.000129
    p0 = 1.29
    h = [0 for i in range(len(t))]
    Vy = [0 for i in range(len(t))]
    Vx = [0 for i in range(len(t))]
    V = [0 for i in range(len(t))]
    changeA = [0 for i in range(len(t))]
    V = [0 for i in range(len(t))]
    changemass = [0 for i in range(len(t))]
    changemass[0] = m0
    changeA[0] = a
    stupen = 0
    steptime = 0
    for i in range(1,len(t)):
        if t[i - 1] > step[stupen][0]:
            steptime = step[stupen][0]
            m0 -= step[stupen][3]
            stupen += 1
            F = step[stupen][1]
            a = (step[stupen][3] - step[stupen][2]) / (step[stupen][0] - step[stupen - 1][0])
        m = m0 - a * (t[i - 1] - steptime)
        changeA[i] = a
        changemass[i] = m
        atime = at(t[i - 1])
        alfax = (F * np.cos(atime) - m * g - (0.5 * c * p0 * exp(-B * h[i - 1]) * s * (V[i - 1] ** 2)) * np.cos(atime)) / m
        alfay = (F * np.sin(atime) - m * g - (0.5 * c * p0 * exp(-B * h[i - 1]) * s * (V[i - 1] ** 2)) * np.sin(atime)) / m
        Vx[i] = Vx[i - 1] + alfax * 0.1
        Vy[i] = Vy[i - 1] + alfay * 0.1
        V[i] = np.sqrt(Vy[i] ** 2 + Vy[i] ** 2)
        h[i] = h[i - 1] + Vy[i] * 0.01
    if (Typegraf == 'H'):
        return h
    elif (Typegraf == 'V'):
        return V
    elif (Typegraf == 'a'):
        return changeA
    else:
        return changemass

t = np.arange(0.0, 937, 0.1)
plt.figure(figsize = (12, 8))
sp = plt.subplot(221)
res = func(t, 'H')
plt.xlabel(r'$t$')
plt.title(r'$H$')
plt.plot(t, func(t, 'H'))
plt.grid(True)
sp = plt.subplot(222)
plt.xlabel(r'$t$')
plt.title(r'$V$')
plt.plot(t, func(t, 'V'))
plt.grid(True)
sp = plt.subplot(223)
plt.xlabel(r'$t$')
plt.title(r'$a$')
plt.plot(t, func(t, 'a'))
plt.grid(True)
sp = plt.subplot(224)
plt.xlabel(r'$t$')
plt.title(r'$m$')
plt.plot(t, func(t, 'm'))
plt.grid(True)
plt.show()
