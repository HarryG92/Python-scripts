# -*- coding: utf-8 -*-
"""
Created on Wed Jun 26 13:24:25 2019

Numerical solvers for ordinary differential equations

@author: gulli
"""

import math
import matplotlib.pyplot as plt
import numpy as np


class ode1:

    def __init__(self, f, x0, stepsize):
        self.x0 = x0
        self.f = f
        self.stepsize = stepsize

    def solution(self, t):
        x, dxdt = self.x0, self.f(self.x0, 0)
        time = 0
        while time <= t:
            x += dxdt * self.stepsize
            time += self.stepsize
            dxdt = self.f(x, time)
        return x

    def plot(self, t_max):
        t = np.arange(0, t_max, self.stepsize)
        x, dxdt = self.x0, self.f(self.x0, 0)
        sol = []

        for time in t:
            x += dxdt * self.stepsize
            dxdt = self.f(x, time + self.stepsize)
            sol.append(x)

        plt.plot(t, sol)
        plt.xlabel('t')
        plt.ylabel('x(t)')

        plt.show()


class RK4:

    def __init__(self, f, x0, stepsize):
        self.x0 = x0
        self.f = f
        self.h = stepsize

    def solution(self, t):
        time = 0
        x = self.x0
        while time <= t - self.h/2:
            k1 = self.f(x, time) * self.h
            k2 = self.f(x + k1/2, time + self.h/2) * self.h
            k3 = self.f(x + k2/2, time + self.h/2) * self.h
            k4 = self.f(x + k3, time + self.h) * self.h
            x += (1/6) * (k1 + 2*k2 + 2*k3 + k4)
            time += self.h
        return x

    def plot(self, t_max):
        t = np.arange(0, t_max, self.h)
        sol = []

        for time in t:
            sol.append(self.solution(time))

        plt.plot(t, sol)
        plt.xlabel('t')
        plt.ylabel('x(t)')

        plt.show()
