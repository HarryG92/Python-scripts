# -*- coding: utf-8 -*-
"""
Created on Tue Jul 16 13:52:23 2019

@author: gulli
"""

from Initial_value_problems import *

import ODE_solvers

import time

import numpy


"""
example comparing Euler and Runge-Kutta methods
the system is dx/dt = x, x(0) = 1, so solution is x=exp(t)
"""

def f(t, x):
    return x[0]

start_time = time.time()
ode = InitialValueProblem("$e^t$", f, 1)
method = EulerMethod(ode)
solver = IVPSolver(method, 0.00001) # tight tolerance on error
print("Euler gives exp(1) = %f" % solver.solutionToTime(1)) # solution should be e
solver.plotToTime(10)
print("Euler takes %s seconds" % (time.time() - start_time))

start_time = time.time()
method = ClassicalRungeKuttaMethod(ode)
solver = IVPSolver(method, 0.1) # much looser tolerance
print("Runge Kutta gives exp(1) = %f" % solver.solutionToTime(1)) # solution as above
solver.plotToTime(10)
print("Runge Kutta takes %s seconds" % (time.time() - start_time))


"""
harmonic motion example. Think of x as the vector [position, velocity]
the single 2nd order ODE x'' + ax' + bx = 0 can be expressed as a system of 2
1st order ODEs: x' = v, v' = -av - bx; so our 1st order solver can solve it
"""

hooke = 1
resistance = 1

def velocity(t, x):
    return x[1] # derivative of position is velocity


def acceleration(t, x):
    return - resistance * x[1] - hooke * x[0] # acceleration in terms of velocity and position

simple_harmonic_motion = InitialValueProblem(["position", "velocity"],
                                             [velocity, acceleration],
                                             [1, 0]) # initial conditions: position 1, velocity 0
method = ClassicalRungeKuttaMethod(simple_harmonic_motion) # default tolerance is 0.0001
solver = IVPSolver(method)
solver.plotToTime(10)
animator = PhaseSpaceAnimator(solver, 30, "Harmonic Motion Phase Plot.mp4")



"""
Lotka-Volterra example. Think of an island with rabbits and foxes, and x is
the vector [rabbits, foxes] of populations
"""

rabbit_growth_rate = 5 # growth rate of rabbit population ignoring predation
predation_rate = 0.1 # number of rabbits a single fox eats per unit time
fox_death_rate = 1 # assuming no rabbits, hence no food
fox_growth_rate = 0.01 # rate at which new foxes are born per rabbit available as food

"""
the Lotka-Volterra model is:
d(rabbits)/dt = rabbit_growth_rate * rabbits - predation_rate*rabbits*foxes,
d(foxes)/dt = fox_growth_rate*rabbits*foxes - fox_death_rate*foxes
"""

def rabbit_derivative(t, x):
    return rabbit_growth_rate * x[0] - predation_rate * x[0] * x[1]

def fox_derivative(t, x):
    return fox_growth_rate * x[0] * x[1] - fox_death_rate * x[1]

lotka_volterra = InitialValueProblem(["rabbits", "foxes"],
                                     [rabbit_derivative, fox_derivative],
                                     [100, 10])
method = ClassicalRungeKuttaMethod(lotka_volterra)
solver = IVPSolver(method)
solver.plotToTime(10)
animator = PhaseSpaceAnimator(solver, 10, "Lotka Volterra Phase Plot.mp4")
