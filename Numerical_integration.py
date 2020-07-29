# -*- coding: utf-8 -*-
"""
Created on Fri Jul  5 13:54:50 2019

@author: gulli
"""

import numpy

import types


class NumericalIntegrator:
    def __init__(self, f):
        if type(f) == types.FunctionType:
            self.integrand = f
        else:
            raise TypeError("argument must be a function")

    def simpson_integrator(self, lower, upper, n):
        if lower < upper:
            A, B = lower, upper
            sign = 1
        else:
            A, B = upper, lower
            sign = -1

        if n % 2 == 0:
            m = n + 1
        else:
            m = n

        partition = numpy.linspace(A, B, m)
        values = self.integrand(partition)
        width = (B - A)/(m - 1)

        Sum = values[0]
        i = 1
        while i < m - 2:
            Sum += 4*values[i] + 2*values[i+1]
            i += 2
        Sum += 4*values[m-2] + values[m-1]
        return sign * (width/3) * Sum

    def integrate(self, lower, upper, tolerance=0.001, max_step_number=1024):
        step_number = 2  # initially try to integrate in 2 steps
        error = float("inf")
        while error > tolerance:
            first_estimate = self.simpson_integrator(lower, upper, step_number)
            second_estimate = self.simpson_integrator(
                                                      lower,
                                                      upper,
                                                      2 * step_number
                                                      )
            error = abs(second_estimate - first_estimate)
            step_number *= 2
            if step_number > max_step_number:
                print("""Step number exceeded maximum. Try a looser tolerance
                      or a higher maximum step number""")
                return None
        return second_estimate
