# -*- coding: utf-8 -*-
"""
Created on Thu Jan 23 14:05:29 2020

A convolution calculator I made to help explain convolution to a student

@author: gulli
"""

import types

from Numerical_integration import NumericalIntegrator


class Convolver:
    def __init__(self, f, g, lower_limit, upper_limit):
        if type(f) == types.FunctionType:
            self.first_function = f
        else:
            raise TypeError("first argument must be a function")
        if type(g) == types.FunctionType:
            self.second_function = g
        else:
            raise TypeError("second argument must be a function")
        self.lower_limit = lower_limit
        self.upper_limit = upper_limit

    def shiftedProduct(self, t, shift):
        return self.first_function(t) * self.second_function(shift - t)

    def convolution(self, tau, tolerance=0.001, max_step_number=1024):
        def tau_shifted_product(t):
            return self.shiftedProduct(t, tau)
        convolution = NumericalIntegrator(tau_shifted_product)
        return convolution.integrate(
                                     self.lower_limit,
                                     self.upper_limit,
                                     tolerance,
                                     max_step_number)

