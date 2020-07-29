# -*- coding: utf-8 -*-
"""
Created on Wed Oct 30 11:02:42 2019

Everyone loves pretty pictures! This code generates the Mandelbrot set;
simplePlot() gives a black-and-white rendering, whereas fullPlot() and
fullPlot2() use black for convergent points, and colour for divergent points;
in fullPlot(), 6 colours are used, whereas in fullPlot2(), a gradient of colour
from blue-green to pink is used; in both cases, the colour indicates how long
it took the iteration to diverge

@author: gulli
"""

import numpy

import pandas

import matplotlib.pyplot as pyplot


class Mandelbrot:
    def __init__(self, iteration_depth, number_of_points):
        self.iteration_depth = iteration_depth
        self.number_of_points = number_of_points
        self.converging_points = []

        diverging_points = []
        divergence_depths = []

        for x in numpy.linspace(-2, 0.6, self.number_of_points):
            for y in numpy.linspace(-1.3, 1.3, self.number_of_points):
                c = complex(x, y)
                current_iterate = c
                step_number = 0
                while (abs(current_iterate) < 2
                       and step_number < self.iteration_depth):
                    current_iterate = current_iterate**2 + c
                    step_number += 1
                if abs(current_iterate) < 2:
                    self.converging_points.append(c)
                else:
                    diverging_points.append(c)
                    divergence_depths.append(step_number)
        self.diverging_points = pandas.DataFrame({'points': diverging_points,
                                                  'depths': divergence_depths})
        self.converging_points = numpy.array(self.converging_points)

    def simplePlot(self):
        pyplot.scatter(self.converging_points.real,
                       self.converging_points.imag, color='black', s=1)
        pyplot.savefig("Mandelbrot set.png")
        pyplot.show()

    def fullPlot(self):
        NUMBER_OF_COLOURS = 6
        divergence_bands = numpy.linspace(-1, self.iteration_depth/100,
                                          NUMBER_OF_COLOURS)
        divergence_bands = numpy.append(divergence_bands,
                                        [self.iteration_depth])
        COLOURS = ['yellow', 'cyan', 'blue', 'green', 'red',
                   'magenta']
        pyplot.axes(aspect='equal')
        for colour in range(NUMBER_OF_COLOURS):
            points_of_colour = self.diverging_points[((
                    self.diverging_points['depths'] > divergence_bands[colour])
                    & (self.diverging_points['depths']
                       <= divergence_bands[colour + 1]))]
            pyplot.scatter(points_of_colour['points'].apply(lambda x: x.real),
                           points_of_colour['points'].apply(lambda x: x.imag),
                           color=COLOURS[colour], s=1)
        self.simplePlot()

    def fullPlot2(self):
        pyplot.figure(num=1, figsize=(5, 5))
        for depth in range(self.iteration_depth + 1):
            red_value = 1 - depth/self.iteration_depth
            points_of_colour = self.diverging_points[(
                    self.diverging_points['depths'] == depth)]
            pyplot.scatter(points_of_colour['points'].apply(lambda x: x.real),
                           points_of_colour['points'].apply(lambda x: x.imag),
                           color=(red_value, 0.5, 0.5), s=1)
        self.simplePlot()


mand = Mandelbrot(100, 1000)
mand.fullPlot2()
