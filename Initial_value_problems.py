# -*- coding: utf-8 -*-
"""
Created on Tue Jul  9 18:05:41 2019

Tools for numerical solutions of systems of explicit first order ODEs. Note
that a higher order ODE can be cast as a system of first order ODEs by adding#
extra variables for the higher derivatives. E.g., the single 2nd order equation
x'' + f(t,x)x' + g(t, x)x = h(t, x) is equivalent to the system of 1st order
equations y=x', y' + f(t,x)y + g(t,x)x = h(t,x)

@author: gulli
"""

import numpy

import matplotlib.pyplot as pyplot

from matplotlib.lines import Line2D

import matplotlib.animation as animation


class ODESystem:
    """
    Represents a system of first order ODEs. Initialised with a list of
    variable names and a matching list of functions f(t,x). If the lists have
    length one, the arguments can be passed as simply a variable name, and a
    function, and they will be listified automatically.
    """
    def __init__(self, list_of_variables, list_of_functions):
        if isinstance(list_of_functions, list):
            if len(list_of_functions) != len(list_of_variables):
                print("Error. Length mismatch!")
            else:
                self.function_list = list_of_functions
                self.number_of_variables = len(list_of_variables)
                self.variable_names = list_of_variables
        elif callable(list_of_functions):
            self.function_list = [list_of_functions]
            self.number_of_variables = 1
            self.variable_names = [list_of_variables]
        else:
            print("Error. Invalid arguments!")


class InitialValueProblem(ODESystem):
    """
    Represents an initial value problem. Initialised with an ODESystem and a
    list of initial values for each variable
    """
    def __init__(self, ode_system, list_of_initial_values):
        ODESystem.__init__(self, ode_system.variable_names,
                           ode_system.function_list)
        self.ode_system = ode_system
        if (ode_system.number_of_variables == 1
                and isinstance(list_of_initial_values, (int, float))):
            self.initial_values = [float(list_of_initial_values)]
        elif ode_system.number_of_variables == len(list_of_initial_values):
            self.initial_values = list_of_initial_values
        else:
            print("Error! Invalid arguments.")


class EulerMethod(InitialValueProblem):
    """
    Allows an estimated current solution to an IVP to be advanced by one
    timestep using Euler's method (first order Runge-Kutta).
    Takes an InitialValueProblem when initialised.
    """
    def __init__(self, initial_value_problem):
        InitialValueProblem.__init__(self,
                                     initial_value_problem.ode_system,
                                     initial_value_problem.initial_values)

    def solutionStep(self, current_time, current_values, step_size):
        if (isinstance(current_values, (list, numpy.ndarray))
                and len(current_values) == self.number_of_variables):
            change_in_values = []
            for variable in range(self.number_of_variables):
                derivative = self.function_list[variable](current_time,
                                                          current_values)
                change_in_values.append(derivative * step_size)
            change_in_values = numpy.array(change_in_values)
            new_values = numpy.array(current_values) + change_in_values
            return new_values
        else:
            print("Error! Invalid arguments.")


class ClassicalRungeKuttaMethod(InitialValueProblem):
    """
    Allows an estimated current solution to an IVP to be advanced by one
    timestep using the classical (fourth order) Runge-Kutta method.
    Takes an InitialValueProblem when initialised.
    """
    def __init__(self, initial_value_problem):
        InitialValueProblem.__init__(self,
                                     initial_value_problem.ode_system,
                                     initial_value_problem.initial_values)

    def solutionStep(self, current_time, current_values, step_size):
        if (isinstance(current_values, (list, numpy.ndarray))
                and len(current_values) == self.number_of_variables):
            derivatives_first_estimate = []
            for variable in range(self.number_of_variables):
                derivatives_first_estimate.append(
                        self.function_list[variable](current_time,
                                                     current_values))
            change_first_estimate = (step_size
                                     * numpy.array(derivatives_first_estimate))

            derivatives_second_estimate = []
            for variable in range(self.number_of_variables):
                derivatives_second_estimate.append(
                    self.function_list[variable](current_time + step_size/2,
                                                 (current_values
                                                  + change_first_estimate/2)))
            change_second_estimate = (step_size
                                      * numpy.array(
                                              derivatives_second_estimate))

            derivatives_third_estimate = []
            for variable in range(self.number_of_variables):
                derivatives_third_estimate.append(
                    self.function_list[variable](current_time + step_size/2,
                                                 (current_values
                                                  + change_second_estimate/2)))
            change_third_estimate = (step_size
                                     * numpy.array(derivatives_third_estimate))

            derivatives_fourth_estimate = []
            for variable in range(self.number_of_variables):
                derivatives_fourth_estimate.append(
                    self.function_list[variable](current_time + step_size,
                                                 (current_values
                                                  + change_third_estimate)))
            change_fourth_estimate = (step_size
                                      * numpy.array(
                                              derivatives_fourth_estimate))

            change_average_estimate = (change_first_estimate
                                       + 2 * change_second_estimate
                                       + 2 * change_third_estimate
                                       + change_fourth_estimate) / 6
            new_values = numpy.array(current_values) + change_average_estimate
            return new_values
        else:
            print("Error! Invalid arguments.")


class IVPSolver(InitialValueProblem):
    """
    Takes a solution method (e.g., Euler's, Classical R-K) and an optional
    tolerance and gives methods to solve the IVP contained in the solver by
    the given method, using the step-doubling approach, and to plot the
    solution.
    """
    def __init__(self, solution_method, tolerance=0.0001):
        InitialValueProblem.__init__(self,
                                     solution_method.ode_system,
                                     solution_method.initial_values)
        self.method = solution_method
        self.values_list = numpy.array(self.initial_values)
        self.current_values = self.initial_values
        self.times_list = [0]
        self.time = 0
        self.tolerance = tolerance

    def solutionToTime(self, end_time):
        """
        Solves using the specified method and step-doubling. This means that to
        advance from time t to time t+h, it does so by two methods: directly in
        one step of size h, and in two steps of size h/2. If these give answers
        within tolerance, then h is good as a step size and so we can try using
        2h at the next step. If the two estimates differ by more than
        tolerance, we must repeat the step with half the step size. In this
        way, we modify our step size as we solve the ODE to find a good
        compromise between speed (large step size) and accuracy (small step
        size).
        """
        step_size = end_time / 1000
        while self.time <= end_time - step_size:
            new_value_full_step = self.method.solutionStep(
                    self.time, self.current_values, step_size)
            new_value_half_step = self.method.solutionStep(
                    self.time, self.current_values, step_size/2)
            new_value_two_half_steps = self.method.solutionStep(
                    self.time + step_size/2, new_value_half_step, step_size/2)
            error_estimate = max(abs(new_value_full_step
                                     - new_value_two_half_steps)/2)
            if error_estimate < self.tolerance:
                self.current_values = new_value_two_half_steps
                self.values_list = numpy.append(self.values_list,
                                                self.current_values)
                self.time += step_size
                self.times_list.append(self.time)
                if 2*step_size < end_time - self.time:
                    step_size = min(2*step_size, end_time/100)
                elif self.time != end_time:
                    step_size = end_time - self.time
            else:
                step_size *= 0.5
        return self.current_values

    def plotToTime(self, end_time):
        self.solutionToTime(end_time)
        figure = pyplot.figure()
        subplots = []
        for variable in range(self.number_of_variables):
            subplots.append(figure.add_subplot(self.number_of_variables,
                                               1, variable + 1))
        figure.subplots_adjust(hspace=.7, wspace=0.5)

        for variable in range(self.number_of_variables):
            subplots[variable].plot(self.times_list,
                                    self.values_list[variable::
                                                     self.number_of_variables])
            subplots[variable].set_xlabel("time")
            subplots[variable].set_ylabel("%s" % self.variable_names[variable])

        pyplot.show()


class SolutionPhaseSpaceAnimator(IVPSolver, animation.TimedAnimation):
    """
    Produces an animation of the solutions of two specified variables of the
    ODE system against time and against each other.
    """
    def __init__(self, ivp_solver, end_time, file_name,
                 xlabel="time", indices=[0, 1], step_number=400):
        if not (isinstance(indices, list) and len(indices) == 2):
            print("Error! Must specify two variable indices to plot together")
            return None
        IVPSolver.__init__(self, ivp_solver.method,
                           ivp_solver.tolerance)

        self.solutionToTime(end_time)
        self.first_index = indices[0]
        self.second_index = indices[1]

        # set up figure and subplots
        fig = pyplot.figure()
        ax1 = fig.add_subplot(2, 2, 1)
        ax2 = fig.add_subplot(2, 2, 2)
        ax3 = fig.add_subplot(2, 1, 2)
        fig.subplots_adjust(hspace=.5, wspace=0.5)

        self.first_variable_values = self.values_list[self.first_index::
                                                      self.number_of_variables]
        first_variable_range = [numpy.amin(self.first_variable_values),
                                numpy.amax(self.first_variable_values)]
        self.second_variable_values = self.values_list[
                self.second_index::self.number_of_variables]
        second_variable_range = [numpy.amin(self.second_variable_values),
                                 numpy.amax(self.second_variable_values)]

        self.first_graph = Line2D([], [], color='black')
        self.first_point = Line2D(
                [], [], color='red', marker='o', markeredgecolor='r')
        ax1.add_line(self.first_graph)
        ax1.add_line(self.first_point)
        ax1.set_xlabel(xlabel)
        ax1.set_xlim(0, self.time)
        ax1.set_ylabel("%s" % self.variable_names[self.first_index])
        ax1.set_ylim(
                first_variable_range[0] - 0.1*abs(first_variable_range[0]),
                first_variable_range[1] + 0.1*abs(first_variable_range[1]))

        self.second_graph = Line2D([], [], color='black')
        self.second_point = Line2D(
                [], [], color='red', marker='o', markeredgecolor='r')
        ax2.add_line(self.second_graph)
        ax2.add_line(self.second_point)
        ax2.set_xlabel(xlabel)
        ax2.set_xlim(0, self.time)
        ax2.set_ylabel("%s" % self.variable_names[self.second_index])
        ax2.set_ylim(
                second_variable_range[0] - 0.1*abs(second_variable_range[0]),
                second_variable_range[1] + 0.1*abs(second_variable_range[1]))

        self.third_graph = Line2D([], [], color='black')
        self.third_point = Line2D(
                [], [], color='red', marker='o', markeredgecolor='r')
        ax3.add_line(self.third_graph)
        ax3.add_line(self.third_point)
        ax3.set_xlabel("%s" % self.variable_names[self.first_index])
        ax3.set_xlim(
                first_variable_range[0] - 0.1*abs(first_variable_range[0]),
                first_variable_range[1] + 0.1*abs(first_variable_range[1]))
        ax3.set_ylabel("%s" % self.variable_names[self.second_index])
        ax3.set_ylim(
                second_variable_range[0] - 0.1*abs(second_variable_range[0]),
                second_variable_range[1] + 0.1*abs(second_variable_range[1]))
        ax3.set_aspect(aspect=1)

        animation.TimedAnimation.__init__(self, fig, interval=25, blit=True)
        self.save(file_name, writer='ffmpeg', fps=30)

    def _draw_frame(self, frame_index):
        head = frame_index - 1  # frame_index counts from 1, but lists are
                                # 0-indexed, so need to offset
        self.first_graph.set_data(self.times_list, self.first_variable_values)
        self.first_point.set_data(self.times_list[head],
                                  self.first_variable_values[head])
        self.second_graph.set_data(self.times_list,
                                   self.second_variable_values)
        self.second_point.set_data(self.times_list[head],
                                   self.second_variable_values[head])
        self.third_graph.set_data(self.first_variable_values[0:head],
                                  self.second_variable_values[0:head])
        self.third_point.set_data(self.first_variable_values[head],
                                  self.second_variable_values[head])

        self._drawn_artists = [self.first_graph, self.first_point,
                               self.second_graph, self.second_point,
                               self.third_graph, self.third_point]

    # sequence of frames
    def new_frame_seq(self):
        return iter(range(len(self.times_list)))

    # start with empty data to plot
    def _init_draw(self):
        lines = [self.first_graph, self.first_point, self.second_graph,
                 self.second_point, self.third_graph, self.third_point]
        for l in lines:
            l.set_data([], [])
