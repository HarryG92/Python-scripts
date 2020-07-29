# -*- coding: utf-8 -*-
"""
Created on Tue Jun 25 11:18:21 2019

One of my first Python projects; a script to produce animations of harmonic
oscillators with various forcing, resistance, etc.

@author: gulli
"""

# importing required modules
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
from matplotlib.lines import Line2D


# animate x/t graph, position on spring, v/t graph, and F/t graph for SHM
class unforced_harmonic_oscillator(animation.TimedAnimation):
    def __init__(self, mass, resistance, hooke, x0, v0, endtime):
        # all inputs are ints or floats. mass, resistance, and hooke should be
        # >=0 for physical meaning. endtime must be >0
        # x0, v0 are initial position, velocity of particle; animation runs
        # from t=0 to t=endtime

        # define a,b to get x''+ax'+cx=0
        a = resistance/mass
        b = hooke/mass

        # define discriminant
        disc = a**2 - 4*b

        # compute position and velocity solutions for each different
        # possible sign of disc
        if disc > 0:  # overdamped; exponential solutions
            rootdisc = np.sqrt(disc)

            # roots of auxiliary equation
            mplus = (-a + rootdisc)/2
            mneg = (-a - rootdisc)/2

            # coefficients of exp terms
            cplus = (v0 - mneg*x0)/rootdisc
            cneg = (mplus*x0 - v0)/rootdisc

            def position(t):
                return cplus*np.exp(mplus*t) + cneg*np.exp(mneg*t)

            def velocity(t):
                return cplus*mplus*np.exp(mplus*t) + cneg*mneg*np.exp(mneg*t)

        elif disc == 0:  # critical damping; exp*linear solution
            def position(t):
                return (x0 + (v0 + a*x0/2)*t)*np.exp(-a*t/2)

            def velocity(t):
                return (v0 - (a/2)*(v0 + a*x0/2)*t)*np.exp(-a*t/2)

        else:  # underdamped; exp*(sin + cos) solutions
            rootdisc = np.sqrt(4*b - a**2)

            # coefficients of sin and cosine terms
            csin = (2*v0 + a)/rootdisc
            ccos = x0

            def position(t):
                return np.exp(-a*t/2) * (ccos*np.cos(rootdisc*t/2) +
                                         csin*np.sin(rootdisc*t/2))

            def velocity(t):
                return (np.exp(-a*t/2)
                        * (
                          (-a*x0 + csin*rootdisc)*np.cos(rootdisc*t/2)
                          + (-a*csin - x0*rootdisc)*np.sin(rootdisc*t/2))/2)

        # create figure and subplots
        fig = plt.figure()
        ax1 = fig.add_subplot(2, 2, 1)
        ax2 = fig.add_subplot(2, 2, 2)
        ax3 = fig.add_subplot(2, 2, 3)
        ax4 = fig.add_subplot(2, 2, 4)
        fig.subplots_adjust(hspace=.5, wspace=0.5)

        # set time, position, velocity, and force data
        self.t = np.linspace(0, endtime, 400)
        self.x = position(self.t)
        self.v = velocity(self.t)
        self.f = -resistance*self.v - hooke*self.x

        # set axes etc.
        ax1.set_xlabel('time')
        ax1.set_ylabel('position')
        self.xgraph = Line2D([], [], color='black')
        self.xpoint = Line2D(
            [], [], color='red', marker='o', markeredgecolor='r')
        ax1.add_line(self.xgraph)
        ax1.add_line(self.xpoint)
        ax1.set_xlim(0, endtime)
        ax1.set_ylim(np.amin(self.x) - 0.1*abs(np.amin(self.x)),
                     np.amax(self.x) + 0.1*abs(np.amin(self.x)))

        ax2.set_ylabel('position')
        self.particle = Line2D(
            [], [], color='red', marker='o', markeredgecolor='r')
        self.arrow = Line2D(
            [], [], color='blue')
        ax2.add_line(self.particle)
        ax2.add_line(self.arrow)
        ax2.set_xlim(-1, 1)
        ax2.set_ylim(np.amin(self.x) - 0.1*abs(np.amin(self.x)),
                     np.amax(self.x) + 0.1*abs(np.amin(self.x)))

        ax3.set_xlabel('time')
        ax3.set_ylabel('velocity')
        self.vgraph = Line2D([], [], color='black')
        self.vpoint = Line2D(
            [], [], color='red', marker='o', markeredgecolor='r')
        ax3.add_line(self.vgraph)
        ax3.add_line(self.vpoint)
        ax3.set_xlim(0, endtime)
        ax3.set_ylim(np.amin(self.v) - 0.1*abs(np.amin(self.v)),
                     np.amax(self.v) + 0.1*abs(np.amin(self.v)))

        ax4.set_xlabel('time')
        ax4.set_ylabel('force')
        self.fgraph = Line2D([], [], color='black')
        self.fpoint = Line2D(
            [], [], color='red', marker='o', markeredgecolor='r')
        ax4.add_line(self.fgraph)
        ax4.add_line(self.fpoint)
        ax4.set_xlim(0, endtime)
        ax4.set_ylim(np.amin(self.f) - 0.1*abs(np.amin(self.f)),
                     np.amax(self.f) + 0.1*abs(np.amin(self.f)))

        # call animation
        animation.TimedAnimation.__init__(self, fig, interval=20, blit=True)

    # how to draw each frame
    def _draw_frame(self, framedata):
        i = framedata  # index of frame
        head = i - 1

        # plot entire position graph, and moving point at current position
        self.xgraph.set_data(self.t, self.x)
        self.xpoint.set_data(self.t[head], self.x[head])

        # plot bouncing ball and force line
        self.particle.set_data(0, self.x[head])
        self.arrow.set_data([0, 0], [self.x[head],
                                     self.x[head] + self.f[head]])

        # plot velocity graph and moving point at current position
        self.vgraph.set_data(self.t, self.v)
        self.vpoint.set_data(self.t[head], self.v[head])

        # force graph and moving point
        self.fgraph.set_data(self.t, self.f)
        self.fpoint.set_data(self.t[head], self.f[head])

        # black box
        self._drawn_artists = [self.xgraph, self.xpoint,
                               self.particle, self.arrow,
                               self.vgraph, self.vpoint,
                               self.fgraph, self.fpoint]

    # sequence of frames
    def new_frame_seq(self):
        return iter(range(self.t.size))

    # initially, set all data to be plotted to empty
    def _init_draw(self):
        lines = [self.xgraph, self.xpoint, self.particle, self.arrow,
                 self.vgraph, self.vpoint, self.fgraph, self.fpoint]
        for l in lines:
            l.set_data([], [])


# as above, but with a constant -9.8m force applied to particle
class gravity_forced_harmonic_oscillator(animation.TimedAnimation):
    def __init__(self, mass, resistance, hooke, x0, v0, endtime):
        # inputs as above

        # define a,b to get x''+ax'+cx=0
        a = resistance/mass
        b = hooke/mass
        g = 9.8

        # discriminant
        disc = a**2 - 4*b

        # define solutions
        if disc > 0:  # overdamped
            rootdisc = np.sqrt(a**2 - 4*b)

            # roots of auxiliary equation
            mplus = (-a + rootdisc)/2
            mneg = (-a - rootdisc)/2

            # coefficients of exp terms
            cplus = (v0 - mneg*(x0 + g/b))/rootdisc
            cneg = (mplus*(x0 + g/b) - v0)/rootdisc

            def position(t):
                return cplus*np.exp(mplus*t) + cneg*np.exp(mneg*t) - g/b

            def velocity(t):
                return cplus*mplus*np.exp(mplus*t) + cneg*mneg*np.exp(mneg*t)

        elif disc == 0:  # critical damping
            def position(t):
                return ((x0 + g/b)*np.exp(-a*t/2)
                        + (v0 + a*(x0 + g/b)/2)*t*np.exp(-a*t/2)
                        - g/b)

            def velocity(t):
                return (v0 - (a/2)*(v0 + a*(x0 + g/b)/2)*t)*np.exp(-a*t/2)

        else:  # underdamped
            rootdisc = np.sqrt(4*b - a**2)

            # trig term coefficients
            ccos = x0 + g/b
            csin = (2*v0 + a*(x0 + g/b))/rootdisc

            def position(t):
                return np.exp(-a*t/2) * (ccos*np.cos(rootdisc*t/2) +
                                         csin*np.sin(rootdisc*t/2)) - g/b

            def velocity(t):
                return (np.exp(-a*t/2)
                        * (
                          (-a*ccos + csin*rootdisc)*np.cos(rootdisc*t/2)
                          + (-a*csin - ccos*rootdisc)*np.sin(rootdisc*t/2))/2)

        # create figure and subplots
        fig = plt.figure()
        ax1 = fig.add_subplot(2, 2, 1)
        ax2 = fig.add_subplot(2, 2, 2)
        ax3 = fig.add_subplot(2, 2, 3)
        ax4 = fig.add_subplot(2, 2, 4)
        fig.subplots_adjust(hspace=.5, wspace=0.5)

        # set time, position, velocity, and force data
        self.t = np.linspace(0, endtime, 400)
        self.x = position(self.t)
        self.v = velocity(self.t)
        self.f = -resistance*self.v - hooke*self.x - mass*g

        # format axes
        ax1.set_xlabel('time')
        ax1.set_ylabel('position')
        self.xgraph = Line2D([], [], color='black')
        self.xpoint = Line2D(
            [], [], color='red', marker='o', markeredgecolor='r')
        ax1.add_line(self.xgraph)
        ax1.add_line(self.xpoint)
        ax1.set_xlim(0, endtime)
        ax1.set_ylim(np.amin(self.x) - 0.1*abs(np.amin(self.x)),
                     np.amax(self.x) + 0.1*abs(np.amin(self.x)))

        ax2.set_ylabel('position')
        self.particle = Line2D(
            [], [], color='red', marker='o', markeredgecolor='r')
        self.arrow = Line2D(
            [], [], color='blue')
        ax2.add_line(self.particle)
        ax2.add_line(self.arrow)
        ax2.set_xlim(-1, 1)
        ax2.set_ylim(np.amin(self.x) - 0.1*abs(np.amin(self.x)),
                     np.amax(self.x) + 0.1*abs(np.amin(self.x)) + 0.1)

        ax3.set_xlabel('time')
        ax3.set_ylabel('velocity')
        self.vgraph = Line2D([], [], color='black')
        self.vpoint = Line2D(
            [], [], color='red', marker='o', markeredgecolor='r')
        ax3.add_line(self.vgraph)
        ax3.add_line(self.vpoint)
        ax3.set_xlim(0, endtime)
        ax3.set_ylim(np.amin(self.v) - 0.1*abs(np.amin(self.v)),
                     np.amax(self.v) + 0.1*abs(np.amin(self.v)))

        ax4.set_xlabel('time')
        ax4.set_ylabel('force')
        self.fgraph = Line2D([], [], color='black')
        self.fpoint = Line2D(
            [], [], color='red', marker='o', markeredgecolor='r')
        ax4.add_line(self.fgraph)
        ax4.add_line(self.fpoint)
        ax4.set_xlim(0, endtime)
        ax4.set_ylim(np.amin(self.f) - 0.1*abs(np.amin(self.f)),
                     np.amax(self.f) + 0.1*abs(np.amin(self.f)))

        # call animation
        animation.TimedAnimation.__init__(self, fig, interval=25, blit=True)

    # how to draw each frame
    def _draw_frame(self, framedata):
        i = framedata  # index of frame
        head = i - 1

        # set data to be plotted
        self.xgraph.set_data(self.t, self.x)
        self.xpoint.set_data(self.t[head], self.x[head])

        self.particle.set_data(0, self.x[head])
        self.arrow.set_data([0, 0], [0, self.x[head]])

        self.vgraph.set_data(self.t, self.v)
        self.vpoint.set_data(self.t[head], self.v[head])

        self.fgraph.set_data(self.t, self.f)
        self.fpoint.set_data(self.t[head], self.f[head])

        self._drawn_artists = [self.xgraph, self.xpoint,
                               self.particle, self.arrow,
                               self.vgraph, self.vpoint,
                               self.fgraph, self.fpoint]

    # sequence of frames
    def new_frame_seq(self):
        return iter(range(self.t.size))

    # initially plot empty data
    def _init_draw(self):
        lines = [self.xgraph, self.xpoint, self.particle, self.arrow,
                 self.vgraph, self.vpoint, self.fgraph, self.fpoint]
        for l in lines:
            l.set_data([], [])


# unforced SHM with varying resistance and spring constant
class harmonic_oscillator_varying_coefficients(animation.TimedAnimation):
    def __init__(self, mass, resistance, hooke, x0, v0, endtime, start_param,
                 end_param):
        # resistance and hooke are functions of a parameter returning ints or
        # floats; all other inputs are ints or floats. Must have endtime > 0,
        # start_param < end_param. For physical solutions, need mass > 0 and
        # resistance, hooke > 0 for all parameters in range. Plots a single
        # x/t graph for t=0 to t=endtime, starting with resistance and hooke
        # determined by start_param, then evolves them through to end_param

        # define a,b to get x''+ax'+cx=0
        def a(param):
            return resistance(param)/mass

        def b(param):
            return hooke(param)/mass

        # define discriminant
        def disc(param):
            return a(param)**2 - 4*b(param)

        # define solution functions for positive, zero, and negative disc
        def position_pos_disc(a, rootdisc, t):
            # roots of auxiliary equation
            mplus = (-a + rootdisc)/2
            mneg = (-a - rootdisc)/2

            # coefficients of exp terms
            cplus = (v0 - mneg*x0)/rootdisc
            cneg = (mplus*x0 - v0)/rootdisc

            return cplus*np.exp(mplus*t) + cneg*np.exp(mneg*t)

        def position_zero_disc(a, t):
            return (x0 + (v0 + a*x0/2)*t)*np.exp(-a*t/2)

        def position_neg_disc(a, rootdisc, t):
            # coefficients of trig terms
            csin = (2*v0 + a)/rootdisc
            ccos = x0

            return np.exp(-a*t/2) * (ccos*np.cos(rootdisc*t/2) +
                                     csin*np.sin(rootdisc*t/2))

        # create a figure, axis and plot element
        fig = plt.figure()
        ax = plt.axes(xlim=(0, endtime), ylim=(-2, 2))

        # set time data and initial position data
        tdata = np.arange(0, endtime, 0.001)
        A = a(start_param)
        if disc(start_param) > 0:
            rootdisc = np.sqrt(disc(start_param))
            xdata = position_pos_disc(A, rootdisc, tdata)
        elif disc(start_param) == 0:
            xdata = position_zero_disc(A, tdata)
        else:
            rootdisc = np.sqrt(-disc(start_param))
            xdata = position_neg_disc(A, rootdisc, tdata)

        # scale y-axis to match data
        ax.set_ylim(np.amin(xdata) - 0.1*abs(np.amin(xdata)),
                    np.amax(xdata) + 0.1*abs(np.amax(xdata)))
        line, = ax.plot(tdata, xdata, lw=2)

        # initialization function
        def init():
            return line,

        frame_number = 500
        parameters = np.linspace(start_param, end_param, frame_number)

        # animation function
        def animate(i):
            # s is the parameter value in frame i
            s = parameters[i]

            # calculate position data for this parameter value
            A = a(s)
            if disc(s) > 0:  # overdamped
                rootdisc = np.sqrt(disc(s))
                xdata = position_pos_disc(A, rootdisc, tdata)
            elif disc(s) == 0:  # critically damped
                xdata = position_zero_disc(A, tdata)
            else:  # underdamped
                rootdisc = np.sqrt(-disc(s))
                xdata = position_neg_disc(A, rootdisc, tdata)

            # rescale y-axis as needed and update position data
            ax.set_ylim(np.amin(xdata) - 0.1*abs(np.amin(xdata)),
                        np.amax(xdata) + 0.1*abs(np.amax(xdata)))
            line.set_ydata(xdata)

            # return line object
            return line,

        # call the animator
        anim = animation.FuncAnimation(fig, animate, init_func=init,
                                       frames=frame_number, interval=20,
                                       blit=True)

        # save the animation as mp4 video file
        anim.save('harmonic oscillator with varying coefficients.mp4',
                  writer='ffmpeg', fps=30)

        # show the final frame
        plt.show()


# plots x/t, v/t, F/t graphs and particle on spring for  harmonic motion with
# a cosine forcing term. I couldn't be bothered to add a phase variable to
# the forcing; feel free to add that if you feel so inclined!
class forced_harmonic_oscillator(animation.TimedAnimation):
    def __init__(self, mass, resistance, hooke, x0, v0, endtime,
                 forcing_frequency, forcing_amplitude):
        # all inputs are ints or floats. Must have endtime>0. Physical
        # solutions for mass, resistance, hooke >0

        # define a,b to get x''+ax'+cx=0
        a = resistance/mass
        b = hooke/mass

        # name forcing coefficients for convenience
        omega = forcing_frequency
        R = forcing_amplitude

        # compute discriminant
        disc = a**2 - 4*b

        # compute coefficients of particular integral (except special case)
        if not (disc < 0 and a == 0 and omega**2 == b):
            det = b**2 - b*(omega**2) + (a**2)*(omega**2) + a*(omega**3)
            cos_term = b*R/det
            sin_term = ((omega + a)*omega*R)/det

        # compute solution for different values of discriminant
        if disc > 0:  # overdamped
            rootdisc = np.sqrt(a**2 - 4*b)

            # roots of auxiliary equation
            mplus = (-a + rootdisc)/2
            mneg = (-a - rootdisc)/2

            # coefficients of exp terms
            cplus = ((v0 - sin_term*omega) - mneg*(x0 - cos_term))/rootdisc
            cneg = (mplus*(x0 - cos_term) - (v0 - sin_term*omega))/rootdisc

            def position(t):
                return (cplus*np.exp(mplus*t) + cneg*np.exp(mneg*t)
                        + cos_term*np.cos(omega*t) + sin_term*np.sin(omega*t))

            def velocity(t):
                return (cplus*mplus*np.exp(mplus*t) + cneg*mneg*np.exp(mneg*t)
                        + omega * (
                                   sin_term * np.cos(omega*t)
                                   - cos_term * np.sin(omega*t)))

        elif disc == 0:  # critically damped
            exp_term = x0 - cos_term  # exp(t) term
            texp_term = v0 - sin_term*omega + a*(x0 - cos_term)/2  # t*exp(t)

            def position(t):
                return ((exp_term + texp_term*t)*np.exp(-a*t/2)
                        + cos_term*np.cos(omega*t) + sin_term*np.sin(omega*t))

            def velocity(t):
                return ((texp_term - a*exp_term/2 - a*texp_term*t/2)
                        * np.exp(-a*t/2)
                        + omega * (
                                   sin_term * np.cos(omega*t)
                                   - cos_term * np.sin(omega*t)))

        else:  # underdamped
            rootdisc = np.sqrt(4*b - a**2)

            # coefficients of trig terms
            csin = (2*v0 + a*(x0 - cos_term) - 2*sin_term*omega)/rootdisc
            ccos = x0 - cos_term
            if not (a == 0 and omega**2 == b):  # not in special case
                def position(t):
                    return (np.exp(-a*t/2) * (ccos*np.cos(rootdisc*t/2) +
                                              csin*np.sin(rootdisc*t/2))
                            + cos_term * np.cos(omega*t)
                            + sin_term * np.sin(omega*t))

                def velocity(t):
                    return ((np.exp(-a*t/2)/2
                             * (
                                (-a*ccos + csin*rootdisc)*np.cos(rootdisc*t/2)
                                 - ((a*csin + ccos*rootdisc)
                                    * np.sin(rootdisc*t/2))))
                            + omega * (
                                       sin_term * np.cos(omega*t)
                                       - cos_term * np.sin(omega*t)))
            # special case where particular integral clashes with general soln;
            # i.e., both are non-decaying trig terms with same frequency, so
            # need to take linear*trig for particular integral
            # physically, this occurs with resonance and zero damping
            else:
                def position(t):
                    return (
                            x0 * np.cos(rootdisc*t/2)
                            + (v0/omega) * np.sin(rootdisc*t/2)
                            + (R/(2*omega)) * t * np.sin(omega*t))

                def velocity(t):
                    return ((v0*rootdisc/(2*omega)) * np.cos(rootdisc*t/2)
                            - (x0*rootdisc/2) * np.sin(rootdisc*t/2)
                            + (
                               (R/(2*omega)) * np.sin(omega*t)
                               + (R/2) * t * np.cos(omega*t)))

        # set up figure and subplots
        fig = plt.figure()
        ax1 = fig.add_subplot(2, 2, 1)
        ax2 = fig.add_subplot(2, 2, 2)
        ax3 = fig.add_subplot(2, 2, 3)
        ax4 = fig.add_subplot(2, 2, 4)
        fig.subplots_adjust(hspace=.5, wspace=0.5)

        # set up data to be plotted
        self.t = np.linspace(0, endtime, 400)
        self.x = position(self.t)
        self.v = velocity(self.t)
        self.f = -resistance*self.v - hooke*self.x + R*np.cos(omega*self.t)

        # set up axes and labels
        ax1.set_xlabel('time')
        ax1.set_ylabel('position')
        self.xgraph = Line2D([], [], color='black')
        self.xpoint = Line2D(
            [], [], color='red', marker='o', markeredgecolor='r')
        ax1.add_line(self.xgraph)
        ax1.add_line(self.xpoint)
        ax1.set_xlim(0, endtime)
        ax1.set_ylim(np.amin(self.x) - 0.1*abs(np.amin(self.x)),
                     np.amax(self.x) + 0.1*abs(np.amin(self.x)))

        ax2.set_ylabel('position')
        self.particle = Line2D(
            [], [], color='red', marker='o', markeredgecolor='r')
        self.arrow = Line2D(
            [], [], color='blue')
        ax2.add_line(self.particle)
        ax2.add_line(self.arrow)
        ax2.set_xlim(-1, 1)
        ax2.set_ylim(np.amin(self.x) - 0.1*abs(np.amin(self.x)),
                     np.amax(self.x) + 0.1*abs(np.amin(self.x)))

        ax3.set_xlabel('time')
        ax3.set_ylabel('velocity')
        self.vgraph = Line2D([], [], color='black')
        self.vpoint = Line2D(
            [], [], color='red', marker='o', markeredgecolor='r')
        ax3.add_line(self.vgraph)
        ax3.add_line(self.vpoint)
        ax3.set_xlim(0, endtime)
        ax3.set_ylim(np.amin(self.v) - 0.1*abs(np.amin(self.v)),
                     np.amax(self.v) + 0.1*abs(np.amin(self.v)))

        ax4.set_xlabel('time')
        ax4.set_ylabel('force')
        self.fgraph = Line2D([], [], color='black')
        self.fpoint = Line2D(
            [], [], color='red', marker='o', markeredgecolor='r')
        ax4.add_line(self.fgraph)
        ax4.add_line(self.fpoint)
        ax4.set_xlim(0, endtime)
        ax4.set_ylim(np.amin(self.f) - 0.1*abs(np.amin(self.f)),
                     np.amax(self.f) + 0.1*abs(np.amin(self.f)))

        animation.TimedAnimation.__init__(self, fig, interval=25, blit=True)

    # how to draw each frame
    def _draw_frame(self, framedata):
        i = framedata  # frame index
        head = i - 1

        self.xgraph.set_data(self.t, self.x)
        self.xpoint.set_data(self.t[head], self.x[head])

        self.particle.set_data(0, self.x[head])
        self.arrow.set_data([0, 0], [0,
                                     self.x[head]])

        self.vgraph.set_data(self.t, self.v)
        self.vpoint.set_data(self.t[head], self.v[head])

        self.fgraph.set_data(self.t, self.f)
        self.fpoint.set_data(self.t[head], self.f[head])

        self._drawn_artists = [self.xgraph, self.xpoint,
                               self.particle, self.arrow,
                               self.vgraph, self.vpoint,
                               self.fgraph, self.fpoint]

    # sequence of frames
    def new_frame_seq(self):
        return iter(range(self.t.size))

    # start with empty data to plot
    def _init_draw(self):
        lines = [self.xgraph, self.xpoint, self.particle, self.arrow,
                 self.vgraph, self.vpoint, self.fgraph, self.fpoint]
        for l in lines:
            l.set_data([], [])


# function to compute resonant frequency, and both print and return it
def resonance_finder(mass, resistance, hooke):
    # all inputs should be ints or floats
    a = resistance/mass
    b = hooke/mass
    if a**2 < 2*b:
        res_freq = np.sqrt(4*b - a**2)/2
        print('The resonant frequency is %f' % res_freq)
        return res_freq
    else:
        print('Too much damping for resonance to occur')


# plot x/t graph for a harmonic oscillator with fixed parameters but varying
# frequency and amplitude of cosine forcing term
class harmonic_oscillator_varying_forcing(animation.TimedAnimation):
    def __init__(self, mass, resistance, hooke, amplitude, frequency, x0, v0,
                 endtime, start_param, end_param):
        # amplitude and frequency are functions of a single parameter returning
        # floats or ints. The rest are floats or ints. Must have endtime>0,
        # start_param < end_param

        # define a,b to get x''+ax'+cx=0
        a, b = resistance/mass, hooke/mass

        # define discriminant
        disc = a**2 - 4*b
        if disc >= 0:
            rootdisc = np.sqrt(disc)
        else:
            rootdisc = np.sqrt(-disc)

        # define solution functions for positive, zero, and negative disc
        def position_pos_disc(param, t):
            # first find particular integral
            omega = frequency(param)
            A = amplitude(param)
            denom = (b - omega**2)**2 + (a*omega)**2  # denominator
            cos_coeff = A*(b - omega**2)/denom
            sin_coeff = (a*omega*A)/denom
            particular_integral = (cos_coeff*np.cos(omega*t)
                                   + sin_coeff*np.sin(omega*t))

            # account for effect of particular integral on coefficients of
            # complementary function
            x_new = x0 - cos_coeff
            v_new = v0 - sin_coeff*omega

            # roots of auxiliary equation
            mplus = (-a + rootdisc)/2
            mneg = (-a - rootdisc)/2

            # coefficients of exp terms
            cplus = (v_new - mneg*x_new)/rootdisc
            cneg = (mplus*x_new - v_new)/rootdisc

            complementary_function = (cplus*np.exp(mplus*t)
                                      + cneg*np.exp(mneg*t)
                                      + particular_integral(param, t))

            return complementary_function + particular_integral

        def position_zero_disc(param, t):
            # first find particular integral
            omega = frequency(param)
            A = amplitude(param)
            denom = (b - omega**2)**2 + (a*omega)**2
            cos_coeff = A*(b - omega**2)/denom
            sin_coeff = (a*omega*A)/denom
            particular_integral = (cos_coeff*np.cos(omega*t)
                                   + sin_coeff*np.sin(omega*t))

            # account for effect of particular integral on coefficients of
            # complementary function
            x_new = x0 - cos_coeff
            v_new = v0 - sin_coeff*omega

            complementary_function = ((x_new + (v_new + a*x_new/2)*t)
                                      * np.exp(-a*t/2))

            return complementary_function + particular_integral

        def position_neg_disc(param, t):
            if a != 0:  # not in special case
                # first find particular integral
                omega = frequency(param)
                A = amplitude(param)
                denom = (b - omega**2)**2 + (a*omega)**2
                cos_coeff = A*(b - omega**2)/denom
                sin_coeff = (a*omega*A)/denom
                particular_integral = (cos_coeff*np.cos(omega*t)
                                       + sin_coeff*np.sin(omega*t))

                # account for effect of particular integral on coefficients of
                # complementary function
                x_new = x0 - cos_coeff
                v_new = v0 - sin_coeff*omega

                # find complementary function
                csin = (2*v_new + a)/rootdisc
                complementary_function = (np.exp(-a*t/2)
                                          * (x_new*np.cos(rootdisc*t/2)
                                             + csin*np.sin(rootdisc*t/2)))
                return complementary_function + particular_integral
            else:  # special case of no resistance and resonance
                # first find particular integral of linear*trig form
                omega = frequency(param)
                A = amplitude(param)
                particular_integral = (A/(2*omega))*t*np.sin(omega*t)

                # then complementary function
                csin = (2*v0 + a)/rootdisc
                complementary_function = (np.exp(-a*t/2)
                                          * (x0*np.cos(rootdisc*t/2)
                                             + csin*np.sin(rootdisc*t/2)))
                return complementary_function + particular_integral

        # create a figure, axis and plot element
        fig = plt.figure()
        ax = plt.axes(xlim=(0, endtime), ylim=(-2, 2))

        # set time data and initial position data
        tdata = np.arange(0, endtime, 0.001)
        if disc > 0:
            xdata = position_pos_disc(start_param, tdata)
        elif disc == 0:
            xdata = position_zero_disc(start_param, tdata)
        else:
            xdata = position_neg_disc(start_param, tdata)

        # scale y-axis to match data
        ax.set_ylim(np.amin(xdata) - 0.1*abs(np.amin(xdata)),
                    np.amax(xdata) + 0.1*abs(np.amax(xdata)))
        line, = ax.plot(tdata, xdata, lw=2)

        # initialization function
        def init():
            return line,

        frame_number = 500
        parameters = np.linspace(start_param, end_param, frame_number)

        # animation function
        if disc > 0:  # overdamped
            def animate(i):
                # s is the value of the parameter in frame i
                s = parameters[i]

                # calculate position data for this parameter value
                xdata = position_pos_disc(s, tdata)

                # rescale y-axis as needed and update position data
                ax.set_ylim(np.amin(xdata) - 0.1*abs(np.amin(xdata)),
                            np.amax(xdata) + 0.1*abs(np.amax(xdata)))
                line.set_ydata(xdata)

            # return line object
            return line,
        elif disc == 0:  # critically damped
            def animate(i):
                # s is a parameter
                s = parameters[i]

                # calculate position data for this parameter value
                xdata = position_zero_disc(s, tdata)

                # rescale y-axis as needed and update position data
                ax.set_ylim(np.amin(xdata) - 0.1*abs(np.amin(xdata)),
                            np.amax(xdata) + 0.1*abs(np.amax(xdata)))
                line.set_ydata(xdata)

            # return line object
            return line,
        else:  # underdamped
            def animate(i):
                # s is a parameter
                s = parameters[i]

                # calculate position data for this parameter value
                xdata = position_neg_disc(s, tdata)

                # rescale y-axis as needed and update position data
                ax.set_ylim(np.amin(xdata) - 0.1*abs(np.amin(xdata)),
                            np.amax(xdata) + 0.1*abs(np.amax(xdata)))
                line.set_ydata(xdata)

                # return line object
                return line,

        # call the animator
        anim = animation.FuncAnimation(fig, animate, init_func=init,
                                       frames=frame_number, interval=20,
                                       blit=True)

        # save the animation as mp4 video file
        anim.save('harmonic oscillator with varying forcing.mp4',
                  writer='ffmpeg', fps=30)

        # show the final frame
        plt.show()
