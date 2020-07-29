"""
When working on Fourier analysis and filtering with a student, I made this
to illustrate using Fourier series to study the effects of low and high pass
filters; this was just hacked together to make the graphs, so is uncommented...
sorry!
"""


import time

import matplotlib.pyplot as pyplot

import numpy

start_time = time.time()

number_of_terms = 12
samples = 10000
a = [2/3]
b = [0]

for n in range(1, number_of_terms):
    if n % 3 == 0:
        a.append(0)
        b.append(0)
    elif n % 3 == 1:
        frequency = 2 * numpy.pi * n
        a.append(-numpy.sqrt(3) / frequency)
        b.append(-3 / frequency)
    else:
        frequency = 2*numpy.pi * n
        a.append(numpy.sqrt(3) / frequency)
        b.append(-3 / frequency)

t = numpy.linspace(0, 12, samples)


unfiltered = numpy.zeros(samples)
for n in range(number_of_terms):
    frequency = 2 * numpy.pi * n / 3
    unfiltered += ((a[n] * numpy.cos(frequency * t))
                   + (b[n] * numpy.sin(frequency * t)))


low_pass = numpy.zeros(samples)
for n in range(number_of_terms):
    frequency = 2 * numpy.pi * n / 3
    trig_argument = (frequency * t) - numpy.arctan(n / 3)
    low_pass += ((3 * a[n] * numpy.cos(trig_argument) / numpy.sqrt(9 + n**2))
                 + (3 * b[n] * numpy.sin(trig_argument)/numpy.sqrt(9 + n**2)))

high_pass = numpy.zeros(samples)
for n in range(number_of_terms):
    frequency = 2 * numpy.pi * n / 3
    trig_argument = (frequency * t) + (numpy.pi / 2) - numpy.arctan(n / 3)
    high_pass += ((n * a[n] * numpy.cos(trig_argument) / numpy.sqrt(9 + n**2))
                  + (n * b[n] * numpy.sin(trig_argument)/numpy.sqrt(9 + n**2)))


pyplot.plot(t, unfiltered)
pyplot.savefig("Unfiltered12.png")
pyplot.show()
pyplot.plot(t, low_pass)
pyplot.savefig("LowPass12.png")
pyplot.show()
pyplot.plot(t, high_pass)
pyplot.savefig("HighPass12.png")
pyplot.show()

print("--- %s seconds ---" % (time.time() - start_time))
