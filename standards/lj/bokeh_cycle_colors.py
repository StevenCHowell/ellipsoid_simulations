
import numpy as np
from bokeh.plotting import figure, output_file, show

# select a palette
from bokeh.palettes import Dark2_5 as palette
# itertools handles the cycling
import itertools

output_file('bokeh_cycle_colors.html')

p = figure(width=400, height=400)
x = np.linspace(0, 10)

# create a color iterator
colors = itertools.cycle(palette)

for m, color in itertools.izip(xrange(10), colors):
    y = m * x
    p.line(x, y, legend='m = {}'.format(m), color=color)

p.legend.location='top_left'
show(p)