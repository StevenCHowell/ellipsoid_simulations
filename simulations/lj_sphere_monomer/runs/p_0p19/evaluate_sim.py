#!/usr/bin/env python
#coding:utf-8
"""
  Author:  Steven C. Howell --<schowell@gwmail.gwu.edu>
  Purpose: evaluate the Lennard Jones NPT results
  Created: 09/16/2016

00000000011111111112222222222333333333344444444445555555555666666666677777777778
12345678901234567890123456789012345678901234567890123456789012345678901234567890
"""

import logging
import os.path as op
import numpy as np
import pandas as pd
import sasmol.sasmol as sasmol

from bokeh.layouts import gridplot
from bokeh.plotting import figure, output_file, show

logging.basicConfig(format=':', level=logging.DEBUG)

""" INPUT HERE """
run_dir = 'run_1'
pdb_file = "run_2.pdb"
goal_density = 0.0213 # atoms/A^3
epsilon = 119.8 * 1.3806505e-23  # J
""" INPUT HERE """

box_mol = sasmol.SasMol(0)
box_mol.read_pdb(op.join(run_dir, pdb_file))
n_atoms = box_mol.natoms()

# step, box_length, density, pressure
box_data = np.loadtxt(op.join(run_dir, 'box_length.txt'))
steps = np.array(box_data[:, 0], dtype=int)
sigma_Ar = 3.405
sigma_Ar_m = sigma_Ar * 1e-10  # convert to meters
n_mean = 10

box_length = box_data[:, 1] * sigma_Ar
box_length_rolling_mean = pd.rolling_mean(box_length, n_mean)

density = n_atoms / box_length**3
density_rolling_mean = pd.rolling_mean(density, n_mean)
density_mc = box_data[:, 2] / sigma_Ar ** 3

pressure_mc = box_data[:, 3] * epsilon / sigma_Ar ** 3  # in Pascals
pressure_mc_rolling_mean = pd.rolling_mean(pressure_mc, n_mean)


# setup the plots
output_file("box_size.html")
size_plot = figure(title="Box Size", x_axis_label='MC step', width=300,
                   height=300, y_axis_label='Side Length (A)')
den_plot = figure(title="Density", x_axis_label='MC step', width=300,
                  height=300, y_axis_label='Density (atoms/A^3)',
                  x_range=size_plot.x_range)
pres_plot = figure(title="Pressure", x_axis_label='MC step', width=300,
                   height=300, y_axis_label='Pressure (Pa)',
                   x_range=size_plot.x_range)
unitless_pres_plot = figure(title="Pressure", x_axis_label='MC step', width=300,
                            height=300, y_axis_label='Pressure (unitless)',
                            x_range=size_plot.x_range)


# plot the data
size_plot.line(steps, box_length, legend="raw", line_width=2)
size_plot.line(steps, box_length.mean(), legend="overall average",
               line_width=2, line_color="green")
size_plot.line(steps, box_length_rolling_mean, legend="rolling average", line_width=2,
               line_color="orange")

den_plot.line(steps, density, legend="density", line_width=2)
den_plot.line(steps, goal_density, legend="goal", line_width=2,
              line_color="green")
den_plot.line(steps, density_rolling_mean, legend="rolling average", line_width=2,
              line_color="orange")
den_plot.line(steps, density_mc, legend="MC density", line_width= 1,
              line_color="firebrick")

pres_plot.line(steps, pressure_mc, legend="MC pressure", line_width= 2,
               line_color="firebrick")
pres_plot.line(steps, pressure_mc_rolling_mean, legend="Mean MC pressure",
               line_width= 2, line_color="orange")

unitless_pres_plot.line(steps, box_data[:, 3], legend="MC pressure", line_width= 2,
                        line_color="firebrick")
unitless_pres_plot.line(steps, pd.rolling_mean(box_data[:, 3], n_mean),
                        legend="MC pressure", line_width= 2, line_color="orange")

# show the results
den_plot.legend.location = "bottom_right"
pres_plot.legend.location = "bottom_right"
unitless_pres_plot.legend.location = "bottom_right"
p = gridplot([[size_plot], [den_plot], [pres_plot, unitless_pres_plot]])
show(p)

logging.info('\m/ >.< \m/')
