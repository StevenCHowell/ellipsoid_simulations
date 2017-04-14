#!/usr/bin/env python
# coding:utf-8
"""
    Author:  Steven C. Howell --<steven.howell@nist.gov>
    Purpose: optimizing a torus to match the scattering an atomic molecule
    Created: 12/05/2016

00000000011111111112222222222333333333344444444445555555555666666666677777777778
12345678901234567890123456789012345678901234567890123456789012345678901234567890

This file contains modules for affine minimization

Hailiang Zhang
Jan 2015
"""

import logging
import os
import sys
import time

import numpy as np

from scipy import optimize

sys.path.append('/home/schowell/data/code/ellipsoid_simulations/hz_affine/')
from affine.shapes import hollowCylinder_equation as shape_eq
from affine import fileIO

sys.path.append('/home/schowell/data/scratch/sascalc_demo/')
import sascalc

logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.DEBUG)

# if len(sys.argv) == 1:
    # NGPU = 1
# else:
    # import locale
    # NGPU = locale.atoi(sys.argv[1])


class AffineParameters:

    def __init__(self, sx, sy, sz, hxy=0.0, hxz=0.0, hyx=0.0, hyz=0.0, hzx=0.0,
                 hzy=0.0, dx=0.0, dy=0.0, dz=0.0):
        self.sx = sx
        self.sy = sy
        self.sz = sz
        self.hxy = hxy
        self.hxz = hxz
        self.hyx = hyx
        self.hyz = hyz
        self.hzx = hzx
        self.hzy = hzy
        self.dx = dx
        self.dy = dy
        self.dz = dz


class DummyMol:

    def __init__(self, coor):
        self._natoms = len(coor)
        self._coor = np.empty((1, self._natoms, 3))
        self._coor[0] = coor

    def coor(self):
        return self._coor

    def natoms(self):
        return self._natoms


def fill_shape(n_coor, full_affine_matrix, shape_eq, coor, extra_args):
    '''
    Purpose:
        Fill a shape in the affine space and do the bondary check
        by reverse affine transformation.

    Return:
        the affine coordinate matrix after transformation
    '''
    affine_matrix = full_affine_matrix[:3, :3]

    # get the boundaries
    # ZHL note: only check for 4 cube vertex because the other 4 is just the
    # reflections; this needs change if translation is included in the affine
    # transformation matrix'''
    cube_vertex_coor_matrix = np.matrix([[-1.0, 1.0, -1.0, 1.0],
                                         [-1.0, -1.0, 1.0, 1.0],
                                         [-1.0, -1.0, -1.0, -1.0]])

    affine_vertex_coor_matrix = affine_matrix * cube_vertex_coor_matrix
    affine_vertex_coor_matrix = np.concatenate((affine_vertex_coor_matrix,
                                                -affine_vertex_coor_matrix),
                                               axis=1)
    affine_vertex_coor_min = affine_vertex_coor_matrix.min(axis=1)
    affine_edge_max = (affine_vertex_coor_matrix.max(axis=1) -
                       affine_vertex_coor_min).max()
    # check for points
    reverse_trans_matrix = affine_matrix.I
    n_coor_found = 0
    while n_coor_found < n_coor:
        points_remaining = n_coor - n_coor_found

        # create random points
        trial_coor = np.random.rand(3, points_remaining)
        trial_coor *= affine_edge_max  # rescale the points
        trial_coor += affine_vertex_coor_min  # offset the points

        trial_coor_orth = reverse_trans_matrix * trial_coor

        good_coor = [shape_eq(*tuple(np.array(coor_i).ravel()), extra=extra_args)
                     for coor_i in trial_coor_orth.T]
        good_coor = np.array(good_coor)

        n_new_coor_found = n_coor_found + sum(good_coor)
        coor[:, n_coor_found:n_new_coor_found] = trial_coor.T[good_coor].T
        n_coor_found = n_new_coor_found

        assert n_coor_found, 'ERROR - check inputs:\ {}'.format(locals())


def calc_x2(parameters, iq_goal, n_atoms, coor=None, iq_model=None):

    sx, sy, sz, hxy, hxz, hyx, hyz, hzx, hzy, dx, dy, dz, inner_radius = list(
        parameters)

    parameters = np.array([sx, sy, sz, hxy, hxz, hyx, hyz, hzx, hzy, dx, dy,
                           dz, inner_radius])
    if coor is None:
        coor = np.matrix(np.zeros((3, n_atoms)))

    n_q = len(iq_goal)
    bs = np.ones(n_atoms)

    if (sx <= 0.0 or sy <= 0.0 or sz <= 0.0 or inner_radius <= 0.0 or
        inner_radius >= 1.0):
        return np.inf

    affine_matrix = np.matrix([[sx, hxy, hxz, dx],
                               [hyx, sy, hyz, dy],
                               [hzx, hzy, sz, dz],
                               [0.0, 0.0, 0.0, 1.0]])

    fill_shape(n_atoms, affine_matrix, shape_eq, coor, [inner_radius])

    mol = DummyMol(coor.T)
    tic = time.time()
    sascalc_instance = sascalc.SasCalc()
    sascalc_instance.calc_GV_fix_n(mol, option='vacuum', B=bs)
    toc = time.time() - tic
    logging.debug('calculation time: {} s'.format(toc))

    if iq_model is None:
        iq_model = np.copy(iq_goal[:, :2])
    iq_model[:, 1] = sascalc_instance.get_I()
    iq_model[:, 1] /= iq_model[0, 1] * iq_goal[0, 1]  # scale to match at I(0)

    x2 = 0
    diff = (iq_model[:, 1] - iq_goal[:, 1]) / iq_goal[:, 2]
    diff2 = diff * diff
    x2 = np.sum(diff2) / (n_q - 1)

    logging.debug('parameters: \n{}'.format(affine_matrix))
    logging.debug('inner_radius: {}'.format(inner_radius))
    logging.debug('reduced X2: {}'.format(x2))

    return x2


if __name__ == '__main__':

    # get the match data
    iq_goal_fname = 'data/Iq_exp_txt_q++.txt'
    assert os.path.exists(iq_goal_fname), 'no such file: {}'.format(
        iq_goal_fname)
    iq_goal = np.loadtxt(iq_goal_fname)
    # logging.debug('Exp data:\n{}'.format(iq_goal))

    n_atoms = 10000
    coor = np.matrix(np.zeros((3, n_atoms)))
    iq_model = np.copy(iq_goal[:, :2])

    # initial guesses
    sx = 40.0
    sy = 40.0
    sz = 30.0
    hxy = 0.0
    hxz = 0.0
    hyx = 0.0
    hyz = 0.0
    hzx = 0.0
    hzy = 0.0
    dx = 0.0
    dy = 0.0
    dz = 0.0
    inner_radius = 0.2

    parameters = np.array([sx, sy, sz, hxy, hxz, hyx, hyz, hzx, hzy, dx, dy,
                           dz, inner_radius])

    calc_x2(parameters, iq_goal, n_atoms, coor, iq_model)

    output_dir = 'output3'
    fileIO.make_and_write_pdb(coor, output_dir, 'model_start.pdb')
    np.savetxt(os.path.join(output_dir, 'model_start.iq'), iq_model)

    # # run minimization
    # methods = ['Powell', 'Nelder-Mead', 'BFGS']
    # args = (iq_goal, n_atoms, coor, iq_model)
    # results = optimize.minimize(calc_x2, parameters, args, methods[0])

    # logging.info('results: {}'.format(results))

    # fileIO.make_and_write_pdb(coor, 'output', 'model_final.pdb')
    # np.savetxt(os.path.join(output_dir, 'model_final.iq'), iq_model)
