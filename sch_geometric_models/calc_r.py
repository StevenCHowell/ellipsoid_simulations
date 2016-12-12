#!/usr/bin/env python
#coding:utf-8
"""
    Author:  Steven C. Howell --<steven.howell@nist.gov>
    Purpose: calculating the R-factor from a SasCalc run
    Created: 12/09/2016

00000000011111111112222222222333333333344444444445555555555666666666677777777778
12345678901234567890123456789012345678901234567890123456789012345678901234567890
"""

import glob
import logging
import os
import pandas as pd
import numpy as np

logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.DEBUG)


def get_r(rf_data, mt_data):
    r, _ = get_r_components(rf_data, mt_data)

    return r


def get_r_components(rf_data, mt_data):
    # R value as defined by doi: 10.1042/bj2670203
    diff = np.abs(mt_data[:, 1] - rf_data[:, 1])
    norm = np.abs(rf_data[:, 1]).sum()
    components = diff / norm
    r = components.sum()

    return r, components

def sascalc_r_factors(run_dir, iq_goal, ext='*.iq', run_name=''):
    r_fname = 'r_factors.txt'
    if run_name:
        r_fname = '{}_{}'.format(run_name, r_fname)
    r_full_fname = os.path.join(run_dir, r_fname)

    if not os.path.exists(r_full_fname):
        goal_iq = np.loadtxt(goal_fname)

        file_search = os.path.join(run_dir, ext)
        iq_fnames = glob.glob(file_search)
        assert iq_fnames, 'ERROR: no files found using {}'.format(file_search)
        iq_fnames.sort()

        index = []
        r_factors = []
        for fname in iq_fnames:
            iq_data = np.asfortranarray(np.loadtxt(fname))
            iq_data[:, 1] *= iq_goal[0, 1] / iq_data[0, 1]
            r = get_r(iq_goal, iq_data)

            r_factors.append(r)
            index.append(int(fname.replace(ext[1:], '')[-5:]))

        data = {'r': r_factors, 'fnames': iq_fnames}
        r_results = pd.DataFrame(data=data, index=index)
        r_results.to_csv(r_full_fname, sep='\t', float_format='%5.10f',
                         columns=['r', 'fnames'])
    else:
        logging.info('output already exists, remove using:\nrm {0}'.format(
            r_full_fname))
        r_results = pd.DataFrame.from_csv(r_full_fname, sep='\t')

    return r_results


def get_p(rf_data, mt_data):
    p, _ = get_p_components(rf_data, mt_data)

    return p


def get_p_components(rf_data, mt_data):
    # percent error
    diff = mt_data[:, 1] - rf_data[:, 1]
    norm = rf_data[:, 1]
    components = np.abs(diff / norm) * 100 / len(norm)
    p = components.sum()

    return p, components


def sascalc_percent_error(run_dir, iq_goal, ext='*.iq', run_name=''):
    p_fname = 'percent_error.txt'
    if run_name:
        p_fname = '{}_{}'.format(run_name, p_fname)
    p_full_fname = os.path.join(run_dir, p_fname)

    if not os.path.exists(p_full_fname):
        goal_iq = np.loadtxt(goal_fname)

        file_search = os.path.join(run_dir, ext)
        iq_fnames = glob.glob(file_search)
        assert iq_fnames, 'ERROR: no files found using {}'.format(file_search)
        iq_fnames.sort()

        index = []
        per_error = []
        for fname in iq_fnames:
            iq_data = np.asfortranarray(np.loadtxt(fname))  # load data
            iq_data[:, 1] *= iq_goal[0, 1] / iq_data[0, 1]  # rescale data
            p = get_p(iq_goal, iq_data)

            per_error.append(p)
            index.append(int(fname.replace(ext[1:], '')[-5:]))

        data = {'p': per_error, 'fnames': iq_fnames}
        p_error = pd.DataFrame(data=data, index=index)
        p_error.to_csv(p_full_fname, sep='\t', float_format='%5.10f',
                         columns=['p', 'fnames'])
    else:
        logging.info('output already exists, remove using:\nrm {0}'.format(
            p_full_fname))

    return p_error


if __name__ == '__main__':
    run_dir = 'lys_ellipsoid_scan1/sascalc/neutron_D2Op_100/'
    assert os.path.exists(run_dir), 'ERROR, bad path: {}'.format(run_dir)

    # goal_fname = 'data/lysozyme_00001.iq'
    goal_fname = 'data/exp_lys/rebinned_exp_data_lysozyme.dat'
    assert os.path.exists(goal_fname), 'ERROR, bad fname: {}'.format(goal_fname)
    goal_data = np.asfortranarray(np.loadtxt(goal_fname)[:, :2])

    r_factors = sascalc_r_factors(run_dir, goal_data, run_name='pdb')
    p_errors = sascalc_percent_error(run_dir, goal_data, run_name='pdb')

    logging.debug('\m/ >.< \m/')