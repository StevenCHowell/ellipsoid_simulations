import logging
import os
import sys
import numpy as np
from scipy.spatial.distance import pdist
import sasmol.sasmol as sasmol

sys.path.append('./')
import gr as fortran_gr

FORMAT = "%(asctime)-15s: %(message)s"
logging.basicConfig(format=FORMAT, level=logging.DEBUG)


def slow_update_gr(xcoor, ycoor, zcoor, box_length, n_bins, n_atoms, delta_r):

    gr = np.zeros(n_bins, np.float)

    for i in xrange(n_atoms - 1):
        for j in xrange(i + 1, n_atoms):

            dx = xcoor[i] - xcoor[j]
            dy = ycoor[i] - ycoor[j]
            dz = zcoor[i] - zcoor[j]

            dx -= box_length * ((dx / box_length).round())
            dy -= box_length * ((dy / box_length).round())
            dz -= box_length * ((dz / box_length).round())

            r = np.sqrt((dx * dx) + (dy * dy) + (dz * dz))

            if (r < box_length / 2.0):
                ig = int(r / delta_r)  # round down
                gr[ig] += 2

    return gr


def calc_gr(coor, box_length, gr, dr):

    dx = pdist(coor[:, :1])
    dy = pdist(coor[:, 0:2])
    dz = pdist(coor[:, 1:])

    dx -= box_length * ((dx / box_length).round())
    dy -= box_length * ((dy / box_length).round())
    dz -= box_length * ((dz / box_length).round())

    r = np.sqrt((dx * dx) + (dy * dy) + (dz * dz))

    r_i = (r[r < box_length / 2.0] / dr).astype(int)  # round down
    r_i_unique = np.unique(r_i, return_counts=True)
    gr[r_i_unique[0]] += 2 * r_i_unique[1]


def main(pdb_fname, run_log_fname, stride=1, sigma=1, dcd_fname=None,
         n_skip=0, n_bins=333):

    box_mol = sasmol.SasMol(0)
    box_mol.read_pdb(pdb_fname)

    if dcd_fname:
        dcd_file = box_mol.open_dcd_read(dcd_fname)
        n_frames = dcd_file[2]
    else:
        n_frames = 1

    run_log = np.loadtxt(run_log_fname)

    if len(run_log) != n_frames:
        if len(run_log) == n_frames + 1:
            run_log = run_log[:-1]
            logging.warning('dcd file  had one more frame than the log file, '
                         'discarding last line\ndcd_fname:\t{}\n'
                         'run_log_fname:\t{}'.format(dcd_fname, run_log_fname))
        else:
            logging.error('mismatch between dcd and log file \n'
                          'dcd_fname:\t{}\nrun_log_fname:\t{}'.format(
                              dcd_fname, run_log_fname))

    n_atoms = box_mol.natoms()
    n_gr = n_frames - n_skip  # number of frames, or g(r) curves, to averages

    box_length = run_log[:, 1] * sigma
    print('box_length: (min, max) = ({}, {})'.format(box_length.min(),
                                                     box_length.max()))

    gr_all = np.zeros((n_gr, n_bins))  # one g(r) for ecah dcd frame
    gr = np.zeros((n_bins, 2))

    # using the same r_grid for each frame
    dr = box_length[n_skip:].max() / (2.0 * n_bins)  # delg in F&S
    bin_index = np.arange(n_bins) + 1
    rho = (n_atoms / (box_length ** 3.0)).reshape(-1, 1)  # frame specific density

    r = (bin_index - 0.5) * dr
    bin_volume = 4.0 / 3.0 * np.pi * ((r + dr / 2) ** 3 - (r - dr / 2) ** 3)
    n_ideal = bin_volume * rho  # expected n for ideal gas

    if n_skip:
        for i in xrange(n_skip):
            # read and throw away these coordinates
            box_mol.read_dcd_step(dcd_file, i)

    for i in xrange(n_skip, n_frames):
        sys.stdout.flush()

        try:
            box_mol.read_dcd_step(dcd_file, i)
        except NameError:
            print('calculating g(r) for {}'.format(pdb_fname))

        coor = box_mol.coor()[0] * sigma
        calc_gr(coor, box_length[i], gr_all[i-n_skip], dr)
        gr_all[i-n_skip] /= n_ideal[i]  # normalize expected n for ideal gas

    box_mol.close_dcd_read(dcd_file)
    gr[:, 0] = r
    gr[:, 1] = np.mean(gr_all, axis=0) / n_atoms  # normalize by frames and atoms

    np.savetxt('gr_all_stride{}.dat'.format(stride), gr, fmt='%.14f')
    gr_cutoff = gr[gr[:, 0] < box_length.min()/2]
    np.savetxt('gr_cutoff_stride{}.dat'.format(stride), gr_cutoff, fmt='%.14f')

    return gr_cutoff


def plot_gr(gr, stride=1, show=False):
    import matplotlib.pyplot as plt

    gr_dat = np.loadtxt('argon_85K_gr.dat')

    fig = plt.figure()
    ax1 = fig.add_subplot(111)

    ax1.set_ylabel('g(r)')
    ax1.set_xlabel('r')
    scale_factor = 1.0
    ax1.plot(gr[:, 0], scale_factor * gr[:, 1], color='red', lw=2)
    ax1.plot(gr_dat[:, 0], gr_dat[:, 1], color='blue', lw=2)

    plt.savefig('steve_gr.png')

    if show:
        plt.show()
    else:
        plt.close('all')


if __name__ == '__main__':
    import os.path as op
    sigma = 3.405

    run_path = '../../simulations/lj_sphere_monomer/runs/p_0p14/output'
    pdb_fname = 'run2.pdb'
    dcd_fname = 'run2.dcd'
    xst_fname = 'box_length.txt'
    pdb_fname = op.join(run_path, pdb_fname)
    dcd_fname = op.join(run_path, dcd_fname)
    xst_fname = op.join(run_path, xst_fname)

    gr = main(pdb_fname, xst_fname, sigma=sigma, dcd_fname=dcd_fname,
              n_skip=1000)

    plot_gr(gr)
