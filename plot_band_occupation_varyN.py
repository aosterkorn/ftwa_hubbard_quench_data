#!/usr/bin/env python3

import sys
import h5py
import argparse

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

# col_oranges = [matplotlib.cm.Oranges(127), matplotlib.cm.Oranges(159), matplotlib.cm.Oranges(191), matplotlib.cm.Oranges(223), matplotlib.cm.Oranges(255)]
col_oranges = [matplotlib.cm.Oranges(155), matplotlib.cm.Oranges(180), matplotlib.cm.Oranges(205), matplotlib.cm.Oranges(230), matplotlib.cm.Oranges(255)]
col = col_oranges

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("path", type=str, nargs=1, help="where to find the .h5 file")
    parser.add_argument("-L", "--length", type=int, default=20, help="lattice sidelength")
    parser.add_argument("-p", "--particles", type=int, default=25, help="number of particles")
    # parser.add_argument("-n", "--n_flavor", type=int, default=100, help="degeneracy of fermions")
    # parser.add_argument("-U", "--U_hub", type=float, default=0.5, help="Hubbard interaction")
    parser.add_argument("-N", "--num_tsteps", type=int, default=201, help="number of time steps")
    parser.add_argument("--T_start", type=float, default=0.0, help="start time")
    parser.add_argument("--T_end", type=float, default=10.0, help="end time")
    parser.add_argument("-s", "--save", nargs="?", const="dump", default="", help="save to pdf file")
    parser.add_argument("--noview", action="store_true", default=False, help="no visual output")
    parser.add_argument("--notitle", action="store_true", default=False, help="no title on plot")
    args = parser.parse_args()
    
    path = args.path[0]
    L = args.length
    N = args.particles
    # U = args.U_hub
    # nc = args.n_flavor
    
    ncs = (2, 10, 100, 1000)
    i_ts = (20, 40, 100, 200)
    U = "5e-1"
    
    pf = [1.0, 1.0, 1.0/4.0, 1.0/10.0]

    font = {'size': 16}
    matplotlib.rc('font', **font)
    matplotlib.rc('text', usetex=True)
    matplotlib.rc('axes', labelsize=23)
    plt.rcParams.update({'text.latex.preamble' : r'\usepackage{amsmath}'})
    
    fig, axs = plt.subplots(4, 1, sharex='col', sharey='row', gridspec_kw={'hspace': 0, 'wspace': 0}, figsize=(8, 8.0), dpi=150)
    
    # fig.suptitle(r'$N U^{-2} \langle \Delta n(\epsilon) \rangle(t)$', x=0.175, fontsize=26)
    fig.suptitle(r'$N U^{-2} \big( n_\epsilon(t) - n_\epsilon(0) \big)$', x=0.215, fontsize=26)
    for ii, i_t in enumerate(i_ts):
        time = i_t * (args.T_end - args.T_start) / (args.num_tsteps - 1)
        axs[0].set_title(rf'$U = {float(U):.1f} t_\text{{h}}$', fontsize=26)

        for i_n, nc in enumerate(ncs):
            label = f'$N = {nc}$'
            data_occ = np.loadtxt(path + f"/n_energy_U{U}_t{int(time)}")
            
            axs[ii].errorbar(data_occ[:,0], pf[ii]*data_occ[:,3+i_n], yerr = pf[ii]*data_occ[:,3+len(ncs)+i_n], label=label)
            if not pf[ii] == 1.0:
                axs[ii].text(0.75, 0.65, rf'$\times {1.0/pf[ii]:.1f}$', transform=axs[ii].transAxes)
        
        axs[ii].errorbar(data_occ[:,0], pf[ii]*data_occ[:,2], color='black', linestyle='dashed', label='p.th.')
    
        axs[ii].set_xticks([-4.0, -2.0, 0.0, 2.0, 4.0])
        axs[ii].set(xlabel = r'band energy $\epsilon/t_\text{h}$', ylabel = rf'$t = {int(time)}$', xlim = (-4.0, 4.0), ylim=(-0.075, 0.075))
        axs[ii].grid(True, 'major', 'y')
        axs[ii].label_outer()

    handles, labels = axs[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='upper center', bbox_to_anchor=(0.7, 1.01), ncol=3, frameon=False, fontsize=14)
    plt.subplots_adjust(left=0.14, right=0.98, top=0.86, bottom=0.1)

    if not args.save == "":
        plt.savefig(args.save + ".pdf", format='pdf')
        plt.savefig(args.save + ".eps", format='eps')
        
    if not args.noview:
        plt.show()

