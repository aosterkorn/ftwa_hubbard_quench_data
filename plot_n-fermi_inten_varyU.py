#!/usr/bin/env python3

import sys
import argparse

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("path", type=str, nargs=1, help="base directory where to find the files")
    parser.add_argument("-L", "--length", type=int, default=10, help="lattice sidelength")
    parser.add_argument("-p", "--particles", type=int, default=25, help="number of particles")
    # parser.add_argument("-n", "--n_flavor", type=int, default=100, help="degeneracy of fermions")
    # parser.add_argument("-U", "--U_hub", type=float, default=0.5, help="Hubbard interaction")
    parser.add_argument("-N", "--num_tsteps", type=int, default=101, help="number of time steps")
    parser.add_argument("--T_start", type=float, default=0.0, help="start time")
    parser.add_argument("--T_end", type=float, default=5.0, help="end time")
    parser.add_argument("-s", "--save", nargs="?", const="dump", default="", help="save to pdf file")
    parser.add_argument("--noview", action="store_true", default=False, help="no visual output")
    parser.add_argument("--notitle", action="store_true", default=False, help="no title on plot")
    args = parser.parse_args()
    
    path = args.path[0]
    L = args.length
    N = args.particles
    V = L*L
    
    ncs = (2, 100)
    Us = ("1e-1", "25e-2", "5e-1", "75e-2", "1")

    dt = (args.T_end - args.T_start) / (args.num_tsteps-1)
    s_i_t = round(args.T_start / dt)
    e_i_t = round(args.T_end / dt) + 1
    
    data_fermi = [ np.loadtxt(path + f"/n_fermi_N{nc}") for nc in ncs ]
    data_inten = [ np.loadtxt(path + f"/n_inten_N{nc}") for nc in ncs ]

    # how the data is presented
    # colors = [matplotlib.cm.Oranges(155), matplotlib.cm.Oranges(180), matplotlib.cm.Oranges(205), matplotlib.cm.Oranges(230), matplotlib.cm.Oranges(255)]
    # colors = [ 'blue', 'orange', 'green', 'pink', 'brown', 'purple', 'gray', 'red', 'yellow' ]
    pfmt_pth = {'color': 'black', 'linewidth': 2.0, 'linestyle': 'dashed', 'label': 'p.th.'}
    # pfmt_num = [ {'color': colors[i_U], 'linewidth': 2.0} for i_U in range(len(Us)) ]
    pfmt_num = [ { 'linewidth': 2.0 } for i_U in range(len(Us)) ]
    
    font = {'size': 16}
    matplotlib.rc('font', **font)
    matplotlib.rc('text', usetex=True)
    matplotlib.rc('axes', labelsize=20)
    plt.rcParams.update({'text.latex.preamble' : r'\usepackage{amsmath}'})
    
    fig, axs = plt.subplots(2, len(ncs), sharex='col', sharey='row', gridspec_kw={'hspace': 0, 'wspace': 0}, figsize=(10, 6), dpi=150)
    
    for i_n, nc in enumerate(ncs):
        # numerical results       
        for i_U, U in enumerate(Us):
            U = float(U)
            label = (f'$U = {U}$' if i_n == 0 else None)
            axs[0,i_n].plot(data_fermi[i_n][s_i_t:e_i_t,0], data_fermi[i_n][s_i_t:e_i_t,2+i_U], **pfmt_num[i_U], label=label)
            axs[1,i_n].plot(data_inten[i_n][s_i_t:e_i_t,0], data_inten[i_n][s_i_t:e_i_t,2+i_U], **pfmt_num[i_U], label=label)

        # perturbation theory
        axs[0,i_n].plot(data_fermi[i_n][:,0], data_fermi[i_n][:,1], **pfmt_pth)
        # axs[1,i_n].plot(data_inten[:,0], -(1/(V * 0.5))*(data_energy[:,2] - pt_data_energy[0,2]), **pfmt)
        axs[1,i_n].plot(data_inten[i_n][:,0], data_inten[i_n][:,1], **pfmt_pth)
 
        axs[0,i_n].set_title(f'$N = {nc}$', fontsize=26)
        axs[0,i_n].tick_params(direction="in")
        axs[1,i_n].tick_params(axis='y', direction="in")
        axs[0,i_n].label_outer()
        axs[1,i_n].label_outer()
    
        axs[1,i_n].set(xlabel = r'time $t/t_\text{h}^{-1}$', xlim = (args.T_start, args.T_end))
        plt.text(0.325, 0.9, r'$N U^{-2} \big( n_{k_F}(t) - n_{k_F}(0) \big)$', horizontalalignment='center', verticalalignment='center', transform=axs[0,i_n].transAxes)
        plt.text(0.325, 0.9, r'$N U^{-1} \big( \langle n_i^2 \rangle(t) - \langle n_i^2 \rangle(0) \big)$', horizontalalignment='center', verticalalignment='center', transform=axs[1,i_n].transAxes)
    
    axs[0,0].tick_params(axis='y', direction="out")
    axs[1,0].tick_params(axis='y', direction="out")
    axs[1,0].set_xticks((0, 1, 2, 3, 4))
    axs[1,1].set_xticks((0, 1, 2, 3, 4, 5))
    axs[1,1].set_xticklabels(('$5/0$', '$1$', '$2$', '$3$', '$4$', '$5$'))
    
    axs[0,0].set(ylabel = 'Fermi surface')
    axs[1,0].set(ylabel = 'interaction energy')

    handles, labels = axs[0,0].get_legend_handles_labels()
    # fig.legend(handles, labels, loc='upper center', ncol=6, fontsize=12)
    axs[0,1].legend(handles, labels, loc='lower right', frameon=False, ncol=2, fontsize=16)

    plt.subplots_adjust(left=0.1, right=0.99, top=0.93, bottom=0.12)

    if not args.save == "":
        plt.savefig(args.save + ".pdf", format='pdf')
        plt.savefig(args.save + ".eps", format='eps')
    
    if not args.noview:
        plt.show()

