#!/bin/bash

figdir="figures/"

mkdir -p $figdir

python3 plot_n-fermi_inten_varyU.py fig2_L10_p25_varyU -s $figdir/fermi_inten_N2_N100_L10
python3 plot_n-fermi_inten_varyN.py fig3_L20_p101_U5e-1_varyN -s $figdir/fermi_inten_U5e-1_L20
python3 plot_band_occupation_varyU.py fig4_L10_p25_varyU -s $figdir/band_occupation_N2_N100_L10
python3 plot_band_occupation_varyN.py fig5_L20_p101_U5e-1_varyN -s $figdir/band_occupation_U5e-1_L20
