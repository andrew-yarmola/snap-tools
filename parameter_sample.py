#!/usr/bin/env python3

import numpy as np
from margulis_futer import get_box_codes
from itertools import product

if __name__ == "__main__" :
    coshmu_vals = np.linspace(1.00542, 1.3174, 10)
    sinhdx_vals = np.linspace(0.0001, 1.0001, 10)
    sinhdy_vals = np.linspace(0.0001, 1.0001, 10)
    cosf_vals = np.linspace(-0.99, 0.99, 10)
    sintx2_vals = np.linspace(-0.99, 0.99, 10)
    sinty2_vals = np.linspace(-0.99, 0.99, 10)
    for coshmu, sinhdx, sinhdy, cosf, sintx2, sinty2 in product(coshmu_vals,
        sinhdx_vals, sinhdy_vals, cosf_vals, sintx2_vals, sinty2_vals):
        box_codes = get_box_codes({'manifold' : "random", 'coshmu' : coshmu,
            'sinhdx' : sinhdx, 'sinhdy' : sinhdy, 'cosf' : cosf,
            'sintx2' : sintx2, 'sinty2' : sinty2 }, depth=90)
        print(box_codes[0])
