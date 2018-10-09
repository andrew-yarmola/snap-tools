#!/usr/bin/env python3

import sys
import csv
from margulis import *

if __name__ == "__main__":
    margulis_bound = 1.4
    print('"name","margulis upper","left geod","left power","right geod","right power","ortho"')
    for values in csv.reader(sys.stdin) :
        # skip leading line
        if values[0] == 'name' : continue        
        name = values[0]
        # print(values, file = sys.stderr)
        l1, l2, d = map(eval, values[1:])
        L1 = l1
        l1_pow = 1
        data = []
        while L1.real < margulis_bound :
            L2 = l2
            l2_pow = 1
            while L2.real < margulis_bound :
                data.append({ 'margulis' : get_margulis_bound(L1,L2,d.real), 'l1_pow' : l1_pow, 'l2_pow' : l2_pow })
                L2 += l2
                l2_pow += 1
            L1 += l1
            l1_pow += 1
        if len(data) > 0 :
            sorted_data = sorted(data, key = lambda x : x['margulis'])
            best = sorted_data[0]
            print('"{}",{},{},{},{},{},{}'.format(name, best['margulis'], repr(l1)[1:-1], best['l1_pow'], repr(l2)[1:-1], best['l2_pow'], repr(d)[1:-1]))
