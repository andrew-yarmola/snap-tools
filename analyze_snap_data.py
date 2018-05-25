#!/usr/bin/env python3

import sys
from scipy.optimize import brentq, minimize
from functools import partial
from numpy import sinh, cosh, arccosh, cos


def cosh_lox_move_dist(l,d) :
    """ Let l be the complex length of a loxodromic element g and
        let x be a point distance d form the axis of g. This function
        returns the cosh of the distance from x to g(x). """
    return cosh(d)**2 * cosh(l.real) - cos(l.imag) * sinh(d)**2

def cosh_lox_dist_diff(l1,l2,d,t) :
    """ We want to find t such that 
            lox_move_dist(l1,t) == lox_move_dist(l2,d-t)
        This function is the one we will optimize to do this. """
    return cosh_lox_move_dist(l1,t) - cosh_lox_move_dist(l2, d-t)

def max_of_two(l1,l2,d,t) :
    """ Maximum displament along an ortholine."""
    return max(cosh_lox_move_dist(l1,t), cosh_lox_move_dist(l2, d-t))

def get_margulis_bound(l1,l2,d) :
    f = partial(cosh_lox_dist_diff,l1,l2,d)
    try :
        t_min = brentq(f,0,d)
    except :
        try :
            g = partial(max_of_two,l1,l2,d)
            res = minimize(g, d/2., bounds = ((0,d),))
            if res.success :
                return res.fun[0]
            else :
                print(res, file = sys.stderr)
                return -1
        except :
            print("Completely failed for {}, {}, and {}".format(l1,l2,d), file = sys.stderr)
            return -1
    return arccosh(cosh_lox_move_dist(l1,t_min))

if __name__ == "__main__":
    for line in sys.stdin :
        values = line.rstrip().split(',')
        if len(values) == 4 :
            name = values[0]
            l1, l2, d = map(eval, values[1:])
            loop_count = 0
            # 0.8 is our upper bound for Margulis number
            L1 = l1
            while L1.real < 0.8 :
                L2 = l2
                while L2.real < 0.8 :
                    print("{} : {}".format(name, get_margulis_bound(L1,L2,d.real)))
                    L2 += L2
                    loop_count += 1
                L1 += L1
            if loop_count > 1 :
                print("{} count {}".format(name, loop_count), file = sys.stderr)
