#!/usr/bin/env python3

import re, sys, time
from scipy.optimize import brentq, minimize
from functools import partial
from numpy import pi, sinh, cosh, arccosh, sin, cos, sqrt
from subprocess import check_output # switch from check_output to run in the future
from copy import deepcopy

twopi = 2*pi
scale_factor = 8
scale = list(map(lambda x : scale_factor * pow(2, x / 6.), range(0,-6,-1)))
COMP_ERR = pow(2,-100)

def shift_imag_to_zero(x) :
    print(x)
    while x.imag > pi :
      x -= twopi * 1j
    print(x.imag, -pi + COMP_ERR, COMP_ERR, -pi + COMP_ERR)
    while x.imag < -pi + COMP_ERR :
      x += twopi * 1j
    print(x)
    return x  

def get_box_codes(validated_params, depth=240) :
    params_printed = False
    sinhdx = validated_params['sinhdx']
    sinhdy = validated_params['sinhdy']
    coshmu = validated_params['coshmu']
    cosf = validated_params['cosf']
    sintx2 = validated_params['sintx2']
    sinty2 = validated_params['sinty2']
    coord = [0]*6
    coord[0] = sinhdx / (scale[0] * 16.) # these are special
    coord[1] = sinhdy / (scale[1] * 16.) # these are special
    coord[2] = coshmu / scale[2]
    coord[3] = cosf / scale[3]
    coord[4] = sintx2 / scale[4]
    coord[5] = sinty2 / scale[5]
    codes_list = [{ 'code' : [], 'coord' : coord }]
    validated_params['possibly_on_box_edge'] = False
    for i in range(0, depth) :
        for idx in range(0,len(codes_list)) :
            code = codes_list[idx]['code'] # Object that will be modified
            coord = codes_list[idx]['coord'] # Object that will be modified
            n = i % 6
            if 2 * coord[n] > COMP_ERR :
                code.append('1')
                coord[n] = 2 * coord[n] - 1
            elif 2 * coord[n] < -COMP_ERR :
                code.append('0')
                coord[n] = 2 * coord[n] + 1
            else :
                assert abs(coord[n]) < COMP_ERR
                print('Warning: Edge condition for manifold {0} with coord {1}. Will generate both children.'.format(validated_params['manifold'], coord[n]), file = sys.stderr)
                if not params_printed:
                  print(validated_params, file = sys.stderr)
                  params_printed = True
                validated_params['possibly_on_box_edge'] = True
                new_code_dict = deepcopy(codes_list[idx])
                new_code = new_code_dict['code']
                new_coord = new_code_dict['coord']
                # Old will go to the right
                code.append('1')
                coord[n] = 2 * coord[n] - 1
                # New will go to the left
                new_code.append('0')
                new_coord[n] = 2 * new_coord[n] + 1
                codes_list.append(new_code_dict)
    box_codes = []
    for code_dict in codes_list :
        box_codes.append(''.join(code_dict['code']))
    return box_codes

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
        # print('Used Brent algo for l1 = {}, l2 = {}, d = {}'.format(l1,l2,d), file = sys.stderr)
        return { 'cosh_margulis' : cosh_lox_move_dist(l1,t_min), 'margulis' : arccosh(cosh_lox_move_dist(l1,t_min)), 'point' : t_min }
    except :
        try :
            g = partial(max_of_two,l1,l2,d)
            res = minimize(g, d/2., bounds = ((0,d),))
            if res.success :
                return { 'cosh_margulis' : res.fun[0], 'margulis' : arccosh(res.fun[0]), 'point' : res.x[0] }
            else :
                print(res, file = sys.stderr)
                return -2
        except :
            print('Completely failed for {}, {}, and {}'.format(l1,l2,d), file = sys.stderr)
            return -1

def run_snap(cmds, timeout) :
    out = None
    try :
        out = check_output('snap', input = cmds, timeout = timeout).decode()
    except :
        print('Failed or timed out.', file = sys.stderr)
    return out

def get_geods(t, n, l, timeout = 240) :
    cmds = ''.join(('read {} {}\n'.format(t,n),
                    'print geodesics {}\n'.format(l))).encode()
    return run_snap(cmds, timeout)

def get_orthos(t, n, l, d, geods, timeout = 600) :
    geod_str =  ' '.join(map(str,geods))
    cmds = ''.join(('read {} {}\n'.format(t,n),
                    'print geodesics {}\n'.format(l),
                    'ortholines {} {}\n'.format(d,geod_str))).encode()
    return run_snap(cmds, timeout)

def get_mfld_info(t, n, timeout = 5) :
    cmds = ''.join(('read {} {}\n'.format(t,n),
                    'print name\n',
                    'print solution_type\n',
                    'print volume\n')).encode()
    return run_snap(cmds, timeout)

def get_ortho_bound(length_list, margulis_bound) :
    """ Given a list of lengths and a Margulis bound, this function returns
    a bound on othroline lengths needed to consider to achieve this Margulis bound. """
    ortho_bound = 0.5
    for l in length_list :
        z = l
        while z.real < margulis_bound :
            d = 2 * arccosh(sqrt( (cosh(margulis_bound) - cos(z.imag)) / (cosh(z.real) - cos(z.imag)) ))
            if ortho_bound < d :
                ortho_bound = d
            z += l
    return ortho_bound + 0.0001 # just up by a little bit so it's really a bound

if __name__ == "__main__" :
    def print_usage() :
            print('Usage : margulis type mfld [margulis_bound] [m1] [l1] [m2] [l2] ...\n'
                  '   type : snap census type, either 5,6,7, closed, manifold, or file.\n'
                  '   mfld : snap census manfiold index or file name\n'
                  '   margulis_bound : optional starting guess, but required if giving slope\n'
                  '   m1,l1,m2,l2... : the meridian and longitude filling slopes, respectivly.')

    census_counts = {'5' : 414, '6' : 961, '7' : 3551, 'closed' : 11031}
    census_letters = {'5' : 'm', '6' : 's', '7' : 'v'}
    other_types = ['file','manifold','link']
    if ( len(sys.argv) < 3 or
         (sys.argv[1] in census_counts.keys() and
         int(sys.argv[2]) > census_counts[sys.argv[1]]) or 
         (sys.argv[1] not in census_counts.keys() and sys.argv[1] not in other_types) ) :
        print_usage() 
        sys.exit(2)

    mfld_type = sys.argv[1]
    mfld_id = sys.argv[2]

    if mfld_type in census_counts.keys() and mfld_type != 'closed' :
        mfld_type = 'census ' + mfld_type

    # cusped and surgery cases
    if len(sys.argv) > 4 :
        mfld_id += ' surgery ' + ' '.join(sys.argv[4:])
 
    margulis_bound = 1.5
    if len(sys.argv) > 3 :
        try :
            margulis_bound =  float(sys.argv[3])
        except :
            print('Using initial margulis guess of {}'.format(margulis_bound))
 
    name = '{} {}'.format(mfld_type, mfld_id)
    volume = None
    info_out = get_mfld_info(mfld_type, mfld_id)
    if not info_out : sys.exit(3)
    for line in info_out.splitlines() :
        name_match = re.match('name : (.*)', line)    
        if name_match :
            name = name_match.group(1)
            continue
        vol_match = re.match('Volume is: (.*)', line)
        if vol_match :
            volume = vol_match.group(1)
            continue
        # fail_match = re.match('.*(other|flat|degenerate|nongeometric)', line)
        fail_match = re.match('.*(other|flat|degenerate)', line)
        if fail_match :
            print('Failed gluing solution of type {} for {}'.format(fail_match.group(1), name),
                  file = sys.stderr)
            sys.exit(4)
    
    # main loop that will update margulis constant if needed
    best = {'margulis' : 0, 'left' : 0, 'right' : 0, 'ortho' : 0}
    while margulis_bound < 3 : # arbitrary bound on max guess
        geods = {}
        geods_out = get_geods(mfld_type, mfld_id, margulis_bound)
        if not geods_out :
            print('Computing geodesics failed or timed out. Try smaller margulis guess.', file = sys.stderr) 
            sys.exit(5)
        for line in geods_out.splitlines() :
            geod_match = re.match('\s*\[(\d+)\]\s+([0-9\.]+)([0-9\.\+\-]+)\*i',line)
            if geod_match :
                idx, real, imag = geod_match.groups()
                geods[idx] = float(real) + float(imag) * 1j
                continue
            fail_match = re.match('.*Problem computing a Dirichlet domain', line)
            if fail_match :
                print('Dirichlet domain failed for {} {}.'.format(mfld_type, mfld_id),
                      file = sys.stderr)
                sys.exit(4)
        # increase guess if no geods are found
        if not geods :
            margulis_bound += 0.125
            continue

        ortho_bound = get_ortho_bound(geods.values(), margulis_bound)
       
        # get all needed ortholines
        orthos = {}       
        orthos_out = get_orthos(mfld_type, mfld_id, margulis_bound, ortho_bound, geods.keys())
        if not orthos_out :
            print('Computing ortholines failed or timed out. Try smaller margulis guess.', file = sys.stderr) 
            sys.exit(6)
        for line in orthos_out.splitlines() :
            ortho_match = re.match('\s*([0-9\.]+)([0-9\.\+\-]+)\*i\s+(\d+):.+\s+(\d+):.*',line)
            if ortho_match :
                real, imag, left, right = ortho_match.groups()
                key = left + ':' + right
                if key not in orthos :
                    orthos[key] = float(real) + float(imag) * 1j
                continue
        # increase guess if no orthos are found
        if not orthos :
            margulis_bound += 0.125
            continue

        margulis_numbers = []
        for key, ortho in orthos.items() :
            left, right = key.split(':')
            l, r = geods[left], geods[right]
            D = ortho.real 
            L = l
            l_pow = 1
            data = []
            while L.real < margulis_bound :
                R = r
                r_pow = 1
                while R.real < margulis_bound :
                    margulis_info = get_margulis_bound(L, R, D)
                    if isinstance(margulis_info, dict) :
                      margulis_info['l_pow'] = l_pow
                      margulis_info['r_pow'] = r_pow
                      data.append(margulis_info)
                    R += r
                    r_pow += 1
                L += l
                l_pow += 1
            if len(data) > 0 :
                data.sort(key = lambda x : x['margulis'])
                margulis_info = data[0]
                margulis_info['left'] = l
                margulis_info['right'] = r
                margulis_info['ortho'] = ortho
                margulis_numbers.append(margulis_info)

        sorted_margulis = sorted(margulis_numbers, key = lambda x : x['margulis'])

        best = sorted_margulis[0]
        if best['margulis'] < margulis_bound :
            break
        else :
            margulis_bound = best['margulis'] + 0.0001 # just up by a little bit so it's really a bound and we termiante above
    # get the box code
    ortho = best['ortho']
    t = best['point']
    if ortho.real == 0 or t == 0 or t == ortho.real : 
        print("Not an internal parameter point!")
        print(margulis_info)
    else :
        coshmu = best['cosh_margulis']
        sinhdx = sinh(t)
        sinhdy = sinh(ortho.real - t)
        # corresponds to applying complex conj to the group to reduce parameter space
        do_conj = True if ortho.imag < 0 else False
        print(do_conj)
        cosf = cos(ortho.imag)
        l = best['left']
        r = best['right']
        if do_conj :
          l = l.conjugate() 
          r = r.conjugate() 
          ortho = ortho.conjugate()
        l_pow = best['l_pow']
        r_pow = best['r_pow']
        x = shift_imag_to_zero(l * l_pow)
        y = shift_imag_to_zero(r * r_pow)
        # complex conj will change the sign on the imag parts
        sintx2 = sin(x.imag / 2.)
        sinty2 = sin(y.imag / 2.)
        box_codes = get_box_codes({'manifold' : name, 'coshmu' : coshmu, 'sinhdx' : sinhdx, 'sinhdy' : sinhdy, 'cosf' : cosf, 'sintx2' : sintx2, 'sinty2' : sinty2 })
        line = '"{}",' + '{},' * 17 +'"{}"'
        print(line.format(name, volume, best['margulis'], l, r, ortho,
                          t, ortho.real - t, l_pow, r_pow, x,y, sinhdx, sinhdy, coshmu, cosf, sintx2, sinty2, box_codes))
