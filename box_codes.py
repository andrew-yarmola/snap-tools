# Computes box_codes from parameters

import sys
from copy import deepcopy

g_scale_factor = 8
g_scale = list(map(lambda x : g_scale_factor * pow(2, x / 6.), range(0,-6,-1)))

def get_box_codes(validated_params, depth=240, scale=None) :
    COMP_ERR = pow(2,-100)
    if scale is None:
        scale = g_scale
    else: # using Arb
        COMP_ERR = 0
    params_printed = False
    sinhdx = validated_params['sinhdx']
    sinhdy = validated_params['sinhdy']
    coshmu = validated_params['coshmu']
    cosf = validated_params['cosf']
    sintx2 = validated_params['sintx2']
    sinty2 = validated_params['sinty2']
    coord = [0]*6
    coord[0] = sinhdx / scale[0]
    coord[1] = sinhdy / scale[1]
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
