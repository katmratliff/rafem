#! /usr/local/bin/python

import math

# this function uses a linear diffusion equation (e.g. Paola 2000, Jerolmack
# and Paola 2007) to compute elevation change along the river course
def smooth_rc(dx, dy, nu, dt, riv_x, riv_y, n, nslope):

    # elevation change along river course due to diffusional smoothing
    dn_rc = [0 for i in range(len(riv_x))]

    #dn_rc[1] = -(dx * nslope)

    for c in range(1, len(riv_x)-1):

        # determine downstream change
        if (((riv_x[c+1]/dx) - (riv_x[c]/dx) == 0) and
                ((riv_y[c+1]/dy) - (riv_y[c]/dy) == -1)):
                k1 = 1
        elif (((riv_x[c+1]/dx) - (riv_x[c]/dx) == 0) and
                ((riv_y[c+1]/dy) - (riv_y[c]/dy) == 1)):
                k1 = 1
        elif (((riv_x[c+1]/dx) - (riv_x[c]/dx) == 1) and
                ((riv_y[c+1]/dy) - (riv_y[c]/dy) == 0)):
                k1 = 1
        elif (((riv_x[c+1]/dx) - (riv_x[c]/dx) == 1) and
                ((riv_y[c+1]/dy) - (riv_y[c]/dy) == -1)):
                k1 = math.sqrt(2)
        elif (((riv_x[c+1]/dx) - (riv_x[c]/dx) == 1) and
                ((riv_y[c+1]/dy) - (riv_y[c]/dy) == 1)):
                k1 = math.sqrt(2)

        dwnst_dn = ((n[riv_x[c+1]/dx][riv_y[c+1]/dy]
                    - n[riv_x[c]/dx][riv_y[c]/dy]) / (dx*k1))

        # determine upstream change
        if (((riv_x[c]/dx) - (riv_x[c-1]/dx) == 0) and
                ((riv_y[c]/dy) - (riv_y[c-1]/dy) == -1)):
                k2 = 1
        elif (((riv_x[c]/dx) - (riv_x[c-1]/dx) == 0) and
                ((riv_y[c]/dy) - (riv_y[c-1]/dy) == 1)):
                k2 = 1
        elif (((riv_x[c]/dx) - (riv_x[c-1]/dx) == 1) and
                ((riv_y[c]/dx) - (riv_y[c-1]/dy) == 0)):
                k2 = 1      # CUS
        elif (((riv_x[c]/dx) - (riv_x[c-1]/dx) == 1) and
                ((riv_y[c]/dy) - (riv_y[c-1]/dy) == -1)):
                k2 = math.sqrt(2)    # LUS
        elif (((riv_x[c]/dx) - (riv_x[c-1]/dx) == 1) and
                ((riv_y[c]/dy) - (riv_y[c-1]/dy) == 1)):
                k2 = math.sqrt(2)    # RUS

        upst_dn = ((n[riv_x[c]/dx][riv_y[c]/dy]
                   - n[riv_x[c-1]/dx][riv_y[c-1]/dy]) / (dx*k2))

        # determine elevation change at cell
        dn_rc[c] = (nu*dt)/(dx**2) * (dwnst_dn - upst_dn)   # eqn 1 J&P 2007

        n[riv_x[c]/dx][riv_y[c]/dy] = n[riv_x[c]/dx][riv_y[c]/dy] + dn_rc[c]

    return (n, dn_rc)
