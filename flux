#! /usr/local/bin/python

def calc_qs(nu, riv_x, riv_y, n, dx, dy, dt):

	sed_flux = 0
	dist = 0

	if (((riv_x[-2]/dx) - (riv_x[-1]/dx) == 0) and
        (riv_y[-2]/dy) - (riv_y[-1]/dy) == -2):
            
            dist = 1
    
    elif (((riv_x[-2]/dx) - (riv_x[-1]/dx) == 0)
        and (riv_y[-2]/dy) - (riv_y[-1]/dy) == 1):
    
            dist = 1
            
    elif (((riv_x[-2]/dx) - (riv_x[-1]/dx) == 1)
        and (riv_y[-2]/dy) - (riv_y[-1]/dy) == 0):

            dist = 1

    elif (((riv_x[-2]/dx) - (riv_x[-1]/dx) == 1)
        and (riv_y[-2]/dy) - (riv_y[-1]/dy) == -2):
            
            dist = math.sqrt(2)

    elif (((riv_x[-2]/dx) - (riv_x[-1]/dx) == 1)
        and (riv_y[-2]/dy) - (riv_y[-1]/dy) == 1):
            
            dist = math.sqrt(2)

	sed_flux = - (nu*dt) * ((n[riv_x[-1]/dx][riv_y[-1]/dy] - 
				n[riv_x[-2]/dx][riv_y[-2]/dy]) / dist)

	return (sed_flux)