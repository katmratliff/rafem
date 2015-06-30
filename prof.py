#! /usr/local/bin/python"""

def make_profile(dx, dy, n, riv_x, riv_y, profile):

    profile = [0 for i in range(len(riv_x))]

    for p in range(len(riv_x)):
        
        profile[p] = n[riv_x[p]/dx][riv_y[p]/dy]
        #profile.append(n[riv_x[p]][riv_y[p]])
        
        p += 1
    
    return profile
