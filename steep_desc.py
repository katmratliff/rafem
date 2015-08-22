#! /usr/local/bin/python


def find_course(dx, dy, imax, jmax, n, riv_x, riv_y):
    # function to find the steepest descent route
    # note: this needs to be improved to remove potential bias that may occur
    # if two or more cells have the steepest descent elevation

    for r in range(imax*3):

        if riv_x[-1]/dx < imax-1:

            # if against right boundary
            if riv_y[-1]/dy + dy == jmax-1:
            # compute surrounding cell elevation differences
                LHS = (n[riv_x[-1]/dx][riv_y[-1]/dy] -
                       n[riv_x[-1]/dx][(riv_y[-1]/dy)-1])
                LDS = (n[riv_x[-1]/dx][riv_y[-1]/dy] -
                       n[(riv_x[-1]/dx)+1][(riv_y[-1]/dy)-1])
                CDS = (n[riv_x[-1]/dx][riv_y[-1]/dy] -
                       n[(riv_x[-1]/dx)+1][riv_y[-1]/dy])
            
                SteepDesc = max(LHS, LDS, CDS)
                
                if CDS == SteepDesc:
                    riv_x = riv_x + [riv_x[-1]+dx]
                    riv_y = riv_y + [riv_y[-1]]
                
                elif LDS == SteepDesc:
                    riv_x = riv_x + [riv_x[-1]+dx]
                    riv_y = riv_y + [riv_y[-1]-dy]
                    
                elif LHS == SteepDesc:
                    riv_x = riv_x + [riv_x[-1]]
                    riv_y = riv_y + [riv_y[-1]-dy]
            
            # if against left boundary
            elif riv_y[-1]/dy - dy == 0:
                RHS = (n[riv_x[-1]/dx][riv_y[-1]/dy] -
                       n[riv_x[-1]/dx][(riv_y[-1]/dy)+1])
                CDS = (n[riv_x[-1]/dx][riv_y[-1]/dy] -
                       n[(riv_x[-1]/dx)+1][riv_y[-1]/dy])
                RDS = (n[riv_x[-1]/dx][riv_y[-1]/dy] -
                   n[(riv_x[-1]/dx)+1][(riv_y[-1]/dy)+1])
                   
                SteepDesc = max(RHS, CDS, RDS)
                
                if CDS == SteepDesc:
                    riv_x = riv_x + [riv_x[-1]+dx]
                    riv_y = riv_y + [riv_y[-1]]
                
                elif RDS == SteepDesc:
                    riv_x = riv_x + [riv_x[-1]+dx]
                    riv_y = riv_y + [riv_y[-1]+dy]
                    
                elif RHS == SteepDesc:
                    riv_x = riv_x + [riv_x[-1]]
                    riv_y = riv_y + [riv_y[-1]+dy]       

            # if not a boundary cell, calculate all 5 possibilities
            else:
                # compute surrounding cell elevation differences
                LHS = (n[riv_x[-1]/dx][riv_y[-1]/dy] -
                       n[riv_x[-1]/dx][(riv_y[-1]/dy)-1])
                RHS = (n[riv_x[-1]/dx][riv_y[-1]/dy] -
                       n[riv_x[-1]/dx][(riv_y[-1]/dy)+1])
                LDS = (n[riv_x[-1]/dx][riv_y[-1]/dy] -
                       n[(riv_x[-1]/dx)+1][(riv_y[-1]/dy)-1])
                CDS = (n[riv_x[-1]/dx][riv_y[-1]/dy] -
                       n[(riv_x[-1]/dx)+1][riv_y[-1]/dy])
                RDS = (n[riv_x[-1]/dx][riv_y[-1]/dy] -
                       n[(riv_x[-1]/dx)+1][(riv_y[-1]/dy)+1])
    
                SteepDesc = max(LHS, RHS, LDS, CDS, RDS)
    
                if CDS == SteepDesc:
                    riv_x = riv_x + [riv_x[-1]+dx]
                    riv_y = riv_y + [riv_y[-1]]
    
                elif LHS == SteepDesc:
                    riv_x = riv_x + [riv_x[-1]]
                    riv_y = riv_y + [riv_y[-1]-dy]
    
                elif RHS == SteepDesc:
                    riv_x = riv_x + [riv_x[-1]]
                    riv_y = riv_y + [riv_y[-1]+dy]
    
                elif LDS == SteepDesc:
                    riv_x = riv_x + [riv_x[-1]+dx]
                    riv_y = riv_y + [riv_y[-1]-dy]
    
                elif RDS == SteepDesc:
                    riv_x = riv_x + [riv_x[-1]+dx]
                    riv_y = riv_y + [riv_y[-1]+dy]

        else:
            return riv_x, riv_y

    return riv_x, riv_y


def find_new_course(dx, dy, imax, jmax, n, new_riv_x, new_riv_y, current_SL):
    # function to find the steepest descent route
    # note: this needs to be improved to remove potential bias that may occur
    # if two or more cells have the steepest descent elevation

    for r in range(imax*3):

        if ((new_riv_x[-1]/dx < imax-1) and
            n[new_riv_x[-1]/dx][new_riv_y[-1]/dy] > current_SL):
            # second part of above conditional might cause problems for an 
            # anastomosing river... think about this in future

            # if against right boundary
            if new_riv_y[-1]/dy + dy == jmax-1:
            # compute surrounding cell elevation differences
                LHS = (n[new_riv_x[-1]/dx][new_riv_y[-1]/dy] -
                       n[new_riv_x[-1]/dx][(new_riv_y[-1]/dy)-1])
                LDS = (n[new_riv_x[-1]/dx][new_riv_y[-1]/dy] -
                       n[(new_riv_x[-1]/dx)+1][(new_riv_y[-1]/dy)-1])
                CDS = (n[new_riv_x[-1]/dx][new_riv_y[-1]/dy] -
                       n[(new_riv_x[-1]/dx)+1][new_riv_y[-1]/dy])
            
                SteepDesc = max(LHS, LDS, CDS)
                
                if CDS == SteepDesc:
                    new_riv_x = new_riv_x + [new_riv_x[-1]+dx]
                    new_riv_y = new_riv_y + [new_riv_y[-1]]
                
                elif LDS == SteepDesc:
                    new_riv_x = new_riv_x + [new_riv_x[-1]+dx]
                    new_riv_y = new_riv_y + [new_riv_y[-1]-dy]
                    
                elif LHS == SteepDesc:
                    new_riv_x = new_riv_x + [new_riv_x[-1]]
                    new_riv_y = new_riv_y + [new_riv_y[-1]-dy]
            
            # if against left boundary
            elif new_riv_y[-1]/dy - dy == 0:
                RHS = (n[new_riv_x[-1]/dx][new_riv_y[-1]/dy] -
                       n[new_riv_x[-1]/dx][(new_riv_y[-1]/dy)+1])
                CDS = (n[new_riv_x[-1]/dx][new_riv_y[-1]/dy] -
                       n[(new_riv_x[-1]/dx)+1][new_riv_y[-1]/dy])
                RDS = (n[new_riv_x[-1]/dx][new_riv_y[-1]/dy] -
                   n[(new_riv_x[-1]/dx)+1][(new_riv_y[-1]/dy)+1])
                   
                SteepDesc = max(RHS, CDS, RDS)
                
                if CDS == SteepDesc:
                    new_riv_x = new_riv_x + [new_riv_x[-1]+dx]
                    new_riv_y = new_riv_y + [new_riv_y[-1]]
                
                elif RDS == SteepDesc:
                    new_riv_x = new_riv_x + [new_riv_x[-1]+dx]
                    new_riv_y = new_riv_y + [new_riv_y[-1]+dy]
                    
                elif RHS == SteepDesc:
                    new_riv_x = new_riv_x + [new_riv_x[-1]]
                    new_riv_y = new_riv_y + [new_riv_y[-1]+dy]       

            else:
                # compute surrounding cell elevation differences
                LHS = (n[new_riv_x[-1]/dx][new_riv_y[-1]/dy] -
                       n[new_riv_x[-1]/dx][(new_riv_y[-1]/dy)-1])
                RHS = (n[new_riv_x[-1]/dx][new_riv_y[-1]/dy] -
                       n[new_riv_x[-1]/dx][(new_riv_y[-1]/dy)+1])
                LDS = (n[new_riv_x[-1]/dx][new_riv_y[-1]/dy] -
                       n[(new_riv_x[-1]/dx)+1][(new_riv_y[-1]/dy)-1])
                CDS = (n[new_riv_x[-1]/dx][new_riv_y[-1]/dy] -
                       n[(new_riv_x[-1]/dx)+1][new_riv_y[-1]/dy])
                RDS = (n[new_riv_x[-1]/dx][new_riv_y[-1]/dy] -
                       n[(new_riv_x[-1]/dx)+1][(new_riv_y[-1]/dy)+1])
    
                SteepDesc = max(LHS, RHS, LDS, CDS, RDS)
    
                if CDS == SteepDesc:
                    new_riv_x = new_riv_x + [new_riv_x[-1]+dx]
                    new_riv_y = new_riv_y + [new_riv_y[-1]]
    
                elif LHS == SteepDesc:
                    new_riv_x = new_riv_x + [new_riv_x[-1]]
                    new_riv_y = new_riv_y + [new_riv_y[-1]-dy]
    
                elif RHS == SteepDesc:
                    new_riv_x = new_riv_x + [new_riv_x[-1]]
                    new_riv_y = new_riv_y + [new_riv_y[-1]+dy]
    
                elif LDS == SteepDesc:
                    new_riv_x = new_riv_x + [new_riv_x[-1]+dx]
                    new_riv_y = new_riv_y + [new_riv_y[-1]-dy]
    
                elif RDS == SteepDesc:
                    new_riv_x = new_riv_x + [new_riv_x[-1]+dx]
                    new_riv_y = new_riv_y + [new_riv_y[-1]+dy]

        else:
            return new_riv_x, new_riv_y

    return new_riv_x, new_riv_y
