import numpy as np

## smooth funciton
def smooth(x,n):

    # ensure n is odd
    if (n % 2) == 0:
        n += 1

    # calculate wing size
    wing_size = int((n-1)/2)

    # initiate x_smooth (no need to change first and last value)
    x_smooth = x.copy() # brackets ensures list x is not edited
    #x_smooth = x_smooth.astype(float)
    j = 1
    rest = len(x)-j
    while j < len(x)-1:
        if (j < wing_size) or (rest-1 < wing_size) > j:
            i_max = j
        elif rest-1 < wing_size:
            i_max = rest-1
        else:
            i_max = wing_size
        sum = x[j]
        for i in range(1,i_max+1):
            sum += x[j-i] + x[j+i]
        x_smooth[j] = sum/(2*i_max+1)
        j += 1
        rest = len(x) - j

    return x_smooth
