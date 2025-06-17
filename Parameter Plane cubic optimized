import numpy as npimport matplotlib.pyplot as plt

import matplotlib.cm as cm
import cmath

import time

from numba import jit

#nit = number of iteration
@jit
def compute_nit(multiplier, xmin, xmax, ymin, ymax, im_width, im_height, zabs_max_square, nit_max):
    parameter = np.zeros((im_height, im_width))
    xwidth = xmax - xmin
    yheight = ymax - ymin
    
    for ix in range(im_width):
        for iy in range(im_height):
            nit = 0
            # Map pixel position to a point in the complex plane
            a = cmath.sqrt(complex(ix / im_width * xwidth + xmin,
                        iy / im_height * yheight + ymin))
            
            crit_neg = (-a-cmath.sqrt(a**2-3*multiplier))/3
            crit_pos = (-a+cmath.sqrt(a**2-3*multiplier))/3
            iter_crit_neg = crit_neg            
            iter_crit_pos = crit_pos

            re_neg = iter_crit_neg.real
            im_neg = iter_crit_neg.imag
            re_pos = iter_crit_pos.real
            im_pos = iter_crit_pos.imag

            re_neg_square = re_neg * re_neg
            im_neg_square = im_neg * im_neg
            re_pos_square = re_neg * re_neg
            im_pos_square = im_neg * im_neg
            # Do the iterations
            while (re_neg_square + im_neg_square) <= zabs_max_square and (re_pos_square + im_pos_square) <= zabs_max_square and nit < nit_max:
                re_neg = iter_crit_neg.real
                im_neg = iter_crit_neg.imag
                re_pos = iter_crit_pos.real
                im_pos = iter_crit_pos.imag

                re_neg_square = re_neg * re_neg
                im_neg_square = im_neg * im_neg
                re_im_neg = re_neg * im_neg

                re_iter_neg_square = re_neg_square - im_neg_square
                im_iter_neg_square = 2*re_im_neg

                re_iter_neg_cube = re_neg * ( re_neg_square - 3*im_neg_square )
                im_iter_neg_cube = im_neg * ( 3*re_neg_square - im_neg_square )

                re_pos_square = re_pos * re_pos
                im_pos_square = im_pos * im_pos
                re_im_pos = re_pos * im_pos

                re_iter_pos_square = re_pos_square - im_pos_square
                im_iter_pos_square = 2*re_im_pos

                re_iter_pos_cube = re_pos * ( re_pos_square - 3*im_pos_square )
                im_iter_pos_cube = im_pos * ( 3*re_pos_square - im_pos_square )

                iter_crit_neg = (re_iter_neg_cube + 1j*im_iter_neg_cube) + a*(re_iter_neg_square + 1j*im_iter_neg_square) + multiplier*iter_crit_neg                            #THIS LINE IS THE FORMULA FOR MY ITERATION
                iter_crit_pos = (re_iter_pos_cube + 1j*im_iter_pos_cube) + a*(re_iter_pos_square + 1j*im_iter_pos_square) + multiplier*iter_crit_pos
                nit += 1
            ratio = np.floor( nit / nit_max )
             # number of row is "y" position in complex plane, vice versa for column = "x"
             # thus the unfortunate situation iy,ix is interchanged
            parameter[iy,ix] = ratio
            
    return parameter

def parameter_plane_cubic(multiplier):
    im_width, im_height = 500, 500

    zabs_max = 4                      # bailout criteria
    zabs_max_square = zabs_max**2
    nit_max =   200                  # max number of iteration

    xmin, xmax = -6, 6
    xwidth = xmax - xmin
    ymin, ymax = -6, 6
    yheight = ymax - ymin
    
    parameter = compute_nit(multiplier, xmin, xmax, ymin, ymax, im_width, im_height, zabs_max_square, nit_max)
    
    fig, ax = plt.subplots(figsize=(6, 5*im_height/im_width))

    # the most important part, this is the coloring scheme. cmap = "color map". I use pre-made coloring called "tab20c"
    ax.pcolormesh(parameter, cmap=cm.tab20c, vmin=0, vmax=1)

    # this is to display defintion of color map "tab20c"
    # fig.colorbar( ax.pcolormesh(parameter, cmap=cm.tab20c, vmin=0, vmax=1), ax=ax)

    print("Below is parameter plane for multiplier " + str(multiplier))

    # Set the labels to the coordinates of z0 in the complex plane
    xtick_labels = np.linspace(xmin, xmax, int(xwidth / 0.5))
    ax.set_xticks([(x-xmin) / xwidth * im_width for x in xtick_labels])
    ax.set_xticklabels(['{:.2f}'.format(xtick) for xtick in xtick_labels])

    ytick_labels = np.linspace(ymin, ymax, int(yheight / 0.5))
    ax.set_yticks([(y-ymin) / yheight * im_height for y in ytick_labels])
    ax.set_yticklabels(['{:.2f}'.format(ytick) for ytick in ytick_labels])
    plt.show()

# Here is an example how to use the code

parameter_plane_cubic(0.9*cmath.exp(2*cmath.pi*1j*0.2))
