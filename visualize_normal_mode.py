#!/usr/bin/env python
"""Usage: visualize_normal_mode.py g03_freq.fchk mode"""

import sys
from numpy import *
from gaussian_file_handling import get_gaussian_results
from normal_coordinates import compute_normal_coord, print_input_files
from conversion import au2cm,au2kcal,au2ang

test

if len(sys.argv)<=2:
    print __doc__
    sys.exit(0)

filename=sys.argv[1]
mode=int(sys.argv[2])
print 'Gaussian formatted checkpoint file:',filename
print 'XYZ structures created for normal mode:',mode

def generate_mode_sequence(atom_l,cart_l,Q,mode):
    nsteps=15
    nr_atoms=len(atom_l)
    (N,nr_modes)=shape(Q)

    out=open('normal_mode_sequence_m'+str(mode)+'.xyz','w')

    for i in range(nsteps):
        factor = 10.0*sin( pi* (float(i)/float(nsteps)-0.5) )
        xyz = cart_l[:] + factor*Q[:,mode-1]

        out.write(str(nr_atoms)+2*'\n')
        for k in range(nr_atoms):
            line='%s%20.8f%20.8f%20.8f\n' %  \
                (atom_l[k],xyz[k*3]*au2ang,xyz[k*3+1]*au2ang,xyz[k*3+2]*au2ang)
            out.write(line)

h,atom_l,mass_l,cart_l=get_gaussian_results(filename)
Q,f=compute_normal_coord(h,mass_l)
print "freq=",f[mode-1]*au2cm
generate_mode_sequence(atom_l,cart_l,Q,mode)

print '*** End of program ***'
