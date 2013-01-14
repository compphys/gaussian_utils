#!/usr/bin/env python
"""Usage: create_input_files.py g03_freq.fchk [template_file] [modes]"""

import sys
from numpy import *

from gaussian_utils.get_hessian_fchk import get_gaussian_results
from gaussian_utils.normal_coordinates import compute_normal_coord
from utils.conversion import au2cm,au2kcal,au2ang

if len(sys.argv)<=1:
    print __doc__
    sys.exit(0)

def print_input_files(filename,atom_l,cart_l,Q,mode_l):
    """Generate input files from template file for finite displacements
    along the normal modes"""

    from utils.utils import from_file, to_file
    import sys

    step_size = 0.01
#    step_l = ['-5','-2','-1','+0','+1','+2','+5']
    step_l = ['-1','+1']

    d=['x','y','z']

    nr_atoms=len(atom_l)
    (N,nr_modes)=shape(Q)
    if (nr_modes!=3*nr_atoms-6 and nr_modes!=3*nr_atoms-5) or N!=3*nr_atoms:
        print 'Incorrect data in print_input_files'
        sys.exit()

    for mode in mode_l:
        i = mode-1 # array index for mode m is m-1, etc.
        for step in step_l:
            factor = step_size*float(step)
            xyz = cart_l[:] + factor*Q[:,i]

            s = from_file(filename)

            for k in range(nr_atoms):
                for l in range(3):
                    replace_by = '%17.14f' % xyz[k*3 + l]
                    if s.find(atom_l[k] + d[l])>0:
                        s = s.replace(atom_l[k] + d[l], str(replace_by), 1)
                    else:
                        print 'print_input_files: error string not found: ',\
                            atom_l[k] + d[l]
                        sys.exit(2)
            if i<9:
                mode = 'm0' + str(i+1)
            else:
                mode = 'm' + str(i+1)
            file = filename.replace('template', mode + step)

            to_file(s, file)
#
# start of main program
#

# read input from Gaussian formatted checkpoint file, name of
# file is the first argument to this script.
print get_gaussian_results.__doc__
h,atom_l,mass_l,cart_l=get_gaussian_results(sys.argv[1])
print 'Hessian=',h
print 'List of atoms=',atom_l
print 'Atomic masses=',mass_l
print 'Atomic Cartesian coordinates (in Angstrom):'
i=0
for atom in atom_l:
    print "%s%18.8f%14.8f%14.8f" % (atom,cart_l[i:i+3][0]*au2ang,
                                    cart_l[i:i+3][1]*au2ang,
                                    cart_l[i:i+3][2]*au2ang)
    i+=3

print compute_normal_coord.__doc__
Q,f=compute_normal_coord(h,mass_l)
print "freq=",f*au2cm
print 'zpve=',sum(f)/2*au2kcal,'kcal/mole'

nmodes=len(f)
if len(sys.argv)==3:
    mode_l = range(1,nmodes+1)
elif len(sys.argv)>3:
    mode_l=map(int,sys.argv[3:])

if len(sys.argv)<=2:
    print 'Q=',Q

if len(sys.argv)>=3:
    print 'Input files created for normal modes:',mode_l
    print 'Mode Freq'
    for mode in mode_l:
        print '%4i%14.8f' % (mode,f[mode-1]*au2cm)
    print print_input_files.__doc__
    print_input_files(sys.argv[2],atom_l,cart_l,Q,mode_l)

print '*** End of program ***'
