#!/usr/bin/env python

#import ase
#from general_surface import oriented_surface
from ase.visualize import view
from ase.spacegroup import crystal
from ase.build.tools import cut
from ase.build import surface, graphene_nanoribbon

### input data ###
Surface1='Au'
###Surface2=graphene_nanoribbon( 3, 4, type='armchair', saturated=True, vacuum=3.5 )
###Surface2.edit()
###exit()
# Miller index of the surface
hkl0=[1,1,1] # initial values
hklF=[1,1,1] # final values
hklD=[1,1,1] # step size (NYI)
# Mixture of the standardize surface vectors
mn0=[1,0] # initial values
mnF=[1,2] # final values
mnD=[1,1] # step size (NYI)
nLayersS=3
#nLayersA=3
fVacuum=10.0
#nDimA=1 # dimension of the "adsorbate": 1D, 2D or 3D
fNumEPS=1.0E-8

### variables ###
hkl=[0,0,0]
mn=[0,0]

for hkl[0] in range(hkl0[0],hklF[0]+1):
 for hkl[1] in range(hkl0[1],hklF[1]+1):
  for hkl[2] in range(hkl0[2],hklF[1]+1):
   print( hkl )
   s1 = surface(Surface1, (hkl[0],hkl[1],hkl[2]), nLayersS, fVacuum, fNumEPS)
   #s1 = oriented_surface(Surface1, (hkl[0],hkl[1],hkl[2]), nLayersS, fVacuum, fNumEPS)
   #s1.edit()
   for mn[0] in range(mn0[0],mnF[0]+1):
    for mn[1] in range(mn0[1],mnF[1]+1):
     print( mn )
     orient1 = cut(s1, a=(mn[0],mn[1],0), b=(-mn[1],mn[0],0), c=(0,0,1))
     orient1.edit()

# Create an aluminium (111) slab with three layers
#a = 4.05
#aluminium = crystal('Al', [(0,0,0)], spacegroup=225, cellpar=[a, a, a, 90, 90, 90])
## Then cut out the slab
#al111 = cut(aluminium, (2,0,0), (0,2,0), nlayers=3)
#al111.edit()
#
