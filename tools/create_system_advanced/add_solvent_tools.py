#Last updated 6/29/12, 1:10 pm

import random
import sys
import string
from math import *
from progressbar import ProgressBar

def add_ghost(ibox, atom, b_ghosts, r0, ndaxis, xlo, ylo, zlo, xhi, yhi, zhi):
  g=False
  i0=0; i1=0; j0=0; j1=0; k0=0; k1=0
  x = atom[0]
  y = atom[1]
  z = atom[2]
  if x-xlo[ibox]<r0:
    i0=-1
    g=True
  if xhi[ibox]-x<r0:
    i1=2
    g=True
  if y-ylo[ibox]<r0:
    j0=-1
    g=True
  if ylo[ibox]-y<r0:
    j1=2
    g=True
  if z-zlo[ibox]<r0:
    k0=-1
    g=True
  if zlo[ibox]-z<r0:
    k1=2
    g=True

  if g==True:
    pbox=[ibox/(ndaxis*ndaxis), ibox%(ndaxis*ndaxis)/ndaxis, ibox%(ndaxis*ndaxis)%ndaxis]
    for i in range(i0,i1):
      for j in range(j0,j1):
        for k in range(k0,k1):
          if i==0 and j==0 and k==0: continue
          ib=pbox[0]+i
          jb=pbox[1]+j
          kb=pbox[2]+k

          if ib<0: ib+=ndaxis
          if ib>=ndaxis: ib-=ndaxis
          if jb<0: jb+=ndaxis
          if jb>=ndaxis: jb-=ndaxis
          if kb<0: kb+=ndaxis
          if kb>=ndaxis: kb-=ndaxis
          
          b = (ib*ndaxis+jb)*ndaxis+kb
          b_ghosts[b].append(atom)

def add_solvent(nsol, box, r0, ndiv):
  # Box bounderies
  xstart = box[0][0]
  xend = box[0][1]
  ystart = box[1][0]
  yend = box[1][1]
  zstart = box[2][0]
  zend = box[2][1]

  # Box dimensions
  Lx = xend - xstart
  Ly = yend - ystart
  Lz = zend - zstart

  nboxes=8**ndiv # total number of boxes
  ndaxis=2**ndiv # number of divitions along each axes
  b_atoms=[] # per box atom arrays
  b_ghosts=[] # lists of ghosts atoms for each box
  xlo=[]; ylo=[]; zlo=[]
  xhi=[]; yhi=[]; zhi=[]

  for i in range(nboxes):
    b_atoms.append([])
    b_ghosts.append([])
    xlo.append(0.0)
    ylo.append(0.0)
    zlo.append(0.0)
    xhi.append(0.0)
    yhi.append(0.0)
    zhi.append(0.0)
  for i in range(ndaxis):
    for j in range(ndaxis):
      for k in range(ndaxis):
        ibox=(i*ndaxis+j)*ndaxis+k
        xlo[ibox] = xstart + i*Lx/ndaxis
        xhi[ibox] = xstart + (i+1)*Lx/ndaxis
        ylo[ibox] = ystart + j*Ly/ndaxis
        yhi[ibox] = ystart + (j+1)*Ly/ndaxis
        zlo[ibox] = zstart + k*Lz/ndaxis
        zhi[ibox] = zstart + (k+1)*Lz/ndaxis

  j=0
  natoms = 0
  fr=nsol/100
  if fr<=0: fr=1
  prb= ProgressBar('green', width=50, block='*', empty='o')
  for isol in range(nsol):
    ibox=random.randint(0,nboxes-1) # randomly choose a box
    g = True
    for i in range(100000):
      x = random.uniform(xlo[ibox],xhi[ibox])
      y = random.uniform(ylo[ibox],yhi[ibox])
      z = random.uniform(zlo[ibox],zhi[ibox])
      g = True
      #check if overlaps with anything
      for ia in b_atoms[ibox]:
        xd = x-ia[0]
        yd = y-ia[1]
        zd = z-ia[2]
        r = sqrt(xd*xd+yd*yd+zd*zd)
        if r < r0: g = False
      if g: break
    if not g:
      print "Warning: Cannot create solvent atom\n\n"
    #add  atoms to array  
    natoms += 1
    atom = [x, y, z]
    b_atoms[ibox].append(atom)

    #add atom as a ghost if needed
    add_ghost(ibox, atom, b_ghosts, r0, ndaxis, xlo, ylo, zlo, xhi, yhi, zhi)

    j += 1
    if j%fr==0:
      p = round(100*float(j)/nsol,2)
      prb.render(int(p), 'Adding solvent atoms\nNumber of tries %d' % i)
  prb.render(100, 'Adding solvent atoms\nNumber of tries %d' % i)

  # Gather atoms
  atoms = []
  for ib in range(nboxes):
    for ia in b_atoms[ib]:
      atoms.append(ia)

  return atoms
