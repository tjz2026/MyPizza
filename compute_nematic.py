#!/usr/local/bin/python
'''
Script:  compute_nematic.py
Purpose: convert a LAMMPS dump file( configurations of polymer brush,currently)
         to XYZ format and then compute the nematic ordering of the whole chains.
Syntax:  compute_nematic.py dumpfile Nid Ntype Nx Ny Nz xyzfile 
         dumpfile = LAMMPS dump file in native LAMMPS format
         Nid,Ntype,Nx,Ny,Nz = columns #s for ID,type,x,y,z
                              (usually 1,2,3,4,5)
         xyzfile = new XYZ file
Author:  Steve Plimpton (Sandia), sjplimp at sandia.gov
Modified by Jiuzhou Tang ( tangjiuzhou@iccas.ac.cn) June, 2015
'''

import sys,os
path = os.environ["LAMMPS_PYTHON_TOOLS"]
sys.path.append(path)
from dump import dump
from Myxyz import xyz
import numpy as NP

#    atoms = snap.atoms
#    ids = atoms[:,id]
#    ordering = np.argsort(ids)
#    for i in xrange(len(atoms[0])):
#      atoms[:,i] = np.take(atoms[:,i],ordering)
#



def cmp_nematic(XYZ):
  n = flag = 0
  while 1:
    which,time,flag = XYZ.data.iterator(flag)
    if flag == -1: break
    time,box,atoms,bonds,tris,lines = XYZ.data.viz(which)
    print atoms[0]
    atom_sorted = sorted(atoms)
    print atom_sorted[0]
    #atom_sorted = sorted(atoms,key=lambda x: (x[0], -x[1]))
    atom_array=NP.array(atoms) 
    print "atom shape is", atom_array.shape 
    print "atom dim %s" % atom_array.ndim
    #ids=atoms[:,1]
    #ordering = NP.argsort(ids)
    #print ids
    #print ordering

    for atom in atom_sorted:
      itype = int(atom[0])
      print itype,atom[2],atom[3],atom[4]
    
    print time,
    n += 1
    
    return None
# --------------------------------------------------------------------



if len(sys.argv) != 8:
  raise StandardError, "Syntax: Compute_nematic.py dumpfile Nid Ntype Nx Ny Nz xyzfile"

dumpfile = sys.argv[1]
nid = int(sys.argv[2])
ntype = int(sys.argv[3])
nx = int(sys.argv[4])
ny = int(sys.argv[5])
nz = int(sys.argv[6])
xyzfile = sys.argv[7]

d = dump(dumpfile)
d.map(nid,"id",ntype,"type",nx,"x",ny,"y",nz,"z")
x = xyz(d)
x.one(xyzfile)
cmp_nematic(x)


# --------------------------------------------------------------------
