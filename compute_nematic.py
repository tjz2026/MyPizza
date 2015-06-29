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



def cmp_nematic(XYZ, chain_num):
  n = flag = 0
  while 1:
    which,time,flag = XYZ.data.iterator(flag)
    if flag == -1: break
    time,box,atoms,bonds,tris,lines = XYZ.data.viz(which)
    print atoms[0]
    atom_sorted = sorted(atoms)
    print atom_sorted[0]
    #atom_sorted = sorted(atoms,key=lambda x: (x[0], -x[1]))
    atom_array=NP.array(atom_sorted) 
    print "atom shape is", atom_array.shape 
    print "atom dim %s" % atom_array.ndim
    atom_num=int(atom_array.shape[0]) 
    print "total atom number =%s" % atom_num
    atom_pos=NP.zeros( (atom_num,4) )
    atom_pos[:,0]=atom_array[:,0]   
    atom_pos[:,1]=atom_array[:,2]   
    atom_pos[:,2]=atom_array[:,3]   
    atom_pos[:,3]=atom_array[:,4]   
    #ids=atoms[:,1]
    #ordering = NP.argsort(ids)
    #print ids
#   the Ith bond is defined as the displacement of atom[I+1]-atom[I]
    bond_vec=NP.array([0.,0.,0.])
    n_bond=0
    s1=0.0
    for atom in atom_pos:
      atom_id = int(atom[0])
      print atom_id,atom[1],atom[2],atom[3]
      if atom_id%chain_num == 0:
        pass
      else:
        bond_vec[0]=atom_pos[atom_id][1]-atom_pos[atom_id-1][1]
        bond_vec[1]=atom_pos[atom_id][2]-atom_pos[atom_id-1][2]
        bond_vec[2]=atom_pos[atom_id][3]-atom_pos[atom_id-1][3]
        bond_vec=bond_vec/NP.linalg.norm(bond_vec) 
        s1=s1+bond_vec[2]
        n_bond=n_bond+1
    s1=s1/n_bond 
    print "total amount of bond=%s" % n_bond
    print "ava of s1=%s" % s1
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
cmp_nematic(x,30)


# --------------------------------------------------------------------
