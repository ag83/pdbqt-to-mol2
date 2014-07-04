import argparse
import sys
from cinfony import pybel

def usage():
   print 'This script get coordinates from docked ligand pdbqt file and after renumbering pick it to original mol2 file. This get mol2 file with connectivity and properly added hydrogens'
   print 'Usage: convert_to_mol2'
   print '-i origanal mol2 file' 
   print '-p undocked pdbqt file'
   print '-r docked pdbqt file'
   print '-o output mol2 file'
   print '-h Print this help'


def get_coord_dict(format, inputmol):
   molH = pybel.readfile(format, inputmol).next()
   molH.OBMol.DeleteHydrogens()
   return {atom.idx: atom.coords for atom in molH}

par_opt = argparse.ArgumentParser()
par_opt.add_argument('-i', action="store", help='origanal mol2 file')
par_opt.add_argument('-p', action="store", help='undocked pdbqt file')
par_opt.add_argument('-r', action="store", help='docked pdbqt file')
par_opt.add_argument('-o', action="store", help='output mol2')
args = vars(par_opt.parse_args())

if None in args.values():
   usage()
   sys.exit(2)
else:
   undocked_pdbqt_mol = args['p']
   docked_pdbqt_mol = args['r']
   origanal_mol2_mol = args['i']
   output_mol2 = args['o']

undocked_pdbqt = get_coord_dict('pdbqt', undocked_pdbqt_mol)
docked_pdbqt = get_coord_dict('pdbqt', docked_pdbqt_mol)
original_mol2 = get_coord_dict('mol2', origanal_mol2_mol)

if len(undocked_pdbqt) == len(docked_pdbqt) == len(original_mol2):
   original_coord = {}
   for key in original_mol2:
      coord_update = [round(x, 3) for x in original_mol2[key]]
      coord_update = tuple(coord_update)
      original_coord.update({key: coord_update})

   coord_map = {}
   for idx, coord in original_coord.items():
      #potential bottleneck for large molecules
      for ind, coordinates in undocked_pdbqt.items():
         n=0
         if coord[0] == coordinates[0]:
            n=n+1
            if coord[1] == coordinates[1]:
               n=n+1
               if coord[2] == coordinates[2]:
                  n=n+1
            else:
               if coord[2] == coordinates[2]:
                  n=n+1       
         else: 
            if coord[1] == coordinates[1]:
               n=n+1
               if coord[2] == coordinates[2]:
                  n=n+1
            else: 
               if coord[2] == coordinates[2]:
                  n=n+1
         if n == 3:
           coord_map.update({idx:ind})
         elif n == 2:
            coord_map.update({idx:ind})
         elif n == 1:
            coord_map.update({idx:ind})
         else:
            pass

   if len(coord_map) == len(original_mol2):
      coord_conform = {}
      for index1, index2 in coord_map.items():
         coord_conform.update({index1:docked_pdbqt.get(index2)})
      mol2 = pybel.readfile('mol2', origanal_mol2_mol).next()
      mol2.OBMol.DeleteHydrogens()
      for atom in mol2:
         atom.OBAtom.SetVector(coord_conform.get(atom.idx)[0], coord_conform.get(atom.idx)[1], coord_conform.get(atom.idx)[2])
      mol2.write('mol2', output_mol2, overwrite=True)
   else:
      print 'Lost coordinates in mapping'

else: 
   print 'Not equal number of atoms in molecules'
print docked_pdbqt_mol
  
