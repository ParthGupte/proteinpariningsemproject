from proteinhelper import *
from PDBtoStructure import *
import time
from hungarian import *

#oneema to itself
one_emapdb = PDB("PDBfiles/1ema.pdb")
print(one_emapdb.file_name)
one_ema_struc = PDBtoStructure(one_emapdb)
one_ema_struc.display('r')
one_ema_struc_copy = one_ema_struc.copy()
one_ema_struc_copy.tumble()
t1 = time.time()
best_neibourhood_find(one_ema_struc,one_ema_struc_copy,5,display = True)
t2 = time.time()
print("That took ",t2-t1,"secs")

#oneqyo to itself
one_qyopdb = PDB("PDBfiles/1qyo.pdb")
print(one_qyopdb.file_name)
one_qyo_struc = PDBtoStructure(one_qyopdb)
one_qyo_struc_copy = one_qyo_struc.copy()
one_qyo_struc_copy.tumble()
t1 = time.time()
best_neibourhood_find(one_qyo_struc,one_qyo_struc_copy,5,display = True)
t2 = time.time()
print("That took ",t2-t1,"secs")

#oneema to oneqyo
t1 = time.time()
best_neibourhood_find(one_ema_struc,one_qyo_struc,5,display = True)
t2 = time.time()
print("That took ",t2-t1,"secs")
