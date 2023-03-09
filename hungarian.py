from proteinhelper import *
from PDBtoStructure import *
import time

def best_neibourhood_find(protA,protB,m):
    neibour_scoredict = {}

    for i in range(protA.points.shape[0]):
        for j in range(protB.points.shape[0]):
            # NA = make_neibourhood(i,protA,m)[0]
            # NB = make_neibourhood(j,protB,m)[0]
            # NB.translate(-NB.points[0])
            # NA.translate(-NA.points[0])
            # if np.all(NA.points[0] != NB.points[0]):
            #     print(NA.points[0],NB.points[0])
            #     print("Not True for",(i,j))
            # vectorsA = Structure(np.concatenate((NA.vectors[0],-NA.vectors[0]))) #so that the translation component of kersley get cancelled and we can work with only rotations
            # vectorsB = Structure(np.concatenate((NB.vectors[0],-NB.vectors[0])))
            # pair_vecAB = Pairing([0,1,2,3,4,5],vectorsA,vectorsB)
            # rot = pair_vecAB.kearsley()[2]
            # NB.rotate(rot)
            # pair_match = Pairing([x for x in range(m+1)],NA,NB)
            # pair_match.hungarian()
            # rmsd = pair_match.rmsd()
            # if i==j:
            #     print(rmsd)
            #     pair_match.display('red','green')
            neibour_scoredict[(i,j)] = neibourhood_match(i,j,protA,protB,m)[0]

    best_match = sorted(neibour_scoredict.items(),key=lambda x:x[1])[0]
    print("Best match is:",best_match)
    print(neibourhood_match(best_match[0][0],best_match[0][1],protA,protB,m))


def neibourhood_match(i,j,protA,protB,m):
    NA, pointsA = make_neibourhood(i,protA,m)
    NB, pointsB = make_neibourhood(j,protB,m)
    NB.translate(-NB.points[0])
    NA.translate(-NA.points[0])
    vectorsA = Structure(np.concatenate((NA.vectors[0],-NA.vectors[0]))) #so that the translation component of kersley get cancelled and we can work with only rotations
    vectorsB = Structure(np.concatenate((NB.vectors[0],-NB.vectors[0])))
    pair_vecAB = Pairing([0,1,2,3,4,5],vectorsA,vectorsB)
    rot = pair_vecAB.kearsley()[2]
    NB.rotate(rot)
    pair_match = Pairing([x for x in range(m+1)],NA,NB)
    new_perm = pair_match.hungarian()
    rmsd = pair_match.rmsd()
    equivs = {} #key = index of points in A, value = index of points in B
    # if i==j:
    #         print(rmsd)
    #         pair_match.display('red','green')
    for i in range(len(new_perm)):
        equivs[pointsA[i]] = pointsB[new_perm[i]] 
    return rmsd , equivs


# cube_A = Structure(cube = True)
# cube_B = cube_A.copy()
# cube_B.tumble()
# best_neibourhood_find(cube_A,cube_B,7)

# anishpdb = PDB("PDBfiles/anish.pdb")
# print(anishpdb.file_name)
# anishstruc = PDBtoStructure(anishpdb)
# anishstruc.display('r')
# anishstruc_copy = anishstruc.copy()
# anishstruc_copy.tumble()
# best_neibourhood_find(anishstruc,anishstruc_copy,16)

one_emapdb = PDB("PDBfiles/1ema.pdb")
print(one_emapdb.file_name)
one_ema_struc = PDBtoStructure(one_emapdb)
one_ema_struc.display('r')
one_ema_struc_copy = one_ema_struc.copy()
one_ema_struc_copy.tumble()
t1 = time.time()
best_neibourhood_find(one_ema_struc,one_ema_struc_copy,5)
t2 = time.time()
print("That took ",t2-t1,"secs")

one_qyopdb = PDB("PDBfiles/1qyo.pdb")
print(one_qyopdb.file_name)
one_qyo_struc = PDBtoStructure(one_qyopdb)
one_qyo_struc_copy = one_qyo_struc.copy()
one_qyo_struc_copy.tumble()
t1 = time.time()
best_neibourhood_find(one_qyo_struc,one_qyo_struc_copy,5)
t2 = time.time()
print("That took ",t2-t1,"secs")

t1 = time.time()
best_neibourhood_find(one_ema_struc,one_qyo_struc,5)
t2 = time.time()
print("That took ",t2-t1,"secs")
