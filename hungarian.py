from proteinhelper import *
from PDBtoStructure import *
import time

def best_neibourhood_find(protA,protB,m,display=False):
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
    final_out = neibourhood_match(best_match[0][0],best_match[0][1],protA,protB,m)
    print(final_out)
    if display:
        display_equiv_fit(final_out[1],protA,protB)
    return final_out


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

def display_equiv_fit(equiv,protA: Structure,protB: Structure):
    '''
    Outputs Pairing Class with the equivalances kearsley fitted and the structures transformed accordingly
    '''
    keys = []
    values = []
    for key, val in equiv.items():
        keys.append(key)
        values.append(val)
    equiv_pointsA = protA.points[keys].copy()
    equiv_pointsB = protB.points[values].copy()
    equiv_vecsA = protA.vectors[keys].copy()
    equiv_vecsB = protB.vectors[values].copy()
    equiv_strucA = Structure(equiv_pointsA,equiv_vecsA)
    equiv_strucB = Structure(equiv_pointsB,equiv_vecsB)
    equiv_pairing = Pairing(list(range(len(keys))),equiv_strucA,equiv_strucB)
    trans, rot = equiv_pairing.kearsley()[1:]
    protB.translate(-trans)
    protB.rotate(rot)
    superimpStruc = Pairing([],protA,protB)
    superimpStruc.display('g','r')


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


# display_equiv_fit(equiv,one_ema_struc,one_ema_struc_copy)

# display_equiv_fit(equiv,one_qyo_struc,one_qyo_struc_copy)

# display_equiv_fit(equiv,one_ema_struc,one_qyo_struc)