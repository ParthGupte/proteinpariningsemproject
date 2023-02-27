from proteinhelper import *

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


def neibourhood_match(i,j,protA,protB,m):
    NA, pointsA = make_neibourhood(i,protA,m)
    NB, pointsB = make_neibourhood(j,protB,m)[0]
    NB.translate(-NB.points[0])
    NA.translate(-NA.points[0])
    vectorsA = Structure(np.concatenate((NA.vectors[0],-NA.vectors[0]))) #so that the translation component of kersley get cancelled and we can work with only rotations
    vectorsB = Structure(np.concatenate((NB.vectors[0],-NB.vectors[0])))
    pair_vecAB = Pairing([0,1,2,3,4,5],vectorsA,vectorsB)
    rot = pair_vecAB.kearsley()[2]
    NB.rotate(rot)
    pair_match = Pairing([x for x in range(m+1)],NA,NB)
    pair_match.hungarian()
    rmsd = pair_match.rmsd()
    equivs = {} #key = index of points in A, value = index of points in B
    for i in range(len(pair_match.perm)):
        equivs[pointsA[i]] = pointsB[pair_match.perm[i]] 
    if i==j:
            print(rmsd)
            pair_match.display('red','green')
    return rmsd , equivs


cube_protA = Structure(cube = True)
pointsB = list(cube_protA.points.copy())
endpt = pointsB.pop(-1)
endpt += np.array([-1,-1,-1])
pointsB.append(endpt)
one_distored_point_cube_protB = Structure(np.array(pointsB),cube_protA.vectors.copy())
cube_protA.display('g')
one_distored_point_cube_protB.display('r')
#translating and rotating randomly
one_distored_point_cube_protB.tumble()
#size of niebourghood

best_neibourhood_find(cube_protA,one_distored_point_cube_protB,7)

cube_protB = Structure(cube_protA.points.copy(),cube_protA.vectors.copy())
cube_protB.tumble()

best_neibourhood_find(cube_protA,cube_protB,7)
