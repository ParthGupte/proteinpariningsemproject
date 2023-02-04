from proteinhelper import *

protA = Structure(random = True, size = 100)
protB = Structure(protA.points,protA.vectors)

#translating and rotating randomly
protB.translate([1,2,3])
protB.rotate(R.from_rotvec(np.random.randint(-5,5,size=3)))
#size of niebourghood
m = 10

neibour_scoredict = {}

for i in range(protA.points.shape[0]):
    for j in range(protB.points.shape[0]):
        NA = make_neibourhood(i,protA,m)[0]
        NB = make_neibourhood(j,protB,m)[0]
        NB.translate(NA.points[0]-NB.points[0])
        if np.all(NA.points[0] != NB.points[0]):
            print(NA.points[0],NB.points[0])
            print("Not True for",(i,j))
        vectorsA = Structure(np.concatenate((NA.vectors[0],-NA.vectors[0]))) #so that the tranlation component of kersley get cancelled and we can work with only rotations
        vectorsB = Structure(np.concatenate((NB.vectors[0],-NB.vectors[0])))
        pair_vecAB = Pairing([0,1,2,3,4,5],vectorsA,vectorsB)
        rot = pair_vecAB.kearsley()[2]
        NB.rotate(rot)
        pair_match = Pairing([x for x in range(m+1)],NA,NB)
        pair_match.hungarian()
        rmsd = pair_match.rmsd()
        if i==j:
            print(pair_match.A.points)
            print("--------------------------")
            print(pair_match.B.points)
            print("----------------------------------------")
            print(pair_match.perm)
        neibour_scoredict[(i,j)] = rmsd

best_match = sorted(neibour_scoredict.items(),key=lambda x:x[1])[0]
print(best_match)



