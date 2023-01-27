from proteinhelper import *
L = []
for i in range(10):
    L.append(np.random.randint(-5,5,size=3))

protA = Structure(np.array(L))
protB = Structure(np.array(L))
protB.translate([1,2,3])
protB.rotate(np.random.randint(-5,5,size=3))

pairing = Pairing([3,4,1,2,5,6,7,8,9,0],protA,protB)
pairing.display('g','r')
pairing.kearsley()
pairing.display('g','r')
# protA.display('r')
# protA.translate([10,10,10])
# protA.display('r')
# protA.rotate(np.array([0,0,1])*mt.pi)
# protA.display(col='r')