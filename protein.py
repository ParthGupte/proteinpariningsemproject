from proteinhelper import *
L = []
for i in range(10):
    L.append(Point(np.random.randint(-5,5,size=3)))

protA = Protein(L)

protA.display('g')
protA.translate(Vec([10,10,10]))
protA.display('g')
protA.rotate(Vec([0,0,1]),mt.pi)
protA.display(col='g')