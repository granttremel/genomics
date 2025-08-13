
import ggene
gm = ggene.genemap()

testpos = (221837334, 221837437)

tp = (testpos[0] - 10000, testpos[1] + 10000)

for g in gm.fetch(1, 0, features = ''):
    
    print(g['feature'])

# l,c,n = gm.neighbors(1, (tp[0] + tp[1])//2)

# print(l)
# print(c)
# print(n)

