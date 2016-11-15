import numpy as np
m = np.array([[[1,2,3],[1,2,3]],[[1,2,3],[1,2,3]]])
n = np.array([[[1,1,1]],[[1,1,1]]])
for i in xrange(4):
    for j in xrange(4):
            try:
                    print np.tensordot(m,n,axes=[i-1,j-1])
                    print i,j
            except Exception as e:print e.message

