import sys

def randomPolys(n):
    """
    Given n, checks rays of trop prevariety vs variety for
    n-1 polys in n unknowns.
    """
    return

def projectedRays(I,n):
    """
    Compute the n-th elimination ideal, then
    find the rays of its tropical prevariety.
    """
    R = I.ring()
    nthvar = R.gens()[n]
    projI = R.remove_var(nthvar)*I.elimination_ideal(nthvar)
    return projI.groebner_fan().tropical_intersection().rays()

def randomCoeffs(I,height_bound=300):
    """
    Given an ideal I, checks rays of trop prevariety vs
    variety for random choices of the coefficients of the
    gens of I.
    """
    d_i = []
    for p in I.gens():
        d_i.append({a:QQ.random_element(height_bound,height_bound) for a in p.dict().keys()})
    newI = R*d_i
    prevariety = newI.groebner_fan().tropical_intersection()
    variety = findTropVar(newI)
    if not prevariety.rays()==variety.rays():
        return newI
    return True


def findTropVar(J):
    """
    Uses gfan. Potentially really slow.
    """
    F = J.groebner_fan().tropical_basis()
    J = R*F
    return J.groebner_fan().tropical_intersection()

def compareCones(n):
    sys.stdout.flush()
    n_i = [ZZ.random_element(2,5) for i in xrange(n)]
    n1,n2 = ZZ.random_element(2,5),ZZ.random_element(2,5)

    #R = PolynomialRing(
    p1,p2 = R.random_element(n1),R.random_element(n2)
    print '-'*40
    if len(p1.dict())==1 or len(p1.dict())==1:return
    print p1
    print p2
    I = R*(p1,p2)
    F = I.groebner_fan()
    f1 = F.tropical_intersection()
    printEm = False
    try:f2 = (I.ring()*F.tropical_basis()).groebner_fan().tropical_intersection()
    except Exception as e:
        print e
        return
    if f1.cones()[1]!=f2.cones()[1]:
        print n1,n2
        print f1.rays()
        print f2.rays()
        print
        print f1.cones()
        print f2.cones()
        print '-'*40

def timeout(func, args=(), kwargs={}, timeout_duration=1, default=None):
    import signal

    class TimeoutError(Exception):
        pass

    def handler(signum, frame):
        raise TimeoutError()

    # set the timeout handler
    signal.signal(signal.SIGALRM, handler) 
    signal.alarm(timeout_duration)
    try:
        result = func(*args, **kwargs)
    except TimeoutError as exc:
        result = default
    finally:
        signal.alarm(0)

    return result

"""
R.<x,y,z> = QQ[]
for i in xrange(100):
    timeout(compareCones,args={R},timeout_duration=5)
"""
