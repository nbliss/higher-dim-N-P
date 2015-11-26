import sys

def compareCones(R):
    sys.stdout.flush()
    n1,n2 = ZZ.random_element(2,5),ZZ.random_element(2,5)
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


R.<x,y,z> = QQ[]
for i in xrange(100):
    timeout(compareCones,args={R},timeout_duration=5)

