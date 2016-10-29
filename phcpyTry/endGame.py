def setupSystem(polys,v,expandVar = 'x1'):
    """
    Given a space curve system polys, returns
    f,g,gsols where gsols are solutions of g,
    for performing polyhedral end game
    """
    from phcpy.solver import solve
    randy = randUnitCircle()
    variables = ['x'+str(i) for i in xrange(1,4)]
    #startSyst = polys+[expandVar+ '-' + str(randy) + ';'] 
    extraPoly = [str(v[i])+'*'+variables[i] for i in xrange(3)]
    extraPoly = ' + '.join(extraPoly) 
    targetPoly = [extraPoly + ';']
    extraPoly += ' + 1;'
    extraPoly = [extraPoly]
    print extraPoly
    startSyst = polys+extraPoly
    startSols = solve(startSyst,silent=True)
    #return (polys+[expandVar+';'],startSyst,startSols)
    print (polys+targetPoly,startSyst)
    return (polys+targetPoly,startSyst,startSols)

def performEndGame(polys,v,expandVar = 'x1'):
    """
    Performs a polyhedral end game on the system in polys,
    which should be n-1 equations in n unknowns, with expandVar
    the expanding variable. Returns the winding numbers, coordinates
    of the direction vectors, and the errors.
    """
    f,g,gsols = setupSystem(polys,v,expandVar)
    #from phcpy.solver import mixed_volume as mv
    #from phcpy.solver import random_coefficient_system as rcs
    from phcpy.tuning import order_endgame_extrapolator_set as setOrder
    setOrder(8)
    #from phcpy.trackers import standard_double_track as track
    from phcpy.trackers import double_double_track as track
    #from phcpy.trackers import track
    sols = track(f, g, gsols)
    #from phcpy.tropisms import standard_retrieve as retrieve
    from phcpy.tropisms import dobldobl_retrieve as retrieve
    (w, d, e) = retrieve(len(sols), len(f))
    #return (w,d,e)
    print e
    #d = [roundVect(v) for v in d]
    #n = len(d[0])/2
    #return [[v[i]+v[i+n] for i in xrange(n)] for v in d]
    #return d
    return [roundVect(v) for v in d]

def randUnitCircle():
    """
    Computes a random point on the complex unit circle.
    Returns as a string with 'i' for complex number i.
    """
    from random import random
    from math import cos,sin
    theta = 6.283185307179586*random()
    return str(cos(theta) + complex(0,1)*sin(theta)).replace('j','*i')

def roundToInt(i):
    # Round a float to nearest int
    try:return int(i)
    except Exception:True
    truncated = i.trunc()
    if i - truncated - 0.5 > 0:toReturn = int(truncated)+1
    else:toReturn = int(truncated)
    return toReturn

def roundVect(v):
    return [roundToInt(i) for i in v]


#l = ['x1^2 - x1 + x2 + x3;', 'x2^2 + x1 + x2 + x3;']
l = ['x2^2 - x1 + x2 + x3;', 'x3^2 + x1 + x2 + x3;']
#print setupSystem(f)
#l = ['x1^2 - x1 + x2 + x3 + x4 + x5;', 'x2^2 + x1 + x2 + x3 + x4 + x5;', 'x3^2 + x1 + x2 + x3 + x4 + x5;', 'x4^2 + x1 + x2 + x3 + x4 + x5;']
v=[2,1,1]
print performEndGame(l,v,'x1')
