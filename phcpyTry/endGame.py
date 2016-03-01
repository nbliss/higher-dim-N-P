def setupSystem(polys,expandVar = 'x1'):
    """
    Given a space curve system polys, returns
    f,g,gsols where gsols are solutions of g,
    for performing polyhedral end game
    """
    from phcpy.solver import solve
    randy = randUnitCircle()
    startSyst = polys+[expandVar+ '-' + str(randy) + ';'] 
    startSols = solve(startSyst,silent=True)
    return (polys+[expandVar+';'],startSyst,startSols)

def performEndGame(polys,expandVar = 'x1'):
    """
    Performs a polyhedral end game on the system in polys,
    which should be n-1 equations in n unknowns, with expandVar
    the expanding variable. Returns the winding numbers, coordinates
    of the direction vectors, and the errrors.
    """
    f,g,gsols = setupSystem(polys,expandVar)
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


#f = ['x1^2 - x1 + x2 + x3;', 'x2^2 + x1 + x2 + x3;']
#print setupSystem(f)
l = ['x1^2 - x1 + x2 + x3 + x4 + x5;', 'x2^2 + x1 + x2 + x3 + x4 + x5;', 'x3^2 + x1 + x2 + x3 + x4 + x5;', 'x4^2 + x1 + x2 + x3 + x4 + x5;']
print performEndGame(l)
