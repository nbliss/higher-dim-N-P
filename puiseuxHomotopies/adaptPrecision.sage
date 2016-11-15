R.<t,x,y> = QQ[]
F = R*(3*x-10,3*y-10)
Fstar = R*(3*x-10 + 9.1 + 10*t*9.1,3*y-1 + 9.1 + 10*t*9.1)
d1 = 0.3
print F.subs({x:d1,y:d1})
