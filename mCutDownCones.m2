needsPackage "gfanInterface";
monInIdeal = I -> (
    m := fold((a,b)->a*b,1,gens ring I);
    varProd := m;
    if saturate(I,m) != 1 then return false;
    while m % I != 0 do (m = m*varProd);
    return m;
);
getRay = vects -> (
    randScale := l -> random(1,1000)*l;
    innerRay := sum apply(vects,randScale);
    innerRay = {0}|innerRay;
    maxy := max(innerRay) + 1;
    --innerRay = append(innerRay,0);
    innerRay = for i in innerRay list maxy-i;
    return innerRay;
);
cutDownCone = (I,vects) -> (
    innerRay := getRay(vects);
    weightedR = newRing(ring I,Variables => flatten({h,gens ring I}),MonomialOrder=>{Weights=>innerRay},Global=>true);
    Ihomog := homogenize(substitute(I,weightedR),weightedR_0);
    xm := monInIdeal(initialIdeal(innerRay,Ihomog));
    if xm==0 then return 0;
    moddy := (xm % Ihomog);
    f := xm - moddy;
    return sub(sub(f,{weightedR_0=>1}),ring I);
);

--R = QQ[x1,x2,x3,x4];
--I = ideal{x1^2 - x1 + x2 + x3 + x4, x2^2 + x1 + x2 + x3 + x4, x3^2 + x1 + x2 + x3 + x4};
--cutDownCone(I,{{1,1,1,1},{0,1,0,0}})
