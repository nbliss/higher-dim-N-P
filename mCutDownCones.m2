--needsPackage "gfanInterface";
innerProd = (v,w) -> (
    return sum(for i from 0 to (#v-1) list v_i*w_i);
);
initialForm = (f,w) -> (
    fTerms := terms f;
    exps := for t in fTerms list (exponents t)_0;
    maxVal := innerProd(w,exps_0);
    maxList := new MutableList from {fTerms_0};
    for i from 1 to (#exps-1) do (
        prodVal := innerProd(exps_i,w);
        if prodVal > maxVal then (
            maxVal = prodVal;
            maxList = {fTerms_i};
        )
        else if prodVal == maxVal then maxList = append(maxList,fTerms_i);
    );
    return sum maxList;
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
monInIdeal = I -> (
    m := fold((a,b)->a*b,1,gens ring I);
    varProd := m;
    if saturate(I,m) != 1 then return 0;
    while m % I != 0 do (m = m*varProd);
    return m;
);
cutDownCone = (I,vects) -> (
    innerRay := getRay(vects);
    varList := flatten({h,gens ring I});
    weightedR = newRing(ring I,Variables => varList,MonomialOrder=>{Weights=>innerRay});
    Ihomog := homogenize(substitute(I,weightedR),weightedR_0);
    initialIhomog = ideal leadTerm(1,Ihomog);
    xm := monInIdeal(initialIhomog);
    if xm==0 then return 0;
    moddy := (xm % Ihomog);
    f := xm - moddy;
    return sub(sub(f,{weightedR_0=>1}),ring I);
);

--R = QQ[x1,x2,x3,x4];
--I = ideal{x1^2 - x1 + x2 + x3 + x4, x2^2 + x1 + x2 + x3 + x4, x3^2 + x1 + x2 + x3 + x4};
--cutDownCone(I,{{1,1,1,1},{0,1,0,0}})
