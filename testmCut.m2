load "mCutDownCones.m2"
badPrevar = n -> (
    R = QQ[x_1..x_n];
    varSum := sum gens R;
    polyList = for i from 1 to n-1 list (
        if i==1 then x_1^2 + varSum - 2*x_1
        else x_i^2 + varSum
    );
    return ideal polyList;
);
for n from 3 to 12 do (
    vects = {for i from 1 to n list 1,{1}|for i from 2 to n list 0};
    I = badPrevar n;
    print n;
    for i from 1 to 3 do (
        time cutDownCone(I,vects);
    );
);
