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
vects = n -> (
    return {for i in 1..n list 1,{1}|for i in 2..n list 0};
);
noetherBadPrevar = n -> (
    R = QQ[x_1..x_n];
    varSum := sum gens R;
    polyList = for i from 1 to n-1 list (
        if i==1 then x_2^2 + varSum - 2*x_1
        else x_(i+1)^2 + varSum
    );
    return ideal polyList;
);
test = n -> (
    I = badPrevar n;
    thisCone = vects(n);
    for i from 1 to 3 do (
        time cutDownCone(I,thisCone);
    );
);
