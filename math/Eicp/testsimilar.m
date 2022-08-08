shiftLambda = 5;
n = 6;
G = [0 1 0 1 0 0;1 0 0 1 0 0;0 0 0 1 0 0;1 1 1 0 1 1;0 0 0 1 0 1;0 0 0 1 1 0] + shiftLambda * eye(n);
H = [0 1 0 0 0 0;1 0 1 1 0 0;0 1 0 1 1 0;0 1 1 0 1 0;0 0 1 1 0 1;0 0 0 0 1 0] + shiftLambda * eye(n);
[flex, fed, fisd] = similarity(G, H, 0.01, shiftLambda)