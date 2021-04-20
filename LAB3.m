% Solve optimization problem
LAB2
% Find optimal LQ controller
Qlqr = diag([6 2 0.2 0.1]);
Rlqr = 1;
Klqr = dlqr(Ad,Bd,Qlqr,Rlqr);
% For simulink:
x = [t' x1 x2 x3 x4];