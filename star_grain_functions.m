
% Geometric relationship between primary angles, epsilon and theta
theta = @(Rp, Ri, e, N) 2*atan((Rp*sin(pi*e/N)*tan(pi*e/N)) ...
    / (Rp * sin(pi*e/N) - Ri * tan(pi*e/N)));

% Three minor arcs burning at a given time
S = @(N, S1, S2, S3) 2*N*(S1+S2+S3); %total burning perimeter
S1 = @(N, Rp, e, o, y, f) Rp * sin(pi*e/N) / sin(o/2) - (y+f)*cot(o/2);
S2 = @(N, e, o, y, f) (y+f)*(pi/2 - o/2 + pi*e/N)
S3 = @(N, Rp, e, y, f) (Rp+y+f)*(pi/N - pi*e/N)

A_phase1 = @(N, Rp, e, th, y, f) 2*N*(1/2*Rp*sin(pi*e/N)*(Rp*cos(pi*e/N) + Rp*sin(pi*e/N)*tan(pi*e/N)) ...
    -1/2*(Rp*sin(pi*e/N)/sin(th/2) - (y+f)*cot(th/2))^2 * tan(th/2) ...
    +1/2*(y+f)^2 * (pi/2-th/2+pi*e/N) + 1/2*(Rp+y+f)^2 * (pi/N - pi*e/N))

phase1_check = @(N, Rp, e, th, y, f) (y+f) > Rp * sin(pi*e/N) / cos(th/2);