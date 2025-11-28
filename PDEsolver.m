function [wk] = solver(H, Nz, kstar, fk)
L = zeros(Nz);
F = zeros(Nz);
for i = 2:Nz+1
    L(i-1:i-1) = (Nz/H)^2;
    L(i-1:i) = -kstar^2 - 2*(Nz/H)^2;
    L(i-1:i+1) = (Nz/H)^2;
    F(i-1:i-1) = (Nz/H)^2;
    F(i-1:i+1) = (Nz/H)^2;
end
wk = L\(f*fk);
end

fk = ones(10000,1);
Nz = 10000;
H = 800;
