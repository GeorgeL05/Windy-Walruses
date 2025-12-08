clear; clc;

%% Parameters
U0 = 20;
Lx = 100000;
H  = 800;
Nx = 4096;
Nz = 101;

x  = linspace(0, Lx, Nx);
z  = linspace(0, H,  Nz);
dz = z(2) - z(1);

%% Turbine forcing f_x(x,z)
A    = 0.5;
x0   = Lx/4;
z0   = 150;
sigx = 50;
sigz = 40;

[X,Z] = meshgrid(x, z);

fx = -A .* exp(-(X-x0).^2/(2*sigx^2)) .* ...
         exp(-(Z-z0).^2/(2*sigz^2));

fz = fx .* (-(Z - z0)/sigz^2);

%% FFT in x
fx_hat = fft(fx, [], 2);
fz_hat = fft(fz, [], 2);

% Wavenumbers corresponding to FFT ordering
kvec = (2*pi/Lx) * [0:(Nx/2-1), -Nx/2:-1];

%% Preallocate Fourier coefficients of w
w_hat = zeros(Nz, Nx);

%% Loop over Fourier modes
for ik = 1:Nx
    k = kvec(ik);
    
    rhs = -1i * k * fz_hat(:,ik);   % Nz x 1 (complex)

    % Build matrix L for this k
    L = zeros(Nz, Nz);

    % Bottom BC: w(0) = 0  -> w_1 = 0
    L(1,1) = 1;
    rhs(1) = 0;

    % Interior points j = 2..Nz-1
    for j = 2:Nz-1
        L(j,j-1) =  U0/dz^2;
        L(j,j)   = -2*U0/dz^2 - U0*k^2;
        L(j,j+1) =  U0/dz^2;
    end

    % Top BC: w_z(H) = 0 -> (w_N - w_{N-1})/dz = 0
    L(Nz, :)    = 0;
    L(Nz, Nz)   =  1;
    L(Nz, Nz-1) = -1;
    rhs(Nz)     = 0;

    % Solve linear system for this Fourier mode
    w_hat(:,ik) = L \ rhs;
end

%% Inverse FFT in x to get w(x,z)
w = ifft(w_hat, [], 2);
w = real(w);

%% Recover u(x,z) from w_hat and plot

%Compute d w_hat / dz for each Fourier mode
Wz_hat = zeros(Nz, Nx);

for ik = 1:Nx
    wmode = w_hat(:,ik);
    dwdz  = zeros(Nz,1);

    % interior points: central difference
    dwdz(2:Nz-1) = ( wmode(3:Nz) - wmode(1:Nz-2) ) / (2*dz);

    % bottom: forward difference
    dwdz(1) = (wmode(2) - wmode(1)) / dz;

    % top: Neumann BC w_z(H)=0
    dwdz(Nz) = 0;

    Wz_hat(:,ik) = dwdz;
end

%Use continuity in Fourier space: u_hat = - (1/(i k)) * d w_hat/dz
u_hat = zeros(Nz, Nx);

for ik = 1:Nx
    k = kvec(ik);
    if abs(k) < 1e-12
        u_hat(:,ik) = 0;
    else
        u_hat(:,ik) = - Wz_hat(:,ik) / (1i * k);
    end
end

%Inverse FFT in x to get u(x,z)
u = ifft(u_hat, [], 2);   % Nz x Nx
u = real(u);              % physical field is real

%% Plot u(x,z) and w(x,z)

% Heatmap of fz
figure;

imagesc(x/1000, z, fz);
set(gca,'YDir','normal');
xlim([23,27]);
xlabel('x (km)');
ylabel('z (m)');
title("fz");
colorbar;


% Heatmap of u
figure;
title("u")
imagesc(x/1000, z, u);
set(gca,'YDir','normal');
xlim([23,27]);
xlabel('x (km)');
ylabel('z (m)');
colorbar;

% Heatmap of w
figure;
title("w")
imagesc(x/1000, z, w);
set(gca,'YDir','normal');
xlim([23,27]);
xlabel('x (km)');
ylabel('z (m)');

colorbar;
