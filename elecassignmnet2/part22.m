clearvars
clearvars -GLOBAL
close all
set(0,'DefaultFigureWindowStyle', 'docked')

global C
global CuCond NoCond
global nx ny

% C.q_0 = 1.60217653e-19;             % electron charge
% C.hb = 1.054571596e-34;             % Dirac constant
% C.h = C.hb * 2 * pi;                    % Planck constant
% C.m_0 = 9.10938215e-31;             % electron mass
% C.kb = 1.3806504e-23;               % Boltzmann constant
% C.eps_0 = 8.854187817e-12;          % vacuum permittivity
% C.mu_0 = 1.2566370614e-6;           % vacuum permeability
% C.c = 299792458;                    % speed of light
% C.g = 9.80665;                      % metres (32.1740 ft) per s²
% 
nx = 30;
ny = 30;

Lb = floor(nx/3);
Wb = floor(ny/3);


CuCond = 100;
NoCond = 10e-9;

%Conductivity map


conductivityMap = zeros(nx,ny);

for i = 1:nx
    for j = 1: ny
        conductivityMap(i,j) = CuCond;
    end
end
for i = 1:nx
    for j = 1:ny
        if (i>=1 && i<=Wb && j>Lb && j<=(2*Lb))
            conductivityMap(i,j) = NoCond;
        end
%         if(i > 1 && i < wb && ((j < Lb || (j > 2*Lb)))
%             conductivity_mapping(i,j) = Nocond;
%         end
        if (i<=ny && i>=(ny-Wb) && j>Lb && j<=(2*Lb))
            conductivityMap(i,j) = NoCond;
        end
    end
end 

G = sparse(nx*ny, nx*ny);
F = zeros(1, nx*ny);

% 


for i = 1:nx
    for j = 1:ny
        n = j + (i - 1) * ny;
        %n1 = j + ((i+1) - 1) * ny;
        %n2 = (j-1) + (i - 1) * ny;
       % n3 = (j+1) + (i - 1) * ny;
        
        if i == 1
            G(n,:) = 0;
            G(n,n) = 1;
             F(n) = 1;
            
        elseif i == nx %(i == 1 && i > 1 && i < nx)
             G(n,:) = 0;
             G(n,n) = 1;
        
        elseif j == 1 
            nxm = j + (i-2) * ny;
            nxp = j + (i) * ny;
            nyp = (j+1) + (i-1) * ny;
            
            rxm = (conductivityMap(i,j) + conductivityMap(i-1,j))/2.0;
            rxp = (conductivityMap(i,j) + conductivityMap(i+1,j))/2.0;
            ryp = (conductivityMap(i,j) + conductivityMap(i,j+1))/2.0;
            
            G(n,n) = -(rxm + rxp + ryp); 
            G(n,nxm) = rxm;
            G(n,nxp) = rxp;
            G(n, nyp) = ryp;
              
        elseif j == ny %(j == ny && i > 1 && i < nx)
            nxm = j + (i-2) * ny;
            nxp = j + (i) * ny;
            nym = (j-1) + (i-1) * ny;
            
            rxm = (conductivityMap(i,j) + conductivityMap(i-1,j))/2.0;
            rxp = (conductivityMap(i,j) + conductivityMap(i+1,j))/2.0;
            rym = (conductivityMap(i,j) + conductivityMap(i,j-1))/2.0;
            
            G(n,n) = -(rxm + rxp + rym);
            G(n,nxm) = rxm;
            G(n,nxp) = rxp;
            G(n, nym) = rym;
      
        else
            nxm = j + (i-2) * ny;
            nxp = j + (i) * ny;
            nym = (j-1) + (i-1) * ny;
            nyp = (j+1) + (i-1) * ny;
            
            rxm = (conductivityMap(i,j) + conductivityMap(i-1,j))/2.0;
            rxp = (conductivityMap(i,j) + conductivityMap(i+1,j))/2.0;
            rym = (conductivityMap(i,j) + conductivityMap(i,j-1))/2.0;
            ryp = (conductivityMap(i,j) + conductivityMap(i,j+1))/2.0;
            
            G(n,n) = -(rxm + rxp + rym + ryp);
            G(n,nxm) = rxm;
            G(n,nxp) = rxp;
            G(n, nym) = rym;
            G(n, nyp) = ryp;
        end
    end
end
V = G\F';
vmap= zeros (nx,ny);
for i = 1:nx %Converting V to matrix to plot
    for j = 1:ny
        n = j + (i-1) * ny;
        VG(i,j) = V(n);
    end   
end



N =150;

a = ny;
b = nx/2;
new = zeros(nx,ny);

for i = 1:nx %Electric Field calculation
    for j = 1:ny
        if i == 1
            Ex(i, j) = (VG(++i, j) - VG(i, j));
        elseif i == nx
            Ex(i, j) = (VG(i, j) - VG(--i, j));
        else
            Ex(i, j) = (VG(i + 1, j) - VG(--i, j)) * 0.5;
        end
        if j == 1
            Ey(i, j) = (VG(i, j + 1) - VG(i, j));
        elseif j == ny
            Ey(i, j) = (VG(i, j) - VG(i, j - 1));
        else
            Ey(i, j) = (VG(i, j + 1) - VG(i, j - 1)) * 0.5;
        end
    end
end

Ex = -Ex;
Ey = -Ey;

Densitymapx = conductivityMap .* Ex;
Densitymapy = conductivityMap .* Ey;

figure

H = surf(conductivityMap');
title('The Conductivity Map')
%set(H, 'linestyle', 'none');
  set(gcf,'DefaultTextColor','r')
view(0, 90) % 2D view

figure 
H = surf(VG');
title('The Vmap with bottle neck')
set(H, 'linestyle', 'none');
 %set(gcf,'DefaultTextColor','r')
view(0, 90) %2D view

figure 
quiver(Ex', Ey');
title('The Electric field Map')
axis([0 nx 0 ny]);

figure
quiver(Densitymapx', Densitymapy');
title(' The Current Density Map')
axis([0 nx 0 ny]);

Current1 = sum(Densitymapx(1, :));
Current2 = sum(Densitymapx(nx, :));
total = (Current1 + Current2) * 1;

%