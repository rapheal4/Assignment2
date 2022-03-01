
clearvars
clearvars -GLOBAL
close all
%set(0,'DefaultFigureWindowStyle', 'docked')

global C
global CuCond
global nx ny

C.q_0 = 1.60217653e-19;             % electron charge
C.hb = 1.054571596e-34;             % Dirac constant
C.h = C.hb * 2 * pi;                    % Planck constant
C.m_0 = 9.10938215e-31;             % electron mass
C.kb = 1.3806504e-23;               % Boltzmann constant
C.eps_0 = 8.854187817e-12;          % vacuum permittivity
C.mu_0 = 1.2566370614e-6;           % vacuum permeability
C.c = 299792458;                    % speed of light
C.g = 9.80665;                      % metres (32.1740 ft) per s²

nx = 90;
ny = 65;

CuCond = 1.7e-8;

%Conductivity Map

cMap = zeros(nx,ny);

for i = 1:nx
    for j = 1: ny
        cMap(i,j) = CuCond;
    end
end

G = sparse(nx*ny, nx*ny);
F = zeros(1, nx*ny);

% G - Matrix Formulation

for i = 1:nx
    for j = 1:ny
        n = j + (i - 1) * ny;
        
        if i == 1
            G(n,:) = 0;
            G(n,n) = 1;
             F(n) = 1;
            
        elseif i == nx
             G(n,:) = 0;
             G(n,n) = 1;
        
        elseif j == 1
            nxm = j + (i-2) * ny;
            nxp = j + (i) * ny;
            nyp = (j+1) + (i-1) * ny;
            
            rxm = (cMap(i,j) + cMap(i-1,j))/2.0;
            rxp = (cMap(i,j) + cMap(i+1,j))/2.0;
            ryp = (cMap(i,j) + cMap(i,j+1))/2.0;
            
            G(n,n) = -(rxm + rxp + ryp); 
            G(n,nxm) = rxm;
            G(n,nxp) = rxp;
            G(n, nyp) = ryp;
              
        elseif j == ny
            nxm = j + (i-2) * ny;
            nxp = j + (i) * ny;
            nym = (j-1) + (i-1) * ny;
            
            rxm = (cMap(i,j) + cMap(i-1,j))/2.0;
            rxp = (cMap(i,j) + cMap(i+1,j))/2.0;
            rym = (cMap(i,j) + cMap(i,j-1))/2.0;
            
            G(n,n) = -(rxm + rxp + rym);
            G(n,nxm) = rxm;
            G(n,nxp) = rxp;
            G(n, nym) = rym;
      
        else
            nxm = j + (i-2) * ny;
            nxp = j + (i) * ny;
            nym = (j-1) + (i-1) * ny;
            nyp = (j+1) + (i-1) * ny;
            
            rxm = (cMap(i,j) + cMap(i-1,j))/2.0;
            rxp = (cMap(i,j) + cMap(i+1,j))/2.0;
            rym = (cMap(i,j) + cMap(i,j-1))/2.0;
            ryp = (cMap(i,j) + cMap(i,j+1))/2.0;
            
            G(n,n) = -(rxm + rxp + rym + ryp);
            G(n,nxm) = rxm;
            G(n,nxp) = rxp;
            G(n, nym) = rym;
            G(n, nyp) = ryp;
        end
    end
end
% for i = 1:nx
%    
%     for j = 1:length_y
%         n = j + (i - 1) * length_y;
%         nxm = j + ((i-1) - 1) * length_y;
%         nxp = j + ((i+1) - 1) * length_y;
%         nym = (j-1) + (i - 1) * length_y;
%         nyp = (j+1) + (i - 1) * length_y;
%        
%         if (i == 1)
%             G(n, :) = 0;
%             G(n, n) = 1;
%             V(1, n) = 1;
%         elseif (i == nx)
%             G(n, :) = 0;
%             G(n, n) = 1;
%            
%         elseif (j == 1 && i > 1 && i < nx)
%             G(n, n) = -1;
%             G(n, nyp) = 1;
%            
%         elseif (j == length_y && i > 1 && i < nx)
%            
%             G_matrix(n, n) = -1;
%             G_matrix(n, nym) = 1;
%            
%         else
%             G(n, n) = -4;
%             G(n, nxm) = 1;
%             G(n, nxp) = 1;
%             G(n, nym) = 1;
%             G(n, nyp) = 1;
%            
%         end
%     end
%    
%    
% end
% Using Finite Difference

V = G\F';

for i = 1:nx
    for j = 1:ny
        n = j + (i-1) * ny;
        VG(i,j) = V(n);
    end   
end

figure
set(surf(VG),'linestyle', 'none');
