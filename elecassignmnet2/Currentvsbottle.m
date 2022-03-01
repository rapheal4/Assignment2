clearvars
clearvars -GLOBAL
close all
set(0,'DefaultFigureWindowStyle', 'docked')


global CuCond NoCond
global nx ny


nx = 80;
ny = 80;

currentset = 150;

Current = zeros(1,currentset);
thenecksizes = zeros(1,currentset);

Lb = floor(nx/3);
Wb = floor(ny/3);

for K = 1:currentset
    
    thenecksizes(K) = Wb;
    
    CuCond = 1;
    NoCond = 10e-2;
    
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
            
            if (i<=ny && i>=(ny-Wb) && j>Lb && j<=(2*Lb))
                conductivityMap(i,j) = NoCond;
            end
        end
    end
    
    G = sparse(nx*ny, nx*ny);
    F = zeros(1, nx*ny);
    
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
                
                rxm = (conductivityMap(i,j) + conductivityMap(i-1,j))/2.0;
                rxp = (conductivityMap(i,j) + conductivityMap(i+1,j))/2.0;
                ryp = (conductivityMap(i,j) + conductivityMap(i,j+1))/2.0;
                
                G(n,n) = -(rxm + rxp + ryp);
                G(n,nxm) = rxm;
                G(n,nxp) = rxp;
                G(n, nyp) = ryp;
                
            elseif j == ny
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
    
    for i = 1:nx
        for j = 1:ny
            n = j + (i-1) * ny;
            VG(i,j) = V(n);
        end
    end
    
    for i = 1:nx
        for j = 1:ny
            if i == 1
                Ex(i, j) = (VG(i + 1, j) - VG(i, j));
            elseif i == nx
                Ex(i, j) = (VG(i, j) - VG(i - 1, j));
            else
                Ex(i, j) = (VG(i + 1, j) - VG(i - 1, j)) * 0.5;
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
    
    Jx = conductivityMap .* Ex;
    Jy = conductivityMap .* Ey;
    
    C0 = sum(Jx(1, :));
    Cnx = sum(Jx(nx, :));
    
    
    Current(K) = (C0 + Cnx) * 0.5;
    
   
    thenecksizes(K) = ny - (2*Wb);
    
   
    Wb = Wb+1;
end

%Current vs Neck Size
figure
plot(thenecksizes, Current,'g');
title('The Current vs bottle-necks')
xlabel('The Bottle neck width')
ylabel('The Current (A)')