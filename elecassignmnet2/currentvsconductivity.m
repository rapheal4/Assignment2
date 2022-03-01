clearvars
clearvars -GLOBAL
close all
set(0,'DefaultFigureWindowStyle', 'docked')


global CuCond NoCond
global nx ny



nx = 80;
ny = 80;

Conductivityset = 150;

newCurrent = zeros(1,Conductivityset);
CuCondVal = zeros(1,Conductivityset);

Lb = floor(nx/3);
Wb = floor(ny/3);

for i = 1:Conductivityset
    CuCondVal(i) = i;
end

for u = 1:Conductivityset
    
    CuCond = CuCondVal(u);
    NoCond = 10e-2;
    
    
    
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
    newset = G\F';
    
    for i = 1:nx
        for j = 1:ny
            n = j + (i-1) * ny;
            VG(i,j) = newset(n);
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
    
    
    newCurrent(u) = (sum(Jx(1, :)) + sum(Jx(nx, :))) * 0.5;

end

%Current vs Neck Size
figure
plot(CuCondVal, newCurrent, 'g');
title('The Current vs Conductivity')
xlabel('Conductivity')
ylabel('Current (A)')




