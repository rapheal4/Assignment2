clc
clear
set(0,'DefaultFigureWindowStyle','docked')


nx=90;
ny=65;


G = sparse((nx * ny), (nx * ny));
V = zeros(1, (nx * ny));
vo = 1;

for i = 1:nx

    for j = 1:ny
        n = j + (i - 1) * ny;
        nxm = j + ((i-1) - 1) * ny;
        nxp = j + ((i+1) - 1) * ny;
        nym = (j-1) + (i - 1) * ny;
        nyp = (j+1) + (i - 1) * ny;
       
        if i == 1
            G(n, :) = 0;
            G(n, n) = 1;
            V(1, n) = vo;
        elseif i == nx
            G(n, :) = 0;
            G(n, n) = 1;
            V(1, n) = vo;
        elseif j == 1
            G(n, :) = 0;
            G(n, n) = 1;
        elseif j == ny
            G(n, :) = 0;
            G(n, n) = 1;
        else
            G(n, :) = 0;
            G(n, n) = -4;
            G(n, nxm) = 1;
            G(n, nxp) = 1;
            G(n, nym) = 1;
            G(n, nyp) = 1;
        end
    end
               

end

NewV = G\V';
figure (3);

surface2 = zeros(nx, ny);


for i = 1:nx

    for j = 1:ny
        n = j + (i - 1) * ny;
        nxm = j + ((i-1) - 1) * ny;
        nxp = j + ((i+1) - 1) * ny;
        nym = (j-1) + (i - 1) * ny;
        nyp = (j+1) + (i - 1) * ny;
        surface2(i, j) = NewV(n);
    end
end
   
surf(surface2);
title('The Numerical Approach');


%newV =zeros(90,65);
newV = zeros(65, 30); 
a = 65;
b = 25;

x = linspace(-25,25,30);
y = linspace(0,65,65);

[x_mesh, y_mesh] = meshgrid(x, y);

for n = 1:2:200
   
    newV =  (newV + (4 * vo/pi).*(cosh((n * pi * x_mesh)/a) .* sin((n * pi * y_mesh)/a)) ./ (n * cosh((n * pi * b)/a)));
   
    figure(4);
    surf(x, y, newV);
    title('The Analytical Approach Solution');
    pause(0.05);
   
end