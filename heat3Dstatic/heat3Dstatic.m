%-------------------------------------------------------------------------%
%
% FEM code for the 3D static heat equation solution
%
% Last update: 06/05/2017
%
% MATLAB version 2016
%
%-------------------------------------------------------------------------%

clc
clear variables
close all

%-------------------------------------------------------------------------%
%
% All length units are micrometers
%
% Parameters imputs
%
T = 300; % Room temperature in Kelvin
%
% Material description Silicon
k_si = 148e-6; %1.45e-4; % Whatt/(micrometers*Kelvin) thermal conductivity
%
% trench thermal conductivity in Whatt/(micrometers*Kelvin)
k_fill = 1.4e-6; %5e-6; % range [1e-6; 5e-6]
%
%-------------------------------------------------------------------------%
%
% Heat source description
%
% x
% heat source width
W = 0.1; %1; % range [0.1, 2]
%
% z
% heat source thickness = 0.12; % fixed value
heat_source_thickness = 0.12; % fixed value
% heat_source_thickness = t
%
% y
% heat source height = 5; % range [2.5, 10]
L = 2.5; %5; % range [2.5, 10]
%
%
% heat source center coordinates
X_center = 250;  Y_center = 250; Z_center = 1.06; %1;
%
%-------------------------------------------------------------------------%
%
% trench geometry definition
%
% x and y
W_t = 1; % fixed value 1
%
% x
S_W = 4; % range [0.5; 4]
%
% y
S_L = 2; % fixed value 2
%
% z
d_t = 5; % range [2.5, 5]
%
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
%
% heat source geometry and units
%
heat_source_volume = W*heat_source_thickness*L; % in micrometers^3
%
volumetric_heat_flux = 1 / heat_source_volume; % in Watt/micrometers^3
%
source_intensity = volumetric_heat_flux;
%
Xmin = X_center - 0.5*W;
Xmax = X_center + 0.5*W;
Ymin = Y_center - 0.5*L;
Ymax = Y_center + 0.5*L;
Zmin = Z_center - 0.5*heat_source_thickness;
Zmax = Z_center + 0.5*heat_source_thickness;
%
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
%
% trench geometry
%
% trench limits
X1_trench = Xmin - S_W - W_t;
X2_trench = Xmin - S_W;
X3_trench = Xmax + S_W;
X4_trench = Xmax + S_W + W_t;
Y1_trench = Ymin - S_L - W_t;
Y2_trench = Ymin - S_L;
Y3_trench = Ymax + S_L;
Y4_trench = Ymax + S_L + W_t;
%
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
%
% Grid generation
%
x_grid = [linspace(0, 200, 5), linspace(200, X1_trench, 5), linspace(X1_trench, X2_trench, 5), ...
      linspace(X2_trench, Xmin, 5), linspace(Xmin, Xmax, 5), linspace(Xmax, X3_trench, 5), ...
      linspace(X3_trench, X4_trench, 5), linspace(X4_trench, 300, 5), linspace(300, 500, 5)];
x_grid = unique(x_grid);
y_grid = [linspace(0, 200, 5), linspace(200, Y1_trench, 5), linspace(Y1_trench, Y2_trench, 5), ...
      linspace(Y2_trench, Ymin, 5), linspace(Ymin, Ymax, 5), linspace(Ymax, Y3_trench, 5), ...
      linspace(Y3_trench, Y4_trench, 5), linspace(Y4_trench, 300, 5), linspace(300, 500, 5)];
y_grid = unique(y_grid);
z_grid = [linspace(0, Zmin, 5), linspace(Zmin, Zmax, 5), linspace(Zmax, d_t, 5), ...
     linspace(d_t, 2*d_t , 5), linspace(2*d_t, 3*d_t, 5), linspace(3*d_t, 100, 5), linspace(100, 500, 5)];
z_grid = unique(z_grid);


%
z_coordinate = z_grid;
%
% Heat source index center node
idx_x_source = find((x_grid-X_center == 0));
idx_y_source = find((y_grid-Y_center == 0));
idx_z_source = find((z_grid-Z_center == 0));
%
m = length(x_grid);
n = length(y_grid);
l = length(z_grid);

[x, y, z] = ndgrid(x_grid, y_grid, z_grid);

% Grid visualization options
r = plot_grid_2d(x, y, z);
r = plot_grid_3d(x, y, z);

p = [x(:), y(:), z(:)]; % Matrix listing x,y,z coordinates of all n*m*l nodes

N_nodes = m*n*l; % number of nodes m*n*l
N_cubes = (size(x,1)-1)*(size(y,2)-1)*(size(z,3)-1); % Total number of cubes

% t is a matrix with the index of each cube nodes
t = cube_nodes(m-1, n-1, l-1);
t = t(:,2:9);

% heat source location
ix_source = (Xmin<=x(t)) .* (x(t)<=Xmax);
iy_source = (Ymin<=y(t)) .* (y(t)<=Ymax);
iz_source = (Zmin<=z(t)) .* (z(t)<=Zmax);
is_source = ix_source.*iy_source.*iz_source;
is_source = is_source(:,1).*is_source(:,2).*is_source(:,3).*is_source(:,4).* ...
            is_source(:,5).*is_source(:,6).*is_source(:,7).*is_source(:,8);

% trench location
%ix_trench = or( (X1_trench<=x(t)).*(x(t)<=X2_trench) , (X3_trench<=x(t)).*(x(t)<=X4_trench) );
%iy_trench = or( (Y1_trench<=y(t)).*(y(t)<=Y2_trench) , (Y3_trench<=y(t)).*(y(t)<=Y4_trench) );

% trench 4
% ix_trench = ( (X2_trench<=x(t)).*(x(t)<=X3_trench) );
% iy_trench = ( (Y3_trench<=y(t)).*(y(t)<=Y4_trench) );

% trench 1
ix_trench = ( (X1_trench<=x(t)).*(x(t)<=X2_trench) );
iy_trench = ( (Y1_trench<=y(t)).*(y(t)<=Y4_trench) );

is_trench = ix_trench.*iy_trench;

% trench 2
ix_trench = ( (X2_trench<=x(t)).*(x(t)<=X3_trench) );
iy_trench = ( (Y1_trench<=y(t)).*(y(t)<=Y2_trench) );

is_trench = or(is_trench, ix_trench.*iy_trench);

% trench 3
ix_trench = ( (X3_trench<=x(t)).*(x(t)<=X4_trench) );
iy_trench = ( (Y1_trench<=y(t)).*(y(t)<=Y4_trench) );

is_trench = or(is_trench, ix_trench.*iy_trench);

% trench 4
ix_trench = ( (X2_trench<=x(t)).*(x(t)<=X3_trench) );
iy_trench = ( (Y3_trench<=y(t)).*(y(t)<=Y4_trench) );

is_trench = or(is_trench, ix_trench.*iy_trench);


iz_trench = (z(t)<=d_t);
is_trench = is_trench.*iz_trench;

is_trench = is_trench(:,1).*is_trench(:,2).*is_trench(:,3).*is_trench(:,4).* ...
            is_trench(:,5).*is_trench(:,6).*is_trench(:,7).*is_trench(:,8);

thermal_cond = repmat(k_si, N_cubes, 1);
thermal_cond(find(is_trench)) = k_fill;

% Boundary nodes selection
plane_x_0 = [];
plane_x_m = [];
plane_y_0 = [];
plane_y_n = [];

for j = 1:l
    plane_x_0 = [plane_x_0, (1:m:n*m) + (j-1)*m*n];
    plane_x_m = [plane_x_m, (m:m:n*m) + (j-1)*m*n];

    plane_y_0 = [plane_y_0, (1:m) + (j-1)*m*n];
    plane_y_n = [plane_y_n, (((n-1)*m+1):m*n) + (j-1)*m*n];
end

plane_z_0 = 1:m*n;
plane_z_l = m*n*(l-1)+1:m*n*l;
% b is a vector listing each boundary node index
b = [plane_x_0, plane_x_m, plane_y_0, plane_y_n, plane_z_0, plane_z_l];

tic()
% Matrix calculations
% K is the matrix of the integral of the basis functions gradients.
K = sparse(N_nodes, N_nodes); % zero matrix in sparse format: zeros(N) would be "dense"
% F is inhomogeneity of the linear system which contains boundary conditions
F = zeros(N_nodes, 1); %

Kx = matrix_Kx(0, 0, 0, 1, 1, 1);
Ky = matrix_Ky(0, 0, 0, 1, 1, 1);
Kz = matrix_Kz(0, 0, 0, 1, 1, 1);

% gradients matrix
x0 = p(t(:,1),1); y0 = p(t(:,1),2); z0 = p(t(:,1),3);
x1 = p(t(:,8),1); y1 = p(t(:,8),2); z1 = p(t(:,8),3);

Kge = (y1-y0).*(z1-z0)./(x1-x0)*Kx + ...
      (x1-x0).*(z1-z0)./(y1-y0)*Ky + ...
      (x1-x0).*(y1-y0)./(z1-z0)*Kz;

Kge = repmat(thermal_cond, 1, 64) .* Kge;

Ige = t(:, kron(1:8, ones(1,8)));
Jge = t(:, kron(ones(1,8), 1:8));

K = sparse(Ige(:), Jge(:), Kge(:), N_nodes, N_nodes);

% inhomogeneity vector
S = vector(0, 0, 0, 1, 1, 1);
volume = (x1-x0).*(y1-y0).*(z1-z0);
Fge = (is_source .* volume) * S;

F = accumarray([t(:), ones(length(t(:)),1)], Fge(:), [N_nodes 1]);

% Source with adecuate units
F = F*source_intensity;
toc()

display(' ')
display('Problem:')
display('   -k(x,y,z) laplacian(U) = source(x,y,z)')
display('   source(x,y,z) with rectangular shape in [Xmin,Xmax]x[Ymin,Ymax]x[Zmin,Zmax] coordinates')
display('   the boundary conditions are adiabatic over all cubes faces but U(z=l) = T.')
display(' ')

% Setting the boundary conditions in K and F
K(plane_z_l,:) = 0;
K(plane_z_l,plane_z_l) = speye(length(plane_z_l));
F(plane_z_l) = T;

%
% Linear equation system resolution
U = K\F;
% U = agmg(K,F);
%
% normalization
U = U - T*ones(size(U)); % Normalized temperature rise in Kelvin/Watt

% %
% % Volumetric average of the normalized temperature inside the heat source
% %
% % Normalized temperature (NT) rize at the source nodes in Kelvin/Watt
% source_te = U(t(find(is_source),:));
% %
% % Average NT at each parallelepiped inside the heat source
% mean_cube_te = mean(source_te')'; % in Kelvin/Watt
% %
% % Volumetric average NT at each parallelepiped inside the heat source
% mean_vol_cube_te = mean_cube_te./volume(find(is_source)); % in Kelvin/(Watt*micrometer^3)
% %
% % Volumetric average of the NT inside the heat source (all source
% % parallelepiped are averaged)
% vol_average_temp = mean(mean_vol_cube_te); % in Kelvin/(Watt*micrometer^3)

%% comparacion con codigo del Alessandro
% Volumetric average of the normalized temperature inside the heat source
%
% Normalized temperature (NT) rize at the source nodes in Kelvin/Watt
source_te = U(t(find(is_source),:));
%
% Average NT at each parallelepiped inside the heat source
mean_cube_te = mean(source_te')'; % in Kelvin/Watt
%
% mean teamperature rize
avarage_temp = mean(mean_cube_te)
%%

% disp('Volumentric average of the normalized temperature rise inside the')
% disp('heat source in Kelvin/(Watt*micrometer^3):')
% disp(vol_average_temp)
% disp(' ')

% re-arrenge U in a mxnxl matrix
U = reshape(U, m, n, l);

% set tempurature limits for visualization
Tmin = min(U(:));
Tmax = max(U(:));
%
figure
plot(squeeze(x(:,idx_y_source,idx_z_source)), squeeze(U(:,idx_y_source,idx_z_source)),'-o')
xlabel('x axis [micrometers]');
ylabel('Normalized temperature rise [K/W]');
title('Normalized temperature rise along x')

figure
plot(squeeze(z(idx_x_source,idx_y_source,:)), squeeze(U(idx_x_source,idx_y_source,:)),'-o')
xlabel('z axis [micrometers]');
ylabel('Normalized temperature rise [K/W]');
title('Normalized temperature rise along z');

% % uncomment the following lines for result visualization
% figure
% for c = 1:size(U,3);
%   mesh(x(:,:,1), y(:,:,1), squeeze(U(:,:,c)) )
%   caxis([Tmin Tmax])
%   % colorbar;
%   zlim([Tmin Tmax]);
%   xlim([0 500]);
%   ylim([0 500]);
%   xlabel('x axis [micrometers]');
%   ylabel('y axis [micrometers]');
%   zlabel('Normalized temperature rise [K/W]');
%
%   graph_title = sprintf('Plane z = %6.4f micrometers', z_coordinate(c));
%   title(graph_title);
%   pause;
% end
% close

%display(' ')
display(' FEM resolution finished OK.')
display(' ')
