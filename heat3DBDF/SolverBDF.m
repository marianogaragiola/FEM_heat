tic

addpath(genpath('./src/'));

TemplateBuild

%Inputs
%in.media_properties

%in.media_positions
IND_MEDIAPOS_x=2;
IND_MEDIAPOS_y=3;
IND_MEDIAPOS_z=4;
IND_MEDIAPOS_Dx=5;
IND_MEDIAPOS_Dy=6;
IND_MEDIAPOS_Dz=7;

%in.sources_position
IND_SOURCEPOS_x=2;
IND_SOURCEPOS_y=3;
IND_SOURCEPOS_z=4;
IND_SOURCEPOS_Dx=5;
IND_SOURCEPOS_Dy=6;
IND_SOURCEPOS_Dz=7;
IND_SOURCEPOS_Val=8;


%BC baseplate
T = 300; %[K]
%Would that be better to scale it? Or is it ok as it is? It is a general
%preconditioning/scaling problem...

% Boundary nodes selection: keep only the top (needed) one!
plane_z_l = m*n*(l-1)+1:m*n*l;

%% Material selection

%Initialize everything to first block (substrate)
thermal_cond = repmat(in.media_properties(in.media_positions(1,end),2)/ ...
                     (in.media_properties(in.media_positions(1,end),3)* ...
                      in.media_properties(in.media_positions(1,end),4)), N_cubes, 1);

%Then edit for following blocks
for k=2:size(in.media_positions,1)
   ix_condition = (in.media_positions(k,IND_MEDIAPOS_x)<=x(t)) &...
                  (x(t)<=(in.media_positions(k,IND_MEDIAPOS_x)+in.media_positions(k,IND_MEDIAPOS_Dx)));
   iy_condition = (in.media_positions(k,IND_MEDIAPOS_y)<=y(t)) &...
                  (y(t)<=(in.media_positions(k,IND_MEDIAPOS_y)+in.media_positions(k,IND_MEDIAPOS_Dy)));
   iz_condition = (in.media_positions(k,IND_MEDIAPOS_z)<=z(t)) &...
                  (z(t)<=(in.media_positions(k,IND_MEDIAPOS_z)+in.media_positions(k,IND_MEDIAPOS_Dz)));
   is_condition_full = ix_condition & iy_condition & iz_condition;
   is_condition=is_condition_full(:,1) & is_condition_full(:,2) & is_condition_full(:,3) & ...
                is_condition_full(:,4) & is_condition_full(:,5) & is_condition_full(:,6) & ...
                is_condition_full(:,7) & is_condition_full(:,8);
   thermal_cond(is_condition) = in.media_properties(in.media_positions(k,end),2)/ ...
                               (in.media_properties(in.media_positions(k,end),3)* ...
                                in.media_properties(in.media_positions(k,end),4));
end

%% Matrix assembly

%Cubes - left bottom corner
x0 = p(t(:,1),1);
y0 = p(t(:,1),2);
z0 = p(t(:,1),3);
%Cubes - top right corner
x1 = p(t(:,8),1);
y1 = p(t(:,8),2);
z1 = p(t(:,8),3);
%Cubes - volume
volume = (x1-x0).*(y1-y0).*(z1-z0);

% Matrix calculations
% K is the matrix of the integral of the basis functions gradients.
K = sparse(N_nodes, N_nodes); % zero matrix in sparse format: zeros(N) would be "dense"
% F is inhomogeneity of the linear system which contains boundary conditions
F = zeros(N_nodes, 1); %

%Construct gradient matrix Kx,Ky,Kz for unit cube ...
Kx = matrix_Kx(0, 0, 0, 1, 1, 1);
Ky = matrix_Ky(0, 0, 0, 1, 1, 1);
Kz = matrix_Kz(0, 0, 0, 1, 1, 1);
%that it will be scaled with Kge
Kge = (y1-y0).*(z1-z0)./(x1-x0)*Kx + ...
      (x1-x0).*(z1-z0)./(y1-y0)*Ky + ...
      (x1-x0).*(y1-y0)./(z1-z0)*Kz;
%And is then scaled again according to the material thermal properties
Kge = repmat(thermal_cond, 1, 64) .* Kge;

Ige = t(:, kron(1:8, ones(1,8)));
Jge = t(:, kron(ones(1,8), 1:8));

K = sparse(Ige(:), Jge(:), Kge(:), N_nodes, N_nodes);

%Construct the superposition matrix A for unit cubre ...
A_sup = superposition_matrix(0, 0, 0, 1, 1, 1);
%that it will be scaled with A
Age = (x1-x0).*(y1-y0).*(z1-z0)*A_sup;

A = sparse(Ige(:), Jge(:), Age(:), N_nodes, N_nodes);

% inhomogeneity vector
S = vector(0, 0, 0, 1, 1, 1);
for k=1:size(in.sources_position,1)
   ix_condition = (in.sources_position(k,IND_SOURCEPOS_x)<=x(t)) &...
                  (x(t)<=(in.sources_position(k,IND_SOURCEPOS_x)+in.sources_position(k,IND_SOURCEPOS_Dx)));
   iy_condition = (in.sources_position(k,IND_SOURCEPOS_y)<=y(t)) &...
                  (y(t)<=(in.sources_position(k,IND_SOURCEPOS_y)+in.sources_position(k,IND_SOURCEPOS_Dy)));
   iz_condition = (in.sources_position(k,IND_SOURCEPOS_z)<=z(t)) &...
                  (z(t)<=(in.sources_position(k,IND_SOURCEPOS_z)+in.sources_position(k,IND_SOURCEPOS_Dz)));
   is_condition_full = ix_condition & iy_condition & iz_condition;
   is_condition=is_condition_full(:,1) & is_condition_full(:,2) & is_condition_full(:,3) & ...
                is_condition_full(:,4) & is_condition_full(:,5) & is_condition_full(:,6) & ...
                is_condition_full(:,7) & is_condition_full(:,8);
   source_intensity=in.sources_position(k,IND_SOURCEPOS_Val)/...
                    (in.sources_position(k,IND_SOURCEPOS_Dx)*...
                     in.sources_position(k,IND_SOURCEPOS_Dy)*...
                     in.sources_position(k,IND_SOURCEPOS_Dz));
  % Assume uniform source in the volume
  Fge = (is_condition .* volume) * S * source_intensity/...
        (in.media_properties(1,3) * in.media_properties(1,4));

  %Temporary - to be kept for evaluation
  is_source=is_condition;
  dissipated_power=in.sources_position(k,end);
end
F = accumarray([t(:), ones(length(t(:)),1)], Fge(:), [N_nodes 1]);

%Call for the BDF routine
BDF
%
% % Setting the boundary conditions in K and F
% K(plane_z_l,:) = 0;
% %Beware: k is set to 0. Is that the correct way to set the isothermal BC?
% %It sould seem to me that the equations has to be removed, otherwise the
% %system will be worse...
% K(plane_z_l,plane_z_l) = speye(length(plane_z_l));
% F(plane_z_l) = T;
%
% % Linear equation system resolution
% U = K\F;
% % U = agmg(K,F);
%
% %% Display result
%
% % normalization
% U = (U-T)/dissipated_power; % Normalized temperature rise [K/W]
%
% % Volumetric average of the normalized temperature inside the heat source
% % Normalized temperature (NT) rize at the source nodes in Kelvin/Watt
% % Average NT at each parallelepiped inside the heat source
% vol_average_temp= mean(mean(U(t(is_source))')') % [K/W]
%
% toc
%
% elapsed_time_write_vtk=tic();
% vtkwrite('test_binary.vtk','structured_grid',x,y,z,'scalars','NormalizedTemperatureRise',U,'BINARY')
% toc(elapsed_time_write_vtk);
