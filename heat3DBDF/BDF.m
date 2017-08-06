%-------------------------------------------------------------------------%
%
% time discretization
%
initial_time = 0.0;
final_time   = 0.01;
%
time_steps = 30;
%
% uniform time steps
% time_vec = linspace(initial_time, final_time, time_steps);
%
% exponential time steps
exp_parameter = 0.08;
%
aux = 1:1:time_steps;
time_vec = exp(exp_parameter*(aux-1)) - 1;
time_vec = (final_time - initial_time)/(time_vec(time_steps) - time_vec(1))*time_vec;
% figure, plot(time_vec, 'o')
%
disp(final_time)
%
q = 5 ; % is the order for de BDF for time aproximation
%
%-------------------------------------------------------------------------%

average_temp = [];
%
% Initial condition for U (isothermal: all grid points at T = 300 K)
C0 = T*ones(N_nodes, 1); % C for Coeficient
%
%
% BDF 1st order
bdf_0 = bdf1o(time_vec(1), time_vec(2), 1, 0);
bdf_1 = bdf1o(time_vec(1), time_vec(2), 0, 1);
%
% Set the left hand side of the equation
LHS = bdf_1*A + K;
% and the boundary conditions
LHS(plane_z_l,:) = 0;
LHS(plane_z_l,plane_z_l) = speye(length(plane_z_l));
%
% Set the right hand side of the equation
RHS = F - bdf_0*A*C0;
%
% Boudary conditions for RHS
RHS(plane_z_l) = T;
%
% Linear equation system resolution
% C1 = agmg(LHS, RHS);
C1 = LHS\RHS;
%
% normalization
U = C1 - T*ones(size(C1)); % Normalized temperature rise in Kelvin/Watt
average_temp = [average_temp; time_vec(2), mean(mean(U(t(is_source))')')]; % [K/W]

%
% BDF 2nd order
bdf_0 = bdf2o(time_vec(1), time_vec(2), time_vec(3), 1, 0, 0);
bdf_1 = bdf2o(time_vec(1), time_vec(2), time_vec(3), 0, 1, 0);
bdf_2 = bdf2o(time_vec(1), time_vec(2), time_vec(3), 0, 0, 1);
% Set the left hand side of the equation
LHS = bdf_2*A + K;
% and the boundary conditions
LHS(plane_z_l,:) = 0;
LHS(plane_z_l,plane_z_l) = speye(length(plane_z_l));
%
% Set the right hand side of the equation
RHS = F - bdf_1*A*C1 - bdf_0*A*C0 ;
%
% Boudary conditions for RHS
RHS(plane_z_l) = T;
%
% Linear equation system resolution
% C2 = agmg(LHS, RHS);
C2 = LHS\RHS;
%
% normalization
U = C2 - T*ones(size(C2)); % Normalized temperature rise in Kelvin/Watt
average_temp = [average_temp; time_vec(3), mean(mean(U(t(is_source))')')]; % [K/W]

%
% BDF 3rd order
bdf_0 = bdf3o(time_vec(1), time_vec(2), time_vec(3), time_vec(4), 1, 0 , 0, 0);
bdf_1 = bdf3o(time_vec(1), time_vec(2), time_vec(3), time_vec(4), 0, 1 , 0, 0);
bdf_2 = bdf3o(time_vec(1), time_vec(2), time_vec(3), time_vec(4), 0, 0 , 1, 0);
bdf_3 = bdf3o(time_vec(1), time_vec(2), time_vec(3), time_vec(4), 0, 0 , 0, 1);
% Set the left hand side of the equation
LHS = bdf_3*A + K;
% and the boundary conditions
LHS(plane_z_l,:) = 0;
LHS(plane_z_l,plane_z_l) = speye(length(plane_z_l));
%
%
% Set the right hand side of the equation
RHS = F - bdf_2*A*C2 - bdf_1*A*C1 - bdf_0*A*C0 ;
%
% Boudary conditions for RHS
RHS(plane_z_l) = T;
%
% Linear equation system resolution
% C3 = agmg(LHS, RHS);
C3 = LHS\RHS;
% normalization
U = C3 - T*ones(size(C3)); % Normalized temperature rise in Kelvin/Watt
average_temp = [average_temp; time_vec(4), mean(mean(U(t(is_source))')')]; % [K/W]

%
% BDF 4th order
bdf_0 = bdf4o(time_vec(1), time_vec(2), time_vec(3), time_vec(4), time_vec(5), 1, 0 , 0, 0, 0);
bdf_1 = bdf4o(time_vec(1), time_vec(2), time_vec(3), time_vec(4), time_vec(5), 0, 1 , 0, 0, 0);
bdf_2 = bdf4o(time_vec(1), time_vec(2), time_vec(3), time_vec(4), time_vec(5), 0, 0 , 1, 0, 0);
bdf_3 = bdf4o(time_vec(1), time_vec(2), time_vec(3), time_vec(4), time_vec(5), 0, 0 , 0, 1, 0);
bdf_4 = bdf4o(time_vec(1), time_vec(2), time_vec(3), time_vec(4), time_vec(5), 0, 0 , 0, 0, 1);
% Set the left hand side of the equation
LHS = bdf_4*A + K;
% and the boundary conditions
LHS(plane_z_l,:) = 0;
LHS(plane_z_l,plane_z_l) = speye(length(plane_z_l));
%
% Set the right hand side of the equation
RHS = F - bdf_3*A*C3 - bdf_2*A*C2 - bdf_1*A*C1 - bdf_0*A*C0 ;
%
% Boudary conditions for RHS
RHS(plane_z_l) = T;
%
% Linear equation system resolution
% C4 = agmg(LHS, RHS);
C4 = LHS\RHS;
%
% normalization
U = C4 - T*ones(size(C4)); % Normalized temperature rise in Kelvin/Watt
average_temp = [average_temp; time_vec(5), mean(mean(U(t(is_source))')')]; % [K/W]
%
% Time loop
%
for time_idx = 6:time_steps
  disp(time_idx)
  toc()
  tn = time_vec(time_idx);
  %
  bdf_0 = bdf5o(time_vec(time_idx-5), time_vec(time_idx-4), time_vec(time_idx-3), time_vec(time_idx-2), time_vec(time_idx-1), time_vec(time_idx), 1, 0 , 0, 0, 0, 0);
  bdf_1 = bdf5o(time_vec(time_idx-5), time_vec(time_idx-4), time_vec(time_idx-3), time_vec(time_idx-2), time_vec(time_idx-1), time_vec(time_idx), 0, 1 , 0, 0, 0, 0);
  bdf_2 = bdf5o(time_vec(time_idx-5), time_vec(time_idx-4), time_vec(time_idx-3), time_vec(time_idx-2), time_vec(time_idx-1), time_vec(time_idx), 0, 0 , 1, 0, 0, 0);
  bdf_3 = bdf5o(time_vec(time_idx-5), time_vec(time_idx-4), time_vec(time_idx-3), time_vec(time_idx-2), time_vec(time_idx-1), time_vec(time_idx), 0, 0 , 0, 1, 0, 0);
  bdf_4 = bdf5o(time_vec(time_idx-5), time_vec(time_idx-4), time_vec(time_idx-3), time_vec(time_idx-2), time_vec(time_idx-1), time_vec(time_idx), 0, 0 , 0, 0, 1, 0);
  bdf_5 = bdf5o(time_vec(time_idx-5), time_vec(time_idx-4), time_vec(time_idx-3), time_vec(time_idx-2), time_vec(time_idx-1), time_vec(time_idx), 0, 0 , 0, 0, 0, 1);
  %
  % Set the left hand side of the equation
  LHS = bdf_5*A + K;
  % and the boundary conditions
  LHS(plane_z_l,:) = 0;
  LHS(plane_z_l,plane_z_l) = speye(length(plane_z_l));
  %
  % Set the right hand side of the equation
  % RHS = (dt*F + A*C);
  RHS = F - bdf_4*A*C4 - bdf_3*A*C3 - bdf_2*A*C2 - bdf_1*A*C1 - bdf_0*A*C0 ;
  %
  % Boudary conditions for RHS
  RHS(plane_z_l) = T;
  %
  % Linear equation system resolution
  % C = agmg(LHS, RHS);
  C = LHS\RHS;
  %
  % normalization
  U = C - T*ones(size(C)); % Normalized temperature rise in Kelvin/Watt
  average_temp = [average_temp; time_vec(time_idx), mean(mean(U(t(is_source))')')]; % [K/W]
  %
  % C matrix update
  C0 = C1;
  C1 = C2;
  C2 = C3;
  C3 = C4;
  C4 = C;
  %
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h = figure('visible', 'on'); hold on;
xlabel('Time [seconds]')
ylabel('Average temperature rise [K/W]')
title('Heat source temperature rise')
plot(average_temp(:,1), average_temp(:,2), '-o')
grid on;
print('average_temp', '-depsc')
print('average_temp', '-dpng')

display(' BDF-FEM resolution finished OK.')
display(' ')
