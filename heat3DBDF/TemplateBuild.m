%-------------------------------------------------------------------------%
% FEM code for 3D heat equation solution
%-------------------------------------------------------------------------%

clc
clear all
close all
format long

%%%Substrate
W_SUB=500; %[um]
L_SUB=500; %[um]
t_SUB=500; %[um]

x_SUB=0;%[um]
y_SUB=0;%[um]
%Isothermal bottom
z_SUB=0;%[um]

%%%Heat source
W_HS=0.1; %[um], range [0.1,2]
L_HS=2.5; %[um], range [2.5,10]
t_HS=0.12; %[um], range [0.01,0.5]

%Centered along x and y
x_HS_mid=0.5*W_SUB; %[um]
y_HS_mid=0.5*L_SUB; %[um]
%At a given distance from the top
z_HS_fromtop=1; %[um]
z_HS_mid=z_HS_fromtop+0.5*t_HS; %[um]

%Evaluate corners
[x_HS,y_HS,z_HS]=mid2corners(x_HS_mid,...
                             y_HS_mid,...
                             z_HS_mid,...
                             W_HS,...
                             L_HS,...
                             t_HS);

%%%Trench
W_T=1;%[um] range [0.3,2] - along x and y as a "frame"
t_T=5;%[um] range [2.5,5] - along z
%Spacing between trench and HS
S_W=4;%[um] range [0.5,4] - spacing along x, i.e., between long side of HS and trench
S_L=2;%[um] range [0.5,4] - spacing along y, i.e., between short side of HS and trench
%All are at the same bottom height
z_T=0;%[um]

%Please see figure for definition of T1...T4
x_T1=x_HS-S_W-W_T;
y_T1=y_HS-S_L-W_T;
W_T1=W_T;
L_T1=2*W_T+2*S_L+L_HS;

x_T2=x_HS-S_W;
y_T2=y_T1;
W_T2=2*S_W+W_HS;
L_T2=W_T;

x_T3=x_HS+S_W+W_HS;
y_T3=y_T1;
W_T3=W_T;
L_T3=L_T1;

x_T4=x_T2;
y_T4=y_HS+L_HS+S_L;
W_T4=W_T2;
L_T4=W_T;

%Additional positions for gridding
x_bufferSUB_low=200;  %[um]
y_bufferSUB_low=200;  %[um]
x_bufferSUB_high=300; %[um]
y_bufferSUB_high=300; %[um]
z_bufferSUB=100;  %[um]

%Beware of the use of *unique* to avoid point repetition!
factor_xy=1;
x_grid_num1=5; %2*factor_xy;%[adim]
x_grid_num2=5; %16*factor_xy;%[adim]
x_grid_num3=5; %7*factor_xy;%[adim]
x_grid_num4=5; %21*factor_xy;%[adim]
x_grid_num5=5; %10*factor_xy;%[adim]
x_grid_num6=5; %21*factor_xy;%[adim]
x_grid_num7=5; %7*factor_xy;%[adim]
x_grid_num8=5; %16*factor_xy;%[adim]
x_grid_num9=5; %11*factor_xy;%[adim]

x_grid=unique([linspace(x_SUB,x_bufferSUB_low,x_grid_num1), ...
               linspace(x_bufferSUB_low,x_T1,x_grid_num2), ...
               linspace(x_T1,x_T2,x_grid_num3), ...
               linspace(x_T2,x_HS,x_grid_num4), ...
               linspace(x_HS,x_HS+W_HS,x_grid_num5), ...
               linspace(x_HS+W_HS,x_T3,x_grid_num6), ...
               linspace(x_T3,x_T3+W_T,x_grid_num7), ...
               linspace(x_T3+W_T,x_bufferSUB_high,x_grid_num8), ...
               linspace(x_bufferSUB_high,W_SUB,x_grid_num9)]);


y_grid_num1=x_grid_num1;
y_grid_num2=x_grid_num2;
y_grid_num3=x_grid_num3;
y_grid_num4=x_grid_num4;
y_grid_num5=x_grid_num5;
y_grid_num6=x_grid_num6;
y_grid_num7=x_grid_num7;
y_grid_num8=x_grid_num8;
y_grid_num9=x_grid_num9;

y_grid=unique([linspace(y_SUB,y_bufferSUB_low,y_grid_num1) ...
               linspace(y_bufferSUB_low,y_T1,y_grid_num2) ...
               linspace(y_T1,y_T1+W_T,y_grid_num3) ...
               linspace(y_T1+W_T,y_HS,y_grid_num4) ...
               linspace(y_HS,y_HS+L_HS,y_grid_num5) ...
               linspace(y_HS+L_HS,y_T4,y_grid_num6) ...
               linspace(y_T4,y_T4+W_T,y_grid_num7) ...
               linspace(y_T4+W_T,y_bufferSUB_high,y_grid_num8) ...
               linspace(y_bufferSUB_high,L_SUB,y_grid_num9)]);

factor_z=1;
z_grid_num1=5; %10*factor_z;%[adim]
z_grid_num2=5; %10*factor_z;%[adim]
z_grid_num3=5; %30*factor_z;%[adim]
z_grid_num4=5; %15*factor_z;%[adim]
z_grid_num5=5; %15*factor_z;%[adim]
z_grid_num6=5; %15*factor_z;%[adim]
z_grid_num7=5; %10*factor_z;%[adim]

z_grid=unique([linspace(z_SUB,z_HS,z_grid_num1) ...
               linspace(z_HS,z_HS+t_HS,z_grid_num2) ...
               linspace(z_HS+t_HS,t_T,z_grid_num3) ...
               linspace(t_T,2*t_T,z_grid_num4) ...
               linspace(2*t_T,3*t_T,z_grid_num5) ...
               linspace(3*t_T,z_bufferSUB,z_grid_num6) ...
               linspace(z_bufferSUB,t_SUB,z_grid_num7)]);

%Number of planes and nodes
m=length(x_grid); %#planes along x
n=length(y_grid); %#planes along y
l=length(z_grid); %#planes along z
N_nodes=m*n*l;    %total number of nodes
%Build full grid from planes
[x,y,z]=ndgrid(x_grid,y_grid,z_grid);
%Total number of cubes
N_cubes=(m-1)*(n-1)*(l-1);
% List of all n*m*l nodes coordinates x,y,z
p = [x(:), y(:), z(:)];
% t is a matrix with the index of each cube nodes
t = cube_nodes_A(m-1, n-1, l-1);

% Grid visualization options
% r = plot_grid_2d(x, y, z);
% r = plot_grid_3d(x, y, z);

%FANTASTIC INPUT DEFINITION (for compatibility and double checking)
%in.media_properties = [#material | k | rho | c]
k_Si = 148e-6; % [W/umK]
rho_Si = 2.3290e-15; % [kg/um^3]
c_Si = 0.712e3; % [J/(kg*K)]
k_SiO2 = 1.4e-6; % fill range [1e-6; 5e-6]  [W/umK]
rho_SiO2 = 2.634e-15; % [kg/um^3]
c_SiO2 = 0.705e3; % [J/(kg*K)] min value 680 max value 730

in.media_properties = [1 k_Si   rho_Si   c_Si
                       2 k_SiO2 rho_SiO2 c_SiO2];

%in.media_positions  = [#id | x | y | z | Dx | Dy | Dz | #material]
in.media_positions  = [1 x_SUB y_SUB z_SUB W_SUB L_SUB t_SUB 1
                       2 x_T1  y_T1  z_T   W_T1  L_T1  t_T   2
                       3 x_T2  y_T2  z_T   W_T2  L_T2  t_T   2
                       4 x_T3  y_T3  z_T   W_T3  L_T3  t_T   2
                       5 x_T4  y_T4  z_T   W_T4  L_T4  t_T   2];

%in.sources_position = [#id | x | y | z | Dx | Dy | Dz | value]
in.sources_position = [1 x_HS y_HS z_HS W_HS L_HS t_HS 1];

%discrete.axes{1}.P=[column vector of points] along x
%discrete.axes{2}.P=[column vector of points] along y
%discrete.axes{3}.P=[column vector of points] along z
