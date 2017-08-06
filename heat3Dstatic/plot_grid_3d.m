function r = plot_grid_3d(x, y, z)

figure
X = x(:); Y = y(:); Z = z(:);
plot3(X, Y, Z, '.')
xlabel('x axis [micrometers]')
ylabel('y axis [micrometers]')
zlabel('z axis [micrometers]')
title('3D Grid')

r = 0;

end