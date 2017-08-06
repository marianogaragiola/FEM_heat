function r = plot_grid_2d(x, y, z)
close all

F = zeros(size(x));

%%% Grid in xy planes
figure
surf(x(:,:,1), y(:,:,1), squeeze(F(:,:,1)));
view(2);
daspect([1 1 1])
xlim([0 500]);
ylim([0 500]);
xlabel('x axis [micrometers]');
ylabel('y axis [micrometers]');
title('Grid in xy planes');

%%% Grid in xz planes
figure
surf(squeeze(x(:,1,:)), squeeze(z(:,1,:)), squeeze(F(:,1,:)));
view(2);
daspect([1 1 1])
xlim([0 500]);
ylim([0 500]);
xlabel('x axis [micrometers]');
ylabel('z axis [micrometers]');
title('Grid in xz planes');

%%% Grid in yz planes
figure
surf(squeeze(y(1,:,:)), squeeze(z(1,:,:)), squeeze(F(1,:,:)));
view(2);
daspect([1 1 1])
xlim([0 500]);
ylim([0 500]);
xlabel('y axis [micrometers]');
ylabel('z axis [micrometers]');
title('Grid in yz planes');

r = 0;

end