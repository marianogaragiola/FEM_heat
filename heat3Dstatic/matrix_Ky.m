function Ky = matrix_Ky(x0, y0, z0, x1, y1, z1);

Ky = [(-1/9).*(x0+(-1).*x1).*(y0+(-1).*y1).^(-1).*(z0+(-1).*z1),(-1/18) ...
  .*(x0+(-1).*x1).*(y0+(-1).*y1).^(-1).*(z0+(-1).*z1),(1/9).*(x0+( ...
  -1).*x1).*(y0+(-1).*y1).^(-1).*(z0+(-1).*z1),(1/18).*(x0+(-1).*x1) ...
  .*(y0+(-1).*y1).^(-1).*(z0+(-1).*z1),(-1/18).*(x0+(-1).*x1).*(y0+( ...
  -1).*y1).^(-1).*(z0+(-1).*z1),(-1/36).*(x0+(-1).*x1).*(y0+(-1).* ...
  y1).^(-1).*(z0+(-1).*z1),(1/18).*(x0+(-1).*x1).*(y0+(-1).*y1).^( ...
  -1).*(z0+(-1).*z1),(1/36).*(x0+(-1).*x1).*(y0+(-1).*y1).^(-1).*( ...
  z0+(-1).*z1),(-1/18).*(x0+(-1).*x1).*(y0+(-1).*y1).^(-1).*(z0+(-1) ...
  .*z1),(-1/9).*(x0+(-1).*x1).*(y0+(-1).*y1).^(-1).*(z0+(-1).*z1),( ...
  1/18).*(x0+(-1).*x1).*(y0+(-1).*y1).^(-1).*(z0+(-1).*z1),(1/9).*( ...
  x0+(-1).*x1).*(y0+(-1).*y1).^(-1).*(z0+(-1).*z1),(-1/36).*(x0+(-1) ...
  .*x1).*(y0+(-1).*y1).^(-1).*(z0+(-1).*z1),(-1/18).*(x0+(-1).*x1).* ...
  (y0+(-1).*y1).^(-1).*(z0+(-1).*z1),(1/36).*(x0+(-1).*x1).*(y0+(-1) ...
  .*y1).^(-1).*(z0+(-1).*z1),(1/18).*(x0+(-1).*x1).*(y0+(-1).*y1).^( ...
  -1).*(z0+(-1).*z1),(1/9).*(x0+(-1).*x1).*(y0+(-1).*y1).^(-1).*(z0+ ...
  (-1).*z1),(1/18).*(x0+(-1).*x1).*(y0+(-1).*y1).^(-1).*(z0+(-1).* ...
  z1),(-1/9).*(x0+(-1).*x1).*(y0+(-1).*y1).^(-1).*(z0+(-1).*z1),( ...
  -1/18).*(x0+(-1).*x1).*(y0+(-1).*y1).^(-1).*(z0+(-1).*z1),(1/18).* ...
  (x0+(-1).*x1).*(y0+(-1).*y1).^(-1).*(z0+(-1).*z1),(1/36).*(x0+(-1) ...
  .*x1).*(y0+(-1).*y1).^(-1).*(z0+(-1).*z1),(-1/18).*(x0+(-1).*x1).* ...
  (y0+(-1).*y1).^(-1).*(z0+(-1).*z1),(-1/36).*(x0+(-1).*x1).*(y0+( ...
  -1).*y1).^(-1).*(z0+(-1).*z1),(1/18).*(x0+(-1).*x1).*(y0+(-1).*y1) ...
  .^(-1).*(z0+(-1).*z1),(1/9).*(x0+(-1).*x1).*(y0+(-1).*y1).^(-1).*( ...
  z0+(-1).*z1),(-1/18).*(x0+(-1).*x1).*(y0+(-1).*y1).^(-1).*(z0+(-1) ...
  .*z1),(-1/9).*(x0+(-1).*x1).*(y0+(-1).*y1).^(-1).*(z0+(-1).*z1),( ...
  1/36).*(x0+(-1).*x1).*(y0+(-1).*y1).^(-1).*(z0+(-1).*z1),(1/18).*( ...
  x0+(-1).*x1).*(y0+(-1).*y1).^(-1).*(z0+(-1).*z1),(-1/36).*(x0+(-1) ...
  .*x1).*(y0+(-1).*y1).^(-1).*(z0+(-1).*z1),(-1/18).*(x0+(-1).*x1).* ...
  (y0+(-1).*y1).^(-1).*(z0+(-1).*z1),(-1/18).*(x0+(-1).*x1).*(y0+( ...
  -1).*y1).^(-1).*(z0+(-1).*z1),(-1/36).*(x0+(-1).*x1).*(y0+(-1).* ...
  y1).^(-1).*(z0+(-1).*z1),(1/18).*(x0+(-1).*x1).*(y0+(-1).*y1).^( ...
  -1).*(z0+(-1).*z1),(1/36).*(x0+(-1).*x1).*(y0+(-1).*y1).^(-1).*( ...
  z0+(-1).*z1),(-1/9).*(x0+(-1).*x1).*(y0+(-1).*y1).^(-1).*(z0+(-1) ...
  .*z1),(-1/18).*(x0+(-1).*x1).*(y0+(-1).*y1).^(-1).*(z0+(-1).*z1),( ...
  1/9).*(x0+(-1).*x1).*(y0+(-1).*y1).^(-1).*(z0+(-1).*z1),(1/18).*( ...
  x0+(-1).*x1).*(y0+(-1).*y1).^(-1).*(z0+(-1).*z1),(-1/36).*(x0+(-1) ...
  .*x1).*(y0+(-1).*y1).^(-1).*(z0+(-1).*z1),(-1/18).*(x0+(-1).*x1).* ...
  (y0+(-1).*y1).^(-1).*(z0+(-1).*z1),(1/36).*(x0+(-1).*x1).*(y0+(-1) ...
  .*y1).^(-1).*(z0+(-1).*z1),(1/18).*(x0+(-1).*x1).*(y0+(-1).*y1).^( ...
  -1).*(z0+(-1).*z1),(-1/18).*(x0+(-1).*x1).*(y0+(-1).*y1).^(-1).*( ...
  z0+(-1).*z1),(-1/9).*(x0+(-1).*x1).*(y0+(-1).*y1).^(-1).*(z0+(-1) ...
  .*z1),(1/18).*(x0+(-1).*x1).*(y0+(-1).*y1).^(-1).*(z0+(-1).*z1),( ...
  1/9).*(x0+(-1).*x1).*(y0+(-1).*y1).^(-1).*(z0+(-1).*z1),(1/18).*( ...
  x0+(-1).*x1).*(y0+(-1).*y1).^(-1).*(z0+(-1).*z1),(1/36).*(x0+(-1) ...
  .*x1).*(y0+(-1).*y1).^(-1).*(z0+(-1).*z1),(-1/18).*(x0+(-1).*x1).* ...
  (y0+(-1).*y1).^(-1).*(z0+(-1).*z1),(-1/36).*(x0+(-1).*x1).*(y0+( ...
  -1).*y1).^(-1).*(z0+(-1).*z1),(1/9).*(x0+(-1).*x1).*(y0+(-1).*y1) ...
  .^(-1).*(z0+(-1).*z1),(1/18).*(x0+(-1).*x1).*(y0+(-1).*y1).^(-1).* ...
  (z0+(-1).*z1),(-1/9).*(x0+(-1).*x1).*(y0+(-1).*y1).^(-1).*(z0+(-1) ...
  .*z1),(-1/18).*(x0+(-1).*x1).*(y0+(-1).*y1).^(-1).*(z0+(-1).*z1),( ...
  1/36).*(x0+(-1).*x1).*(y0+(-1).*y1).^(-1).*(z0+(-1).*z1),(1/18).*( ...
  x0+(-1).*x1).*(y0+(-1).*y1).^(-1).*(z0+(-1).*z1),(-1/36).*(x0+(-1) ...
  .*x1).*(y0+(-1).*y1).^(-1).*(z0+(-1).*z1),(-1/18).*(x0+(-1).*x1).* ...
  (y0+(-1).*y1).^(-1).*(z0+(-1).*z1),(1/18).*(x0+(-1).*x1).*(y0+(-1) ...
  .*y1).^(-1).*(z0+(-1).*z1),(1/9).*(x0+(-1).*x1).*(y0+(-1).*y1).^( ...
  -1).*(z0+(-1).*z1),(-1/18).*(x0+(-1).*x1).*(y0+(-1).*y1).^(-1).*( ...
  z0+(-1).*z1),(-1/9).*(x0+(-1).*x1).*(y0+(-1).*y1).^(-1).*(z0+(-1) ...
  .*z1)];

end
