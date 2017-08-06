function Kz = matrix_Kz(x0, y0, z0, x1, y1, z1);

Kz = [(-1/9).*(x0+(-1).*x1).*(y0+(-1).*y1).*(z0+(-1).*z1).^(-1),(-1/18) ...
  .*(x0+(-1).*x1).*(y0+(-1).*y1).*(z0+(-1).*z1).^(-1),(-1/18).*(x0+( ...
  -1).*x1).*(y0+(-1).*y1).*(z0+(-1).*z1).^(-1),(-1/36).*(x0+(-1).* ...
  x1).*(y0+(-1).*y1).*(z0+(-1).*z1).^(-1),(1/9).*(x0+(-1).*x1).*(y0+ ...
  (-1).*y1).*(z0+(-1).*z1).^(-1),(1/18).*(x0+(-1).*x1).*(y0+(-1).* ...
  y1).*(z0+(-1).*z1).^(-1),(1/18).*(x0+(-1).*x1).*(y0+(-1).*y1).*( ...
  z0+(-1).*z1).^(-1),(1/36).*(x0+(-1).*x1).*(y0+(-1).*y1).*(z0+(-1) ...
  .*z1).^(-1),(-1/18).*(x0+(-1).*x1).*(y0+(-1).*y1).*(z0+(-1).*z1) ...
  .^(-1),(-1/9).*(x0+(-1).*x1).*(y0+(-1).*y1).*(z0+(-1).*z1).^(-1),( ...
  -1/36).*(x0+(-1).*x1).*(y0+(-1).*y1).*(z0+(-1).*z1).^(-1),(-1/18) ...
  .*(x0+(-1).*x1).*(y0+(-1).*y1).*(z0+(-1).*z1).^(-1),(1/18).*(x0+( ...
  -1).*x1).*(y0+(-1).*y1).*(z0+(-1).*z1).^(-1),(1/9).*(x0+(-1).*x1) ...
  .*(y0+(-1).*y1).*(z0+(-1).*z1).^(-1),(1/36).*(x0+(-1).*x1).*(y0+( ...
  -1).*y1).*(z0+(-1).*z1).^(-1),(1/18).*(x0+(-1).*x1).*(y0+(-1).*y1) ...
  .*(z0+(-1).*z1).^(-1),(-1/18).*(x0+(-1).*x1).*(y0+(-1).*y1).*(z0+( ...
  -1).*z1).^(-1),(-1/36).*(x0+(-1).*x1).*(y0+(-1).*y1).*(z0+(-1).* ...
  z1).^(-1),(-1/9).*(x0+(-1).*x1).*(y0+(-1).*y1).*(z0+(-1).*z1).^( ...
  -1),(-1/18).*(x0+(-1).*x1).*(y0+(-1).*y1).*(z0+(-1).*z1).^(-1),( ...
  1/18).*(x0+(-1).*x1).*(y0+(-1).*y1).*(z0+(-1).*z1).^(-1),(1/36).*( ...
  x0+(-1).*x1).*(y0+(-1).*y1).*(z0+(-1).*z1).^(-1),(1/9).*(x0+(-1).* ...
  x1).*(y0+(-1).*y1).*(z0+(-1).*z1).^(-1),(1/18).*(x0+(-1).*x1).*( ...
  y0+(-1).*y1).*(z0+(-1).*z1).^(-1),(-1/36).*(x0+(-1).*x1).*(y0+(-1) ...
  .*y1).*(z0+(-1).*z1).^(-1),(-1/18).*(x0+(-1).*x1).*(y0+(-1).*y1).* ...
  (z0+(-1).*z1).^(-1),(-1/18).*(x0+(-1).*x1).*(y0+(-1).*y1).*(z0+( ...
  -1).*z1).^(-1),(-1/9).*(x0+(-1).*x1).*(y0+(-1).*y1).*(z0+(-1).*z1) ...
  .^(-1),(1/36).*(x0+(-1).*x1).*(y0+(-1).*y1).*(z0+(-1).*z1).^(-1),( ...
  1/18).*(x0+(-1).*x1).*(y0+(-1).*y1).*(z0+(-1).*z1).^(-1),(1/18).*( ...
  x0+(-1).*x1).*(y0+(-1).*y1).*(z0+(-1).*z1).^(-1),(1/9).*(x0+(-1).* ...
  x1).*(y0+(-1).*y1).*(z0+(-1).*z1).^(-1),(1/9).*(x0+(-1).*x1).*(y0+ ...
  (-1).*y1).*(z0+(-1).*z1).^(-1),(1/18).*(x0+(-1).*x1).*(y0+(-1).* ...
  y1).*(z0+(-1).*z1).^(-1),(1/18).*(x0+(-1).*x1).*(y0+(-1).*y1).*( ...
  z0+(-1).*z1).^(-1),(1/36).*(x0+(-1).*x1).*(y0+(-1).*y1).*(z0+(-1) ...
  .*z1).^(-1),(-1/9).*(x0+(-1).*x1).*(y0+(-1).*y1).*(z0+(-1).*z1).^( ...
  -1),(-1/18).*(x0+(-1).*x1).*(y0+(-1).*y1).*(z0+(-1).*z1).^(-1),( ...
  -1/18).*(x0+(-1).*x1).*(y0+(-1).*y1).*(z0+(-1).*z1).^(-1),(-1/36) ...
  .*(x0+(-1).*x1).*(y0+(-1).*y1).*(z0+(-1).*z1).^(-1),(1/18).*(x0+( ...
  -1).*x1).*(y0+(-1).*y1).*(z0+(-1).*z1).^(-1),(1/9).*(x0+(-1).*x1) ...
  .*(y0+(-1).*y1).*(z0+(-1).*z1).^(-1),(1/36).*(x0+(-1).*x1).*(y0+( ...
  -1).*y1).*(z0+(-1).*z1).^(-1),(1/18).*(x0+(-1).*x1).*(y0+(-1).*y1) ...
  .*(z0+(-1).*z1).^(-1),(-1/18).*(x0+(-1).*x1).*(y0+(-1).*y1).*(z0+( ...
  -1).*z1).^(-1),(-1/9).*(x0+(-1).*x1).*(y0+(-1).*y1).*(z0+(-1).*z1) ...
  .^(-1),(-1/36).*(x0+(-1).*x1).*(y0+(-1).*y1).*(z0+(-1).*z1).^(-1), ...
  (-1/18).*(x0+(-1).*x1).*(y0+(-1).*y1).*(z0+(-1).*z1).^(-1),(1/18) ...
  .*(x0+(-1).*x1).*(y0+(-1).*y1).*(z0+(-1).*z1).^(-1),(1/36).*(x0+( ...
  -1).*x1).*(y0+(-1).*y1).*(z0+(-1).*z1).^(-1),(1/9).*(x0+(-1).*x1) ...
  .*(y0+(-1).*y1).*(z0+(-1).*z1).^(-1),(1/18).*(x0+(-1).*x1).*(y0+( ...
  -1).*y1).*(z0+(-1).*z1).^(-1),(-1/18).*(x0+(-1).*x1).*(y0+(-1).* ...
  y1).*(z0+(-1).*z1).^(-1),(-1/36).*(x0+(-1).*x1).*(y0+(-1).*y1).*( ...
  z0+(-1).*z1).^(-1),(-1/9).*(x0+(-1).*x1).*(y0+(-1).*y1).*(z0+(-1) ...
  .*z1).^(-1),(-1/18).*(x0+(-1).*x1).*(y0+(-1).*y1).*(z0+(-1).*z1) ...
  .^(-1),(1/36).*(x0+(-1).*x1).*(y0+(-1).*y1).*(z0+(-1).*z1).^(-1),( ...
  1/18).*(x0+(-1).*x1).*(y0+(-1).*y1).*(z0+(-1).*z1).^(-1),(1/18).*( ...
  x0+(-1).*x1).*(y0+(-1).*y1).*(z0+(-1).*z1).^(-1),(1/9).*(x0+(-1).* ...
  x1).*(y0+(-1).*y1).*(z0+(-1).*z1).^(-1),(-1/36).*(x0+(-1).*x1).*( ...
  y0+(-1).*y1).*(z0+(-1).*z1).^(-1),(-1/18).*(x0+(-1).*x1).*(y0+(-1) ...
  .*y1).*(z0+(-1).*z1).^(-1),(-1/18).*(x0+(-1).*x1).*(y0+(-1).*y1).* ...
  (z0+(-1).*z1).^(-1),(-1/9).*(x0+(-1).*x1).*(y0+(-1).*y1).*(z0+(-1) ...
  .*z1).^(-1)];
  
end
