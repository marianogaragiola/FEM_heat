function du = bdf1o(x0, x1, u0, u1);

  du = ((-1).*u0+u1).*((-1).*x0+x1).^(-1);

end
