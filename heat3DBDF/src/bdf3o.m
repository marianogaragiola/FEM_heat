function du = bdf3o(x0, x1, x2, x3, u0, u1, u2, u3);

  du = (x0+(-1).*x1).^(-1).*(x0+(-1).*x2).^(-1).*(x1+(-1).*x2).^(-1).*( ...
        x0+(-1).*x3).^(-1).*(x1+(-1).*x3).^(-1).*(x2+(-1).*x3).^(-1).*( ...
        u2.*(x0+(-1).*x1).*(x0+(-1).*x3).^2.*(x1+(-1).*x3).^2+((-1).*u1.*( ...
        x0+(-1).*x2).*(x0+(-1).*x3).^2+u0.*(x1+(-1).*x2).*(x1+(-1).*x3) ...
        .^2).*(x2+(-1).*x3).^2+(-1).*u3.*(x0+(-1).*x1).*(x0+(-1).*x2).*( ...
        x1+(-1).*x2).*(x0.*x1+x0.*x2+x1.*x2+(-2).*(x0+x1+x2).*x3+3.*x3.^2) ...
        );

end