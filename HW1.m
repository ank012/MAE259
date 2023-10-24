syms z h
% Coefficients for AB3 method
beta = [23/12, -16/12, 5/12];

poly = z^3 - z^2 - z*h*(beta(1)*z^2 + beta(2)*z^1 + beta(3));
roots_AB3 = solve(poly,z);

v=sin(z);
sympref('PolynomialDisplayStyle','ascend');
taylor(poly,z,'Order',12)