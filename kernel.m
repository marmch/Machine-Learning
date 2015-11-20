% Gaussian radial basis
function kernel = kernel(xn,xm)
a = 100;
kernel = exp(-a*(abs(xn-xm).^2));
end