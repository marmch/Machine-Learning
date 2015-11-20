%Mini MATLAB #1
%Mark Mchedlishvili

mu = 0;
sigma = 0.1;
N = 10;
x = 0:0.01:1;
noise_pd = makedist('Normal', 'mu', mu, 'sigma', sigma);
noise = random(noise_pd, 1, N);

xd = rand(1,N).';
yd = sin(2*pi*xd) + noise.';

figure;
plot(x,sin(2*pi*x));
title('10 Data Points');
hold on;
scatter(xd, yd);

C = zeros(N,N);
beta = 10;
delta = 1;

for i = 1:N
    for j = 1:N
        C(i,j) = kernel(xd(i),xd(j)) + 1/beta*delta;
    end
end

mux = zeros(1,101);
varx = zeros(1,101);
for i=1:101
    kx = zeros(N,1);
    for j = 1:N
        kx(j) = kernel(x(i),xd(j));
    end
    c = kernel(x(i),x(i)) + 1/beta;
    
    mux(i) = kx.' * inv(C) * yd;
    varx(i) = c - kx.' * inv(C) * kx;
end
%figure;
plot(x,mux);
xx = [x, fliplr(x)];
yy = [mux+varx, fliplr(mux-varx)];
fill(xx, yy, 'r', 'FaceAlpha', 0.1);
hold off;