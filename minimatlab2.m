%Mini MATLAB #1
%Mark Mchedlishvili
%% Part 1

N = 20;
a0 = -0.3;
a1 = 0.5;
mu = 0;
sigma = 0.2;
beta = 25;
alpha = 2;

noise_pd = makedist('Normal', 'mu', mu, 'sigma', sigma);
noise = random(noise_pd, 1, N);

xd = (rand(1,N).*2 - 1).';
yd = a0 + a1*xd + noise.';
xr = -1:0.1:1;
yr = a0 + a1.*xr;
figure;
plot(xr,yr);
title('Actual Plot and Data');
hold on;
scatter(xd, yd);
hold off;

x1 = -1:0.01:1;
x2 = -1:0.01:1;



SN = inv(alpha*eye(2));
mN = zeros(2,1);
X00 = -1:0.01:1;

p = zeros(201);
for i = 1:201
    X01 = (i-101)/100*ones(1,201);
    p(:,i) = mvnpdf([X00;X01].', mN.', SN);
end
figure;
contourf(x1,x2, p);
title('Prior');

r = mvnrnd(mN.', SN, 6);
figure;
title('Data Space (Prior)');
hold on;
for i = 1:6
    plot(x1,r(i,1)+r(i,2)*x1);
end
hold off;

sx1 = size(x1);
sx1 = sx1(2);
sx2 = size(x2);
sx2 = sx2(2);

likelihood = zeros(sx1,sx2);
for i = 1:sx1
    for j = 1:sx2
        distance = abs((-(x1(i).*xd(1)) + yd(1) - x2(j))./sqrt(x1(i)^2+1));
        likelihood(i,j) = 1/(sigma*sqrt(2*pi))*exp(-(distance^2)/(2*sigma^2));
    end
end
figure;
contourf(x1,x2, likelihood);
title('Likelihood (1)');

phi0 = ones(1,1);
phi1 = xd(1);
phi = [phi0 phi1];
SN = inv(alpha*eye(2) + beta*(phi.')*phi);
mN = beta*SN*(phi.')*yd(1); %#ok<MINV>
X00 = -1:0.01:1;

p = zeros(201);
for i = 1:201
    X01 = (i-101)/100*ones(1,201);
    p(:,i) = mvnpdf([X00;X01].', mN.', SN);
end
figure;
contourf(x1,x2, p);
title('Posterior (1)');

r = mvnrnd(mN.', SN, 6);
figure;
title('Data Space (1)');
hold on;
for i = 1:6
    plot(x1,r(i,1)+r(i,2)*x1);
end
hold off;


likelihood = zeros(sx1,sx2);
for k = 1:2
    for i = 1:sx1
        for j = 1:sx2
            distance = abs((-(x1(i).*xd(k)) + yd(k) - x2(j))./sqrt(x1(i)^2+1));
            likelihood(i,j) = 1/(sigma*sqrt(2*pi))*exp(-(distance^2)/(2*sigma^2));
        end
    end
end
figure;
contourf(x1,x2, likelihood);
title('Likelihood (2)');

phi0 = ones(2,1);
phi1 = xd(1:2);
phi = [phi0 phi1];
SN = inv(alpha*eye(2) + beta*(phi.')*phi);
mN = beta*SN*(phi.')*yd(1:2); %#ok<MINV>
X00 = -1:0.01:1;

p = zeros(201);
for i = 1:201
    X01 = (i-101)/100*ones(1,201);
    p(:,i) = mvnpdf([X00;X01].', mN.', SN);
end
figure;
contourf(x1,x2, p);
title('Posterior (2)');

r = mvnrnd(mN.', SN, 6);
figure;
title('Data Space (2)');
hold on;
for i = 1:6
    plot(x1,r(i,1)+r(i,2)*x1);
end
hold off;

likelihood = zeros(sx1,sx2);
for k = 1:N
    for i = 1:sx1
        for j = 1:sx2
            distance = abs((-(x1(i).*xd(k)) + yd(k) - x2(j))./sqrt(x1(i)^2+1));
            likelihood(i,j) = 1/(sigma*sqrt(2*pi))*exp(-(distance^2)/(2*sigma^2));
        end
    end
end
figure;
contourf(x1,x2, likelihood);
title('Likelihood (20)');

phi0 = ones(20,1);
phi1 = xd;
phi = [phi0 phi1];
SN = inv(alpha*eye(2) + beta*(phi.')*phi);
mN = beta*SN*(phi.')*yd; %#ok<MINV>
X00 = -1:0.01:1;

p = zeros(201);
for i = 1:201
    X01 = (i-101)/100*ones(1,201);
    p(:,i) = mvnpdf([X00;X01].', mN.', SN);
end
figure;
contourf(x1,x2, p);
title('Posterior (20)');

r = mvnrnd(mN.', SN, 6);
figure;
title('Data Space (20)');
hold on;
for i = 1:6
    plot(x1,r(i,1)+r(i,2)*x1);
end
hold off;

%% Part 2
mu = 0;
sigma = 0.1;
N = 25;
x = 0:0.01:1;
beta = 5;
alpha = 1;
noise_pd = makedist('Normal', 'mu', mu, 'sigma', sigma);
noise = random(noise_pd, 1, N);

xd = rand(1,N).';
yd = sin(2*pi*xd) + noise.';



figure;
plot(x,sin(2*pi*x));
title('1 Data Point');
hold on;
scatter(xd(1), yd(1));

phi = zeros(1,9);
phi(:,1) = ones(1,1);
for i = 2:9
    phi(:,i) = exp(-((xd(1)-(i-2)/7).^2)/(2*(0.25.^2)));
end

SN = inv(alpha*eye(9) + beta*(phi.')*phi);
mN = beta*SN*(phi.')*yd(1); %#ok<MINV>

stuff = mN(1);
for i = 2:9
    stuff = stuff + mN(i)*exp(-((x-(i-2)/7).^2)/(2*(0.25.^2)));
end

plot(x,stuff);

xphi = zeros(101,9);
xphi(:,1) = ones(101,1);
for i = 2:9
    xphi(:,i) = exp(-((x-(i-2)/7).^2)/(2*(0.25.^2)));
end

var = zeros(1,101);
for i = 1:101
    var(i) = 1/beta + xphi(i,:)*SN*(xphi(i,:).');
end
xx = [x, fliplr(x)];
yy = [stuff+var, fliplr(stuff-var)];
fill(xx, yy, 'r', 'FaceAlpha', 0.1);
hold off;




figure;
plot(x,sin(2*pi*x));
title('2 Data Points');
hold on;
scatter(xd(1:2), yd(1:2));

phi = zeros(2,9);
phi(:,1) = ones(2,1);
for i = 2:9
    phi(:,i) = exp(-((xd(1:2)-(i-2)/7).^2)/(2*(0.25.^2)));
end

SN = inv(alpha*eye(9) + beta*(phi.')*phi);
mN = beta*SN*(phi.')*yd(1:2); %#ok<MINV>

stuff = mN(1);
for i = 2:9
    stuff = stuff + mN(i)*exp(-((x-(i-2)/7).^2)/(2*(0.25.^2)));
end

plot(x,stuff);

xphi = zeros(101,9);
xphi(:,1) = ones(101,1);
for i = 2:9
    xphi(:,i) = exp(-((x-(i-2)/7).^2)/(2*(0.25.^2)));
end

var = zeros(1,101);
for i = 1:101
    var(i) = 1/beta + xphi(i,:)*SN*(xphi(i,:).');
end

xx = [x, fliplr(x)];
yy = [stuff+var, fliplr(stuff-var)];
fill(xx, yy, 'r', 'FaceAlpha', 0.1);




figure;
plot(x,sin(2*pi*x));
title('4 Data Points');
hold on;
scatter(xd(1:4), yd(1:4));

phi = zeros(4,9);
phi(:,1) = ones(4,1);
for i = 2:9
    phi(:,i) = exp(-((xd(1:4)-(i-2)/7).^2)/(2*(0.25.^2)));
end

SN = inv(alpha*eye(9) + beta*(phi.')*phi);
mN = beta*SN*(phi.')*yd(1:4); %#ok<MINV>

stuff = mN(1);
for i = 2:9
    stuff = stuff + mN(i)*exp(-((x-(i-2)/7).^2)/(2*(0.25.^2)));
end

plot(x,stuff);

xphi = zeros(101,9);
xphi(:,1) = ones(101,1);
for i = 2:9
    xphi(:,i) = exp(-((x-(i-2)/7).^2)/(2*(0.25.^2)));
end

var = zeros(1,101);
for i = 1:101
    var(i) = 1/beta + xphi(i,:)*SN*(xphi(i,:).');
end

xx = [x, fliplr(x)];
yy = [stuff+var, fliplr(stuff-var)];
fill(xx, yy, 'r', 'FaceAlpha', 0.1);




figure;
plot(x,sin(2*pi*x));
title('20 Data Points');
hold on;
scatter(xd, yd);

phi = zeros(N,9);
phi(:,1) = ones(N,1);
for i = 2:9
    phi(:,i) = exp(-((xd-(i-2)/7).^2)/(2*(0.25.^2)));
end

SN = inv(alpha*eye(9) + beta*(phi.')*phi);
mN = beta*SN*(phi.')*yd; %#ok<MINV>

stuff = mN(1);
for i = 2:9
    stuff = stuff + mN(i)*exp(-((x-(i-2)/7).^2)/(2*(0.25.^2)));
end

plot(x,stuff);

xphi = zeros(101,9);
xphi(:,1) = ones(101,1);
for i = 2:9
    xphi(:,i) = exp(-((x-(i-2)/7).^2)/(2*(0.25.^2)));
end

var = zeros(1,101);
for i = 1:101
    var(i) = 1/beta + xphi(i,:)*SN*(xphi(i,:).');
end

xx = [x, fliplr(x)];
yy = [stuff+var, fliplr(stuff-var)];
fill(xx, yy, 'r', 'FaceAlpha', 0.1);
