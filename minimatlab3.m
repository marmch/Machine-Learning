%Mark Mchedlishvili
%Mini MATLAB 3

%% Gaussian Generative

% BIMODAL
figure('Position', [0,0, 512, 1024]);
subplot(4,2,1);
scatter(bimodal.x(:,1), bimodal.x(:,2), [], bimodal.y);
axis([-8,6,-8,6]);

xm = bimodal.x(:,1) + bimodal.x(:,2) + 4;
x2 = abs(abs(xm)-4);

t = bimodal.y;
n1 = sum(t);
n2 = length(t)-sum(t);
mu1 = 1/n1*sum(t.*x2);
mu2 = 1/n2*sum((1-t).*x2);
s1 = 1/n1*sum(t.*(x2-mu1).*(x2-mu1));
s2 = 1/n2*sum((1-t).*(x2-mu2).*(x2-mu2));
s = n1/length(t)*s1 + n2/length(t)*s2;
w = (1/s)*mu1;
w0 = -1/2*mu1*(1/s)*mu1;

sigmat = zeros(141,141);
[xr1,xr2] = meshgrid(-8:0.1:6);
for i = 1:141
    for j = 1:141
        xm = xr1(1,i) + xr2(j,1) + 4;
        xx = abs(abs(xm)-4);
        sigmat(i,j) = sigmf(xx*w+w0,[1,0]);
    end
end

subplot(4,2,2);
%scatter(x1, x2, [], t);
contourf(xr2,xr1,sigmat);
axis([-8,6,-8,6]);

% CIRCLES
subplot(4,2,3);
scatter(circles.x(:,1), circles.x(:,2), [], circles.y);
axis([-0.5,2.5,-0.5,2.5]);
x2 = (circles.x(:,1)-1).*(circles.x(:,1)-1)+(circles.x(:,2)-1).*(circles.x(:,2)-1);
t = circles.y;
n1 = sum(t);
n2 = length(t)-sum(t);
mu1 = 1/n1*sum(t.*x2);
mu2 = 1/n2*sum((1-t).*x2);
s1 = 1/n1*sum(t.*(x2-mu1).*(x2-mu1));
s2 = 1/n2*sum((1-t).*(x2-mu2).*(x2-mu2));
s = n1/length(t)*s1 + n2/length(t)*s2;
w = (1/s)*mu1;
w0 = -1/2*mu1*(1/s)*mu1;
sigmat = zeros(31,31);
[xr1,xr2] = meshgrid(-0.5:0.1:2.5);
for i = 1:31
    for j = 1:31
        xx = (xr1(1,i)-1)*(xr1(1,i)-1)+(xr2(j,1)-1)*(xr2(j,1)-1);
        sigmat(i,j) = sigmf(w*xx + w0,[1,0]);
    end
end

subplot(4,2,4);
%scatter(x1, x2, [], bimodal.y);
contourf(xr2,xr1,sigmat);
axis([-0.5,2.5,-0.5,2.5]);

% SPIRAL
subplot(4,2,5);
scatter(spiral.x(:,1), spiral.x(:,2), [], spiral.y);
axis([-6,8,-6,8]);
x1 = (spiral.x(:,1).*spiral.x(:,1))+(spiral.x(:,2).*spiral.x(:,2));
x2 = atan(spiral.x(:,2)./spiral.x(:,1));
for i = 1:400
    if -pi/15*x1(i)+pi/2+pi/15*20 < x2(i)
        x2(i) = x2(i) - 4*pi;
    elseif -pi/10*x1(i) + pi/2 + pi/10*8 < x2(i)
        x2(i) = x2(i) - 3*pi;
    elseif -pi/5*x1(i) + pi/2 + pi/5*1 < x2(i)
        x2(i) = x2(i) - 2*pi;
    elseif -2*x1(i) + 0.5 < x2(i)
        x2(i) = x2(i) - pi;
    end
    x2(i) = x2(i)+2.55*(x1(i).^0.45);
end
t = spiral.y;
n1 = sum(t);
n2 = length(t)-sum(t);
mu1 = 1/n1*sum(t.*x2);
mu2 = 1/n2*sum((1-t).*x2);
s1 = 1/n1*sum(t.*(x2-mu1).*(x2-mu1));
s2 = 1/n2*sum((1-t).*(x2-mu2).*(x2-mu2));
s = n1/length(t)*s1 + n2/length(t)*s2;
w = (1/s)*mu1;
w0 = -1/2*mu1*(1/s)*mu1;
sigmat = zeros(141,141);
[xr1,xr2] = meshgrid(-6:0.1:8);
for i = 1:141
    for j = 1:141
        xx = (xr1(1,i)*xr1(1,i))+(xr2(j,1)*xr2(j,1));
        yy = atan(xr2(j,1)/xr1(1,i));
        if -pi/15*xx+pi/2+pi/15*20 < yy
            yy = yy - 4*pi;
        elseif -pi/10*xx + pi/2 + pi/10*8 < yy
            yy = yy - 3*pi;
        elseif -pi/5*xx + pi/2 + pi/5*1 < yy
            yy = yy - 2*pi;
        elseif -2*xx + 0.5 < yy
            yy = yy - pi;
        end
        yy = yy+2.55*(xx.^0.45);
        sigmat(i,j) = sigmf(w*yy + w0,[1,0]);
    end
end

subplot(4,2,6);
%scatter(x1, x2, [], bimodal.y);
contourf(xr2,xr1,sigmat);
axis([-6,8,-6,8]);

% UNIMODAL
subplot(4,2,7);
scatter(unimodal.x(:,1), unimodal.x(:,2), [], unimodal.y);
axis([-6,8,-6,8]);
x1 = unimodal.x(:,1);
x2 = unimodal.x(:,2);

t = unimodal.y;
n1 = sum(t);
n2 = length(t)-sum(t);
mu1 = 1/n1*[sum(t.*x1),sum(t.*x2)];
mu2 = 1/n2*[sum((1-t).*x1),sum((1-t).*x2)];
s1 = 1/n1*(repmat(t,1,2).*([x1,x2]-repmat(mu1,400,1))).'*([x1,x2]-repmat(mu1,400,1));
s2 = 1/n2*(repmat(1-t,1,2).*([x1,x2]-repmat(mu2,400,1))).'*([x1,x2]-repmat(mu2,400,1));
s = n1/length(t)*s1 + n2/length(t)*s2;
w = inv(s)*(mu1-mu2).';
w0 = -1/2*mu1*inv(s)*mu1.' + 1/2*mu2*inv(s)*mu2.';

sigmat = zeros(141,141);
[xr1,xr2] = meshgrid(-6:0.1:8);
for i = 1:141
    for j = 1:141
        xx = (xr1(1,i)*xr1(1,i))+(xr2(j,1)*xr2(j,1));
        yy = atan(xr2(j,1)/xr1(1,i));
        yy = yy+2.55*(xx.^0.45);
        sigmat(i,j) = sigmf(w.'*[xr1(1,i);xr2(j,1)] + w0,[1,0]);
    end
end

subplot(4,2,8);
%scatter(x1, x2, [], bimodal.y);
contourf(xr2,xr1,sigmat);
axis([-6,8,-6,8]);

%% Logistic Regression

% BIMODAL
figure('Position', [0,0, 512, 1024]);
subplot(4,2,1);
scatter(bimodal.x(:,1), bimodal.x(:,2), [], bimodal.y);
axis([-8,6,-8,6]);
xm = bimodal.x(:,1) + bimodal.x(:,2) + 4;
x2 = abs(abs(xm)-4);
t = bimodal.y;
phi = [ones(length(t),1),x2];

w = inv(phi.' * phi) * (phi.' * t);
y = sigmf(w.'*phi.',[1,0]).';
for i = 1:100000
    y_ = y;
    R = diag(y.*(1-y));
    z = phi*w - inv(R)*(y-t);
    w = inv(phi.'*R*phi)*phi.'*R*z;
    y = sigmf(w.'*phi.',[1,0]).';
    if sum(abs(y_-y)/length(y)) < 0.0001
        break
    end
end
sigmat = zeros(141,141);
[xr1,xr2] = meshgrid(-8:0.1:6);
for i = 1:141
    for j = 1:141
        xxm = xr1(1,i) + xr2(j,1) + 4;
        xx = abs(abs(xxm)-4);
        sigmat(i,j) = sigmf(w.'*[1;xx],[1,0]);
    end
end

subplot(4,2,2);
%scatter(x1, x2, [], bimodal.y);
contourf(xr2,xr1,sigmat);
axis([-8,6,-8,6]);

% CIRCLES
subplot(4,2,3);
scatter(circles.x(:,1), circles.x(:,2), [], circles.y);
axis([-0.5,2.5,-0.5,2.5]);
x2 = (circles.x(:,1)-1).*(circles.x(:,1)-1)+(circles.x(:,2)-1).*(circles.x(:,2)-1);
t = circles.y;
phi = [ones(length(t),1),x2];

w = inv(phi.' * phi) * (phi.' * t);
y = sigmf(w.'*phi.',[1,0]).';
for i = 1:100000
    y_ = y;
    R = diag(y.*(1-y));
    z = phi*w - inv(R)*(y-t);
    w = inv(phi.'*R*phi)*phi.'*R*z;
    y = sigmf(w.'*phi.',[1,0]).';
    if sum(abs(y_-y)/length(y)) < 0.0001
        break
    end
end

sigmat = zeros(31,31);
[xr1,xr2] = meshgrid(-0.5:0.1:2.5);
for i = 1:31
    for j = 1:31
        xx = (xr1(1,i)-1)*(xr1(1,i)-1)+(xr2(j,1)-1)*(xr2(j,1)-1);
        sigmat(i,j) = sigmf(w.'*[1;xx],[1,0]);
    end
end

subplot(4,2,4);
%scatter(x1, x2, [], bimodal.y);
contourf(xr2,xr1,sigmat);
axis([-0.5,2.5,-0.5,2.5]);


% SPIRAL 
subplot(4,2,5);
scatter(spiral.x(:,1), spiral.x(:,2), [], spiral.y);
axis([-6,8,-6,8]);
x1 = (spiral.x(:,1).*spiral.x(:,1))+(spiral.x(:,2).*spiral.x(:,2));
x2 = atan(spiral.x(:,2)./spiral.x(:,1));
for i = 1:400
    if -pi/15*x1(i)+pi/2+pi/15*20 < x2(i)
        x2(i) = x2(i) - 4*pi;
    elseif -pi/10*x1(i) + pi/2 + pi/10*8 < x2(i)
        x2(i) = x2(i) - 3*pi;
    elseif -pi/5*x1(i) + pi/2 + pi/5*1 < x2(i)
        x2(i) = x2(i) - 2*pi;
    elseif -2*x1(i) + 0.5 < x2(i)
        x2(i) = x2(i) - pi;
    end
    x2(i) = x2(i)+2.55*(x1(i).^0.45);
end
t = spiral.y;
phi = [ones(length(t),1),x2];

w = inv(phi.' * phi) * (phi.' * t);
y = sigmf(w.'*phi.',[1,0]).';
for i = 1:100000
    y_ = y;
    R = diag(y.*(1-y));
    z = phi*w - inv(R)*(y-t);
    w = inv(phi.'*R*phi)*phi.'*R*z;
    y = sigmf(w.'*phi.',[1,0]).';
    if sum(abs(y_-y)/length(y)) < 0.0001
        break
    end
end

sigmat = zeros(141,141);
[xr1,xr2] = meshgrid(-6:0.1:8);
for i = 1:141
    for j = 1:141
        xx = (xr1(1,i)*xr1(1,i))+(xr2(j,1)*xr2(j,1));
        yy = atan(xr2(j,1)/xr1(1,i));
        if -pi/15*xx+pi/2+pi/15*20 < yy
            yy = yy - 4*pi;
        elseif -pi/10*xx + pi/2 + pi/10*8 < yy
            yy = yy - 3*pi;
        elseif -pi/5*xx + pi/2 + pi/5*1 < yy
            yy = yy - 2*pi;
        elseif -2*xx + 0.5 < yy
            yy = yy - pi;
        end
        yy = yy+2.55*(xx.^0.45);
        sigmat(i,j) = sigmf(w.'*[1;yy],[1,0]);
    end
end

subplot(4,2,6);
%scatter(x1, x2, [], bimodal.y);
contourf(xr2,xr1,sigmat);
axis([-6,8,-6,8]);


% UNIMODAL
subplot(4,2,7);
scatter(unimodal.x(:,1), unimodal.x(:,2), [], unimodal.y);
axis([-4,5,-4,5]);
x1 = unimodal.x(:,1);
x2 = unimodal.x(:,2);
t = unimodal.y;
phi = [ones(length(t),1),x1,x2];

w = inv(phi.' * phi) * (phi.' * t);
y = sigmf(w.'*phi.',[1,0]).';
for i = 1:100000
    y_ = y;
    R = diag(y.*(1-y));
    z = phi*w - inv(R)*(y-t);
    w = inv(phi.'*R*phi)*phi.'*R*z;
    y = sigmf(w.'*phi.',[1,0]).';
    if sum(abs(y_-y)/length(y)) < 0.0001
        break
    end
end

sigmat = zeros(91,91);
[xr1,xr2] = meshgrid(-4:0.1:5);
for i = 1:91
    for j = 1:91
        xx = xr1(1,i);
        yy = xr2(j,1);
        sigmat(i,j) = sigmf(w.'*[1;xx;yy],[1,0]);
    end
end

subplot(4,2,8);
%scatter(x1, x2, [], bimodal.y);
contourf(xr2,xr1,sigmat);
axis([-4,5,-4,5]);