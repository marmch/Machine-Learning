%Mark Mchedlishvili
%Mini MATLAB 3
%% Gaussian
scatter(bimodal.x(:,1), bimodal.x(:,2), [], bimodal.y);
%scatter(unimodal.x(:,1), unimodal.x(:,2), [], unimodal.y);
x1 = abs(bimodal.x(:,1))-2;
x2 = abs(bimodal.x(:,2)+2)-2;
x = [x1, x2];
t = bimodal.y;
n1 = sum(t);
n2 = length(t)-sum(t);
pi = n1/length(t);

mu1 = 1/n1*[sum(x1.*t),sum(x2.*t)];
mu2 = 1/n2*[sum(x1.*(1-t)),sum(x2.*(1-t))];
s1 = zeros(2,2);
s2 = zeros(2,2);
for i = 1:length(t)
    s1 = s1 + t(i)*(x(i,:)-mu1).'*(x(i,:)-mu1);
    s2 = s2 + (1-t(i))*(x(i,:)-mu2).'*(x(i,:)-mu2);
end
s1 = s1/n1;
s2 = s2/n2;
s = n1/length(t)*s1 + n2/length(t)*s2;