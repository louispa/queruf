function[Z,Abs,Coo]=Q2()
% time period
T = 0.5;
% initial vector
x0 = [1 1 1 1]';
% preallocating
X = zeros(4,200);
Z = zeros(1,200);
% construction of matrix F
F = eye(4,4);
F(1,3) = T;
F(2,4) = T;
% G
G = [T^2/2 0; 0 T^2/2; T 0; 0 T];
% sigma_a and sigma_theta
sigm_a = 0.01;
sigm_th = 10;
% v random gaussian noise of variance sigma_a*I and zero-mean (iid) 
% w random gaussian noise of variance sigma_theta and zero-mean
%for m=1:5000
v = randn(2,200)*sigm_a;
w = randn(1,200)*sigm_th;
for k = 1:200
    eps = G*v(:,k);
    if k == 1
        X(:,k) = F*x0 + eps;
    else
        X(:,k) = F*X(:,k-1) + eps;
    end
    
    if X(1,k)>=0 && X(2,k)>=0
        Z(1,k) = atan(abs(X(1,k)/X(2,k)))+w(1,k);
    elseif X(1,k)>0 && X(2,k)<0
        Z(1,k) = pi - atan(abs(X(1,k)/X(2,k)))+w(1,k);
    elseif X(1,k)<0 && X(2,k)<0
        Z(1,k) = pi + atan(abs(X(1,k)/X(2,k)))+w(1,k);
    else
        Z(1,k) = 2*pi - atan(abs(X(1,k)/X(2,k)))+w(1,k);
    end
end

% representations
subplot(1,2,2)
plot(X(1,:),X(2,:),'r')
xlabel('x(k)');
ylabel('y(k)');
title('trajectory of the relatives positions x and y')

Zpi = Z.*(1/pi);
subplot(1,2,1)
plot(Zpi,'b')
xlabel('time [s]')
ylabel('z(k)/\pi')
title('Evolution of the measured theta angle')
end
