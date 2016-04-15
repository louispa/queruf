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
sigm_a = 0.1;
sigm_th = 0.1;
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
    Z(1,k) = atan(X(1,k)/X(2,k))+w(1,k);
end
%second part of the algortihm
%end
Abs=X(1,:);

Coo=X(2,:);
% representations
figure(1)
plot(X(1,:),X(2,:),'r')
xlabel('x(k)');
ylabel('y(k)');
title('trajectory of the relatives positions x and y')

Zpi = Z.*(1/pi);
figure(2)
plot(Zpi,'b')
xlabel('time [s]')
ylabel('z(k)/\pi')
title('Evolution of the measured theta angle')
end
