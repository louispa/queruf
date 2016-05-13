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
% Gamma
G = [T^2/2 0; 0 T^2/2; T 0; 0 T];
% sigma_a and sigma_theta
% choose various values for sigma_a and sigma_theta
sigm_a = 0.01;
sigm_th = 10;
% v random gaussian noise of variance sigma^2_a*I and zero-mean (iid) 
% w random gaussian noise of variance sigma^2_theta and zero-mean
v = randn(2,200)*sigm_a;
w = randn(1,200)*sigm_th;
for k = 1:200
    eps = G*v(:,k);
    
    % We suppose U to be zero for all k
    % state equation: x_k+1 = F*x_k + epsilon
    if k == 1
        X(:,k) = F*x0 + eps;
    else
        X(:,k) = F*X(:,k-1) + eps;
    end
    
    % state equation: theta_k = h(x_k) + w_k
    if X(1,k)>=0 && X(2,k)>=0
        Z(1,k) = atan(abs(X(1,k)/X(2,k)))+w(k);
    elseif X(1,k)>0 && X(2,k)<0
        Z(1,k) = pi - atan(abs(X(1,k)/X(2,k)))+w(k);
    elseif X(1,k)<0 && X(2,k)<0
        Z(1,k) = pi + atan(abs(X(1,k)/X(2,k)))+w(k);
    else
        Z(1,k) = 2*pi - atan(abs(X(1,k)/X(2,k)))+w(k);
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
