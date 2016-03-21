T = 0.5;
x0 = [1 1 1 1]';
X = zeros(4,200);
Z = zeros(1,200);
F = eye(4,4);
F(1,3) = T;
F(2,4) = T;
G = [T^2/2 0; 0 T^2/2; T 0; 0 T];
sigm_a = 0.01;
sigm_th = 0.001;
v = randn(2,200)*sigm_a;
w = randn(1,200)*sigm_th;

for k = 1:200
    eps = G*v(:,k);
    if k == 1
        X(:,k) = F*x0 + eps;
    else
        X(:,k+1) = F*X(:,k) + eps;
    end
    Z(1,k) = atan(X(1,k)/X(2,k))+w(1,k);
end

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
