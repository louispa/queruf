function [ bound ] = Q6fun(observer,target,bearingMeasurements,mu_r,mu_theta,mu_s,mu_c,X)


t_f=26;
sigma_s = sqrt(0.1);
sigma_theta = sqrt(10^-4);
sigma_r = sqrt(0.1);
sigma_c =sqrt(0.1);

Pxx=mu_r^2*sigma_theta^2*cos(mu_theta)^2+sigma_r^2*sin(mu_theta)^2;
Pyy=mu_r^2*sigma_theta^2*sin(mu_theta)^2+sigma_r^2*cos(mu_theta)^2;
Pxy=sigma_r^2-mu_r^2*sigma_theta^2;

Ppxx=mu_s^2*sigma_c^2*cos(mu_c)^2+sigma_s^2*sin(mu_c)^2;
Ppyy=mu_s^2*sigma_c^2*sin(mu_c)^2+sigma_s^2*cos(mu_c)^2;
Ppxy=(sigma_s^2-mu_s^2*sigma_theta^2)*sin(mu_c)*cos(mu_c);

P1=[ [Pxx Pxy; Pxy Pyy] zeros(2,2) ;
    zeros(2,2) [Ppxx Ppxy; Ppxy Ppyy] ];

T=1;
f=[eye(2) T*eye(2);zeros(2) eye(2)];
f_t=inv(f.');
f_1=inv(f);
%initialisation
Jk=cell(t_f,1);
Jk{1}=inv(P1);
%loop
for i=2:26
    Jk{i}= f_t*Jk{i-1}*f_1+sigma_theta^(-2)*dG(X(:,i))*dG(X(:,i)).';
end

bound=zeros(26,1);
for i=1:26
   bound(i)=sqrt(Jk{i}(1,1)^-1 +Jk{i}(2,2)^-1); 
end



end

function[y_out] = dG(x_in)
    if x_in(1)>=0 && x_in(2)>=0
       y_out = [sign(x_in(1))/(abs(x_in(2))*(abs(x_in(1))^2/abs(x_in(2))^2 + 1));...
           -(abs(x_in(1))*sign(x_in(2)))/(abs(x_in(2))^2*(abs(x_in(1))^2/abs(x_in(2))^2 + 1)); 0 ; 0];
    elseif x_in(1)>0 && x_in(2)<0
       y_out = [-sign(x_in(1))/(abs(x_in(2))*(abs(x_in(1))^2/abs(x_in(2))^2 + 1));...
           (abs(x_in(1))*sign(x_in(2)))/(abs(x_in(2))^2*(abs(x_in(1))^2/abs(x_in(2))^2 + 1)); 0 ; 0];
    elseif x_in(1)<0 && x_in(2)<0
       y_out = [sign(x_in(1))/(abs(x_in(2))*(abs(x_in(1))^2/abs(x_in(2))^2 + 1));...
           -(abs(x_in(1))*sign(x_in(2)))/(abs(x_in(2))^2*(abs(x_in(1))^2/abs(x_in(2))^2 + 1)); 0 ; 0];
    else
       y_out = [-sign(x_in(1))/(abs(x_in(2))*(abs(x_in(1))^2/abs(x_in(2))^2 + 1));...
           (abs(x_in(1))*sign(x_in(2)))/(abs(x_in(2))^2*(abs(x_in(1))^2/abs(x_in(2))^2 + 1)); 0 ; 0];
    end;
end

