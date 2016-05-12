function [ bound ] = Q6fun(observer,target,bearingMeasurements,mu_r,mu_theta,mu_s,mu_c,X)


t_f=26;
sigma_s_squared = 0.1;
sigma_theta_squared = 10^-4;
sigma_r_squared = 0.1;
sigma_c_squared =0.1;

Pxx=mu_r^2*sigma_theta_squared*cos(mu_theta)^2+sigma_r_squared*sin(mu_theta)^2;
Pyy=mu_r^2*sigma_theta_squared*sin(mu_theta)^2+sigma_r_squared*cos(mu_theta)^2;
Pxy=(sigma_r_squared-mu_r^2*sigma_theta_squared)*sin(mu_theta)*cos(mu_theta);

Ppxx=mu_s^2*sigma_c_squared*cos(mu_c)^2+sigma_s_squared*sin(mu_c)^2;
Ppyy=mu_s^2*sigma_c_squared*sin(mu_c)^2+sigma_s_squared*cos(mu_c)^2;
Ppxy=(sigma_s_squared-mu_s^2*sigma_theta_squared)*sin(mu_c)*cos(mu_c);

P1=[ [Pxx Pxy; Pxy Pyy] zeros(2,2) ;
    zeros(2,2) [Ppxx Ppxy; Ppxy Ppyy] ];

T=1;
f=[eye(2) T*eye(2);zeros(2) eye(2)];
f_t=inv(f');
f_1=inv(f);
%initialisation
Jk=cell(t_f,1);
Jk{1}=inv(P1);
%loop
for i=2:26
    Jk{i}= f_t*Jk{i-1}*f_1+sigma_theta_squared^-1*dG(X(:,i))*dG(X(:,i))';
end

bound=zeros(26,1);
for i=1:26
   P=inv(Jk{i});
   bound(i)=sqrt(P(1,1) +P(2,2)); 
end

end

%G= constant +- artan( abs(x/y)
function[y_out] = dG(x_in)
    if x_in(1)>=0 && x_in(2)>=0
       y_out = [1/(x_in(2)*(x_in(1)^2/x_in(2)^2 + 1));...
           -x_in(1)/(x_in(2)^2*(x_in(1)^2/x_in(2)^2 + 1)); 0 ; 0];
    elseif x_in(1)>0 && x_in(2)<0
       y_out = [-1/(x_in(2)*(x_in(1)^2/x_in(2)^2 + 1));...
           x_in(1)/(x_in(2)^2*(x_in(1)^2/x_in(2)^2 + 1)); 0 ; 0];
    elseif x_in(1)<0 && x_in(2)<0
       y_out = [1/(x_in(2)*(x_in(1)^2/x_in(2)^2 + 1));...
           -x_in(1)/(x_in(2)^2*(x_in(1)^2/x_in(2)^2 + 1)); 0 ; 0];
    else
       y_out = [-1/(x_in(2)*(x_in(1)^2/x_in(2)^2 + 1));...
           x_in(1)/(x_in(2)^2*(x_in(1)^2/x_in(2)^2 + 1)); 0 ; 0];
    end;
end
