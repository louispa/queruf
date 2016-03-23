
%Sequential Monte Carlo applied on the problem of the project
function[]=MonteCarlo()
%%
%our parameters

%the systeme
%x_(k+1)=F (x_k)+Gamma * v_k
%z_k = G(x_k)+w_k
t_f=200;
d_x=4;
d_z=1;
d_v=4;
T=0.5;
%F is defined in the function below
Gamma = [T^2/2 * eye(2);T*eye(2)];
%G is defined in the function below
mu_v = 0;
Sigma_v = sqrt(0.01);%sqrt(variance)
mu_w= 0;
Sigma_w = sqrt(0.01); %sqrt(variance)

out_noise_pdf=0; %pour l'isntant n'est-il pas :p

x_true = zeros(d_x,t_f +1);
z_true = zeros(d_z,t_f +1);

x_true(:,0 +1)=ones(4,1);%initialisation
for t= 0:t_f-1
   v_true=Gamma*eye(2)*(mu_v+ Sigma_v.*randn(d_v,1));
   x_true(:,t+1 +1)=F(x_true(:,t +1))+ v_true;
   w_true=mu_w+ Sigma_w*randn(d_z,1);
   z_true(:,t +1)=G(x_true(:,t +1))+w_true;
end
w_true=mu_w+ Sigma_w*randn(d_z,1);
z_true(:,t_f +1)=G(x_true(:,t_f +1))+w_true;


%%
%Start of the algorithm

%ca marche pas pour l'instant ici car on devrait avoir un tableau de
%5000*200*4 
n= 5000;
X=cell(n,t_f +1); %right X
Xtilde=cell(n,t_f +1); %predicitons
t=0;
X{:,t +1} = ones(4,1);

%Beginning of the loop

for t=0:t_f-1
    %Prediction
    v=Gamma*eye(2)*(mu_v+ Sigma_v.*randn(d_v,1));
    Xtilde={i,t+1 +1}=0;
end


end
%%
%Function of the systems
%f(x)
function[x_out]=F(x_in)
%delta T
f=[eye(2) T*eye(2);zeros(2) eye(2)]; 
x_out= f*x_in;

end
%g(x)
function[y_out]=G(x_in)
y_out=atan(x_in(1)/x_in(2));
end