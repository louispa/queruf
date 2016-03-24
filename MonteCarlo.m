
%Sequential Monte Carlo applied on the problem of the project
%based on the received script by the professor Pierre-Antoine Absil
function[X]=MonteCarlo()
%%
%our parameters
%if nargin<1
%   show=0 
%end
%the systeme
%x_(k+1)=F (x_k)+Gamma * v_k
%z_k = G(x_k)+w_k
t_f=200;
d_x=4;
d_z=1;
d_v=2;
T=0.5;
%F is defined in the function below
Gamma = [T^2/2 * eye(2);T * eye(2)];
%G is defined in the function below
mu_v = 0;
Sigma_v = sqrt(0.01);%sqrt(variance)
mu_w= 0;
Sigma_w = sqrt(0.01); %sqrt(variance)

out_noise_pdf= @(w) 1/sqrt((2*pi)^d_z*abs(det(Sigma_w))) * exp(-.5*(w-mu_w)'*inv(Sigma_w)*(w-mu_w));; %normal sigma theta

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
%je pense qu'on devrait monter à 8 tableaux: 2 par valeur
X=cell(n,t_f +1); %right X
Xtilde=cell(n,t_f +1); %predicitons
t=0;
for i=1:n
    X{i,t +1} = ones(4,1);
end
%Beginning of the loop

for t=0:t_f-1
    %Prediction
    
    for i=1:n
        epsilon=Gamma*eye(2)*(mu_v+ Sigma_v.*randn(d_v,1)); %epsilon_k
        Xtilde{i,t+1 +1} = F(X{i,t +1})+ epsilon;
    end
    
    z = z_true(:,t+1 +1);
    %weights
    weights = zeros(1,n);
    for i=1:n
        weights(i) = out_noise_pdf(z-G(Xtilde{i,t+1 +1}));
    end
    ind_sample = randsample(n,n,true,weights);
    for i=1:n
        X{i,t+1 +1} = Xtilde{ind_sample(i),t+1 +1};
    end
    
end




%if show
    %draft(X,n);
%end

end



%%
%Function of the systems
%f(x)
function[x_out]=F(x_in)
T=0.5;
f=[eye(2) T*eye(2);zeros(2) eye(2)];
x_out= f*x_in;

end
%g(x)
function[y_out]=G(x_in)
y_out=atan(x_in(1)/x_in(2));
end

%%
%Function to draw the histograms and the plot of the trajectory
function[]=draft(X,n)
one=zeros(n,2);
for i=1:n
   one(i,:)=[X{i,1}(1) X{i,1}(2)]; 
end
fifty=zeros(n,2);
for i=1:n
   fifty(i,:)=[X{i,50}(1) X{i,1}(2)]; 
end
hundred=zeros(n,2);
for i=1:n
   hundred(i,:)=[X{i,100}(1) X{i,1}(2)]; 
end
twohundred=zeros(n,2);
for i=1:n
   twohundred(i,:)=[X{i,200}(1) X{i,1}(2)]; 
end
figure(1)
hold on;
subplot(2,2,1)
histogram(one(:,1));

end