
%Sequential Monte Carlo applied on the problem of the project
%based on the received script by the professor Pierre-Antoine Absil
function[X_target,Xtilde_target]=Q4fun(observer,target,bearingMeasurements,r,theta,s,c)
%%
%the systeme
%x_(k+1)=F (x_k)+Gamma * v_k
%z_k = G(x_k)+w_k
t_f=26;
d_x=4;
d_z=1;
d_v=2;
T=1;
%F is defined in the function below
Gamma = [T^2/2 * eye(2);T * eye(2)];
%G is defined in the function below
mu_v = 0;
Sigma_v = sqrt(10^(-6) );%sqrt(variance)

Sigma_r = sqrt(0.1); %sqrt(variance) of the relative distance

Sigma_theta = sqrt(10^(-4) ); %sqrt(variance) of the initial bearing

Sigma_s = sqrt(0.1); %sqrt(variance) of the speed target

Sigma_c = sqrt(0.1); %sqrt(variance) of the course of the target

out_noise_pdf= @(w) 1/sqrt((2*pi)^d_z*abs(det(Sigma_theta)))...
    * exp(-.5*(w-theta)'*inv(Sigma_theta)*(w-theta)); %normal sigma theta

x_true = zeros(d_x,t_f +1);%we will work on that( the relative distance)
z_true = bearingMeasurements; %according to the hypothesis, z is known

x_true(:,0 +1)=[c+Sigma_c*randn(2,1);s+Sigma_s*randn(2,1)]-...
    [observer(:,1);zeros(2,1)];%initialisation DOIT ETRE MODIFIE
for t=0:t_f-1
    epsilon_true=Gamma*eye(2)*(mu_v+ Sigma_v.*randn(d_v,1)); %process noise
    x_true(:,t+1 +1)=F(x_true(:,t +1))+ epsilon_true;
end

%%
%Start of the algorithm

%During the loops we only work on the relative state vector X
%We create the target state vector X^t afterwards by adding X^0 to X
n= 5000;
X=cell(n,t_f +1); %right relative X postions
Xtilde=cell(n,t_f +1); %predicitons for relative X positions 
t=0;
for i=1:n
    X{i,t +1} =[c+Sigma_c*randn(2,1);s+Sigma_s*randn(2,1)]-...
        [(r+Sigma_r*randn(2,1));zeros(2,1)];
end
%Beginning of the loop

for t=0:t_f-1
    %Prediction
    
    for i=1:n
        epsilon=Gamma*eye(2)*(mu_v+ Sigma_v.*randn(d_v,1)); %epsilon_k
        Xtilde{i,t+1 +1} = F(X{i,t +1})+ epsilon;
    end
    
    z = z_true(t+1);
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

%%Création of the neeeded vector to plot
Xtilde_target=cell(n,t_f);%from x^t_1 to x^t_n
X_target=cell(n,t_f);
for i=1:n
   for j=1:t_f
       Xtilde_target{i,j}=Xtilde{i,j+1}+[observer(:,t_f);zeros(2,1)];
       X_target{i,j}=X{i,j+1}+[observer(:,t_f);zeros(2,1)];
   end
end

%draft(observer,target,X_target,Xtilde_target)

end



%%
%Function of the systems
%f(x)
function[x_out]=F(x_in)
T=1;
f=[eye(2) T*eye(2);zeros(2) eye(2)];
x_out= f*x_in;

end
%g(x)
function[y_out]=G(x_in)
y_out=atan(x_in(1)/x_in(2));
end

%%
%Function to draw the figure with the two subplots
function[]=draft(obs,real_target,X,Xtilde)
subplot(2,3,1)
plot(obs(1,:),obs(2,:))

end