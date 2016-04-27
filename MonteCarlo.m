
%Sequential Monte Carlo applied on the problem of the project
%based on the received script by the professor Pierre-Antoine Absil

function[X]=MonteCarlo()
%the system
%x_(k+1)=F (x_k)+Gamma * v_k
%z_k = G(x_k)+w_k
% time steps
t_f=200;
% dimensions of the vectors x,z,v
d_x=4;
d_z=1;
d_v=2;
% period
T=0.5;
%F and G are defined in the function below
Gamma = [T^2/2 * eye(2);T * eye(2)];
% v is a zero mean noise of variance sigma_v*I (iid)
mu_v = 0;
Sigma_v = sqrt(0.01);
% w is a zero mean noise of variance sigma_w
mu_w= 0;
Sigma_w = sqrt(0.01);
% probability distribution function of noise w (Gaussian)
out_noise_pdf= @(w) 1/sqrt((2*pi)^d_z*abs(det(Sigma_w^2))) * exp(-.5*(w-mu_w)'*inv(Sigma_w^2)*(w-mu_w)); %normal sigma theta

x_true = zeros(d_x,t_f +1);
z_true = zeros(d_z,t_f +1);

%initialisation (x0 = [1 1 1 1]')
x_true(:,0 +1)=ones(4,1);

for t= 0:t_f-1
    v_true=Gamma*(mu_v+ Sigma_v.*randn(d_v,1));
    x_true(:,t+1 +1)=F(x_true(:,t +1))+ v_true;
    w_true=mu_w+ Sigma_w*randn(d_z,1);
    z_true(:,t +1)=G(x_true(:,t +1))+w_true;
end
w_true=mu_w+ Sigma_w*randn(d_z,1);
z_true(:,t_f +1)=G(x_true(:,t_f +1))+w_true;


%%
%               *** SEQUENTIAL MONTE CARLO METHOD ***

n= 5000;
X=cell(n,t_f +1); %right X
Xtilde=cell(n,t_f +1); %predicitons

% generating initial sample set {x_0^i,...,x_0^n}:
t=0;
for i=1:n
    X{i,t +1} = ones(4,1);
end

%Beginning of the loop on time
for t=0:t_f-1
    
    % PREDICTION
    
    for i=1:n
        epsilon=Gamma*(mu_v+ Sigma_v.*randn(d_v,1)); %epsilon_k
        Xtilde{i,t+1 +1} = F(X{i,t +1})+ epsilon;
    end
    
    % CORRECTION
    
    z = z_true(:,t+1 +1);
    %weights
    weights = zeros(1,n);
    for i=1:n
        % w = z-G(x)
        weights(i) = out_noise_pdf(z-G(Xtilde{i,t+1 +1}));
    end
    % resampling
    ind_sample = randsample(n,n,true,weights);
    for i=1:n
        X{i,t+1 +1} = Xtilde{ind_sample(i),t+1 +1};
    end
    
end

draft(X,n);

end



%%

%Functions of the system
%F(x)
function[x_out]=F(x_in)
T=0.5;
f=[eye(2) T*eye(2);zeros(2) eye(2)];
x_out= f*x_in;
end

%G(x)
function[y_out]=G(x_in)
    if x_in(1)>=0 && x_in(2)>=0
       y_out = atan(abs(x_in(1)/x_in(2)));
    elseif x_in(1)>0 && x_in(2)<0
       y_out = pi - atan(abs(x_in(1)/x_in(2)));
    elseif x_in(1)<0 && x_in(2)<0
       y_out = pi + atan(abs(x_in(1)/x_in(2)));
    else
       y_out = 2*pi - atan(abs(x_in(1)/x_in(2)));
    end;
end

%%

function[]=draft(X,n)
%Draw the histograms at time t=1,50,100 and 200 for x and y
%Draw the trajectory of the target 
one=zeros(n,2);
%our first state is in t=0
for i=1:n
   one(i,:)=[X{i,2}(1) X{i,2}(2)]; 
end
fifty=zeros(n,2);
for i=1:n
   fifty(i,:)=[X{i,51}(1) X{i,51}(2)]; 
end
hundred=zeros(n,2);
for i=1:n
   hundred(i,:)=[X{i,101}(1) X{i,101}(2)]; 
end
twohundred=zeros(n,2);
for i=1:n
   twohundred(i,:)=[X{i,201}(1) X{i,201}(2)]; 
end
figure(1)
hold on;
subplot(2,2,1)
histogram(one(:,1));
title('Histograms for x at t=1');
subplot(2,2,2)
histogram(fifty(:,1));
title('Histograms for x at t=50');
subplot(2,2,3)
histogram(hundred(:,1));
title('Histograms for x at t=100');
subplot(2,2,4)
histogram(twohundred(:,1));
title('Histograms for x at t=200');

hold off;
figure(2)
hold on;
subplot(2,2,1)
histogram(one(:,2));
title('Histograms for y at t=1');
subplot(2,2,2)
histogram(fifty(:,2));
title('Histograms for y at t=50');
subplot(2,2,3)
histogram(hundred(:,2));
title('Histograms for y at t=100');
subplot(2,2,4)
histogram(twohundred(:,2));
title('Histograms for y at t=200');

traj=zeros(200,2);
for t=2:length(X(1,:))
    helper=0;
   for i=1:n
       helper=helper+X{i,t}([1 2]);
   end
    traj(i,:)=helper./n;
end

figure(3)
plot(traj(:,1),traj(:,2));
title('Trajectory of the target');

end