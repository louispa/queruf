
% Sequential Monte Carlo applied on the problem of the project
% based on the script received from professor Pierre-Antoine Absil.

function[X,Xtilde,n,t_f] = Q3fun()
% the system
% x_(k+1)=F (x_k)+Gamma * v_k
% z_k = H(x_k)+w_k

% time steps
t_f=200;
% dimensions of the vectors x,z,v
d_x=4;
d_z=1;
d_v=2;
% period
T=0.5;
% F and G are defined in the function below
Gamma = [T^2/2 * eye(2);T * eye(2)];
% v is a zero mean noise of variance sigma^2_v*I (iid)
mu_v = 0;
Sigma_a = sqrt(0.01);
% w is a zero mean noise of variance sigma^2_w
mu_w= 0;
Sigma_w = sqrt(0.01);

% true positions and true measurements
x_true = zeros(d_x,t_f +1);
z_true = zeros(d_z,t_f +1);

% initialisation (x0 = [1 1 1 1]')
x_true(:,0 +1)=ones(4,1);

% computation of x_true and z_true
for t= 0:t_f-1
    v_true = Gamma*(mu_v + Sigma_a.*randn(d_v,1));
    x_true(:,t+1 +1) = F(x_true(:,t +1)) + v_true;
    w_true = mu_w + Sigma_w*randn(d_z,1);
    z_true(:,t +1)=H(x_true(:,t +1)) + w_true;
end
w_true = mu_w + Sigma_w*randn(d_z,1);
z_true(:,t_f +1) = H(x_true(:,t_f +1)) + w_true;


%%
%               *** SEQUENTIAL MONTE CARLO METHOD ***

n= 5000;
X=cell(n,t_f +1); %right X
Xtilde=cell(n,t_f +1); %predicitons

% generating initial sample set {x_0^i,...,x_0^n}:
% normally we sample from the distribution of x0 but here it is
% deterministic
t=0;
for i=1:n
    X{i,t +1} = ones(4,1);
end

%Beginning of the loop on time
for t=0:t_f-1
    
    % PREDICTION
    
    for i=1:n
        epsilon=Gamma*(mu_v+ Sigma_a.*randn(d_v,1)); %epsilon_k
        Xtilde{i,t+1 +1} = F(X{i,t +1})+ epsilon;
    end
    
    % CORRECTION
    
    z = z_true(:,t+1 +1);
    
    %weights
    weights = zeros(1,n);
    for i=1:n
        % w = z - h(x)
        weights(i) = W(z-H(Xtilde{i,t+1 +1}));
    end
    % resampling
    ind_sample = randsample(n,n,true,weights);
    for i=1:n
        X{i,t+1 +1} = Xtilde{ind_sample(i),t+1 +1};
    end
end
end



%%

% Functions of the system
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%F(x)
function[x_out]=F(x_in)
T=0.5;
f=[eye(2) T*eye(2);zeros(2) eye(2)];
x_out= f*x_in;
end

%H(x)
function[y_out]=H(x_in)
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

%W(w)
function y_out = W(w)
% w is a zero mean noise of variance sigma^2_w
mu_w= 0;
Sigma_w = sqrt(0.01);
d_z = 1;
% probability distribution function of noise w (Gaussian)
y_out = 1/sqrt((2*pi)^d_z*abs(det(Sigma_w^2))) * exp(-.5*(w-mu_w)'*inv(Sigma_w^2)*(w-mu_w)); 
end
