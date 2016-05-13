function[X,Xtilde]=Q5fun(observer,bearingMeasurements,mu_r,mu_theta,mu_s,mu_c)
%%

% the system :
%   state equation: x_(k+1) = F(x_k) - U_(k,k+1) + Gamma*v_k
%   z_k = H(x_k) + w_k
%   F and H are defined in the function below

t_f = 25;
% dimensions of v and x
d_v=2;
d_x=4;
T=1; % time period

% reconstruction of state vector of dimension 4 of the observer by computing 
% the velocities at time k with (x(k+1) - x(k))/T. The velocity at time t_f
% is supposed to be the same as the one at time t_f-1

obs = zeros(4,t_f+1);
for k = 1:t_f
    obs(1:2,k) = observer(:,k);
    obs(3:4,k) = (observer(:,k+1) - observer(:,k))/T;
end
obs(1:2,t_f+1) = observer(1:2,t_f+1);
obs(3:4,t_f+1) = obs(3:4,t_f);
    
Gamma = [T^2/2 * eye(2);T * eye(2)];

% v is a zero-mean Gaussian noise with variance Sigma^2_a.

mu_v = 0;
Sigma_a = 0;
Sigma_r = sqrt(0.1); %sqrt(variance) of the relative distance
Sigma_theta = sqrt(10^(-4)); %sqrt(variance) of the initial bearing
Sigma_s = sqrt(0.1); %sqrt(variance) of the initial speed of the target
Sigma_c = sqrt(0.1); %sqrt(variance) of the initial course of the target

% true measurements (given)
z_true = bearingMeasurements; 

%%
%Start of the algorithm

%During the loops we only work on the relative state vector X
%We create the target state vector X^t afterwards by adding X^0 to X
n = 5000;
X=cell(n,t_f +1); %right relative X postions
Xtilde=cell(n,t_f +1); %predicitons for relative X positions 

% We have to sample from the distribustion of x0:
% c = mu_c+Sigma_c*randn(1,1);
% s = mu_s+Sigma_s*randn(1,1);
% r = mu_r+Sigma_r*randn(1,1);
% theta = mu_theta+Sigma_theta*randn(1,1);
% Afterwards : x0 = [r*sin(theta); r*cos(theta); s*sin(c)-obs(3,1); s*cos(c)-obs(4,1)];
t=0;
for i=1:n
    c = mu_c+Sigma_c*randn(1,1);
    s = mu_s+Sigma_s*randn(1,1);
    r = mu_r+Sigma_r*randn(1,1);
    theta = mu_theta+Sigma_theta*randn(1,1);
    X{i,t +1} = [r*sin(theta); r*cos(theta);s*sin(c) - obs(3,1); s*cos(c) - obs(4,1)];
end

%Beginning of the loop on time
e = zeros(t_f,1);
for t=0:t_f-1
    
    %PREDICTION
    
    for i=1:n
        epsilon = Gamma*(mu_v + Sigma_a.*randn(d_v,1)); %epsilon_k
        % state equation
        Xtilde{i,t+1 +1} = F(X{i,t +1}) - U(obs(:,t +1),obs(:,t+1 +1)) + epsilon;
    end
    
    % CORRECTION
    
    z = z_true(:,t+1 +1);
    
    %weights
    weights = zeros(1,n);
    for i=1:n
       weights(i) = W(z-H(Xtilde{i,t+1 +1}));
    end
    
    % REGULARIZATION
    
    % normalization of the weights
    somme = sum(weights);
    weights = weights/somme;
    
    % computation of Neff, Nth, h
    sumWeight=sum(weights.^2);
    Neff=1/sumWeight;
    Nth = n/3;
    h = (4/(n*(d_x+2)))^(1/(d_x+4));
    
    % creation of a 5000x4 matric for every sample i (i = 1,...,26)
    xI = Xtilde(:,t+1 +1);
    xI = cell2mat(xI);
    xI = reshape(xI,d_x,n);
    xI = xI';
    
    if Neff <= Nth
        e(t+1) = 1;
        S = cov(xI); % empirical covariance matrix
        A = sqrtm(S); % square root of the empirical covariance matrix
        
        % resampling
        ind_sample = randsample(n,n,true,weights);
        for i=1:n
            X{i,t+1 +1} = Xtilde{ind_sample(i),t+1 +1};
        end
        
        % generation of epsilon from the Gaussian Kernel and recursion step
        [~,ca]  = size(A);
        for i = 1:n
            epsilon = randn(ca,1);
            X{i,t+1 +1} = X{i,t+1 +1} + h*A*epsilon;
        end
    else
        % resampling
        ind_sample = randsample(n,n,true,weights);
        for i=1:n
            X{i,t+1 +1} = Xtilde{ind_sample(i),t+1 +1};
        end
    end
end
e
end
%%

%Functions of the system
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%f(x)
function[x_out]=F(x_in)
T=1;
f=[eye(2) T*eye(2);zeros(2) eye(2)];
x_out = f*x_in;
end

% h(x)
function[y_out] = H(x_in)
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

%u(x)
function[u_out] = U(x1,x2)
T = 1;
u_out = zeros(4,1);
u_out(1:2,1) = x2(1:2)-x1(1:2)-T*x1(3:4);
u_out(3:4,1) = x2(3:4)-x1(3:4);
end

% w(x)
% normal distribution zero mean and variance = 10^-4
function[w_out] = W(w)
Sigma_w = 0.01;
mu_w = 0;
w_out = 1/(sqrt(2*pi)*Sigma_w) * exp(-0.5*((w-mu_w)/Sigma_w)^2);
end
