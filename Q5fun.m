function[X,Xtilde]=Q5fun(observer,target,bearingMeasurements,mu_r,mu_theta,mu_s,mu_c)
%%

% 
% On reprend le script de PA absil avec t_f = 25 et non 26


% the system :
%   x_(k+1)=F (x_k)+Gamma * v_k
%   z_k = G(x_k)+w_k
%   F and G are defined in the function below

%t_f=26;
t_f = 25;
d_z=1;
d_v=2;
d_x=4;
T=1;
% reconstruction du state vector de dim 4 en calculant les vitesses par
% x(k+1) - x(k). La vitesse au temps t_f est supposee la meme que celle au
% temps t_f-1


obs = zeros(4,t_f+1);
for k = 1:t_f
    obs(1:2,k) = observer(:,k);
    obs(3:4,k) = (observer(:,k+1) - observer(:,k))/T;
end
obs(1:2,t_f+1) = observer(1:2,t_f+1);
obs(3:4,t_f+1) = obs(3:4,t_f);
    
Gamma = [T^2/2 * eye(2);T * eye(2)];
mu_v = 0;
%Sigma_v = sqrt(10^(-6));%sqrt(variance)
Sigma_v = 0;
Sigma_r = sqrt(0.1); %sqrt(variance) of the relative distance
Sigma_theta = sqrt(10^(-4)); %sqrt(variance) of the initial bearing
Sigma_s = sqrt(0.1); %sqrt(variance) of the initial speed of the target
Sigma_c = sqrt(0.1); %sqrt(variance) of the initial course of the target

%  Sigma_v = 0;%sqrt(variance)
%  Sigma_r = 0; %sqrt(variance) of the relative distance
%  Sigma_theta = 0; %sqrt(variance) of the initial bearing
%  Sigma_s = 0; %sqrt(variance) of the initial speed of the target
%  Sigma_c = 0; %sqrt(variance) of the initial course of the target

%out_noise_pdf= @(w) 1/sqrt((2*pi)^d_z*abs(det(Sigma_theta)))...
%    * exp(-.5*w'*inv(Sigma_theta)*w); %normal sigma theta

x_true = target-observer;
z_true = bearingMeasurements; %according to the hypothesis, z is known

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
% Afterwards : x0 = [r*sin(theta); r*cos(theta); s*sin(c); s*cos(c)];
t=0;
for i=1:n
    c = mu_c+Sigma_c*randn(1,1);
    s = mu_s+Sigma_s*randn(1,1);
    r = mu_r+Sigma_r*randn(1,1);
    theta = mu_theta+Sigma_theta*randn(1,1);
    X{i,t +1} = [r*sin(theta); r*cos(theta);s*sin(c) - obs(3,1); s*cos(c) - obs(4,1)];
end

%Beginning of the loop on time
for t=0:t_f-1
    
    %PREDICTION
    
    for i=1:n
        epsilon = Gamma*(mu_v + Sigma_v.*randn(d_v,1)); %epsilon_k
        Xtilde{i,t+1 +1} = F(X{i,t +1}) - U(obs(:,t +1),obs(:,t+1 +1)) + epsilon;
        %Xtilde{i,t+1 +1} = F(X{i,t +1}) + epsilon;
    end
    % CORRECTION
    z = z_true(:,t+1 +1);
    %weights
    weights = zeros(1,n);
    for i=1:n
       weights(i) = W(z-G(Xtilde{i,t+1 +1}));
    end
    
    % merci Quentin :)
    sumWeight=sum(weights.*weights);
    Neff=1/sumWeight;
    Nth = n/3;
    nx = 4;
    xI = Xtilde(:,t+1 +1);
    xI = cell2mat(xI);
    xI = reshape(xI,d_x,n);
    xk = zeros(d_x, n);
    h = (4/(n*(nx+2)))^(1/(nx+4));
    for i=1:d_x
        if Neff < Nth
            try
            [xk(i,:), ~] = ksdensity(sort(xI(i,:)), 'npoints', n, 'bandwidth', h, 'weights', weights, 'function', 'icdf');
            catch error
            ind_sample = randsample(n,n,true,weights);
            x = xI(i,:);
            xk(i,:) = x(ind_sample);
            %xk(i,:) = xI(ind_sample(i));
            end
        else
            ind_sample = randsample(n,n,true,weights);
            x = xI(i,:);
            xk(i,:) = x(ind_sample);
            %xk(i,:) = xI(ind_sample(i));
        end
    end
    
    for i=1:n
        X{i,t+1 +1} = xk(:,i);
    end
    
    % resampling
    %ind_sample = randsample(n,n,true,weights);
    %for i=1:n
    %    X{i,t+1 +1} = Xtilde{ind_sample(i),t+1 +1};
    %end
end
end
%%

%Functions of the system

%f(x)
function[x_out]=F(x_in)
T=1;
f=[eye(2) T*eye(2);zeros(2) eye(2)];
x_out = f*x_in;
end

%g(x)

% g de Quentin
% function[y_out] = G(x_in)
%     if x_in(1)>0 && x_in(2)>0
%        y_out = atan(x_in(1)/x_in(2));
%     elseif x_in(1)<0 && x_in(2)>0
%        y_out = atan(x_in(2)/x_in(1))+3*pi/2;
%     elseif x_in(1)<0 && x_in(2)<0
%        y_out = pi + atan(x_in(1)/x_in(2));
%     else
%        y_out = atan(x_in(2)/x_in(1))+pi/2;
%     end;
% end

function[y_out] = G(x_in)
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

%  function[y_out]= G(x_in)
%     y_out = atan(x_in(1)/x_in(2));
%  end

% function[y_out]=G(x_in)
%     if x_in(1)>=0 && x_in(2)>=0
%        y_out = atan((x_in(1)/x_in(2)));
%     elseif x_in(1)>0 && x_in(2)<0
%        y_out = pi - atan(x_in(1)/x_in(2));
%     elseif x_in(1)<0 && x_in(2)<0
%        y_out = pi + atan(x_in(1)/x_in(2));
%     else
%        y_out = 2*pi - atan(x_in(1)/x_in(2));
%     end;
% end

%u(x)
% on a besoin de la vitesse de l'obs pour calculer u. Mais on a que sa
% position... quid?
% x1 au temps k et x2 au temps k+1
function[u_out] = U(x1,x2)
T = 1;
u_out = zeros(4,1);
u_out(1:2,1) = x2(1:2)-x1(1:2)-T*x1(3:4);
u_out(3:4,1) = x2(3:4)-x1(3:4);
end

function[u_out]=U2(x1,x2)%version de l'autre
T=1;
u_out=zeros(4,1);
u_out(1:2)=T*(x2(3:4)-x1(3:4));
u_out(3:4) = x2(3:4)-x1(3:4);
end


function[w_out] = W(w)
Sigma_w = 0.01;
mu_w = 0;
w_out = 1/(sqrt(2*pi)*Sigma_w) * exp(-0.5*((w-mu_w)/Sigma_w)^2); %normal sigma theta
%w_out = 1/sqrt((2*pi)^d_y*abs(det(Sigma_w))) * exp(-.5*(w-mu_w)'*inv(Sigma_w)*(w-mu_w)); 
end


% % 
% % On reprend le script de PA absil avec t_f = 25 et non 26
% 
% 
% % the system :
% %   x_(k+1)=F (x_k)+Gamma * v_k
% %   z_k = G(x_k)+w_k
% %   F and G are defined in the function below
% 
% %t_f=26;
% t_f = 25;
% d_z=1;
% d_v=2;
% T=1;
% 
% % reconstruction du state vector de dim 4 en calculant les vitesses par
% % x(k+1) - x(k). La vitesse au temps t_f est supposee la meme que celle au
% % temps t_f-1
% 
% obs = zeros(4,t_f+1);
% for k = 1:t_f
%     obs(1:2,k) = observer(:,k);
%     obs(3:4,k) = (observer(:,k+1) - observer(:,k))/T;
% end
% obs(1:2,t_f+1) = observer(1:2,t_f+1);
% obs(3:4,t_f+1) = obs(3:4,t_f);
%     
% Gamma = [T^2/2 * eye(2);T * eye(2)];
% mu_v = 0;
% Sigma_v = sqrt(10^(-6));%sqrt(variance)
% Sigma_r = sqrt(0.1); %sqrt(variance) of the relative distance
% Sigma_theta = sqrt(10^(-4)); %sqrt(variance) of the initial bearing
% Sigma_s = sqrt(0.1); %sqrt(variance) of the initial speed of the target
% Sigma_c = sqrt(0.1); %sqrt(variance) of the initial course of the target
% 
% %out_noise_pdf= @(w) 1/sqrt((2*pi)^d_z*abs(det(Sigma_theta)))...
% %    * exp(-.5*w'*inv(Sigma_theta)*w); %normal sigma theta
% 
% x_true = target-observer;
% z_true = bearingMeasurements; %according to the hypothesis, z is known
% 
% %%
% %Start of the algorithm
% 
% %During the loops we only work on the relative state vector X
% %We create the target state vector X^t afterwards by adding X^0 to X
% n = 5000;
% Nth=n/3;
% X=cell(n,t_f +1); %right relative X postions
% Xtilde=cell(n,t_f +1); %predicitons for relative X positions 
% 
% % We have to sample from the distribustion of x0:
% % c = mu_c+Sigma_c*randn(1,1);
% % s = mu_s+Sigma_s*randn(1,1);
% % r = mu_r+Sigma_r*randn(1,1);
% % theta = mu_theta+Sigma_theta*randn(1,1);
% % Afterwards : x0 = [r*sin(theta); r*cos(theta); s*sin(c); s*cos(c)];
% t=0;
% for i=1:n
%     c = mu_c+Sigma_c*randn(1,1);
%     s = mu_s+Sigma_s*randn(1,1);
%     r = mu_r+Sigma_r*randn(1,1);
%     theta = mu_theta+Sigma_theta*randn(1,1);
%     X{i,t +1} = [r*sin(theta); r*cos(theta); s*sin(c); s*cos(c)];
% end
% 
% %Beginning of the loop on time
% for t=0:t_f-1
%     
%     %PREDICTION
%     
%     for i=1:n
%         epsilon = Gamma*(mu_v + Sigma_v.*randn(d_v,1)); %epsilon_k
%         Xtilde{i,t+1 +1} = F(X{i,t +1}) - U(obs(:,t +1),obs(:,t+1 +1)) + epsilon;
%         %Xtilde{i,t+1 +1} = F(X{i,t +1}) + epsilon;
%     end    
%     % CORRECTION
%     z = z_true(:,t+1 +1);
%     %weights
%     weights = zeros(1,n);
%     for i=1:n
%        weights(i) = W(z-G(Xtilde{i,t+1 +1}));
%     end
%     
%     sumWeight=sum(weights.*weights);
%     Neff=1/sumWeight;
%     if Neff < Nth
%         % resampling
%         ind_sample = randsample(n,n,true,weights);
%         for i=1:n
%             X{i,t+1 +1} = Xtilde{ind_sample(i),t+1 +1};
%         end
%         A=(2/3)^(1/8);
%         h=A*n^(-1/8);
%         for i=1:n
%            epsilonV=randn(2,1);
%            X{i,t+1 +1}=X{i,t+1 +1} + h*Gamma*epsilonV;
%         end
%     end
% end
% end
% %%
% 
% %Functions of the system
% 
% %f(x)
% function[x_out]=F(x_in)
% T=1;
% f=[eye(2) T*eye(2);zeros(2) eye(2)];
% x_out = f*x_in;
% end
% 
% %g(x)
% 
% % g de Quentin
% % function[y_out] = G(x_in)
% %     if x_in(1)>0 && x_in(2)>0
% %        y_out = atan(x_in(1)/x_in(2));
% %     elseif x_in(1)<0 && x_in(2)>0
% %        y_out = atan(x_in(1)/x_in(2))+3*pi/2;
% %     elseif x_in(1)<0 && x_in(2)<0
% %        y_out = pi + atan(x_in(1)/x_in(2));
% %     else
% %        y_out = atan(x_in(1)/x_in(2))+pi/2;
% %     end;
% % end
% 
% function[y_out] = G(x_in)
%     if x_in(1)>=0 && x_in(2)>=0
%        y_out = atan(abs(x_in(1)/x_in(2)));
%     elseif x_in(1)>0 && x_in(2)<0
%        y_out = pi - atan(abs(x_in(1)/x_in(2)));
%     elseif x_in(1)<0 && x_in(2)<0
%        y_out = pi + atan(abs(x_in(1)/x_in(2)));
%     else
%        y_out = 2*pi - atan(abs(x_in(1)/x_in(2)));
%     end;
% end
% 
% % function[y_out]=G(x_in)
% %     y_out = atan(abs(x_in(1)/x_in(2)));
% % end
% 
% % function[y_out]=G(x_in)
% %     if x_in(1)>=0 && x_in(2)>=0
% %        y_out = atan((x_in(1)/x_in(2)));
% %     elseif x_in(1)>0 && x_in(2)<0
% %        y_out = pi - atan(x_in(1)/x_in(2));
% %     elseif x_in(1)<0 && x_in(2)<0
% %        y_out = pi + atan(x_in(1)/x_in(2));
% %     else
% %        y_out = 2*pi - atan(x_in(1)/x_in(2));
% %     end;
% % end
% 
% %u(x)
% % on a besoin de la vitesse de l'obs pour calculer u. Mais on a que sa
% % position... quid?
% % x1 au temps k et x2 au temps k+1
% function[u_out] = U(x1,x2)
% T = 1;
% u_out = zeros(4,1);
% u_out(1:2,1) = -x2(1:2)+x1(1:2);%+T*x1(3:4);
% u_out(3:4,1) = x2(3:4)-x1(3:4);
% end
% 
% function[w_out] = W(w)
% Sigma_w = 10^(-2);
% mu_w = 0;
% w_out = 1/sqrt((2*pi*Sigma_w)) * exp((-0.5/Sigma_w)*(w-mu_w)^2); %normal sigma theta
% end
% 
% %%
% 
% %Function to draw the figure with the two subplots
% % function[]=draft(obs,real_target,X,Xtilde)
% % subplot(2,3,1)
% % % observer's trajectory
% % plot(obs(1,:),obs(2,:))
% % end