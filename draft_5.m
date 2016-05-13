
function[]=draft_5(observer,real_target,X,Xtilde,X1,Xtilde1)
%%
t_f = 25;

% computation of the observer's and target's velocities at time k with 
% (x(k+1) - x(k))/T. The velocity at time t_f is supposed to be the same as
% the one at time t_f-1
T = 1;
obs = zeros(4,t_f+1);
v_target = zeros(4,t_f+1);
for k = 1:t_f
    obs(1:2,k) = observer(:,k);
    obs(3:4,k) = (observer(:,k+1) - observer(:,k))/T;
    %v_target(1:2,k) = real_target(:,k);
    %v_target(3:4,k) = (real_target(:,k+1) - real_target(:,k))/T;
end
obs(1:2,t_f+1) = observer(1:2,t_f+1);
obs(3:4,t_f+1) = obs(3:4,t_f);
%v_target(1:2,t_f+1) = real_target(1:2,t_f+1);
%v_target(3:4,t_f+1) = v_target(3:4,t_f);

figure(1)
plot(obs(1,:),obs(2,:),'.'); hold on; % real trajectory of the observer
plot(real_target(1,:),real_target(2,:),'.')% real trajectory of the target

%plot the predicted target's trajectory
target=zeros(4,t_f+1);
for t=1:t_f+1
    helper=[0 0 0 0]';
    for i=1:length(X)
        helper=helper+X{i,t}; %size(X) = 5000 26
    end
    target(:,t)=helper/length(X);
end
pred_traj = zeros(2,t_f+1);
pred_traj(1:2,:) = target(1:2,:) + obs(1:2,:);
plot(pred_traj(1,:),pred_traj(2,:),'g.'); % predicted target's trajectory
title('trajectories of the regularized particle filter')
xlabel('x')
ylabel('y')

figure(2);
% number of particles after resampling
numberx = zeros(1,t_f+1);
% number of particles before resampling
numberxt = zeros(1,t_f+1);
% since Xtilde(:,1) is empty
Xtilde(:,1) = X(:,1);
for t = 0:t_f
    % generation of a 4x5000 matric for every sample after resampling
    x = X(:,t+1);
    x = cell2mat(x);
    x = reshape(x,4,5000);
    y = x';
    % we only keep unique elements and put the number in numberx at the
    % right index
    unix = unique(y,'rows');
    numberx(t+1) = length(unix(:,1));
    
    % generation of a 4x5000 matric for every sample before resampling
    xt = Xtilde(:,t+1);
    xt = cell2mat(xt);
    xt = reshape(xt,4,5000);
    yt = xt';
    % we only keep unique elements and put the number in numberx at the
    % right index
    unixt = unique(yt,'rows');
    numberxt(t+1) = length(unixt(:,1));
end
% semilogarithmic plots. Log scale on the y_axis
semilogy(numberx,'-ob'); hold on;
semilogy(numberxt,'-or');
title('number of different particles for the regularized particle filter')
xlabel('k')
ylabel('n')

figure(3)
plot(obs(1,:),obs(2,:),'.'); hold on; % real trajectory of the observer
plot(real_target(1,:),real_target(2,:),'.')% real trajectory of the target

%plot the predicted target's trajectory
target=zeros(4,t_f+1);
for t=1:t_f+1
    helper=[0 0 0 0]';
    for i=1:length(X1)
        helper=helper+X1{i,t}; %size(X) = 5000 26
    end
    target(:,t)=helper/length(X1);
end
pred_traj = zeros(2,t_f+1);
pred_traj(1:2,:) = target(1:2,:) + obs(1:2,:);
plot(pred_traj(1,:),pred_traj(2,:),'g.'); % predicted target's trajectory
title('trajectories for the particle filter')
xlabel('x')
ylabel('y')

figure(4);
numberx = zeros(1,t_f+1);
numberxt = zeros(1,t_f+1);
Xtilde1(:,1) = X1(:,1);
for t = 0:t_f
    x = X1(:,t+1);
    x = cell2mat(x);
    x = reshape(x,4,5000);
    y = x';
    unix = unique(y,'rows');
    numberx(t+1) = length(unix(:,1));
    
    xt = Xtilde1(:,t+1);
    xt = cell2mat(xt);
    xt = reshape(xt,4,5000);
    yt = xt';
    unixt = unique(yt,'rows');
    numberxt(t+1) = length(unixt(:,1));
end
semilogy(numberx,'-ob'); hold on;
semilogy(numberxt,'-or');
title('number of different particles for the particle filter')
xlabel('k')
ylabel('n')
end