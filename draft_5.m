
function[]=draft_5(observer,real_target,X,Xtilde,X1,Xtilde1)
%%
t_f = 25;
% computation of the observer's and target's velocities
T = 1;
obs = zeros(4,t_f+1);
v_target = zeros(4,t_f+1);
for k = 1:t_f
    obs(1:2,k) = observer(:,k);
    obs(3:4,k) = (observer(:,k+1) - observer(:,k))/T;
    v_target(1:2,k) = real_target(:,k);
    v_target(3:4,k) = (real_target(:,k+1) - real_target(:,k))/T;
end
obs(1:2,t_f+1) = observer(1:2,t_f+1);
obs(3:4,t_f+1) = obs(3:4,t_f);
v_target(1:2,t_f+1) = real_target(1:2,t_f+1);
v_target(3:4,t_f+1) = v_target(3:4,t_f);

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
numberx = zeros(1,t_f+1);
numberxt = zeros(1,t_f+1);
Xtilde(:,1) = X(:,1);
for t = 0:t_f
    x = X(:,t+1);
    x = cell2mat(x);
    x = reshape(x,4,5000);
    y = x';
    unix = unique(y,'rows');
    numberx(t+1) = length(unix(:,1));
    
    xt = Xtilde(:,t+1);
    xt = cell2mat(xt);
    xt = reshape(xt,4,5000);
    yt = xt';
    unixt = unique(yt,'rows');
    numberxt(t+1) = length(unixt(:,1));
end
semilogy(numberx,'b'); hold on;
semilogy(numberxt,'r');
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
semilogy(numberx,'b'); hold on;
semilogy(numberxt,'r');
title('number of different particles for the particle filter')
xlabel('k')
ylabel('n')
end