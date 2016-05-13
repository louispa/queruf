function[]=draft_4(observer,real_target,X,Xtilde)
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
    v_target(1:2,k) = real_target(:,k);
    v_target(3:4,k) = (real_target(:,k+1) - real_target(:,k))/T;
end
obs(1:2,t_f+1) = observer(1:2,t_f+1);
obs(3:4,t_f+1) = obs(3:4,t_f);
v_target(1:2,t_f+1) = real_target(1:2,t_f+1);
v_target(3:4,t_f+1) = v_target(3:4,t_f);


% We will comment the first plot since we repeat the almost exact same code
% 4 times. Comments are therefore similar


% k = 1
% ------------------------------------------------------------------------
% plot the true trajectory of the observer and the target
figure(1)

% first subplot for the different trajectories and the position of the 
% particles 
subplot(2,1,1)
plot(obs(1,:),obs(2,:),'.'); hold on; % real trajectory of the observer
plot(real_target(1,:),real_target(2,:),'.')% real trajectory of the target
% plot the predicted target's trajectory
target=zeros(4,1);
for t=1:1
    helper=[0 0 0 0]';
    for i=1:length(X)
        helper=helper+X{i,t}; %size(X) = 5000 26
    end
    target(:,t)=helper/length(X);
end
pred_traj = zeros(2,1);
pred_traj(1:2,:) = target(1:2,:) + obs(1:2,1);
plot(pred_traj(1,:),pred_traj(2,:),'g.'); % predicted target's trajectory

% the particules will be stored by row
particules_before=zeros(5000,4);
particules_after=zeros(5000,4);
for j=1:length(X)
     xtilde = X{j,1};
     x = X{j,1};
     particules_before(j,:)= xtilde'+ obs(:,1)';
     particules_after(j,:)= x'+ obs(:,1)';
end

% particles before resampling are plotted in magenta.
% particles after resampling are plotted in blue
plot(particules_before(:,1),particules_before(:,2),'m.');
plot(particules_after(:,1),particules_after(:,2),'y.');
title('positions');
xlabel('x');
ylabel('y');

% second subplot for the different velocities
% velocities of the particles before resampling are plotted in blue
% velocities of the particles after resampling are plotted in yellow
% true velocity of the target is plotted in red (constant (0.04,0.04))
subplot(2,1,2)
plot(particules_before(:,3),particules_before(:,4),'m.'); hold on;
plot(particules_after(:,3),particules_after(:,4),'y.');
plot(v_target(3,:),v_target(4,:),'r*')
title('velocities');

% k = 2
% ------------------------------------------------------------------------
figure(2)
subplot(2,1,1)
plot(obs(1,:),obs(2,:),'.'); hold on; % real trajectory of the observer
plot(real_target(1,:),real_target(2,:),'.')% real trajectory of the target

% plot the predicted target's trajectory
target=zeros(4,2);
for t=1:2
    helper=[0 0 0 0]';
    for i=1:length(X)
        helper=helper+X{i,t}; %size(X) = 5000 26
    end
    target(:,t)=helper/length(X);
end
pred_traj = zeros(2,2);
pred_traj(1:2,:) = target(1:2,:) + obs(1:2,1:2);
plot(pred_traj(1,:),pred_traj(2,:),'g.'); % predicted target's trajectory

% the particules will be stored by row
particules_before=zeros(5000,4);
particules_after=zeros(5000,4);
for j=1:length(X)
     xtilde = Xtilde{j,2};
     x = X{j,2};
     particules_before(j,:)= xtilde'+ obs(:,2)';
     particules_after(j,:)= x'+ obs(:,2)';
 end
plot(particules_before(:,1),particules_before(:,2),'m.');
plot(particules_after(:,1),particules_after(:,2),'y.');
title('positions'); 
xlabel('x');
ylabel('y');

subplot(2,1,2)
plot(particules_before(:,3),particules_before(:,4),'m.'); hold on;
plot(particules_after(:,3),particules_after(:,4),'y.');
plot(v_target(3,:),v_target(4,:),'r*')
title('velocities');

% k = 3
% -----------------------------------------------------------------------
figure(3)
subplot(2,1,1)
plot(obs(1,:),obs(2,:),'.'); hold on; % real trajectory of the observer
plot(real_target(1,:),real_target(2,:),'.')% real trajectory of the target

% plot the predicted target's trajectory
target=zeros(4,3);
for t=1:3
    helper=[0 0 0 0]';
    for i=1:length(X)
        helper=helper+X{i,t}; %size(X) = 5000 26
    end
    target(:,t)=helper/length(X);
end
pred_traj = zeros(2,3);
pred_traj(1:2,:) = target(1:2,:) + obs(1:2,1:3);
plot(pred_traj(1,:),pred_traj(2,:),'g.'); % predicted target's trajectory

% the particules will be stored by row
particules_before=zeros(5000,4);
particules_after=zeros(5000,4);
for j=1:length(X)
     xtilde = Xtilde{j,3};
     x = X{j,3};
     particules_before(j,:)= xtilde'+ obs(:,3)';
     particules_after(j,:)= x'+ obs(:,3)';
 end
plot(particules_before(:,1),particules_before(:,2),'m.');
plot(particules_after(:,1),particules_after(:,2),'y.');
title('positions'); 
xlabel('x');
ylabel('y');

subplot(2,1,2)
plot(particules_before(:,3),particules_before(:,4),'m.'); hold on;
plot(particules_after(:,3),particules_after(:,4),'y.');
plot(v_target(3,:),v_target(4,:),'r*')
title('velocities');

% k = 15
% ------------------------------------------------------------------------
figure(4)
subplot(2,1,1)
plot(obs(1,:),obs(2,:),'.'); hold on; % real trajectory of the observer
plot(real_target(1,:),real_target(2,:),'.')% real trajectory of the target

% plot the predicted target's trajectory
target=zeros(4,15);
for t=1:15
    helper=[0 0 0 0]';
    for i=1:length(X)
        helper=helper+X{i,t}; %size(X) = 5000 26
    end
    target(:,t)=helper/length(X);
end
pred_traj = zeros(2,15);
pred_traj(1:2,:) = target(1:2,:) + obs(1:2,1:15);
plot(pred_traj(1,:),pred_traj(2,:),'g.'); % predicted target's trajectory

% the particules will be stored by row
particules_before=zeros(5000,4);
particules_after=zeros(5000,4);
for j=1:length(X)
     xtilde = Xtilde{j,15};
     x = X{j,15};
     particules_before(j,:)= xtilde'+ obs(:,15)';
     particules_after(j,:)= x'+ obs(:,15)';
 end
plot(particules_before(:,1),particules_before(:,2),'m.');
plot(particules_after(:,1),particules_after(:,2),'y.');
title('positions');
xlabel('x');
ylabel('y');

subplot(2,1,2)
plot(particules_before(:,3),particules_before(:,4),'m.'); hold on;
plot(particules_after(:,3),particules_after(:,4),'y.');
plot(v_target(3,:),v_target(4,:),'r*')
title('velocities');

% k = 26
% ------------------------------------------------------------------------
figure(5)
subplot(2,1,1)
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
title('trajectories')
xlabel('x')
ylabel('y')

% the particules will be stored by row
particules_before=zeros(5000,4);
particules_after=zeros(5000,4);
for j=1:length(X)
     xtilde = Xtilde{j,t_f+1};
     x = X{j,t_f+1};
     particules_before(j,:)= xtilde'+ obs(:,t_f+1)';
     particules_after(j,:)= x'+ obs(:,t_f+1)';
 end
plot(particules_before(:,1),particules_before(:,2),'m.');
plot(particules_after(:,1),particules_after(:,2),'y.');
title('positions');
xlabel('x');
ylabel('y');

subplot(2,1,2)
plot(particules_before(:,3),particules_before(:,4),'m.'); hold on;
plot(particules_after(:,3),particules_after(:,4),'y.');
plot(v_target(3,:),v_target(4,:),'r*')
title('velocities');
 end