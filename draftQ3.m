% draft
function[]=draftQ3(X,n,t_f)
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

% histogram at k = 1
figure(1)
subplot(1,2,1)
histogram(one(:,1)); hold on;
title('Histograms for x at t=1');
subplot(1,2,2)
histogram(one(:,2));
title('Histograms for y at t=1');

% histogram at k = 50
figure(2)
subplot(1,2,1)
histogram(fifty(:,1)); hold on;
title('Histograms for x at t=50');
subplot(1,2,2)
histogram(fifty(:,2));
title('Histograms for y at t=50');

% histogram at k = 100
figure(3)
subplot(1,2,1)
histogram(hundred(:,1)); hold on;
title('Histograms for x at t=100');
subplot(1,2,2)
histogram(hundred(:,2));
title('Histograms for y at t=100');

% histogram at k = 200
figure(4)
subplot(1,2,1)
histogram(twohundred(:,1)); hold on;
title('Histograms for x at t=200');
subplot(1,2,2)
histogram(twohundred(:,2));
title('Histograms for y at t=200');

% the realtive trajectory is computed by taking the mean of the relative
% trajectory on the sample. We thus work with X
traj=zeros(2,200);
for t=1:t_f
   helper=[0 0]';
   for i=1:n
       % sum of the x and y of all particles in the sample
       helper = helper + X{i,t}(1:2);
   end
   % we take the mean
    traj(:,t) = helper./n;
end

% plot of the relative trajectory
figure(5)
plot(traj(1,:),traj(2,:),'.');
title('Relative trajectory');
end