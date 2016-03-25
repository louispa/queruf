function[one]=draft(X,n)
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