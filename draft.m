function[one]=draft(X,n)
one=zeros(n,2);
for i=1:n
   one(i,:)=[X{i,1}(1) X{i,1}(2)]; 
end
fifty=zeros(n,2);
for i=1:n
   fifty(i,:)=[X{i,50}(1) X{i,50}(2)]; 
end
hundred=zeros(n,2);
for i=1:n
   hundred(i,:)=[X{i,100}(1) X{i,100}(2)]; 
end
twohundred=zeros(n,2);
for i=1:n
   twohundred(i,:)=[X{i,200}(1) X{i,200}(2)]; 
end
figure(1)
hold on;
subplot(2,2,1)
a=histogram(one(:,1));
subplot(2,2,2)
histogram(fifty(:,1));
subplot(2,2,3)
histogram(hundred(:,1));
subplot(2,2,4)
histogram(twohundred(:,1));
title('Histograms for x');
hold off;
figure(2)
hold on;
subplot(2,2,1)
a=histogram(one(:,2));
subplot(2,2,2)
histogram(fifty(:,2));
subplot(2,2,3)
histogram(hundred(:,2));
subplot(2,2,4)
histogram(twohundred(:,2));
title(a,'Histograms for y');
end