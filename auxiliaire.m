Trm=zeros(1,1);
Zm=zeros(5000,200);
Xm=zeros(5000,200);
Ym=zeros(5000,200);
for i=1:5000
  
    [Zm(i,:),Xm(i,:),Ym(i,:)]=Q2();
    
end

figure(3)
%histogram(Zm)
histogram(Xm)
%histogram(Ym)
% figure(1)
% plot(mean(Zm))
% 
% figure(2)
% plot(mean(Xm),mean(Ym))


figure(1)
hold on
histogram(Xm(:,200)')
histogram(Xm(:,100)')
histogram(Xm(:,50)')
histogram(Xm(:,1)')
figure(2)
hold on
histogram(Ym(:,200)')
histogram(Ym(:,100)')
histogram(Ym(:,50)')
histogram(Ym(:,1)')
