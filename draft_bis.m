function[]=draft_bis(obs,real_target,X,Xtilde)
%%

t_f = 26;

%At the first iteration
% plot the true trajectory of the observer and the target
figure(1)
plot(obs(1,:),obs(2,:),'.'); hold on; % real trajectory of the observer
plot(real_target(1,:),real_target(2,:),'.')% real trajectory of the target

% plot the predicted target's trajectory
target=zeros(4,t_f);
for t=1:t_f
    helper=[0 0 0 0]';
    for i=1:length(X)
        helper=helper+X{i,t}; %size(X) = 5000 26
    end
    target(:,t)=helper/length(X);
end
pred_traj = zeros(2,t_f);
pred_traj(1:2,:) = target(1:2,:) + obs;
plot(pred_traj(1,:),pred_traj(2,:),'.'); hold off % predicted target's trajectory


% particules_before=zeros(5000,4);
% particules_after=zeros(5000,4);
% for j=1:length(X)
%     particules_before(j,:)=(Xtilde{j,1}(:))';
%     particules_after(j,:)=(X{j,1}(:))';
% end
% %plot(particules_before(:,1),particules_before(:,2),'m*');
% %plot(particules_before(:,1),particules_before(:,2),'g.');
% hold off;
% 
% subplot(2,1,2)
% hold on;
% plot(particules_before(:,3),particules_before(:,4),'b.');
% plot(particules_after(:,3),particules_after(:,4),'y.');
% plot(target(:,3),target(:,4),'r.')
% hold off;
% 
% % %%
% % figure(2)
% % subplot(2,1,1)
% % hold on;
% % plot(obs(1,:),obs(2,:),'.') %trajectory of the observer
% % 
% % plot(real_target(1,:),real_target(2,:),'.')%trajectory of the target
% % 
% % target=zeros(2,4);
% % for i=1:2%need changes
% %     helper=[0 0 0 0]';
% %     for j=1:length(X)
% %         helper=helper+X{j,i}(:); %size(X) = 5000 26
% %     end
% %     target(i,:)=helper/length(X);
% % end
% % 
% % plot(target(:,1),target(:,2),'b.')
% % 
% % 
% % particules_before=zeros(5000,4);
% % particules_after=zeros(5000,4);
% % for j=1:length(X)
% %     particules_before(j,:)=(Xtilde{j,2}(:))';
% %     particules_after(j,:)=(X{j,2}(:))';
% % end
% % plot(particules_before(:,1),particules_before(:,2),'m*');
% % plot(particules_before(:,1),particules_before(:,2),'g.');
% % hold off;
% % 
% % subplot(2,1,2)
% % hold on;
% % plot(particules_before(:,3),particules_before(:,4),'b.');
% % plot(particules_after(:,3),particules_after(:,4),'y.');
% % plot(target(:,3),target(:,4),'r.')
% % hold off;
% % 
% % %%
% % figure(3)
% % subplot(2,1,1)
% % hold on;
% % plot(obs(1,:),obs(2,:),'.') %trajectory of the observer
% % 
% % plot(real_target(1,:),real_target(2,:),'.')%trajectory of the target
% % 
% % target=zeros(3,4);
% % for i=1:3%need changes
% %     helper=[0 0 0 0]';
% %     for j=1:length(X)
% %         helper=helper+X{j,i}(:); %size(X) = 5000 26
% %     end
% %     target(i,:)=helper/length(X);
% % end
% % 
% % plot(target(:,1),target(:,2),'.')
% % 
% % 
% % particules_before=zeros(5000,4);
% % particules_after=zeros(5000,4);
% % for j=1:length(X)
% %     particules_before(j,:)=(Xtilde{j,3}(:))';
% %     particules_after(j,:)=(X{j,3}(:))';
% % end
% % plot(particules_before(:,1),particules_before(:,2),'m*');
% % plot(particules_before(:,1),particules_before(:,2),'g.');
% % hold off;
% % 
% % subplot(2,1,2)
% % hold on;
% % plot(particules_before(:,3),particules_before(:,4),'b.');
% % plot(particules_after(:,3),particules_after(:,4),'y.');
% % plot(target(:,3),target(:,4),'r.')
% % hold off;
% % 
% % %%
% % figure(4)
% % subplot(2,1,1)
% % hold on;
% % plot(obs(1,:),obs(2,:),'.') %trajectory of the observer
% % 
% % plot(real_target(1,:),real_target(2,:),'.')%trajectory of the target
% % 
% % target=zeros(15,4);
% % for i=1:15%need changes
% %     helper=[0 0 0 0]';
% %     for j=1:length(X)
% %         helper=helper+X{j,i}(:); %size(X) = 5000 26
% %     end
% %     target(i,:)=helper/length(X);
% % end
% % 
% % plot(target(:,1),target(:,2),'.')
% % 
% % 
% % particules_before=zeros(5000,4);
% % particules_after=zeros(5000,4);
% % for j=1:length(X)
% %     particules_before(j,:)=(Xtilde{j,15}(:))';
% %     particules_after(j,:)=(X{j,15}(:))';
% % end
% % plot(particules_before(:,1),particules_before(:,2),'m*');
% % plot(particules_before(:,1),particules_before(:,2),'g.');
% % hold off;
% % 
% % subplot(2,1,2)
% % hold on;
% % plot(particules_before(:,3),particules_before(:,4),'b.');
% % plot(particules_after(:,3),particules_after(:,4),'y.');
% % plot(target(:,3),target(:,4),'r.')
% % hold off;
% 
% %%
% figure(5)
% subplot(2,1,1)
% hold on;
% plot(obs(1,:),obs(2,:),'.') %trajectory of the observer
% 
% plot(real_target(1,:),real_target(2,:),'.')%trajectory of the target
% 
% target=zeros(26,4);
% for i=1:26%need changes
%     helper=[0 0 0 0]';
%     for j=1:length(X)
%         helper=helper+X{j,i}(:); %size(X) = 5000 26
%     end
%     target(i,:)=helper/length(X);
% end
% target
% plot(target(:,1),target(:,2),'.')
% 
% 
% particules_before=zeros(5000,4);
% particules_after=zeros(5000,4);
% for j=1:length(X)
%     particules_before(j,:)=(Xtilde{j,26}(:))';
%     particules_after(j,:)=(X{j,26}(:))';
% end
% plot(particules_before(:,1),particules_before(:,2),'m*');
% plot(particules_before(:,1),particules_before(:,2),'g.');
% hold off;
% 
% subplot(2,1,2)
% hold on;
% plot(particules_before(:,3),particules_before(:,4),'b.');
% plot(particules_after(:,3),particules_after(:,4),'y.');
% plot(target(:,3),target(:,4),'r.')
% hold off;
% end
% 
% 
% 
% % function[]=lol()
% % % function[one]=draft(X,n)
% % % %Draw the histograms at time t=1,50,100 and 200 for x and y
% % % %Draw the trajectory of the target 
% % % one=zeros(n,2);
% % % %our first state is in t=0
% % % for i=1:n
% % %    one(i,:)=[X{i,2}(1) X{i,2}(2)]; 
% % % end
% % % fifty=zeros(n,2);
% % % for i=1:n
% % %    fifty(i,:)=[X{i,51}(1) X{i,51}(2)]; 
% % % end
% % % hundred=zeros(n,2);
% % % for i=1:n
% % %    hundred(i,:)=[X{i,101}(1) X{i,101}(2)]; 
% % % end
% % % twohundred=zeros(n,2);
% % % for i=1:n
% % %    twohundred(i,:)=[X{i,201}(1) X{i,201}(2)]; 
% % % end
% % % figure(1)
% % % hold on;
% % % subplot(2,2,1)
% % % histogram(one(:,1));
% % % title('Histograms for x at t=1');
% % % subplot(2,2,2)
% % % histogram(fifty(:,1));
% % % title('Histograms for x at t=50');
% % % subplot(2,2,3)
% % % histogram(hundred(:,1));
% % % title('Histograms for x at t=100');
% % % subplot(2,2,4)
% % % histogram(twohundred(:,1));
% % % title('Histograms for x at t=200');
% % % 
% % % hold off;
% % % figure(2)
% % % hold on;
% % % subplot(2,2,1)
% % % histogram(one(:,2));
% % % title('Histograms for y at t=1');
% % % subplot(2,2,2)
% % % histogram(fifty(:,2));
% % % title('Histograms for y at t=50');
% % % subplot(2,2,3)
% % % histogram(hundred(:,2));
% % % title('Histograms for y at t=100');
% % % subplot(2,2,4)
% % % histogram(twohundred(:,2));
% % % title('Histograms for y at t=200');
% % % 
% % % traj=zeros(200,2);
% % % for t=2:length(X(1,:))
% % %     helper=0;
% % %    for i=1:n
% % %        helper=helper+X{i,t}([1 2]);
% % %    end
% % %     traj(i,:)=helper./n;
% % % end
% % % 
% % % figure(3)
% % % plot(traj(:,1),traj(:,2));
% % % title('Trajectory of the target');
% % % 
% % % 
% % % 
% % % 
% % % end
 end