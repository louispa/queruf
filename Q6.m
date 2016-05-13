p=load('data.mat');
[X_target,Xtilde_target]=Q4fun(p.observer,p.measurements, p.r, p.theta, p.s, p.c);

t_f=26;

% computation of the realtive trajectory. We work with the mean of the
% sample
estimated_relative=zeros(4,t_f);
for i=1:t_f
    helper=[0 0 0 0]';
    for j=1:length(X_target)
        helper=helper+X_target{j,i}; 
    end
    estimated_relative(:,i)=helper/length(X_target);
end

% computation of the CRLB(RMS) bound
[bound]=Q6fun(p.r, p.theta, p.s, p.c,estimated_relative);

% plot of the estimated relative trajectory
figure(1)
plot(estimated_relative(1,:),estimated_relative(2,:),'.-b');
title('relative trajectory')

true_relative=p.target-p.observer; % true relative trajectory

% computation of the errors between the estimated realtive trajectory and
% the true relative trajectory in order to compute the RMS. Again, we work
% with the mean of the errors on every sample
estimated_relative_error=zeros(26,1);
for i=1:t_f
    helper=0;
    for j=1:length(X_target)
        helper=helper+ (X_target{j,i}(1)-true_relative(1,i))^2 ...
            +(X_target{j,i}(2)-true_relative(2,i))^2; 
    end
    estimated_relative_error(i)=sqrt(helper/length(X_target));
end

% graphs
figure(2)
plot(bound,'-ob');
hold on;
plot(estimated_relative_error,'-or');
hold off;
title('Evolution of the CRLMB of the RMS error and the RMS error in function of the time');
legend('Cramer Rao Lower Bound','RMS postition error');
xlabel('iterations');
ylabel('error');




