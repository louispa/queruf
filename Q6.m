%p=load('data.mat');
%[X_target,Xtilde_target]=Q4fun_bis(p.observer, p.target, p.measurements, p.r, p.theta, p.s, p.c);

t_f=26;

estimated_relative=zeros(4,t_f);
for i=1:t_f
    helper=[0 0 0 0]';
    for j=1:length(X_target)
        helper=helper+X_target{j,i}; %size(X) = 5000 26
    end
    estimated_relative(:,i)=helper/length(X_target);
end
[bound]=Q6fun(p.observer, p.target, p.measurements, p.r, p.theta, p.s, p.c,estimated_relative);


true_relative=p.target-p.observer;
estimated_relative_error=zeros(26,1);
for i=1:t_f
    helper=0;
    for j=1:length(X_target)
        helper=helper+ (X_target{j,i}(1)-true_relative(1,i))^2 ...
            +(X_target{j,i}(2)-true_relative(2,i))^2; %size(X) = 5000 26
    end
    estimated_relative_error(i)=sqrt(helper/length(X_target));
end
%graphs
plot(bound);
hold on;
plot(estimated_relative_error);
hold off;
legend('Cramer Rao Lower Bound','RMS postition error');




