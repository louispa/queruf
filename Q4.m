p=load('data.mat');
[X_target,Xtilde_target]=Q4fun(p.observer, p.target, p.measurements, p.r, p.theta, p.s, p.c)
draft(p.observer,p.target,X_target,Xtilde_target)