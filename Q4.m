p=load('data.mat');
[X,Xtilde]=Q4fun_bis(p.observer, p.target, p.measurements, p.r, p.theta, p.s, p.c);
draft_bis(p.observer,p.target,X,Xtilde)