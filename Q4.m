p=load('data.mat');
[X,Xtilde]=Q4fun(p.observer,p.measurements,p.r,p.theta,p.s,p.c);
draft_4(p.observer,p.target,X,Xtilde);