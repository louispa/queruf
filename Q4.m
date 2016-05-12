p=load('data.mat');
% computation of X and Xtilde
[X,Xtilde]=Q4fun(p.observer,p.measurements,p.r,p.theta,p.s,p.c);
% create the plots asked
draft_4(p.observer,p.target,X,Xtilde);