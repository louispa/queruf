p=load('data.mat');
[X,Xtilde]=Q5fun(p.observer,p.measurements,p.r,p.theta,p.s,p.c);

% don't forget to put sigma_a = 0 in the Q4fun file
[X1,Xtilde1]=Q4fun(p.observer,p.measurements,p.r,p.theta,p.s,p.c);

draft_5(p.observer,p.target,X,Xtilde,X1,Xtilde1);
