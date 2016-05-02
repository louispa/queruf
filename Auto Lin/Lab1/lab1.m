%%
% données
t=250;
%h=19.6;
h = 28.8;
Sr=43;
g=981;
qp=30;
Ss30=qp/sqrt(2*g*h);
kp=2.5;
ki=2;

% fonction de transfert G
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%G=tf([1],[Sr Ss30*sqrt(g/2/h)])
G=(1/Sr)*tf(1,[1 (Ss30/Sr)*sqrt((g/(2*h)))]);
[gy,tg]=step(G,250);
plot(tg,gy);
hold on;
% importation des données
loop=reader1('test1.txt');
plot(loop(:,1),loop(:,3)./h,'r');
%plot(loop(:,1),loop(:,3)./19.6,'r');
title('fonction de transfert de la commande (G)')
legend('Résultat théorique','Résultat expérimental','Location','southeast');
hold off;
%print -djpeg99 10_0.jpeg

%%
% fonction de transfert H
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
H=tf(-sqrt(2*g*h),[Sr Ss30*sqrt(g/2/h)]);
[hy,th]=step(H,250);
figure(2)
plot(th,hy);
title('fonction de transfert de la perturubation (H)')

% fonction de transfert Tr (closed loop)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C = kp*tf([1 ki],[1 0]);
Tr = feedback(C*G,1);
%Tr=tf([kp kp*ki],[Sr 2*Ss30*sqrt(g/2/h)+kp kp*ki]);
[ty,tt] = step(Tr,250);
Tv = H*feedback(1,C*G);
[vy,vt] = step(Tr,250);
figure(3)
plot(tt,ty); hold on
title('fonction de transfert de la commande en boucle fermée (Tr) avec kp = 2.5 et ki = 2')
p = reader1('test2.txt');
plot(p(:,1),p(:,3)./h,'r');
figure(4)
plot(vt,vy);
title('fonction de transfert de la perturbation en boucle fermée (Tv) avec kp = 2.5 et ki = 2')
 

kp=10;
ki=0;
C = kp*tf([1 ki],[1 0]);
Tr = feedback(C*G,1);
%Tr=tf([kp kp*ki],[Sr 2*Ss30*sqrt(g/2/h)+kp kp*ki]);
[ty,tt] = step(Tr,250);
figure(5)
plot(tt,ty);
title('fonction de transfert de la commande en boucle fermée (Tr) avec kp = 10 et ki = 0')
% tr100=tf([kp/Sr kp*ki/Sr],[1 2*Ss30*sqrt(g/2/h)/Sr+kp/Sr kp*ki/Sr]);
% step(tr100,250);
hold on;
p=reader1('test3.txt');
plot(p(:,1),p(:,3)./h,'r');
legend('Résultat théorique','Résultat expérimental','Location','southeast');
hold off;
% print -djpeg99 10_0.jpeg
% 

kp=3;
ki=1;
C = kp*tf([1 ki],[1 0]);
Tr = feedback(C*G,1);
%Tr=tf([kp kp*ki],[Sr 2*Ss30*sqrt(g/2/h)+kp kp*ki]);
[ty,tt] = step(Tr,250);
figure(6)
plot(tt,ty);
title('fonction de transfert de la commande en boucle fermée (Tr) avec kp = 3 et ki = 1')
hold on;
p=reader1('test4.txt');
plot(p(:,1),p(:,3)./h,'r');
legend('Résultat théorique','Résultat expérimental','Location','southeast');
hold off;
% tr100=tf([kp/Sr kp*ki/Sr],[1 2*Ss30*sqrt(g/2/h)/Sr+kp/Sr kp*ki/Sr]);
% step(tr100,250);
% print -djpeg99 3_1.jpeg

kp=10;
ki=0.1;
C = kp*tf([1 ki],[1 0]);
Tr = feedback(C*G,1);
%Tr=tf([kp kp*ki],[Sr 2*Ss30*sqrt(g/2/h)+kp kp*ki]);
[ty,tt] = step(Tr,250);
figure(7)
plot(tt,ty);
title('fonction de transfert de la commande en boucle fermée (Tr) avec kp = 10 et ki = 0.1')
hold on;
p=reader1('test5.txt');
plot(p(:,1),p(:,3)./h,'r');
legend('Résultat théorique','Résultat expérimental','Location','southeast');
hold off;
% print -djpeg99 10_01.jpeg
% 


% Kp et Ki optimaux
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

ki = (Ss30/Sr)*sqrt(g/(2*h));
kp = -3*Ss30*sqrt(g/(2*h));
C = kp*tf([1 ki],[1 0]);
Tr = feedback(C*G,1);
%Tr=tf([kp kp*ki],[Sr 2*Ss30*sqrt(g/2/h)+kp kp*ki]);
[ty,tt] = step(Tr,250);
figure(8)
plot(tt,ty);
title('Tr avec kp et ki optimaux versus G')
hold on;
plot(tg,gy);
hold off;


figure(9)
plot(tt,ty);
title('Tr avec kp et ki optimaux')
hold on;
p=reader1('kp_ki_opti_1.txt');
plot(p(:,1),p(:,3)./h,'r');
legend('Résultat théorique','Résultat expérimental','Location','southeast');
hold off;

%omega=12/t
% kp=(2*omega-2*Ss30/Sr*sqrt(g/2/h))*Sr
% ki=omega^2*Sr/kp
% Tropti=tf([kp/Sr kp*ki/Sr],[1 2*Ss30*sqrt(g/2/h)/Sr+kp/Sr kp*ki/Sr]);
% figure(1)
% step(Tropti,250);
% hold on;
% p=reader('kp_ki_opti_1.txt');
% plot(p(:,1),p(:,3)./22,'r');
% legend('Résultat théorique','Résultat expérimentale','Location','southeast');
% hold off;
% %print -djpeg99 bon_kc_ki.jpeg

