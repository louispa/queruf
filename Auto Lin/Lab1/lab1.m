t=250;
h=19.6;
Sr=43;
g=981;
qp=30
Ss30=qp/sqrt(2*g*h)
G=tf([1],[Sr 2*Ss30*sqrt(g/2/h)]);
[gy,tg]=step(G,250);

H=tf([-sqrt(2*g*h)],[Sr Ss30*sqrt(g/2/h)]);
[hy,th]=step(G,250);
close all;
kp=2.5;ki=2;
Tr=tf([kp kp*ki],[Sr 2*Ss30*sqrt(g/2/h)+kp kp*ki]);
%step(Tr,250)

figure()
kp=10;ki=0;
tr100=tf([kp/Sr kp*ki/Sr],[1 2*Ss30*sqrt(g/2/h)/Sr+kp/Sr kp*ki/Sr]);
step(tr100,250);
hold on;
p=reader('pi_10_0_premier.txt');
plot(p(:,1),p(:,3)./22,'r');
legend('Résultat théorique','Résultat expérimentale','Location','southeast');
hold off;
print -djpeg99 10_0.jpeg

figure()
kp=10;ki=0.1;
tr100=tf([kp/Sr kp*ki/Sr],[1 2*Ss30*sqrt(g/2/h)/Sr+kp/Sr kp*ki/Sr]);
step(tr100,250);
hold on;
p=reader('pi_10_0.1_premier.txt');
plot(p(:,1),p(:,3)./19.6,'r');
legend('Résultat théorique','Résultat expérimentale','Location','southeast');
hold off;
print -djpeg99 10_01.jpeg

figure()
kp=3;ki=1;
tr100=tf([kp/Sr kp*ki/Sr],[1 2*Ss30*sqrt(g/2/h)/Sr+kp/Sr kp*ki/Sr]);
step(tr100,250);
hold on;
p=reader('pi_3_1_premier.txt');
plot(p(:,1),p(:,3)./19.6,'r');
legend('Résultat théorique','Résultat expérimentale','Location','southeast');
hold off;
print -djpeg99 3_1.jpeg





% omega=12/t
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

