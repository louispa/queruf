%Paramètres
v1=4.7;
v2=4.7;
rp=0.47;
c1=2.2;
c2=2.2;
r12=4.7;
%Paramètre du système linearise
a11=1/c1*(1/rp+1/r12);
a12=1/c1/r12;
a21=1/c1/r12;
a22=1/c2/r12;
b=1/5/c1;
d=v1/c1/rp^2;


%%%%%%%%%%%%%%%%%%%%%%%
% OPEN LOOP           %
%%%%%%%%%%%%%%%%%%%%%%%


%A minimum de phase
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ga=zpk(tf(10*b*a21,[1 (a11+a22) a11*a22-a12*a21]));
Ha=zpk(tf(10*d*a21,[1 (a11+a22) a11*a22-a12*a21]));
figure(1)
[ga, tga] = step(Ga);
plot(tga,ga); hold on;
m = reader3('Exp_5_3_1.txt')
plot(m(:,1),(m(:,4)-40)./10,'r'); hold off;
title('réponse indicielle de Ga');

figure(2)
[va,tva] = step(49-200*Ha);
plot(tva,va); hold on;
m = reader3('Exp_5_3_1.txt')
plot(m(:,1),(m(:,4)-40)./10,'r'); hold off;
title('réponse indicielle de Ha');

tau=1/0.087;
tr=tau*log(50);


%B non minimum de phase
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Gb=zpk(tf(-10*b.*[1 -2*a21+a22],[1 (a11+a22) a11*a22-a12*a21]));
Hb=zpk(tf(-10*d.*[1 -2*a21+a22],[1 (a11+a22) a11*a22-a12*a21]));
[g,tg] = step(Gb);
figure(3)
plot(tg,g); hold on;
mb=reader3('linear_2');
plot(mb(:,1),(mb(:,4)-40)./10,'r'); hold off;
title('réponse indicielle de Gb (sytème à non-minimum de phase)');

figure(4)
plot(tga,ga); hold on;
plot(tg,g); hold off;
title('comparaison de Ga (min. phase) et Gb (non min. phase)')

figure(5)
step(-Hb); hold on;
%[gb,tgb] = step(49-200*Hb); hold on;
%plot(tgb,gb); hold on;
mb=reader3('Exp_5_3_2.txt');
plot(mb(9:end,1),(mb(9:end,4)),'r'); hold off;
title('réponse indicielle de Hb (sytème à non-minimum de phase)');

%%%%%%%%%%%%%%%%%%%%%
% CLOSED LOOP       %
%%%%%%%%%%%%%%%%%%%%%

alpha=0.536;
pb=30.84;
ti=11.48;

%Tv
%minimum phase
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tva=tf(10*d*a21,[1 0.087])*tf(1, [1 2*alpha alpha^2]);
figure(6)
[tva,ttva] = step(tva);
plot(ttva,tva); hold on;
mm=reader3('close_loop_minimum_phase_1');
plot(mm(:,1),mm(:,4)./40,'r');
hold off;
title('Tva en boucle fermée')

%non minimum phase
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tvb=zpk(tf(-10*d.*[1 -2*a21+a22 0],[1 a11+a22-10^3/pb*b...
    a11*a22-a21*a12+10^3/pb*b*(2*a21-a22-1/ti)... 
    10^3/pb/ti*b*(2*a21-a22)]))
figure(7)
[tvb,ttvb]=step(tvb);
plot(ttvb,tvb);hold on;
mnm=reader3('non_minimum_phase_equilibrum_1_complet');
plot(mnm(:,1),mnm(:,4),'r');hold off;
title('Tvb en boucle fermée (non min de phase')

%Tr b non minimum phase
pb=156;
trb=zpk(tf(1/pb*[-91 8.23],[1 1.072-91/pb 8.23/pb]))
tvb=zpk(tf(-10*d.*[1 -2*a21+a22 0],[1 a11+a22-10^3/pb*b...
    a11*a22-a21*a12+10^3/pb*b*(2*a21-a22-1/ti)... 
    10^3/pb/ti*b*(2*a21-a22)]));
[trbb,ttrbb] = step(trb);
[tvbb,ttvbb] = step(tvb);
figure(8)
plot(ttrbb,trbb);
title('step response de Trb avec Pb = 156')
figure(9)
plot(ttvbb,tvbb);
title('step response de Tvb avec Pb = 156')


