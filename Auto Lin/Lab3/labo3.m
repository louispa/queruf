%Param�tres
v1=4.7;
v2=4.7;
rp=0.47;
c1=2.2;
c2=2.2;
r12=4.7
%Param�tre du syst�me lin�aris�
a11=1/c1*(1/rp+1/r12)
a12=1/c1/r12
a21=1/c1/r12
a22=1/c2/r12
b=1/5/c1
d=v1/c1/rp^2


%A minimum de phase
Ga=zpk(tf([10*b*a21],[1 (a11+a22) a11*a22-a12*a21]))
Ha=zpk(tf([10*d*a21],[1 (a11+a22) a11*a22-a12*a21]))
% figure(1)
% step(Ga);hold on;
% m=reader('linear_2')
% plot(m(:,1),(m(:,4)-40)./10,'r');hold off;
tau=1/0.087
tr=tau*log(50)


%B non minimum de phase
Gb=zpk(tf(-10*b.*[1 -2*a21+a22],[1 (a11+a22) a11*a22-a12*a21]))
step(Gb);
% hold on;
% mb=reader('linear_2');
% plot(mb(:,1),(mb(:,4)-40)./10,'r');hold off;
Hb=zpk(tf(-10*d.*[1 -2*a21+a22],[1 (a11+a22) a11*a22-a12*a21]))



%Tv
%minimum phase
alpha=0.536;
pb=30.68;
ti=11.5;
tva=zpk(tf([10*d*a21 0],[1 a11+a22...
    a11*a22-a21*a12+10^3/pb*b*a21 10^3/pb/ti*b*a21]));
figure(2)
step(tva)
hold on;
mm=reader('minimum_phase_1');
plot(mm(:,1),mm(:,4)./40,'r');
hold off;

% %non minimum phase
% tvb=zpk(tf(-10*d.*[1 -2*a21+a22 0],[1 a11+a22-10^3/pb*b...
%     a11*a22-a21*a12+10^3/pb*b*(2*a21-a22-1/ti)... 
%     10^3/pb/ti*b*(2*a21-a22)]))
% 
% figure(3)
% step(tvb);hold on;
% mnm=reader('non_minimum_phase_equilibrum_1_complet');
% plot(mnm(:,1),mnm(:,4),'r');hold off;
% %Tr b non minimum phase
% pb=156;
% trb=zpk(tf(1/pb*[-91 8.23],[1 1.072-91/pb 8.23/pb]))
% tvb=zpk(tf(-10*d.*[1 -2*a21+a22 0],[1 a11+a22-10^3/pb*b...
%     a11*a22-a21*a12+10^3/pb*b*(2*a21-a22-1/ti)... 
%     10^3/pb/ti*b*(2*a21-a22)]))
% %step(trb)


