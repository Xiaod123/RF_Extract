%function f=deembedding()




% Ltsv=0.03;
Ltsv=300*1e-6;
freq=xlsread('deembedding','8.2.2a','A2:A202');

S11r=xlsread('deembedding','8.2.2a','B2:B202');
S11i=xlsread('deembedding','8.2.2a','C2:C202');

S12r=xlsread('deembedding','8.2.2a','D2:D202');
S12i=xlsread('deembedding','8.2.2a','E2:E202');

S21r=xlsread('deembedding','8.2.2a','F2:F202');
S21i=xlsread('deembedding','8.2.2a','G2:G202');

S22r=xlsread('deembedding','8.2.2a','H2:H202');
S22i=xlsread('deembedding','8.2.2a','I2:I202');

% S11a=S11r+S11i*j;
% S12a=S12r+S12i*j;
% S21a=S21r+S21i*j;
% S22a=S22r+S22i*j;

S11a=10.^(S11r/20).*cos(S11i/180*pi)+10.^(S11r/20).*sin(S11i/180*pi)*j;
S12a=10.^(S12r/20).*cos(S12i/180*pi)+10.^(S12r/20).*sin(S12i/180*pi)*j;
S21a=10.^(S21r/20).*cos(S21i/180*pi)+10.^(S21r/20).*sin(S21i/180*pi)*j;
S22a=10.^(S22r/20).*cos(S22i/180*pi)+10.^(S22r/20).*sin(S22i/180*pi)*j;
% Z11r=xlsread('deembedding','8.2.2b','B2:B10');
% Z11i=xlsread('deembedding','8.2.2b','C2:C10');
% 
% Z12r=xlsread('deembedding','8.2.2b','D2:D10');
% Z12i=xlsread('deembedding','8.2.2b','E2:E10');
% 
% Z21r=xlsread('deembedding','8.2.2b','F2:F10');
% Z21i=xlsread('deembedding','8.2.2b','G2:G10');
% 
% Z22r=xlsread('deembedding','8.2.2b','H2:H10');
% Z22i=xlsread('deembedding','8.2.2b','I2:I10');
% 
% Z11a=Z11r+Z11i*j;
% Z12a=Z12r+Z12i*j;
% Z21a=Z21r+Z21i*j;
% Z22a=Z22r+Z22i*j;

% S11r=xlsread('deembedding','8.2.2b','B2:B10');
% S11i=xlsread('deembedding','8.2.2b','C2:C10');
% 
% S12r=xlsread('deembedding','8.2.2b','D2:D10');
% S12i=xlsread('deembedding','8.2.2b','E2:E10');
% 
% S21r=xlsread('deembedding','8.2.2b','F2:F10');
% S21i=xlsread('deembedding','8.2.2b','G2:G10');
% 
% S22r=xlsread('deembedding','8.2.2b','H2:H10');
% S22i=xlsread('deembedding','8.2.2b','I2:I10');
% 
% S11b=S11r+S11i*j;
% S12b=S12r+S12i*j;
% S21b=S21r+S21i*j;
% S22b=S22r+S22i*j;
% Zcr=xlsread('deembedding','8.2.2','K2:K10');
% Zci=xlsread('deembedding','8.2.2','L2:L10');
% Zc=Zcr+Zci*j;
% gammar=xlsread('deembedding','8.2.2','M2:M10');
% gammai=xlsread('deembedding','8.2.2','N2:N10');
% gamma=gammar+gammai*j;
for (i=1:201)
%     K1=(((S11a(i)^2-S21a(i)^2+1)^2-(2*S11a(i))^2)/(2*S21a(i))^2)^0.5;
%     K2=(((S11b(i)^2-S21b(i)^2+1)^2-(2*S11b(i))^2)/(2*S21b(i))^2)^0.5;
% % %     Key(i)=K;
%     gamma(i)=1/0.0049*(-log((1-S11b(i)^2+S21b(i)^2)/2/S21b(i)+K2)+log((1-S11a(i)^2+S21a(i)^2)/2/S21a(i)+K1));
% %     gamma(i)=abs(real(gammap(i)))+abs(imag(gammap(i)))*j;
%     Zc(i)=50*(((1+S11(i))^2-S21(i)^2)/((1-S11(i))^2-S21(i)^2))^0.5;
    Sa=[S11a(i) S12a(i); S21a(i) S22a(i)];
    T1=s2abcd(Sa, 50);
    
%     Sb=[S11b(i) S12b(i); S21b(i) S22b(i)];
% %     I=[1 0;0 1];
% %     Z=50*(I+Sa)*(I-Sa)^-1;
% %     A=Z(1)*Z(2)^-1;a
% %     B=Z(1)*Z(2)^-1*Z(4)-Z(3);
% %     C=Z(2)^-1;
% %     D=Z(2)^-1*Z(4);
% %     T1=[A B;C D];
% %     Z=50*(I+Sb)*(I-Sb)^-1;
% 
%     
%     T2=s2abcd(Sb, 50);
    
%     X=T1+T2;
%     Y=T1-T2;
%     M1=(Y*X^-1+X*Y^-1)*(-Y*X^-1+X*Y^-1)^-1;
%     M2=(-Y*X^-1+X*Y^-1)^-1;

%     gamma(i)=1/0.002*acosh(M1(4));
% gamma(i)=log(Sb(2)/Sa(2))/0.005
%     Zc(i)=0.5*M2(2)^-1*sinh(gamma(i)*0.002);
    gamma(i)=1/Ltsv*acosh(T1(4));
%     gammac(i)=1/Ltsv*acosh(T1(1));
%     gamma(i)=real(gammap(i))+abs(imag(gammap(i)))*j
    Zc(i)=T1(2)^-1*sinh(gamma(i)*Ltsv);
    
%     Zcc(i)=T1(3)*sinh(gamma(i)*Ltsv)^-1;
    R(i)=real(gamma(i)*Zc(i));
    L(i)=1/2/pi/freq(i)*imag(gamma(i)*Zc(i));
    G(i)=real(gamma(i)/Zc(i));
    C(i)=1/2/pi/freq(i)*imag(gamma(i)/Zc(i));
    losstan(i)=real(gamma(i)/Zc(i))/imag(gamma(i)/Zc(i));
    
    Attenuation(i) = 20*log10(abs(exp(-gamma(i)*Ltsv)));
%     AA=Z11a(i)*Z21a(i)^-1;
%     BB=Z11a(i)*Z21a(i)^-1*Z22a(i)-Z12a(i);
%     CC=Z21a(i)^-1;
%     DD=Z21a(i)^-1*Z22a(i);
%     T2=[AA BB;CC DD];
%        gamma2(i)=1/Ltsv*acosh(T2(4));
% %     gamma(i)=real(gammap(i))+abs(imag(gammap(i)))*j
%     Zc2(i)=T2(2)^-1*sinh(gamma2(i)*Ltsv);
%     R2(i)=real(gamma2(i)*Zc2(i));
%     L2(i)=1/2/pi/freq(i)/1e9*imag(gamma2(i)*Zc2(i));
%     G2(i)=real(gamma2(i)/Zc2(i));
%     C2(i)=1/2/pi/freq(i)/1e9*imag(gamma2(i)/Zc2(i));
%     losstan2(i)=real(gamma2(i)/Zc2(i))/imag(gamma2(i)/Zc2(i));
%     result(i)=s2rlgc(S,0.005,freq(i),50);
%     R(i)=result(i).R
%     L(i)=result(i).L
%     G(i)=result(i).G
%     C(i)=result(i).C
%    
end
% Sa
xlswrite('Verfication_TVinter',freq,'8.2.2','A2:A202');
xlswrite('Verfication_TVinter',R','8.2.2','B2:B202');
xlswrite('Verfication_TVinter',L','8.2.2','C2:C202');
xlswrite('Verfication_TVinter',G','8.2.2','D2:D202');
xlswrite('Verfication_TVinter',C','8.2.2','E2:E202');
xlswrite('Verfication_TVinter',real(gamma)','8.2.2','I2:I202');
xlswrite('Verfication_TVinter',imag(gamma)','8.2.2','J2:J202');
xlswrite('Verfication_TVinter',real(Zc)','8.2.2','G2:G202');
xlswrite('Verfication_TVinter',imag(Zc)','8.2.2','H2:H202');
% xlswrite('Verfication_TVinter',real(gammac)','8.2.2','N2:N247');
% xlswrite('Verfication_TVinter',imag(gammac)','8.2.2','O2:O247');
% xlswrite('Verfication_TVinter',real(Zcc)','8.2.2','L2:L247');
% xlswrite('Verfication_TVinter',imag(Zcc)','8.2.2','M2:M247');
% % 
% xlswrite('Verfication_TVinter',Attenuation','8.2.2','L2:L247');
% xlswrite('Verfication_TVinter',20*log10(abs(S21a)),'8.2.2','M2:M247');


% xlswrite('Verfication_TVinter',freq,'8.2.2b','A2:A10');
% xlswrite('Verfication_TVinter',R2','8.2.2b','B2:B10');
% xlswrite('Verfication_TVinter',L2','8.2.2b','C2:C10');
% xlswrite('Verfication_TVinter',G2','8.2.2b','D2:D10');
% xlswrite('Verfication_TVinter',C2','8.2.2b','E2:E10');
% xlswrite('Verfication_TVinter',real(gamma2)','8.2.2b','I2:I10');
% xlswrite('Verfication_TVinter',imag(gamma2)','8.2.2b','J2:J10');
% xlswrite('Verfication_TVinter',real(Zc2)','8.2.2b','G2:G10');
% xlswrite('Verfication_TVinter',imag(Zc2)','8.2.2b','H2:H10');

% xlswrite('Verfication_TVinter',20*log10(abs(S21a)),'8.2.2','L2:L10');
%end
