function f=deembedding_HO()




% Ltsv=0.03;
% Ltsv=300*1e-6;
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

    Sa=[S11a(i) S12a(i); S21a(i) S22a(i)];
    ZZ=s2z(Sa, 50);
    YY=s2y(Sa, 50);
    Zdiff = 50*(ZZ(1)-ZZ(2)-ZZ(3)+ZZ(4));
    Ycomm = YY(1)+YY(2)+YY(3)+YY(4);
    R(i)=real(Zdiff);
    L(i)=1/2/pi/freq(i)*imag(Zdiff);
    G(i)=real(Ycomm);
    C(i)=1/2/pi/freq(i)*imag(Ycomm);
    
end
% Sa
xlswrite('Verfication_TVinter',freq,'8.2.2','A2:A202');
xlswrite('Verfication_TVinter',R','8.2.2','B2:B202');
xlswrite('Verfication_TVinter',L','8.2.2','C2:C202');
xlswrite('Verfication_TVinter',G','8.2.2','D2:D202');
xlswrite('Verfication_TVinter',C','8.2.2','E2:E202');

end
