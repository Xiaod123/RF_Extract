function f=twostepdeembedding()




% Ltsv=0.03;
% Ltsv=1000*1e-6;
freq=xlsread('2stepdeembedding','2Lpad','A2:A202');


% %%%% 2L whole %%%%%
% S11r=xlsread('2stepdeembedding','2Lwhole','B2:B251');
% S11i=xlsread('2stepdeembedding','2Lwhole','C2:C251');
% 
% S12r=xlsread('2stepdeembedding','2Lwhole','D2:D251');
% S12i=xlsread('2stepdeembedding','2Lwhole','E2:E251');
% 
% S21r=xlsread('2stepdeembedding','2Lwhole','F2:F251');
% S21i=xlsread('2stepdeembedding','2Lwhole','G2:G251');
% 
% S22r=xlsread('2stepdeembedding','2Lwhole','H2:H251');
% S22i=xlsread('2stepdeembedding','2Lwhole','I2:I251');
% 
% S11a=S11r+S11i*j;
% S12a=S12r+S12i*j;
% S21a=S21r+S21i*j;
% S22a=S22r+S22i*j;
% 
% %%%%% L whole %%%%%
% S11r=xlsread('2stepdeembedding','Lwhole','B2:B251');
% S11i=xlsread('2stepdeembedding','Lwhole','C2:C251');
% 
% S12r=xlsread('2stepdeembedding','Lwhole','D2:D251');
% S12i=xlsread('2stepdeembedding','Lwhole','E2:E251');
% 
% S21r=xlsread('2stepdeembedding','Lwhole','F2:F251');
% S21i=xlsread('2stepdeembedding','Lwhole','G2:G251');
% 
% S22r=xlsread('2stepdeembedding','Lwhole','H2:H251');
% S22i=xlsread('2stepdeembedding','Lwhole','I2:I251');
% 
% S11b=S11r+S11i*j;
% S12b=S12r+S12i*j;
% S21b=S21r+S21i*j;
% S22b=S22r+S22i*j;



%%%%% 2L Pad %%%%%
S11r=xlsread('2stepdeembedding','2Lpad','B2:B202');
S11i=xlsread('2stepdeembedding','2Lpad','C2:C202');

S12r=xlsread('2stepdeembedding','2Lpad','D2:D202');
S12i=xlsread('2stepdeembedding','2Lpad','E2:E202');

S21r=xlsread('2stepdeembedding','2Lpad','F2:F202');
S21i=xlsread('2stepdeembedding','2Lpad','G2:G202');

S22r=xlsread('2stepdeembedding','2Lpad','H2:H202');
S22i=xlsread('2stepdeembedding','2Lpad','I2:I202');

S11c=10.^(S11r/20).*cos(S11i/180*pi)+10.^(S11r/20).*sin(S11i/180*pi)*j;
S12c=10.^(S12r/20).*cos(S12i/180*pi)+10.^(S12r/20).*sin(S12i/180*pi)*j;
S21c=10.^(S21r/20).*cos(S21i/180*pi)+10.^(S21r/20).*sin(S21i/180*pi)*j;
S22c=10.^(S22r/20).*cos(S22i/180*pi)+10.^(S22r/20).*sin(S22i/180*pi)*j;

%%%%% L Pad %%%%%
S11r=xlsread('2stepdeembedding','Lpad','B2:B202');
S11i=xlsread('2stepdeembedding','Lpad','C2:C202');

S12r=xlsread('2stepdeembedding','Lpad','D2:D202');
S12i=xlsread('2stepdeembedding','Lpad','E2:E202');

S21r=xlsread('2stepdeembedding','Lpad','F2:F202');
S21i=xlsread('2stepdeembedding','Lpad','G2:G202');

S22r=xlsread('2stepdeembedding','Lpad','H2:H202');
S22i=xlsread('2stepdeembedding','Lpad','I2:I202');

S11d=10.^(S11r/20).*cos(S11i/180*pi)+10.^(S11r/20).*sin(S11i/180*pi)*j;
S12d=10.^(S12r/20).*cos(S12i/180*pi)+10.^(S12r/20).*sin(S12i/180*pi)*j;
S21d=10.^(S21r/20).*cos(S21i/180*pi)+10.^(S21r/20).*sin(S21i/180*pi)*j;
S22d=10.^(S22r/20).*cos(S22i/180*pi)+10.^(S22r/20).*sin(S22i/180*pi)*j;
%%%%% L TSV %%%%%
S11r=xlsread('2stepdeembedding','LTSV','B2:B202');
S11i=xlsread('2stepdeembedding','LTSV','C2:C202');

S12r=xlsread('2stepdeembedding','LTSV','D2:D202');
S12i=xlsread('2stepdeembedding','LTSV','E2:E202');

S21r=xlsread('2stepdeembedding','LTSV','F2:F202');
S21i=xlsread('2stepdeembedding','LTSV','G2:G202');

S22r=xlsread('2stepdeembedding','LTSV','H2:H202');
S22i=xlsread('2stepdeembedding','LTSV','I2:I202');

S11e=10.^(S11r/20).*cos(S11i/180*pi)+10.^(S11r/20).*sin(S11i/180*pi)*j;
S12e=10.^(S12r/20).*cos(S12i/180*pi)+10.^(S12r/20).*sin(S12i/180*pi)*j;
S21e=10.^(S21r/20).*cos(S21i/180*pi)+10.^(S21r/20).*sin(S21i/180*pi)*j;
S22e=10.^(S22r/20).*cos(S22i/180*pi)+10.^(S22r/20).*sin(S22i/180*pi)*j;

%%%%% 2L TSV %%%%%
S11r=xlsread('2stepdeembedding','2LTSV','B2:B202');
S11i=xlsread('2stepdeembedding','2LTSV','C2:C202');

S12r=xlsread('2stepdeembedding','2LTSV','D2:D202');
S12i=xlsread('2stepdeembedding','2LTSV','E2:E202');

S21r=xlsread('2stepdeembedding','2LTSV','F2:F202');
S21i=xlsread('2stepdeembedding','2LTSV','G2:G202');

S22r=xlsread('2stepdeembedding','2LTSV','H2:H202');
S22i=xlsread('2stepdeembedding','2LTSV','I2:I202');

S11f=10.^(S11r/20).*cos(S11i/180*pi)+10.^(S11r/20).*sin(S11i/180*pi)*j;
S12f=10.^(S12r/20).*cos(S12i/180*pi)+10.^(S12r/20).*sin(S12i/180*pi)*j;
S21f=10.^(S21r/20).*cos(S21i/180*pi)+10.^(S21r/20).*sin(S21i/180*pi)*j;
S22f=10.^(S22r/20).*cos(S22i/180*pi)+10.^(S22r/20).*sin(S22i/180*pi)*j;

for (i=1:201)

%     Sa=[S11a(i) S12a(i); S21a(i) S22a(i)];
%     T2LW=s2abcd(Sa, 50);
%     Sb=[S11b(i) S12b(i); S21b(i) S22b(i)];
%     TLW=s2abcd(Sb, 50);
    Sc=[S11c(i) S12c(i); S21c(i) S22c(i)];
    T2LP=s2abcd(Sc, 50);
    Sd=[S11d(i) S12d(i); S21d(i) S22d(i)];
    TLP=s2abcd(Sd, 50);
    Se=[S11e(i) S12e(i); S21e(i) S22e(i)];
    TLV=s2abcd(Se, 50);
    Sf=[S11f(i) S12f(i); S21f(i) S22f(i)];
    T2LV=s2abcd(Sf, 50);
    
%     TP1sq = ((TLP)^-1*T2LP*(TLP)^-1)^-1
%     [V1sq,D1sq] = eig(TP1sq)
%     D1 = sqrtm(D1sq)
%     TP1a = V1sq*D1/V1sq
%     D1 = -sqrtm(D1sq)
%     TP1b = V1sq*D1/V1sq
%     D1 = sqrtm(D1sq)*[1 0;0 -1]
%     TP1c = V1sq*D1/V1sq
%     D1 = sqrtm(D1sq)*[-1 0;0 1]
%     TP1d = V1sq*D1/V1sq
%     
%     ST1a =  abcd2s(TP1a,50)
%     ST1b =  abcd2s(TP1b,50)
%     ST1c =  abcd2s(TP1c,50)
%     ST1d =  abcd2s(TP1d,50)
    
    
    TP1 = sqrtm(((TLP)^-1*T2LP*(TLP)^-1)^-1);
    
    
%     TL1 = TP1^-1*TLV*TP1^-1;
    TL1 = TP1;
    TL2 = TP1^-1*T2LV*TP1^-1;
    
    Sfinal = abcd2s(TL1,50);
%     Sfinal = Sc
    
    Sf11r(i)=20*log10(abs(Sfinal(1)));
    Sf11i(i)=asin(imag(Sfinal(1))/abs(Sfinal(1)))/pi*180;

    Sf12r(i)=20*log10(abs(Sfinal(3)));
    Sf12i(i)=asin(imag(Sfinal(3))/abs(Sfinal(3)))/pi*180;

    Sf21r(i)=20*log10(abs(Sfinal(2)));
    Sf21i(i)=asin(imag(Sfinal(2))/abs(Sfinal(2)))/pi*180;

    Sf22r(i)=20*log10(abs(Sfinal(4)));
    Sf22i(i)=asin(imag(Sfinal(4))/abs(Sfinal(4)))/pi*180;
    
        Sfinal1 = abcd2s(TL2,50);
    
    Sf11r1(i)=20*log10(abs(Sfinal1(1)));
    Sf11i1(i)=asin(imag(Sfinal1(1))/abs(Sfinal1(1)))/pi*180;

    Sf12r1(i)=20*log10(abs(Sfinal1(3)));
    Sf12i1(i)=asin(imag(Sfinal1(3))/abs(Sfinal1(3)))/pi*180;

    Sf21r1(i)=20*log10(abs(Sfinal1(2)));
    Sf21i1(i)=asin(imag(Sfinal1(2))/abs(Sfinal1(2)))/pi*180;

    Sf22r1(i)=20*log10(abs(Sfinal1(4)));
    Sf22i1(i)=asin(imag(Sfinal1(4))/abs(Sfinal1(4)))/pi*180;

    
end

xlswrite('2stepdeembeddingResult',freq,'TLV','A2:A202');
xlswrite('2stepdeembeddingResult',Sf11r','TLV','B2:B202');
xlswrite('2stepdeembeddingResult',Sf11i','TLV','C2:C202');
xlswrite('2stepdeembeddingResult',Sf12r','TLV','D2:D202');
xlswrite('2stepdeembeddingResult',Sf12i','TLV','E2:E202');
xlswrite('2stepdeembeddingResult',Sf21r','TLV','F2:F202');
xlswrite('2stepdeembeddingResult',Sf21i','TLV','G2:G202');
xlswrite('2stepdeembeddingResult',Sf22r','TLV','H2:H202');
xlswrite('2stepdeembeddingResult',Sf22i','TLV','I2:I202');

xlswrite('2stepdeembeddingResult',freq,'T2LV','A2:A202');
xlswrite('2stepdeembeddingResult',Sf11r1','T2LV','B2:B202');
xlswrite('2stepdeembeddingResult',Sf11i1','T2LV','C2:C202');
xlswrite('2stepdeembeddingResult',Sf12r1','T2LV','D2:D202');
xlswrite('2stepdeembeddingResult',Sf12i1','T2LV','E2:E202');
xlswrite('2stepdeembeddingResult',Sf21r1','T2LV','F2:F202');
xlswrite('2stepdeembeddingResult',Sf21i1','T2LV','G2:G202');
xlswrite('2stepdeembeddingResult',Sf22r1','T2LV','H2:H202');
xlswrite('2stepdeembeddingResult',Sf22i1','T2LV','I2:I202');


end
