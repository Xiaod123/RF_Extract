function abcd = s2abcd(s, Z0)

S11 = s(1,1);
S12 = s(1,2);
S21 = s(2,1);
S22 = s(2,2);

Z01 = Z0;
Z02 = Z0;
R01 = Z0;
R02 = Z0;

denom = 2*S21*sqrt(R01*R02);

A = ( (conj(Z01) + S11*Z01)*(1-S22) + S12*S21*Z01) / denom;
B = ( ( conj(Z01) + S11*Z01)*(conj(Z02) + S22*Z02) - S12*S21*Z01*Z02 ) / denom;
C = ( (1-S11)*(1-S22) - S12*S21 ) / denom;
D = ( (1-S11)*(conj(Z02) + S22*Z02) + S12*S21*Z02 ) / denom;

abcd = [A B ; C D];

%function [abcd] = s2abcd(s,z0)
% z0 = Z0;
% s11 = s(1,1);
% s12 = s(1,2);
% s21 = s(2,1);
% s22 = s(2,2);
% 
% a = ((1+s11)*(1-s22)+s12*s21)/(2*s21);
% b = z0*(((1+s11)*(1+s22)-s12*s21)/(2*s21));
% c = ((1-s11)*(1-s22)-s12*s21)/(2*s21*z0);
% d = ((1-s11)*(1+s22)+s12*s21)/(2*s21); 
% abcd = [a b ;c d];
% % s = [s11 s12;s21 s21];