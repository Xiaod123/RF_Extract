function S = abcd2s(abcd, z0)

R01 = z0;
R02 = z0;
Z01 = z0;
Z02 = z0;

A = abcd(1,1);
B = abcd(1,2);
C = abcd(2,1);
D = abcd(2,2);

denom = (A*Z02 + B + C*Z01*Z02 + D*Z01);

S11 = ( A*Z02 + B - C*conj(Z01)*Z02 - D*conj(Z01) ) / denom;
S12 = ( 2*(A*D - B*C)*sqrt(R01*R02) ) / denom;
S21 = ( 2*sqrt(R01*R02) ) / denom;
S22 = (-A*conj(Z02) + B - C*Z01*conj(Z02) + D*Z01 ) / denom;

S = [ S11 S12 ; S21 S22 ];


% function [s_new] = abcd2s(abcd,z0)
% 
% A = abcd(1,1);
% B = abcd(1,2);
% C = abcd(2,1);
% D = abcd(2,2);
% divider = A+(B/z0)+(C*z0)+D;
% s11new = (A+(B/z0)-(C*z0)-D)/divider;
% s12new = (2*(A*D-B*C))/divider;
% s21new = 2/divider;
% s22new = (-A+(B/z0)-(C*z0)+D)/divider;
% s_new = [s11new s12new; s21new s22new];