clc; clear all; close all;

% 시스템 파라미터 (Example 1.2)
Ap = 0.8;
Bp = 0.1;
Cp = 1;

% 제어 및 예측 구간 설정
Nc = 4;
Np = 10;

% 함수 실행
[Phi_Phi, Phi_F, Phi_R, A_e, B_e, C_e] = mpcgain(Ap, Bp, Cp, Nc, Np);

