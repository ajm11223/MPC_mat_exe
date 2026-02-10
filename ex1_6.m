clc; clear all; close all;
omega=10;
numc=omega^2; %전달 함수 분자
denc=[1 0.1*omega omega^2]; % 전달함수 분모
[Ac,Bc,Cc,Dc]=tf2ss(numc,denc);

Delta_t=0.01; % 샘플링 간격 1초

%from tutorial_1.m
sysc = ss(Ac, Bc, Cc, Dc); % 1. 연속시간 상태공간 모델 생성 [cite: 557]
sysd = c2d(sysc, Delta_t);  % 2. c2d 함수를 사용하여 이산화 수행 
[Ap, Bp, Cp, Dp] = ssdata(sysd); % 3. 이산화된 행렬 추출

Nc = 3;   Np = 20;
rw = 0.5;