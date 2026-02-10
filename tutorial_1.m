clear all; clc; close all;
% 변수 정의
Ac= [0 1 0; 3 0 1; 0 1 0];
Bc= [1; 1; 3];
Cc= [0 1 0];
Dc=zeros(1,1); % 입력이 출력에 동시에 영향을 끼치지 않으므로 D=0mat

Delta_t=1; % 샘플링 간격 1초

sysc = ss(Ac, Bc, Cc, Dc); % 1. 연속시간 상태공간 모델 생성 [cite: 557]
sysd = c2d(sysc, Delta_t);  % 2. c2d 함수를 사용하여 이산화 수행 
[Ad, Bd, Cd, Dd] = ssdata(sysd); % 3. 이산화된 행렬 추출


% 교재의 Am,Bm 등이 현재 코드의 Ac, Cc 등이다..

[m1,n1]=size(Cd); % Cd는 현재 1x3 mat 즉, 출력 y의 개수는 1, 상태변수 Xm의 개수는 3
[n1,n_in]=size(Bd); % Bd는 현재 3x1 mat 즉, 상태 변수의 개수:3, 입력의 개수: 1
A_e=eye(n1+m1,n1+m1); % Ae를 초기에 단위행렬로 설정, 우측 하단의 적분기에 I 삽입.
A_e(1:n1,1:n1)=Ad; % Ae의 좌상단 mat Ad로 채우기
A_e(n1+1:n1+m1,1:n1)=Cd*Ad; % Ae의 좌하단 mat CdAd로 채우기 -> 나머지는 이미 0과 1이므로 Ae mat 끝
B_e=zeros(n1+m1,n_in); % Be 행렬의 크기만큼의 0행렬 만들기
B_e(1:n1,:)=Bd; % 상단엔 Bd
B_e(n1+1:n1+m1,:)=Cd*Bd; % 하단엔 CdBd
C_e=zeros(m1,n1+m1); % Ce 행렬의 크기만큼의 0행렬 만들기
C_e(:,n1+1:n1+m1)=eye(m1,m1); % 위 코드 분석하면 ABCD 매트릭스에 대한 일반적인 성질을 학습 가능 할 듯