clc; clear all; close all;

% 시스템 파라미터 (Example 1.2 & 1.4) [cite: 827]
Ap = 0.8; Bp = 0.1; Cp = 1;
Nc = 4; Np = 10;
rw = 0; % 가중치 (rw = 0 또는 10으로 테스트 가능) 

% 1. MPC Gain 계산 (mpcgain 함수가 같은 폴더에 있어야 합니다)
[Phi_Phi, Phi_F, Phi_R, A_e, B_e, C_e] = mpcgain(Ap, Bp, Cp, Nc, Np);

% 2. 시뮬레이션 설정 
N_sim = 15;             % 시뮬레이션 반복 횟수
x_k = [0.1; 0.2];       % k=10 시점의 초기 확장 상태 [delta_xm(10); y(10)]
u_k_prev = 0;           % u(9) = 0 (초기 제어 입력) 
r_k = 1;                % 목표치 (Set-point)

% 결과 저장용 변수 초기화 (초기 k=10 시점의 값 저장)
y_history = [x_k(2)];   % y(10) = 0.2 
u_history = [u_k_prev]; % u(9) = 0 

% 역행렬 미리 계산 (rw 적용) [cite: 640]
R_bar = rw * eye(Nc);
H_inv = (Phi_Phi + R_bar) \ eye(Nc); 

% 3. Receding Horizon 제어 루프
for k = 1:N_sim
    % (1) 최적 제어 증분 시퀀스 Delta_U 계산 [cite: 649, 1174]
    Delta_U = H_inv * (Phi_R * r_k - Phi_F * x_k);
    
    % (2) 첫 번째 요소만 선택 (Receding Horizon Principle) [cite: 823, 1515]
    delta_u = Delta_U(1, 1);
    
    % (3) 제어 입력 업데이트: u(k) = u(k-1) + delta_u(k) 
    u_k = u_k_prev + delta_u;
    
    % (4) 플랜트 시뮬레이션 (다음 상태 계산) [cite: 831]
    % Cp=1이므로 xm(k) = y(k)입니다. [cite: 653]
    xm_current = x_k(2); 
    xm_next = Ap * xm_current + Bp * u_k; 
    y_next = Cp * xm_next;
    
    % (5) 데이터 저장 (오류 해결 지점: 계산 후 저장 ⭐)
    u_history = [u_history; u_k];
    y_history = [y_history; y_next];
    
    % (6) 다음 루프를 위한 확장 상태 x(k+1) 구성 [cite: 833]
    % x(k+1) = [delta_xm(k+1); y(k+1)] = [xm(k+1)-xm(k); y(k+1)]
    x_k = [(xm_next - xm_current); y_next];
    u_k_prev = u_k; % 현재 u(k)가 다음 루프의 u(k-1)이 됨
end

% 4. 결과 시각화
k_axis = 10 : (10 + N_sim); % k=10부터 시작하는 축 생성 
figure;
subplot(2,1,1); plot(k_axis, y_history, 'b-o', 'LineWidth', 1.5); hold on;
plot(k_axis, ones(size(k_axis))*r_k, 'r--'); % 목표선 (r=1)
ylabel('Output (y)'); title(['MPC Receding Horizon Control (rw = ', num2str(rw), ')']);
grid on;

subplot(2,1,2); stairs(k_axis, u_history, 'm-s', 'LineWidth', 1.5);
ylabel('Control (u)'); xlabel('Sampling Instant (k)');
grid on;