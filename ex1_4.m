clc; clear all; close all;

%% 1. 시스템 파라미터 및 설정 (Ex 1.2 & 1.4) [cite: 827]
Ap = 0.8; Bp = 0.1; Cp = 1;
Nc = 4;   Np = 10;
rw = 0;   % 가중치 (0 또는 10으로 테스트 가능) [cite: 628]

% MPC 이득 및 확장 모델 행렬 계산
[Phi_Phi, Phi_F, Phi_R, A_e, B_e, C_e] = mpcgain(Ap, Bp, Cp, Nc, Np);

%% 2. 시뮬레이션 준비
N_sim = 15;
r_k = 1;                % 목표치 (Set-point) [cite: 620]
x_k = [0.1; 0.2];       % 초기 확장 상태 [delta_xm(10); y(10)] [cite: 828]
u_k_prev = 0;           % u(k-1) 초기값

% 결과 저장용 벡터 미리 할당 (속도 최적화)
y_history = zeros(N_sim + 1, 1);
u_history = zeros(N_sim + 1, 1);

y_history(1) = x_k(2);  % 초기 출력 y(10)
u_history(1) = u_k_prev;

% 제어 이득(K_mpc, K_y) 사전 계산 (루프 내부 연산 감소) [cite: 893, 1179]
% Delta_u(k) = Ky*r(k) - Kmpc*x(k)
H_inv = (Phi_Phi + rw * eye(Nc)) \ eye(Nc);
K_gain = H_inv * Phi_F;
Ky_gain = H_inv * Phi_R;

% 첫 번째 행만 추출 (Receding Horizon Principle) [cite: 823, 1515]
K_mpc = K_gain(1, :);
Kx = K_gain(1, 1);
Ky = Ky_gain(1, :);

% closed-loop eigenvalue 찾기
clo_loop = A_e-B_e*K_mpc;
e = eig(clo_loop);
e_data = [real(e)'; imag(e)'];

%% 3. Receding Horizon 제어 루프
for k = 1:N_sim
    % (1) 최적 제어 증분 계산 (첫 번째 요소만 사용) [cite: 891, 1515]
    delta_u = Ky * r_k - K_mpc * x_k;
    
    % (2) 현재 제어 입력 업데이트 [cite: 830]
    u_k = u_k_prev + delta_u;
    
    % (3) 플랜트 시뮬레이션 (실제 시스템 반영) [cite: 831]
    % 현재 출력 y(k)는 x_k(2)에 저장되어 있음
    y_current = x_k(2);
    xm_next = Ap * y_current + Bp * u_k; % Cp=1 가정
    y_next = Cp * xm_next;
    
    % (4) 데이터 저장
    u_history(k+1) = u_k;
    y_history(k+1) = y_next;
    
    % (5) 다음 루프를 위한 확장 상태 업데이트 [cite: 833]
    x_k = [(y_next - y_current); y_next];
    u_k_prev = u_k;
end

%% 4. 결과 시각화
k_axis = 10 : (10 + N_sim);
figure('Name', 'MPC Optimized Simulation');

subplot(2,1,1);
plot(k_axis, y_history, 'b-o', 'LineWidth', 1.5); hold on;
line([k_axis(1) k_axis(end)], [r_k r_k], 'Color', 'r', 'LineStyle', '--');
ylabel('Output (y)'); grid on;
title(['MPC Receding Horizon Control (rw = ', num2str(rw), ')']);

subplot(2,1,2);
stairs(k_axis, u_history, 'm-s', 'LineWidth', 1.5);
ylabel('Control (u)'); xlabel('Sampling Instant (k)'); grid on;

fprintf('K_mpc is');
disp(K_mpc);

fprintf('\n K_x is');
disp(Kx);

fprintf('\n K_y is');
disp(Ky);

fprintf('\n closed-loop eigenvalue is %f + (%f)i ', e_data);
fprintf('\n');