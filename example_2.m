clc;
clear all;
close all;

% 변수 정의
a = 0.8;
b = 0.1;
N_p = 10;
N_c = 4;
% N_c = 9;
r_ki = 1;
Rs_bar = ones(1,N_p)';
x_ki = [0.1 0.2]';
A_mat = [a 0;a 1];
B_mat = [b;b];
C_mat = [0 1];
F_mat = zeros(N_p,2);
lil_Phi_mat = zeros(N_p,1);
Phi_mat = zeros(N_p,N_c);

% Find F_mat
for i = 1:N_p
    F_mat(i,:) = C_mat*(A_mat)^i;
end
% F 매트릭스는 N_p 행 매트릭스, CA부터 CA^10까지 C와 A의 k제곱근의 매트릭스 곱


% Find lil_Phi_mat
for i = 1:N_p
    lil_Phi_mat(i,:) = C_mat*(A_mat^(i-1))*B_mat;
end
% Phi 매트릭스 만들기 위한 작업, 어차피 Phi 매트릭스는 CB, CAB, C(A^2)B,.....
...C(A^Np-1)B 로 구성된 (NpX1)column matrix로 이루어진 (NpXNc) matrx 이므로..
...하나의 열만 구성하면 그것의 반복적 사용으로 Phi 매트릭스를 만들 수 있다.


% Find Phi_mat
for i = 1:N_c
    Phi_mat(i:end,i) = lil_Phi_mat(1:end-i+1);
end
    
% calculate phiTphi phiTf phiTRs_bar
PTP=Phi_mat'*Phi_mat;
PTF=Phi_mat'*F_mat;
PTR=Phi_mat'*Rs_bar;

% 비교할 r_w 값들을 배열로 정의
rw_values = [0, 10]; 

figure; % 하나의 창 생성

for k = 1:length(rw_values)
    r_w = rw_values(k);
    
    % 1. 해당 r_w에 대한 Delta U 계산
    R_bar = r_w * eye(N_c);
    Del_U = (PTP + R_bar)^(-1) * (PTR - PTF * x_ki);
    
    % 2. 예측 궤적 시뮬레이션 (초기화)
    x_temp = x_ki;
    x_history = zeros(2, N_p+1);
    x_history(:, 1) = x_ki;
    y_history = zeros(1, N_p+1);
    y_history(1) = C_mat * x_ki;
    
    % 3. 미래 상태 예측 루프
    for t = 1:N_p
        if t <= N_c
            u_k = Del_U(t);
        else
            u_k = 0;
        end
        x_temp = A_mat * x_temp + B_mat * u_k;
        x_history(:, t+1) = x_temp;
        y_history(t+1) = C_mat * x_temp;
    end
    
    % 4. Subplot으로 그리기
    subplot(1, 2, k); % 2행 1열의 그리드 중 k번째 위치에 그림
    time_steps = 10 : (10 + N_p);
    
    plot(time_steps, x_history(1,:), '.-', 'DisplayName', '\Delta x_m'); hold on;
    plot(time_steps, y_history, 'o-', 'DisplayName', 'y');
    
    title(['Simulation result with r_w = ', num2str(r_w)]);
    xlabel('Sampling Instant'); 
    ylabel('Response');
    legend('Location', 'best');
    grid on;
    ylim([-0.2 1.2]); % y축 범위를 고정하면 비교하기 더 좋습니다
    fprintf('\n r_w가 %s일 경우 Delta U \n',num2str(r_w));

    fprintf('%f \n',Del_U);
end