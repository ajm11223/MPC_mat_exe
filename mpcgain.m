function [Phi_Phi, Phi_F, Phi_R, A_e, B_e, C_e] = mpcgain(Ap, Bp, Cp, Nc, Np)
% mpcgain 함수는 MPC gain 행렬을 계산합니다.
% 입력: 
%   Ap, Bp, Cp: 플랜트의 상태 공간 모델 행렬
%   Nc: 제어 구간 (Control Horizon)
%   Np: 예측 구간 (Prediction Horizon)
%   r_w: 가중치

    % 1. 시스템 차원 확인
    [m1, n1] = size(Cp); %Cp = Cm (m1xn1) 매트릭스
    [n1, n_in] = size(Bp); %Bp = Bm (n1xn_in) 매트릭스

    % 2. 확장 상태 공간 모델 (Augmented State-space Model) 생성 
    % A_e: (n1+m1) x (n1+m1) % 교과서 상으로 A mtrx
    A_e = eye(n1 + m1, n1 + m1);
    A_e(1:n1, 1:n1) = Ap;
    A_e(n1+1:n1+m1, 1:n1) = Cp * Ap;

    % B_e: (n1+m1) x n_in
    B_e = zeros(n1 + m1, n_in);
    B_e(1:n1, :) = Bp;
    B_e(n1+1:n1+m1, :) = Cp * Bp;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 01-23 하기

    % C_e: m1 x (n1+m1)
    C_e = zeros(m1, n1 + m1);
    C_e(:, n1+1:n1+m1) = eye(m1, m1);

    % 3. F와 Phi 행렬 계산 
    n = n1 + m1;
    % 메모리 할당 및 초기화
    h = zeros(Np, m1 * (n1 + m1)); % 임시 저장용
    F = zeros(Np, n);             % F 행렬 초기화
    
    h(1,:) = C_e;
    F(1,:) = C_e * A_e;

    % 반복문을 통해 F와 h 계산 (A_e의 거듭제곱 이용)
    for kk = 2:Np
        h(kk,:) = h(kk-1,:) * A_e;
        F(kk,:) = F(kk-1,:) * A_e;
    end

    % v 벡터 계산 (Phi 행렬의 첫 번째 열을 구성하는 요소)
    v = h * B_e;
    
    % Phi 행렬 생성 (Toeplitz 구조)
    Phi = zeros(Np, Nc); % Phi 차원 선언
    Phi(:,1) = v;        % 첫 번째 열
    
    for i = 2:Nc
        % 이전 열을 아래로 한 칸씩 이동 (Shift)
        Phi(:,i) = [zeros(i-1,1); v(1:Np-i+1,1)]; 
    end

    % 4. 비용 함수(J) 최적화에 필요한 행렬 계산 
    BarRs = ones(Np, 1);
    Phi_Phi = Phi' * Phi;     % Hessian Matrix의 주요 부분
    Phi_F = Phi' * F;         % 상태 피드백 이득 계산용
    Phi_R = Phi' * BarRs;     % 기준 신호(Reference) 추적용
    % r_w = 10;
    r_w = 0;
    R_bar = r_w * eye(Nc);
    u_k = 0;
    x_k = [0.1; 0.2];


    % 6. 결과 출력
    fprintf('DelU 제외 원소 매트릭스 \n\n')

    disp('Phi_Phi:');
    disp(Phi_Phi);

    disp('Phi_F:');
    disp(Phi_F);

    disp('Phi_R:');
    disp(Phi_R);

    % 5. step by step Delta U 계산

    for i = 1:Nc

    Delta_U = (Phi_Phi+R_bar)^(-1) * (Phi') * (BarRs*1-F*x_k);
    fprintf('Delta_U at k=%i \n',9+i);
    disp(Delta_U);
    %u(k+1) = u(k) + delu(k) 이므로
    u_k = u_k + Delta_U(1,1);
    x_m = x_k(2,1);
    
    x_m_next = Ap * x_m + Bp * u_k;
    y_next = Cp*x_m_next;
    x_k = [x_m_next-x_m; y_next]
        
        if x_k(1,1) < 1e-4 % 완벽하게 0이 될 수는 없어서
            break;
        end
    end




end