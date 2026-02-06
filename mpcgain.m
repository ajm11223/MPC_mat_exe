function [Phi_Phi, Phi_F, Phi_R, A_e, B_e, C_e] = mpcgain(Ap, Bp, Cp, Nc, Np)
% mpcgain: MPC 설계에 필요한 확장 모델 및 예측 행렬 계산 

    %% 1. 시스템 차원 확인 및 확장 모델(Augmented Model) 생성
    [m1, n1] = size(Cp); 
    [~, n_in] = size(Bp); 

    % A_e: 확장 시스템 행렬 [cite: 527]
    A_e = eye(n1 + m1);
    A_e(1:n1, 1:n1) = Ap;
    A_e(n1+1:end, 1:n1) = Cp * Ap;

    % B_e: 확장 입력 행렬 [cite: 527]
    B_e = [Bp; Cp * Bp];

    % C_e: 확장 출력 행렬 [cite: 527]
    C_e = zeros(m1, n1 + m1);
    C_e(:, n1+1:end) = eye(m1);

    %% 2. F와 Phi 행렬 계산 [cite: 617]
    n = n1 + m1;
    F = zeros(Np * m1, n);
    Phi = zeros(Np * m1, Nc * n_in);
    
    % h_temp는 C_e * A_e^i를 계산하기 위한 임시 행렬
    h_temp = C_e;
    
    % 첫 번째 행 계산
    F(1:m1, :) = C_e * A_e;
    v = C_e * B_e; % Phi의 첫 번째 열 블록 요소
    
    % 재귀적 행렬 계산 (연산량 최적화)
    for kk = 2:Np
        F((kk-1)*m1+1 : kk*m1, :) = F((kk-2)*m1+1 : (kk-1)*m1, :) * A_e;
        h_temp = h_temp * A_e;
        v_block = h_temp * B_e;
        
        % Toeplitz 구조 생성 
        % (v_block을 아래로 쌓아서 Phi의 첫 열 구성)
        v((kk-1)*m1+1 : kk*m1, :) = v_block;
    end

    % Phi 행렬 구성 (Toeplitz 구조 적용) [cite: 807]
    for i = 1:Nc
        Phi((i-1)*m1+1 : end, (i-1)*n_in+1 : i*n_in) = v(1 : (Np-i+1)*m1, :);
    end

    %% 3. 최적화용 행렬 계산 [cite: 649, 1174]
    BarRs = ones(Np * m1, 1);
    Phi_Phi = Phi' * Phi;     % Hessian Matrix 부분 [cite: 644]
    Phi_F = Phi' * F;         % 상태 피드백 이득용 [cite: 888]
    Phi_R = Phi' * BarRs;     % 기준 신호 추적용 [cite: 888]

end