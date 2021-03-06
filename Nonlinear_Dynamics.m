% Compute the natural frequency and the associated period for the
% continuous system dynamics

function [period_continuous, frequency_continuous] = Nonlinear_Dynamics(x, N,m_L1,m_L2)   

%test = [0, 1110];
% Simulation Parameters


r1 = zeros(N,1);
r2 = zeros(N,1);

A_prime = zeros(8,8,2);
B_prime = zeros(8,2,2);
bias_prime = zeros(8,1,2);
theta1 = zeros(N,1);
theta2_unshift = zeros(N,1);
theta2_shift = zeros(N,1);
theta_dot1 = zeros(N,1);
theta_dot2 = zeros(N,1);



%% Iterative Simulation
max_eig = 0; 
i = 1;
hold on
for j = 1:N
    %for i = 1:2   
        % x(+,:) = A*x(-,:) + B(:,1)*u(1) + B(:,2)*u(2) +
        %          + B(:,3)*u(3) + B(:,4)*u(4)          ;
        %m_L1 = test(i);
        %m_L2 = test(i);
        
        %% Continuous Time System
        [~,~,~,A_1,B_1] = get_continuous_matrices(m_L1,m_L2);
        %% Matrices for Linearization
        s1 = [0 1 0 1 0 0 0 0];
        s2 = [0 0 0 0 0 1 0 1];
        s = [s1;s2];
        B1_1 = B_1(:,1);
        B1_2 = B_1(:,2);
        B1_3 = B_1(:,3);
        B1_4 = B_1(:,4);
        B1_12 = B_1(:,1:2);
        B1_34 = B_1(:,3:4);
        
        x1_est = x(j,1)+x(j,3); % link 1 displacement
        dx1_est = x(j,2)+x(j,4); % link 1 velocity 
        x2_est = x(j,5)+x(j,7); % link 1 displacement
        dx2_est = x(j,6)+x(j,8); % link 1 velocity          
        % Note that theta2 is with respect to link 1. The angle of link 2
        % with respect to the horizontal axis is theta1+theta2
        theta1(j) = get_theta1(x1_est);
        theta2_unshift(j) = get_theta_2(x2_est);
        theta2_shift(j)= theta2_unshift(j)+pi;

        [r1(j), r2(j)] = getr(theta1(j),theta2_unshift(j));
       
        theta_dot1(j) = dx1_est/r1(j);
        theta_dot2(j) = dx2_est/r2(j);

        [H11, H22, H12, H21, h, G1, G2] = get_dynamic_components(theta1(j),theta2_shift(j));

        [F_x_components, F_u_components, F_b_components] = ...
            get_F_components(H11, H22, H12, H21, h, G1, G2,...
            theta1(j), theta2_unshift(j), r1(j), r2(j), theta_dot1(j), theta_dot2(j),...
            A_1, B1_1, B1_2, B1_3, B1_4, s1, s2,m_L1,m_L2);

        [A_prime(:,:,i),B_prime(:,:,i),bias_prime(:,:,i)] = linearize_matrices(A_1,B1_12,B1_34,F_x_components,...
        F_u_components,F_b_components);
 
        
        eigenvalues = eig( A_prime(:,:,i) );
        eig_test = max(abs(eigenvalues));
        if eig_test > max_eig(i)
            max_eig(i) = eig_test;
        end 
               
    %end 
    
%     if(norm((A_prime(:,:,1)-A_prime(:,:,2)),'fro') > 1e-7)   
%         stop = 0;
%         disp(A_prime(:,:,1))
%         disp(A_prime(:,:,2))
%         
%     end 


end

% Plots 
% period_0 =1/((max_eig(1))/(2*pi))
% period_mass =1/((max_eig(2))/(2*pi))

% Maximum eigenvalue and corresponding period
period_continuous = 1/((max_eig)/(2*pi));
frequency_continuous = max_eig;





end

