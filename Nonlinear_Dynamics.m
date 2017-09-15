% Compute the natural frequency and the associated period for the
% continuous system dynamics

% Notes
% theta1 is an absolute angle
% theta2 is shifted by pi in order to be with respect to link 1
% The absolute angle of link 2 is theta1+theta2

function [period_continuous, frequency_continuous] = Nonlinear_Dynamics(x, N,m_L1,m_L2)   

%test = [0, 1110];
%% System Parameters
 
% moment arm (Jacobian) 
get_r1 = @(theta) (21*sin(theta - atan(1/2)))/(2500*(457/10000 - (21*cos(theta - atan(1/2)))/1250)^(1/2));
get_r2 = @(theta) ((sin(pi/2 + asin(sin(pi/6 - acos(cos(asin((3*sin(theta))/(100*(9/5000 - (9*cos(theta))/5000)^(1/2))) + asin((1 - (390625*((9*cos(theta))/5000 + 7/5000)^2)/4)^(1/2)/(25*(9/5000 - (9*cos(theta))/5000)^(1/2))))))/(20*(17/400 - cos(pi/6 - acos(cos(asin((3*sin(theta))/(100*(9/5000 - (9*cos(theta))/5000)^(1/2))) + asin((1 - (390625*((9*cos(theta))/5000 + 7/5000)^2)/4)^(1/2)/(25*(9/5000 - (9*cos(theta))/5000)^(1/2))))))/50)^(1/2))))*(17/400 - cos(pi/6 - acos(cos(asin((3*sin(theta))/(100*(9/5000 - (9*cos(theta))/5000)^(1/2))) + asin((1 - (390625*((9*cos(theta))/5000 + 7/5000)^2)/4)^(1/2)/(25*(9/5000 - (9*cos(theta))/5000)^(1/2))))))/50)^(1/2)*((sin(asin((3*sin(theta))/(100*(9/5000 - (9*cos(theta))/5000)^(1/2))) + asin((1 - (390625*((9*cos(theta))/5000 + 7/5000)^2)/4)^(1/2)/(25*(9/5000 - (9*cos(theta))/5000)^(1/2))))*cos(pi/6 - acos(cos(asin((3*sin(theta))/(100*(9/5000 - (9*cos(theta))/5000)^(1/2))) + asin((1 - (390625*((9*cos(theta))/5000 + 7/5000)^2)/4)^(1/2)/(25*(9/5000 - (9*cos(theta))/5000)^(1/2))))))*(((9*sin(theta)*(1 - (390625*((9*cos(theta))/5000 + 7/5000)^2)/4)^(1/2))/(250000*(9/5000 - (9*cos(theta))/5000)^(3/2)) - (225*sin(theta)*((9*cos(theta))/5000 + 7/5000))/(32*(9/5000 - (9*cos(theta))/5000)^(1/2)*(1 - (390625*((9*cos(theta))/5000 + 7/5000)^2)/4)^(1/2)))/(1 - ((390625*((9*cos(theta))/5000 + 7/5000)^2)/4 - 1)/((9*cos(theta))/8 - 9/8))^(1/2) - ((3*cos(theta))/(100*(9/5000 - (9*cos(theta))/5000)^(1/2)) - (27*sin(theta)^2)/(1000000*(9/5000 - (9*cos(theta))/5000)^(3/2)))/((9*sin(theta)^2)/(18*cos(theta) - 18) + 1)^(1/2)))/(20*(17/400 - cos(pi/6 - acos(cos(asin((3*sin(theta))/(100*(9/5000 - (9*cos(theta))/5000)^(1/2))) + asin((1 - (390625*((9*cos(theta))/5000 + 7/5000)^2)/4)^(1/2)/(25*(9/5000 - (9*cos(theta))/5000)^(1/2))))))/50)^(1/2)*(1 - cos(asin((3*sin(theta))/(100*(9/5000 - (9*cos(theta))/5000)^(1/2))) + asin((1 - (390625*((9*cos(theta))/5000 + 7/5000)^2)/4)^(1/2)/(25*(9/5000 - (9*cos(theta))/5000)^(1/2))))^2)^(1/2)) - (sin(asin((3*sin(theta))/(100*(9/5000 - (9*cos(theta))/5000)^(1/2))) + asin((1 - (390625*((9*cos(theta))/5000 + 7/5000)^2)/4)^(1/2)/(25*(9/5000 - (9*cos(theta))/5000)^(1/2))))*sin(pi/6 - acos(cos(asin((3*sin(theta))/(100*(9/5000 - (9*cos(theta))/5000)^(1/2))) + asin((1 - (390625*((9*cos(theta))/5000 + 7/5000)^2)/4)^(1/2)/(25*(9/5000 - (9*cos(theta))/5000)^(1/2))))))^2*(((9*sin(theta)*(1 - (390625*((9*cos(theta))/5000 + 7/5000)^2)/4)^(1/2))/(250000*(9/5000 - (9*cos(theta))/5000)^(3/2)) - (225*sin(theta)*((9*cos(theta))/5000 + 7/5000))/(32*(9/5000 - (9*cos(theta))/5000)^(1/2)*(1 - (390625*((9*cos(theta))/5000 + 7/5000)^2)/4)^(1/2)))/(1 - ((390625*((9*cos(theta))/5000 + 7/5000)^2)/4 - 1)/((9*cos(theta))/8 - 9/8))^(1/2) - ((3*cos(theta))/(100*(9/5000 - (9*cos(theta))/5000)^(1/2)) - (27*sin(theta)^2)/(1000000*(9/5000 - (9*cos(theta))/5000)^(3/2)))/((9*sin(theta)^2)/(18*cos(theta) - 18) + 1)^(1/2)))/(2000*(17/400 - cos(pi/6 - acos(cos(asin((3*sin(theta))/(100*(9/5000 - (9*cos(theta))/5000)^(1/2))) + asin((1 - (390625*((9*cos(theta))/5000 + 7/5000)^2)/4)^(1/2)/(25*(9/5000 - (9*cos(theta))/5000)^(1/2))))))/50)^(3/2)*(1 - cos(asin((3*sin(theta))/(100*(9/5000 - (9*cos(theta))/5000)^(1/2))) + asin((1 - (390625*((9*cos(theta))/5000 + 7/5000)^2)/4)^(1/2)/(25*(9/5000 - (9*cos(theta))/5000)^(1/2))))^2)^(1/2))))/(50*(sin(pi/6 - acos(cos(asin((3*sin(theta))/(100*(9/5000 - (9*cos(theta))/5000)^(1/2))) + asin((1 - (390625*((9*cos(theta))/5000 + 7/5000)^2)/4)^(1/2)/(25*(9/5000 - (9*cos(theta))/5000)^(1/2))))))^2/(8*cos(pi/6 - acos(cos(asin((3*sin(theta))/(100*(9/5000 - (9*cos(theta))/5000)^(1/2))) + asin((1 - (390625*((9*cos(theta))/5000 + 7/5000)^2)/4)^(1/2)/(25*(9/5000 - (9*cos(theta))/5000)^(1/2)))))) - 17) + 1)^(1/2)) + (sin(asin((3*sin(theta))/(100*(9/5000 - (9*cos(theta))/5000)^(1/2))) + asin((1 - (390625*((9*cos(theta))/5000 + 7/5000)^2)/4)^(1/2)/(25*(9/5000 - (9*cos(theta))/5000)^(1/2))))*sin(pi/6 - acos(cos(asin((3*sin(theta))/(100*(9/5000 - (9*cos(theta))/5000)^(1/2))) + asin((1 - (390625*((9*cos(theta))/5000 + 7/5000)^2)/4)^(1/2)/(25*(9/5000 - (9*cos(theta))/5000)^(1/2))))))*(((9*sin(theta)*(1 - (390625*((9*cos(theta))/5000 + 7/5000)^2)/4)^(1/2))/(250000*(9/5000 - (9*cos(theta))/5000)^(3/2)) - (225*sin(theta)*((9*cos(theta))/5000 + 7/5000))/(32*(9/5000 - (9*cos(theta))/5000)^(1/2)*(1 - (390625*((9*cos(theta))/5000 + 7/5000)^2)/4)^(1/2)))/(1 - ((390625*((9*cos(theta))/5000 + 7/5000)^2)/4 - 1)/((9*cos(theta))/8 - 9/8))^(1/2) - ((3*cos(theta))/(100*(9/5000 - (9*cos(theta))/5000)^(1/2)) - (27*sin(theta)^2)/(1000000*(9/5000 - (9*cos(theta))/5000)^(3/2)))/((9*sin(theta)^2)/(18*cos(theta) - 18) + 1)^(1/2)))/(50*(1 - cos(asin((3*sin(theta))/(100*(9/5000 - (9*cos(theta))/5000)^(1/2))) + asin((1 - (390625*((9*cos(theta))/5000 + 7/5000)^2)/4)^(1/2)/(25*(9/5000 - (9*cos(theta))/5000)^(1/2))))^2)^(1/2)) - (cos(pi/2 + asin(sin(pi/6 - acos(cos(asin((3*sin(theta))/(100*(9/5000 - (9*cos(theta))/5000)^(1/2))) + asin((1 - (390625*((9*cos(theta))/5000 + 7/5000)^2)/4)^(1/2)/(25*(9/5000 - (9*cos(theta))/5000)^(1/2))))))/(20*(17/400 - cos(pi/6 - acos(cos(asin((3*sin(theta))/(100*(9/5000 - (9*cos(theta))/5000)^(1/2))) + asin((1 - (390625*((9*cos(theta))/5000 + 7/5000)^2)/4)^(1/2)/(25*(9/5000 - (9*cos(theta))/5000)^(1/2))))))/50)^(1/2))))*sin(asin((3*sin(theta))/(100*(9/5000 - (9*cos(theta))/5000)^(1/2))) + asin((1 - (390625*((9*cos(theta))/5000 + 7/5000)^2)/4)^(1/2)/(25*(9/5000 - (9*cos(theta))/5000)^(1/2))))*sin(pi/6 - acos(cos(asin((3*sin(theta))/(100*(9/5000 - (9*cos(theta))/5000)^(1/2))) + asin((1 - (390625*((9*cos(theta))/5000 + 7/5000)^2)/4)^(1/2)/(25*(9/5000 - (9*cos(theta))/5000)^(1/2))))))*(((9*sin(theta)*(1 - (390625*((9*cos(theta))/5000 + 7/5000)^2)/4)^(1/2))/(250000*(9/5000 - (9*cos(theta))/5000)^(3/2)) - (225*sin(theta)*((9*cos(theta))/5000 + 7/5000))/(32*(9/5000 - (9*cos(theta))/5000)^(1/2)*(1 - (390625*((9*cos(theta))/5000 + 7/5000)^2)/4)^(1/2)))/(1 - ((390625*((9*cos(theta))/5000 + 7/5000)^2)/4 - 1)/((9*cos(theta))/8 - 9/8))^(1/2) - ((3*cos(theta))/(100*(9/5000 - (9*cos(theta))/5000)^(1/2)) - (27*sin(theta)^2)/(1000000*(9/5000 - (9*cos(theta))/5000)^(3/2)))/((9*sin(theta)^2)/(18*cos(theta) - 18) + 1)^(1/2)))/(5000*(17/400 - cos(pi/6 - acos(cos(asin((3*sin(theta))/(100*(9/5000 - (9*cos(theta))/5000)^(1/2))) + asin((1 - (390625*((9*cos(theta))/5000 + 7/5000)^2)/4)^(1/2)/(25*(9/5000 - (9*cos(theta))/5000)^(1/2))))))/50)^(1/2)*(1 - cos(asin((3*sin(theta))/(100*(9/5000 - (9*cos(theta))/5000)^(1/2))) + asin((1 - (390625*((9*cos(theta))/5000 + 7/5000)^2)/4)^(1/2)/(25*(9/5000 - (9*cos(theta))/5000)^(1/2))))^2)^(1/2)))/(2*(213/5000 - (cos(pi/2 + asin(sin(pi/6 - acos(cos(asin((3*sin(theta))/(100*(9/5000 - (9*cos(theta))/5000)^(1/2))) + asin((1 - (390625*((9*cos(theta))/5000 + 7/5000)^2)/4)^(1/2)/(25*(9/5000 - (9*cos(theta))/5000)^(1/2))))))/(20*(17/400 - cos(pi/6 - acos(cos(asin((3*sin(theta))/(100*(9/5000 - (9*cos(theta))/5000)^(1/2))) + asin((1 - (390625*((9*cos(theta))/5000 + 7/5000)^2)/4)^(1/2)/(25*(9/5000 - (9*cos(theta))/5000)^(1/2))))))/50)^(1/2))))*(17/400 - cos(pi/6 - acos(cos(asin((3*sin(theta))/(100*(9/5000 - (9*cos(theta))/5000)^(1/2))) + asin((1 - (390625*((9*cos(theta))/5000 + 7/5000)^2)/4)^(1/2)/(25*(9/5000 - (9*cos(theta))/5000)^(1/2))))))/50)^(1/2))/50 - cos(pi/6 - acos(cos(asin((3*sin(theta))/(100*(9/5000 - (9*cos(theta))/5000)^(1/2))) + asin((1 - (390625*((9*cos(theta))/5000 + 7/5000)^2)/4)^(1/2)/(25*(9/5000 - (9*cos(theta))/5000)^(1/2))))))/50)^(1/2));

%% Simulation Parameters


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

        r1(j) = get_r1(theta1(j));
        r2(j) = get_r2(theta2_unshift(j));        
        theta_dot1(j) = dx1_est/r1(j);
        theta_dot2(j) = dx2_est/r2(j);

        [H11, H22, H12, H21, h, G1, G2] = get_dynamic_components(theta1(j),theta2_shift(j));

        [F_x_components, F_u_components, F_b_components] = ...
            get_F_components(H11, H22, H12, H21, h, G1, G2,...
            theta1(j), theta2_unshift(j), r1(j), r2(j), theta_dot1(j), theta_dot2(j),...
            A_1, B1_1, B1_2, B1_3, B1_4, s1, s2,m_L1,m_L2);

        [A_prime(:,:,i),B_prime(:,:,i),bias_prime(:,:,i)] = linearize(A_1,B1_12,B1_34,F_x_components,...
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

period_continuous = 1/((max_eig)/(2*pi));
frequency_continuous = max_eig;





end

