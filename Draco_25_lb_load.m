% Main File for simulation and analysis of the Draco leg 

% Notes
% theta1 is an absolute angle
% theta2 is shifted by pi in order to be with respect to link 1
% The absolute angle of link 2 is theta1+theta2

clear all;
close all;

global m_M b_M k_M Km_M m_S b_S work k_L b_L Km_L...
    J1 J2 m_1 m_2 L_1 L_2 L_c1 L_c2 time_step k_S


% There are different options to run this program:
% 1) Input = 1: Solve the optimization 
% 1A) utilize_compliance = 1: Solve the problem while considering the actuator
% subsystems' dynamics
%   i. Euler_method = 0: 
%   ii. Euler_method = 1: Use Euler's method to linearize problem in the convex
%   optimization sub-problem
% 1B) utilize_compliance = 0: Solve the problem without considering the
% actuator subsystem's dynamics
% 2) Input = 0: Analyze the zero input behavior of the system


%% Parameters for Comparisons and Analysis

% 1 equals true, utilize compliance
% 0 equals false, no actuator dynamics considered
utilize_compliance =1; 

% 1 = true, use Euler method (utilize_compliance should = 1 and input
% should = 1)
% 0 = false
Euler_method = 0;

% Set input to 0 to analyze zero input behavior
% For input = 0, utilize_compliance should = 1 and Euler_method should = 0
% Set input to 1 to run convex optimization program
input = 1;

% 1 = true, compare compliant to rigid behavior
% This sets the cost function such that the penalty is equal to 1e-3 in the
% cost function
compare_compliant_to_rigid = 0;

% 1 = true, analyze convergence rate 
% In this configuration, set:
% utilize_compliance = 1
% Euler_method = 0
% input = 1
analyze_convergence = 0;
if analyze_convergence == 1
    filename = 'nominal_compliant_trajectory.mat';
    load(filename)    
end 

% Mass of M_L+M_P
%mass = [100:100:1500];
mass = [1110];

% Scales for spring constant
%spring_scale = [1/3, 2/3, 1 4/3, 5.5/3, 6/3, 7/3, 8/3, 9/3, 10/3];
spring_scale = 1;

% Penalization in cost function
%penalty_array = [1e-9 1e-8 1e-7 1e-6 1e-5 1e-4 1e-3 1e-2 1e-1 1];
penalty_array = [1e-9];

% save data
% z trajectory (rigid behavior)
save_z_data = 0;
% upward velocity (rigid or compliant)
save_y_velocity_data = 1; 
% save x trajectory (compliant behaivor)
save_x_data = 0;

% use rigid trajectory as nominal trajectory for compliant case
use_rigid_nominal = 0;
if use_rigid_nominal == 1
    % load reference data
    filename = 'nominal_rigid_trajectory.mat';
    load(filename)
end 

%% Eigenvalues
period_a_lin =zeros(1,length(mass));
period_a_1= zeros(1,length(mass));
optimal = zeros(1,length(mass));
frequency_a_lin = zeros(1,length(mass));
frequency_a_1 = zeros(1,length(mass));
period_continuous = zeros(1,length(mass));
frequency_continuous = zeros(1,length(mass));
if length(spring_scale) > 1
    k_values = zeros(1,length(spring_scale));
end 

%% Begin Evaluation
% Different tests can be run to evaluate performance based on varying
% m_L+m_P, spring constant, and penalty on the input. Only one test should
% be run at a time.

if length(mass)>1
    % test differemt m_L+m_P
    end_test = length(mass);
elseif length(spring_scale) > 1
    % test different spring constants
    end_test = length(spring_scale);
elseif length(penalty_array) > 1
    % test different penalties on the inputs u1, u2
    end_test = length(penalty_array);
else
    end_test = 1; % no test is being run
end 

for r = 1: end_test
    if length(mass) > 1 
        disp(mass(r))
    end
    % configure iterations
    if length(penalty_array) > 1
        iters = 20;
    else
        iters = 15;
    end 
    % configure N and step_size
    if utilize_compliance == 1 && input == 1 && Euler_method == 0
        N = 50;
        time_step = .0095;
    elseif utilize_compliance == 1 && input == 1 && Euler_method == 1
        N = 50;
        time_step = .006;
    elseif utilize_compliance == 1 && input == 0
        N = 5000; %5800;%5000;
        time_step = .0001;
    else
        N = 50;
        time_step = .0095;
    end 

    if length(mass)>1
        m_L1 = mass(r);
        m_L2 = mass(r);
    else
        m_L1 = mass(1);
        m_L2 = mass(1);
    end
    if length(spring_scale) > 1
        scale = spring_scale(r);
    else
        scale = spring_scale(1); 
    end
    k_S = 300000*scale; % N/m 
    disp(k_S)
    
    if length(penalty_array) > 1
        penalty = penalty_array(r);
        disp(penalty)
    else 
        penalty = penalty_array(1);
    end 
    if length(spring_scale) > 1
        k_values(r) = k_S;
    end 
    %% Model Parameters
    %Motor Subsystem
    %Spring Subsystem
    %Load Subsystem

    get_parameters(scale)

    % %% Frequency
    % mass_lumped1 = 1/((1/m_L1)+(1/m_M));
    % mass_lumped2 = 1/((1/m_L2)+(1/m_M));
    % omega1 = sqrt(k_S/mass_lumped1);
    % omega2 = sqrt(k_S/mass_lumped2);

    %% Continuous Time System

    [~,~,~,A_1,B_1] = get_continuous_matrices(m_L1,m_L2);
    eigenvalues = eig( A_1 );
    max_eig_a1 = max(abs(eigenvalues));

    %% Discrete time system
    [A,B] = get_discrete_matrices(A_1,B_1,1);
    if utilize_compliance == 0
        % configure discrete model for rigid system
        A_1_no_compliance = [0 1 0 0; 0 0 0 0; 0 0 0 1; 0 0 0 0];
        B_1_no_compliance = [0 0; 1 0; 0 0; 0 1];
        [A_no_compliance,B_no_compliance] = get_discrete_matrices(A_1_no_compliance,B_1_no_compliance,0);
        eigenvalues_no_compliance_a1 = eig( A_1_no_compliance );
        max_eig_a1_no_compliance = max(abs(eigenvalues_no_compliance_a1));
        eigenvalues_no_compliance_a = eig( A_no_compliance );
        max_eig_a_no_compliance = max(abs(eigenvalues_no_compliance_a));

    end 

    %% Matrices for Linearization
    s1v = [0 1 0 1 0 0 0 0];
    s2v = [0 0 0 0 0 1 0 1];
    s_velocity = [s1v;s2v];
    s1z = [1 0 1 0 0 0 0 0];
    s2z = [0 0 0 0 1 0 1 0];
    s_position = [s1z; s2z];

    % components of continuous matrix B1
    B1_1 = B_1(:,1);
    B1_2 = B_1(:,2);
    B1_3 = B_1(:,3);
    B1_4 = B_1(:,4);
    if Euler_method == 1
        B1_12 = B_1(:,1:2);
        B1_34 = B_1(:,3:4);
    end

    % components of discrete matrix B
    B_12 = B(:,1:2);
    B_34 = B(:,3:4);

    %% System Parameters                                        

    % center of mass position in y direction 
    my = 0.0; % meters

    % locations of end of link 1
    get_l1 = @(theta) [cos(theta)*L_1; sin(theta)*L_1];
    % locations of end of link 2
    get_l2 = @(theta1,theta2) [cos(theta1)*L_1+cos(theta1+theta2)*L_2; sin(theta1)*L_1+sin(theta1+theta2)*L_2];
    % COM locations
    get_mass1 = @(theta) [cos(theta)*L_c1; sin(theta)*L_c1];
    get_mass2 = @(theta1,theta2) [cos(theta1)*L_1+cos(theta1+theta2)*L_c2; sin(theta1)*L_1+sin(theta1+theta2)*L_c2];
    % actuator length 
    get_x1 = @(theta) (457/10000 - (21*cos(theta - atan(1/2)))/1250)^(1/2);
    get_x2 = @(theta) (213/5000 - (cos(pi/2 + asin(sin(pi/6 - acos(cos(asin((3*sin(theta))/(100*(9/5000 - (9*cos(theta))/5000)^(1/2))) + asin((1 - (390625*((9*cos(theta))/5000 + 7/5000)^2)/4)^(1/2)/(25*(9/5000 - (9*cos(theta))/5000)^(1/2))))))/(20*(17/400 - cos(pi/6 - acos(cos(asin((3*sin(theta))/(100*(9/5000 - (9*cos(theta))/5000)^(1/2))) + asin((1 - (390625*((9*cos(theta))/5000 + 7/5000)^2)/4)^(1/2)/(25*(9/5000 - (9*cos(theta))/5000)^(1/2))))))/50)^(1/2))))*(17/400 - cos(pi/6 - acos(cos(asin((3*sin(theta))/(100*(9/5000 - (9*cos(theta))/5000)^(1/2))) + asin((1 - (390625*((9*cos(theta))/5000 + 7/5000)^2)/4)^(1/2)/(25*(9/5000 - (9*cos(theta))/5000)^(1/2))))))/50)^(1/2))/50 - cos(pi/6 - acos(cos(asin((3*sin(theta))/(100*(9/5000 - (9*cos(theta))/5000)^(1/2))) + asin((1 - (390625*((9*cos(theta))/5000 + 7/5000)^2)/4)^(1/2)/(25*(9/5000 - (9*cos(theta))/5000)^(1/2))))))/50)^(1/2);
    % moment arm (Jacobian) 
    get_r1 = @(theta) (21*sin(theta - atan(1/2)))/(2500*(457/10000 - (21*cos(theta - atan(1/2)))/1250)^(1/2));
    get_r2 = @(theta) ((sin(pi/2 + asin(sin(pi/6 - acos(cos(asin((3*sin(theta))/(100*(9/5000 - (9*cos(theta))/5000)^(1/2))) + asin((1 - (390625*((9*cos(theta))/5000 + 7/5000)^2)/4)^(1/2)/(25*(9/5000 - (9*cos(theta))/5000)^(1/2))))))/(20*(17/400 - cos(pi/6 - acos(cos(asin((3*sin(theta))/(100*(9/5000 - (9*cos(theta))/5000)^(1/2))) + asin((1 - (390625*((9*cos(theta))/5000 + 7/5000)^2)/4)^(1/2)/(25*(9/5000 - (9*cos(theta))/5000)^(1/2))))))/50)^(1/2))))*(17/400 - cos(pi/6 - acos(cos(asin((3*sin(theta))/(100*(9/5000 - (9*cos(theta))/5000)^(1/2))) + asin((1 - (390625*((9*cos(theta))/5000 + 7/5000)^2)/4)^(1/2)/(25*(9/5000 - (9*cos(theta))/5000)^(1/2))))))/50)^(1/2)*((sin(asin((3*sin(theta))/(100*(9/5000 - (9*cos(theta))/5000)^(1/2))) + asin((1 - (390625*((9*cos(theta))/5000 + 7/5000)^2)/4)^(1/2)/(25*(9/5000 - (9*cos(theta))/5000)^(1/2))))*cos(pi/6 - acos(cos(asin((3*sin(theta))/(100*(9/5000 - (9*cos(theta))/5000)^(1/2))) + asin((1 - (390625*((9*cos(theta))/5000 + 7/5000)^2)/4)^(1/2)/(25*(9/5000 - (9*cos(theta))/5000)^(1/2))))))*(((9*sin(theta)*(1 - (390625*((9*cos(theta))/5000 + 7/5000)^2)/4)^(1/2))/(250000*(9/5000 - (9*cos(theta))/5000)^(3/2)) - (225*sin(theta)*((9*cos(theta))/5000 + 7/5000))/(32*(9/5000 - (9*cos(theta))/5000)^(1/2)*(1 - (390625*((9*cos(theta))/5000 + 7/5000)^2)/4)^(1/2)))/(1 - ((390625*((9*cos(theta))/5000 + 7/5000)^2)/4 - 1)/((9*cos(theta))/8 - 9/8))^(1/2) - ((3*cos(theta))/(100*(9/5000 - (9*cos(theta))/5000)^(1/2)) - (27*sin(theta)^2)/(1000000*(9/5000 - (9*cos(theta))/5000)^(3/2)))/((9*sin(theta)^2)/(18*cos(theta) - 18) + 1)^(1/2)))/(20*(17/400 - cos(pi/6 - acos(cos(asin((3*sin(theta))/(100*(9/5000 - (9*cos(theta))/5000)^(1/2))) + asin((1 - (390625*((9*cos(theta))/5000 + 7/5000)^2)/4)^(1/2)/(25*(9/5000 - (9*cos(theta))/5000)^(1/2))))))/50)^(1/2)*(1 - cos(asin((3*sin(theta))/(100*(9/5000 - (9*cos(theta))/5000)^(1/2))) + asin((1 - (390625*((9*cos(theta))/5000 + 7/5000)^2)/4)^(1/2)/(25*(9/5000 - (9*cos(theta))/5000)^(1/2))))^2)^(1/2)) - (sin(asin((3*sin(theta))/(100*(9/5000 - (9*cos(theta))/5000)^(1/2))) + asin((1 - (390625*((9*cos(theta))/5000 + 7/5000)^2)/4)^(1/2)/(25*(9/5000 - (9*cos(theta))/5000)^(1/2))))*sin(pi/6 - acos(cos(asin((3*sin(theta))/(100*(9/5000 - (9*cos(theta))/5000)^(1/2))) + asin((1 - (390625*((9*cos(theta))/5000 + 7/5000)^2)/4)^(1/2)/(25*(9/5000 - (9*cos(theta))/5000)^(1/2))))))^2*(((9*sin(theta)*(1 - (390625*((9*cos(theta))/5000 + 7/5000)^2)/4)^(1/2))/(250000*(9/5000 - (9*cos(theta))/5000)^(3/2)) - (225*sin(theta)*((9*cos(theta))/5000 + 7/5000))/(32*(9/5000 - (9*cos(theta))/5000)^(1/2)*(1 - (390625*((9*cos(theta))/5000 + 7/5000)^2)/4)^(1/2)))/(1 - ((390625*((9*cos(theta))/5000 + 7/5000)^2)/4 - 1)/((9*cos(theta))/8 - 9/8))^(1/2) - ((3*cos(theta))/(100*(9/5000 - (9*cos(theta))/5000)^(1/2)) - (27*sin(theta)^2)/(1000000*(9/5000 - (9*cos(theta))/5000)^(3/2)))/((9*sin(theta)^2)/(18*cos(theta) - 18) + 1)^(1/2)))/(2000*(17/400 - cos(pi/6 - acos(cos(asin((3*sin(theta))/(100*(9/5000 - (9*cos(theta))/5000)^(1/2))) + asin((1 - (390625*((9*cos(theta))/5000 + 7/5000)^2)/4)^(1/2)/(25*(9/5000 - (9*cos(theta))/5000)^(1/2))))))/50)^(3/2)*(1 - cos(asin((3*sin(theta))/(100*(9/5000 - (9*cos(theta))/5000)^(1/2))) + asin((1 - (390625*((9*cos(theta))/5000 + 7/5000)^2)/4)^(1/2)/(25*(9/5000 - (9*cos(theta))/5000)^(1/2))))^2)^(1/2))))/(50*(sin(pi/6 - acos(cos(asin((3*sin(theta))/(100*(9/5000 - (9*cos(theta))/5000)^(1/2))) + asin((1 - (390625*((9*cos(theta))/5000 + 7/5000)^2)/4)^(1/2)/(25*(9/5000 - (9*cos(theta))/5000)^(1/2))))))^2/(8*cos(pi/6 - acos(cos(asin((3*sin(theta))/(100*(9/5000 - (9*cos(theta))/5000)^(1/2))) + asin((1 - (390625*((9*cos(theta))/5000 + 7/5000)^2)/4)^(1/2)/(25*(9/5000 - (9*cos(theta))/5000)^(1/2)))))) - 17) + 1)^(1/2)) + (sin(asin((3*sin(theta))/(100*(9/5000 - (9*cos(theta))/5000)^(1/2))) + asin((1 - (390625*((9*cos(theta))/5000 + 7/5000)^2)/4)^(1/2)/(25*(9/5000 - (9*cos(theta))/5000)^(1/2))))*sin(pi/6 - acos(cos(asin((3*sin(theta))/(100*(9/5000 - (9*cos(theta))/5000)^(1/2))) + asin((1 - (390625*((9*cos(theta))/5000 + 7/5000)^2)/4)^(1/2)/(25*(9/5000 - (9*cos(theta))/5000)^(1/2))))))*(((9*sin(theta)*(1 - (390625*((9*cos(theta))/5000 + 7/5000)^2)/4)^(1/2))/(250000*(9/5000 - (9*cos(theta))/5000)^(3/2)) - (225*sin(theta)*((9*cos(theta))/5000 + 7/5000))/(32*(9/5000 - (9*cos(theta))/5000)^(1/2)*(1 - (390625*((9*cos(theta))/5000 + 7/5000)^2)/4)^(1/2)))/(1 - ((390625*((9*cos(theta))/5000 + 7/5000)^2)/4 - 1)/((9*cos(theta))/8 - 9/8))^(1/2) - ((3*cos(theta))/(100*(9/5000 - (9*cos(theta))/5000)^(1/2)) - (27*sin(theta)^2)/(1000000*(9/5000 - (9*cos(theta))/5000)^(3/2)))/((9*sin(theta)^2)/(18*cos(theta) - 18) + 1)^(1/2)))/(50*(1 - cos(asin((3*sin(theta))/(100*(9/5000 - (9*cos(theta))/5000)^(1/2))) + asin((1 - (390625*((9*cos(theta))/5000 + 7/5000)^2)/4)^(1/2)/(25*(9/5000 - (9*cos(theta))/5000)^(1/2))))^2)^(1/2)) - (cos(pi/2 + asin(sin(pi/6 - acos(cos(asin((3*sin(theta))/(100*(9/5000 - (9*cos(theta))/5000)^(1/2))) + asin((1 - (390625*((9*cos(theta))/5000 + 7/5000)^2)/4)^(1/2)/(25*(9/5000 - (9*cos(theta))/5000)^(1/2))))))/(20*(17/400 - cos(pi/6 - acos(cos(asin((3*sin(theta))/(100*(9/5000 - (9*cos(theta))/5000)^(1/2))) + asin((1 - (390625*((9*cos(theta))/5000 + 7/5000)^2)/4)^(1/2)/(25*(9/5000 - (9*cos(theta))/5000)^(1/2))))))/50)^(1/2))))*sin(asin((3*sin(theta))/(100*(9/5000 - (9*cos(theta))/5000)^(1/2))) + asin((1 - (390625*((9*cos(theta))/5000 + 7/5000)^2)/4)^(1/2)/(25*(9/5000 - (9*cos(theta))/5000)^(1/2))))*sin(pi/6 - acos(cos(asin((3*sin(theta))/(100*(9/5000 - (9*cos(theta))/5000)^(1/2))) + asin((1 - (390625*((9*cos(theta))/5000 + 7/5000)^2)/4)^(1/2)/(25*(9/5000 - (9*cos(theta))/5000)^(1/2))))))*(((9*sin(theta)*(1 - (390625*((9*cos(theta))/5000 + 7/5000)^2)/4)^(1/2))/(250000*(9/5000 - (9*cos(theta))/5000)^(3/2)) - (225*sin(theta)*((9*cos(theta))/5000 + 7/5000))/(32*(9/5000 - (9*cos(theta))/5000)^(1/2)*(1 - (390625*((9*cos(theta))/5000 + 7/5000)^2)/4)^(1/2)))/(1 - ((390625*((9*cos(theta))/5000 + 7/5000)^2)/4 - 1)/((9*cos(theta))/8 - 9/8))^(1/2) - ((3*cos(theta))/(100*(9/5000 - (9*cos(theta))/5000)^(1/2)) - (27*sin(theta)^2)/(1000000*(9/5000 - (9*cos(theta))/5000)^(3/2)))/((9*sin(theta)^2)/(18*cos(theta) - 18) + 1)^(1/2)))/(5000*(17/400 - cos(pi/6 - acos(cos(asin((3*sin(theta))/(100*(9/5000 - (9*cos(theta))/5000)^(1/2))) + asin((1 - (390625*((9*cos(theta))/5000 + 7/5000)^2)/4)^(1/2)/(25*(9/5000 - (9*cos(theta))/5000)^(1/2))))))/50)^(1/2)*(1 - cos(asin((3*sin(theta))/(100*(9/5000 - (9*cos(theta))/5000)^(1/2))) + asin((1 - (390625*((9*cos(theta))/5000 + 7/5000)^2)/4)^(1/2)/(25*(9/5000 - (9*cos(theta))/5000)^(1/2))))^2)^(1/2)))/(2*(213/5000 - (cos(pi/2 + asin(sin(pi/6 - acos(cos(asin((3*sin(theta))/(100*(9/5000 - (9*cos(theta))/5000)^(1/2))) + asin((1 - (390625*((9*cos(theta))/5000 + 7/5000)^2)/4)^(1/2)/(25*(9/5000 - (9*cos(theta))/5000)^(1/2))))))/(20*(17/400 - cos(pi/6 - acos(cos(asin((3*sin(theta))/(100*(9/5000 - (9*cos(theta))/5000)^(1/2))) + asin((1 - (390625*((9*cos(theta))/5000 + 7/5000)^2)/4)^(1/2)/(25*(9/5000 - (9*cos(theta))/5000)^(1/2))))))/50)^(1/2))))*(17/400 - cos(pi/6 - acos(cos(asin((3*sin(theta))/(100*(9/5000 - (9*cos(theta))/5000)^(1/2))) + asin((1 - (390625*((9*cos(theta))/5000 + 7/5000)^2)/4)^(1/2)/(25*(9/5000 - (9*cos(theta))/5000)^(1/2))))))/50)^(1/2))/50 - cos(pi/6 - acos(cos(asin((3*sin(theta))/(100*(9/5000 - (9*cos(theta))/5000)^(1/2))) + asin((1 - (390625*((9*cos(theta))/5000 + 7/5000)^2)/4)^(1/2)/(25*(9/5000 - (9*cos(theta))/5000)^(1/2))))))/50)^(1/2));
    % Actuator Location, for graphing
    get_A1 = @(theta,delta) [cos(theta)*(L_c1+delta); sin(theta)*(L_c1+delta)];
    get_A2 = @(theta1,theta2,delta) [cos(theta1)*L_1+cos(theta1+theta2)*(L_c2-delta); sin(theta1)*L_1+sin(theta1+theta2)*(L_c2-delta)];

    % actuator length limits, based on geometric limits
    limit1_lower = .171045; % theta = .67
    limit1_upper = .246427; % theta = pi
    limit2_lower = .15098814; % theta = pi
    limit2_upper = .24357; % theta = .01


    %% Simulation Parameters

    % optimal value array, in order to evaluate convergence with tolerance
    % \epsilon
    opt_val = zeros(iters,1);

    % states for all time
    x = zeros(N,8);

    if input == 1
        % starting possition (equals ending position)
        theta_init1 = 5*pi/8;
        theta_unshift_init2 = 3*pi/4; 
    else
        % test no-input behavior from nearly-vertical position
        theta_init1 = pi/2;
        theta_unshift_init2 = pi-.1;    
    end

    % moment arms
    r1 = zeros(N,1);
    r2 = zeros(N,1);

    % Linearization matrices
    A_lin = zeros(8,8,N);
    B_lin = zeros(8,2,N);
    bias_lin = zeros(8,1,N);

    if Euler_method == 1
        A_prime = zeros(8,8,N);
        B_prime = zeros(8,2,N);
        bias_prime = zeros(8,1,N);
    end 

    % generalized joint positions
    theta1 = zeros(N,1);
    theta2_unshift = zeros(N,1);
    theta2_shift = zeros(N,1);
    x(:,:)= ones(N,1)*[0,0,get_x1(theta_init1),0,0,0,get_x2(theta_unshift_init2),0];
    theta_dot1 = zeros(N,1);
    theta_dot2 = zeros(N,1);

    if utilize_compliance == 0
        % Matrices for rigid behavior evaluation
        Linv_M_Linv = zeros(2,2,N-1);
        F_b_components_no_compliance= zeros(2,1,N-1);
    end

    if analyze_convergence == 1
        % Used to evaluate error behavior as iterations increases
        error = zeros(1,iters);
    end 

    foot_length = .254; % meters, = 10 inches
    mu = .8; 
    ddx1_est = zeros(N-1,1);
    ddx2_est= zeros(N-1,1);

    %% Getting nominal trajectory      
    % x(+,:) = A*x(-,:) + B(:,1)*u(1) + B(:,2)*u(2) +
    %          + B(:,3)*u(3) + B(:,4)*u(4)          ;
    if input == 1
        % Nominal trajectory at equilibrium with springs 
        if use_rigid_nominal == 0
        x1_est = x(1,1)+x(1,3); % link 1 displacement
        dx1_est = x(1,2)+x(1,4); % link 1 velocity 
        x2_est = x(1,5)+x(1,7); % link 1 displacement
        dx2_est = x(1,6)+x(1,8); % link 1 velocity          
        % Note that theta2_unshift is with respect to link 1. The angle of link 2
        % with respect to the horizontal axis is theta1+theta2

        theta1 = get_theta1(x1_est);
        theta2_unshift = get_theta_2(x2_est);
        theta2_shift= theta2_unshift+pi;

        r1 = get_r1(theta1);
        r2 = get_r2(theta2_unshift);        
        theta_dot1 = dx1_est/r1;
        theta_dot2 = dx2_est/r2;
        
        % Lagrangian Dynamic
        [H11, H22, H12, H21, h, G1, G2] = get_dynamic_components(theta1,theta2_shift);

        [F_x_components, F_u_components, F_b_components] = ...
            get_F_components(H11, H22, H12, H21, h, G1, G2,...
            theta1, theta2_unshift, r1, r2, theta_dot1, theta_dot2,...
            A_1, B1_1, B1_2, B1_3, B1_4, s1v, s2v,m_L1,m_L2);

        [A_lin1(:,:),B_lin1(:,:),bias_lin1(:,:)] = linearize(A,B_12,B_34,F_x_components,...
            F_u_components,F_b_components);

        S = [(eye(8)-A_lin1), -B_lin1; 1 0 1 0 0 0 0 0 0 0; 0 0 0 0 1 0 1 0 0 0];
        z_init = [bias_lin1; get_x1(theta_init1); get_x2(theta_unshift_init2)];
        x_int = S\z_init;
        x0 = x_int(1:8,1)';
        u0 = x_int(9:10,1)';

        % configure x_baseline or z_baseline from nominal trajectory
        if utilize_compliance == 1 
            x(:,:)= ones(N,1)*x0;
            x_baseline = x;
        else
           z1 = x0(1)+x0(3);
           z1_dot = x0(2)+x0(4);
           z2 = x0(5)+x0(7);
           z2_dot = x0(6)+x0(8); 
           z0 = [z1 z1_dot z2 z2_dot];
           z(:,:)= ones(N,1)*z0;
           z_baseline = z;
        end

        % use the optimal trajectory for the rigid system as the nominal
        % trajectory for the compliant system
        elseif use_rigid_nominal == 1
            x(:,3) = zcvx(:,1);
            x(:,4) = zcvx(:,2);
            x(:,7) = zcvx(:,3);
            x(:,8) = zcvx(:,4); 
            x_baseline = x;
        end

        % Parameters to configure final upward velocity at N
        x1_final = get_x1(theta_init1);
        x2_final = get_x2(theta_unshift_init2);   
        theta_2_unshift_final = get_theta_2(x2_final);
        % Change in joint angle with respect to actuator length for ankle
        Tdot1_final = get_ankle_jacobian(x1_final);
        % Change in joint angle with respect to actuator length for knee
        [Tdot2_final, ~] = solve_knee_jacobians2(theta_2_unshift_final);
    end 
    %% Check: Initial Forces
    if input == 1
        wrench_nominal = zeros(3,N-1);
        checkY = zeros(N-1,1);
        checkT = zeros(N-1,1);
        for j = 1:N-1
                    ddx1_est(j) = (x(j+1,2)+x(j+1,4)-x(j,2)-x(j,4))/time_step;
                    ddx2_est(j) = (x(j+1,6)+x(j+1,8)-x(j,6)-x(j,8))/time_step;
                    [wrench_nominal(:,j),~,checkY(j),checkT(j)] = ContraintComponents(x(j:j+1,:),ddx1_est(j),ddx2_est(j),A_lin(:,:,j),bias_lin(:,:,j));;
        end
    end 
    %% Iterative Simulation
    Constraint_Jacobian = zeros(3,2,N-1);
    H11_graph = zeros(N,1);
    H22_graph = zeros(N,1);
    max_eig = 0; 
    if input == 0
        % no input for zero-input evaluation
        u = zeros(N-1,2);
    end 

    hold on
    for i = 1:iters
        for j = 1:N
            % x(+,:) = A*x(-,:) + B(:,1)*u(1) + B(:,2)*u(2) +
            %          + B(:,3)*u(3) + B(:,4)*u(4)          ;
            if utilize_compliance == 1 
                % states for compliant system
                x1_est = x(j,1)+x(j,3); % link 1 displacement
                dx1_est = x(j,2)+x(j,4); % link 1 velocity 
                x2_est = x(j,5)+x(j,7); % link 1 displacement
                dx2_est = x(j,6)+x(j,8); % link 1 velocity          
                % Note that theta2 is with respect to link 1. The angle of link 2
                % with respect to the horizontal axis is theta1+theta2
            else
                % states for rigid system
                x1_est = z(j,1); % link 1 displacement
                dx1_est = z(j,2); % link 1 velocity 
                x2_est = z(j,3); % link 2 displacement
                dx2_est = z(j,4); % link 2 velocity   
            end
            theta1(j) = get_theta1(x1_est);
            theta2_unshift(j) = get_theta_2(x2_est);
            % formulation for lagrangian dynamics
            theta2_shift(j)= theta2_unshift(j)+pi;

            % moment arm
            r1(j) = get_r1(theta1(j));
            r2(j) = get_r2(theta2_unshift(j));  
            % angular velocity
            theta_dot1(j) = dx1_est/r1(j);
            theta_dot2(j) = dx2_est/r2(j);
    
            % Lagrangian dynamics 
            [H11, H22, H12, H21, h, G1, G2] = get_dynamic_components(theta1(j),theta2_shift(j));
            H11_graph(j) = H11/r1(j)^2;
            H22_graph(j) = H22/(r2(j)^2);

            if j == 30
                stop = 0;
            end

            if utilize_compliance == 1
                % Compute F' for compliant case
                [F_x_components, F_u_components, F_b_components] = ...
                    get_F_components(H11, H22, H12, H21, h, G1, G2,...
                    theta1(j), theta2_unshift(j), r1(j), r2(j), theta_dot1(j), theta_dot2(j),...
                    A_1, B1_1, B1_2, B1_3, B1_4, s1v, s2v,m_L1,m_L2);
            end 
            if utilize_compliance == 1 && Euler_method == 0
                % Create linear approximation of system, using discrete A
                % and B
                [A_lin(:,:,j),B_lin(:,:,j),bias_lin(:,:,j)] = ...
                    linearize(A,B_12,B_34,F_x_components,F_u_components,F_b_components);
            elseif utilize_compliance == 1 && Euler_method == 1
                % Utilize Euler method

                % Notice that the change here compared to above is the use of
                % the A_1 and B_1 matrix as opposed to the A and B discrete
                % matrices 
                [A_prime(:,:,j),B_prime(:,:,j),bias_prime(:,:,j)] = ...
                    linearize(A_1,B1_12,B1_34,F_x_components,F_u_components,F_b_components);        
                A_lin(:,:,j) = eye(8,8)+time_step*A_prime(:,:,j) ;
                B_lin(:,:,j) = time_step*B_prime(:,:,j);
                bias_lin(:,:,j) = time_step*bias_prime(:,:,j);  
            else
                % no actuator dynamics considered
                [Linv_M_Linv(:,:,j), F_b_components_no_compliance(:,:,j)] = ...
                    get_F_components_no_compliance(H11, H22, H12, H21, h, G1, G2,...
                    theta1(j), theta2_unshift(j), r1(j), r2(j), theta_dot1(j), theta_dot2(j));
            end

            if j < N
                % Used to compute reaction forces at foot, bolted to ground
                ddx1_est(j) = (x(j+1,2)+x(j+1,4)-x(j,2)-x(j,4))/time_step;
                ddx2_est(j) = (x(j+1,6)+x(j+1,8)-x(j,6)-x(j,8))/time_step;
                [wrench_nominal(:,j), Constraint_Jacobian(:,:,j),~,~] = ContraintComponents(x(j:j+1,:),ddx1_est(j),ddx2_est(j),A_lin(:,:,j),bias_lin(:,:,j));
            end 
        end 

        % Obtain maximum eigenvalue
        for eig_check = 1:N
            eigenvalues = eig( A_lin(:,:,eig_check) );
            eigenvalues = log(eigenvalues); % take the log because A_lin is discrete 
            test = max(abs(eigenvalues));
            if test > max_eig
                max_eig = test;
            end 
        end 

        if input == 1
        % do convex optimization
            cvx_begin
        %    cvx_precision best
            if length(penalty_array) > 1
                % for penalty test, the default solver is SDPT3
                if penalty == 1e-2 % in the penalty test, the solver is sensitive at penalty = 1e-2
                    cvx_solver sedumi  
                end                          
            elseif length(penalty_array) == 1
                cvx_solver sedumi        
            end 

            if utilize_compliance == 1 
                variable x(N,8)
            else
                variable z(N,4)
                variable z_double_dot(N-1,2)
            end
            variable c_constraint(N-1,4) 
            variable u(N-1,2) 
            variable extracost
            variable F(2,N-1)
            variable xdd(2,N-1)
          %  function below varifies no motion
           %minimize 0.00001*sum( abs(x(:,1)+x(:,3)-x_baseline(:,1)-x_baseline(:,3) ))+0.00001*sum( abs(x(:,5)+x(:,7)-x_baseline(:,5)-x_baseline(:,7) ))
        %    minimize -(L_1*cos(theta_init1)*Tdot1_final*(x(N,2)+x(N,4))+L_2*cos(theta_init1+theta_init2+pi)*Tdot1_final*(x(N,2)+x(N,4))+L_2*cos(theta_init1+theta_init2+pi)*Tdot2_final*(x(N,6)+x(N,8)))+0.00001*sum( abs(x(:,1)+x(:,3)-x_baseline(:,1)-x_baseline(:,3) ))+0.00001*sum( abs(x(:,5)+x(:,7)-x_baseline(:,5)-x_baseline(:,7) ));%+(1e-8)*norm(c_constraint, 1)
            if utilize_compliance == 1 
                if compare_compliant_to_rigid == 1
                    minimize -(L_1*cos(theta_init1)*Tdot1_final*(x(N,2)+x(N,4))+...
                     L_2*cos(theta_init1+theta_unshift_init2+pi)*Tdot1_final*(x(N,2)+x(N,4))...
                     +L_2*cos(theta_init1+theta_unshift_init2+pi)*Tdot2_final*(x(N,6)+x(N,8)))...
                     +(1e-3)*sum(abs(u(:,1)))+(1e-3)*sum(abs(u(:,2)));
                else
                   minimize -(L_1*cos(theta_init1)*Tdot1_final*(x(N,2)+x(N,4))+...
                     L_2*cos(theta_init1+theta_unshift_init2+pi)*Tdot1_final*(x(N,2)+x(N,4))...
                     +L_2*cos(theta_init1+theta_unshift_init2+pi)*Tdot2_final*(x(N,6)+x(N,8)))...
                     +penalty*sum(abs(u(:,1)))+penalty*sum(abs(u(:,2))); 
                end
            else
                    minimize -(L_1*cos(theta_init1)*Tdot1_final*z(N,2)+...
                     L_2*cos(theta_init1+theta_unshift_init2+pi)*Tdot1_final*z(N,2)...
                     +L_2*cos(theta_init1+theta_unshift_init2+pi)*Tdot2_final*z(N,4))...
                     +(1e-3)*sum(abs(u(:,1)))+(1e-3)*sum(abs(u(:,2)));
            end  

            subject to

            if utilize_compliance == 0 % rigid behavior
                % No actuator dynamics 
                for n = 1:N-1
                    F(:,n) == Linv_M_Linv(:,:,n)*z_double_dot(n,:)'+...
                        +F_b_components_no_compliance(:,:,n)
                end
                    % Joint Limit Constraints
                    limit1_lower<= z(:,1) <= limit1_upper
                    limit2_lower <= z(:,3) <= limit2_upper
                    get_x1(0+.5) <= z(:,1) <= get_x1(pi-1.05)
                    get_x2(pi-.5) <= z(:,3) <= get_x2(.01+.5)
 

                    % Starting Point Constraint
                    z(1,:) == z_baseline(1,:);

                    % End Dispacement Constraints
                    z(N,1)==get_x1(theta_init1);
                    z(N,3)== get_x2(theta_unshift_init2);

                    % No x-velocity at last time step
                    (-L_1*sin(theta_init1)*Tdot1_final*(z(N,2))...
                        -L_2*sin(theta_init1+theta_unshift_init2+pi)*Tdot1_final*(z(N,2))...
                        -L_2*sin(theta_init1+theta_unshift_init2+pi)*Tdot2_final*(z(N,4)))...
                        == 0

                        for j = 1: N-1
        %                     z_double_dot(j,:)' == Linv_M_Linv(:,:,j)\F(:,j)...
        %                         -Linv_M_Linv(:,:,j)\F_b_components_no_compliance(:,:,j);
        %                     z(j+1, :) == z(j,:)*A_no_compliance'+...
        %                         z_double_dot(j,:)*B_no_compliance' 
                            z(j+1, :) == z(j,:)*A_no_compliance'+...
                                (B_no_compliance*(Linv_M_Linv(:,:,j)\...
                                (F(:,j)-F_b_components_no_compliance(:,:,j))))';

                        end 
                        
                        %Ballscrew Limit
                        abs(z(:,2)) <= .3; 
                        abs(z(:,4)) <= .3;
                        
                        % Additional u limitation
                         Km_M*u(:,:)' == F(:,:);


            else
                % Utilize the compliance

                % Spring Deflection Constraint
                -0.01*(1/scale)  <= x(:,1) <= 0.01*(1/scale) ; % spring deflection limit (meters) 
                -0.01*(1/scale)  <= x(:,5) <= 0.01*(1/scale); % spring deflection limit (meters)

                % Joint Limit Constraints
                limit1_lower<= x(:,1)+x(:,3) <= limit1_upper
                limit2_lower <= x(:,5)+x(:,7) <= limit2_upper
                if compare_compliant_to_rigid == 1
                    get_x1(0+.5) <= x(:,1)+x(:,3) <= get_x1(pi-1.05)
                    get_x2(pi-.5) <= x(:,5)+x(:,7) <= get_x2(.01+.5)
                elseif compare_compliant_to_rigid == 0
                    get_x1(0+.5) <= x(:,1)+x(:,3) <= get_x1(pi-.5)
                    get_x2(pi-.5) <= x(:,5)+x(:,7) <= get_x2(.01+.5)         
                end
                extracost >= 0;
                % Starting Point Constraint
                x(1,:) == x_baseline(1,:);

                % End Dispacement Constraints
                x(N,1)+x(N,3)==get_x1(theta_init1);
                x(N,5)+x(N,7)== get_x2(theta_unshift_init2);

                % Baseline Deviation Constraints
                abs(x(:,1)+x(:,3)-x_baseline(:,1)-x_baseline(:,3)) <= .1*ones(N,1); 
                abs(x(:,5)+x(:,7)-x_baseline(:,5)-x_baseline(:,7)) <= .1*ones(N,1);


                % No x-velocity at last time step
                (-L_1*sin(theta_init1)*Tdot1_final*(x(N,2)+x(N,4))...
                    -L_2*sin(theta_init1+theta_unshift_init2+pi)*Tdot1_final*(x(N,2)+x(N,4))...
                    -L_2*sin(theta_init1+theta_unshift_init2+pi)*Tdot2_final*(x(N,6)+x(N,8)))...
                    == 0

                % Ballscrew Limit
                abs(x(:,4)) <= .3; 
                abs(x(:,8)) <= .3; 



                % Reaction Force Constraints for foot bolted to ground
                for n = 1:N-1  
                    c_constraint(n,:)*[mu 1 foot_length; -mu 1 foot_length; mu 1 0; -mu 1 0] == (...
                    wrench_nominal(:,n)'...
                    +(Constraint_Jacobian(:,:,n)*[0,1,0,1,0,0,0,0;0,0,0,0,0,1,0,1]*...
                    (B_lin(:,:,n)*u(n,:)')*(1/time_step))');       
        %          c_constraint(n,1)>=0;
        %          c_constraint(n,2)>=0;
        %          c_constraint(n,3)>=0;
        %          c_constraint(n,4)>=0;         

                end 

                % Linearization Constraint
                for j = 1: N-1
                    x(j+1, :) == x(j,:)*A_lin(:,:,j)'+...
                        u(j,:)*B_lin(:,:,j)'+...
                        bias_lin(:,:,j)';
                end 
                

            end

                % Current Constraint, applicable for both situations
                -15.0 <= u(:,:) <= 15.0;

            cvx_end
            % update baseline trajectory
            if utilize_compliance == 1 
                x_baseline = x; 
            else
                z_baseline = z;
            end
            
            % If we are evaluating zero-input behavior, cvx is not
            % performed 
        elseif input == 0
            for j = 1: N-1
                x(j+1, :) = x(j,:)*A_lin(:,:,j)'+...
                    u(j,:)*B_lin(:,:,j)'+...
                    bias_lin(:,:,j)';
            end               

        end 
        var = num2str(i);

        if utilize_compliance == 1
            %plot(real(x(:,1))+real(x(:,3)),'DisplayName',var)
            %plot(real(x(:,2))+real(x(:,4)),'DisplayName',var) 
            %plot(real(x(:,5))+real(x(:,7)),'DisplayName',var) 
            %plot(real(x(:,6))+real(x(:,8)),'DisplayName',var) 
            plot(real(x(:,1))+real(x(:,3)), real(x(:,5))+real(x(:,7)),'DisplayName',var,'Linewidth',1)   
        else
            plot(z(:,1),z(:,3),'DisplayName',var,'Linewidth',1); 
        end
        if input == 1
            opt_val(i) = cvx_optval;
            if isequal(cvx_status, 'Solved') 
                solve = 1; 
            else 
                solve = 0;
            end 
        end 
        
        % Stop iterations if convergence occurs within specified tolerance
        if input == 1
            if i > 2
                if abs(opt_val(i) - opt_val(i-1)) < .001 && solve == 1
                    break
                end
            end 
        end 
        
        % Error analysis
        if analyze_convergence == 1
            error(i) = norm(x_compare_convergence-x)/norm(x_compare_convergence);
        end

    end

    %% Plots 

    % Used to track eigenvalues
    if utilize_compliance == 1
        % Evauluate eigenvalues
        maximum = max_eig;
        period_a_lin(r) =1/((max_eig/time_step)/(2*pi));
        period_a_1(r)=1/((max_eig_a1)/(2*pi)); % recall that A1 is a continuous matrix,
                                               % and hence time_step is not
                                               % inclued 
        frequency_a_lin(r) = max_eig/time_step;
        frequency_a_1(r) = max_eig_a1;
    end
    if input == 1 && length(penalty_array) == 1
        optimal(r) = cvx_optval;   
%      if length(spring_scale) > 1 && input == 1 && length(penalty_array) == 1
%         optimal(r) = cvx_optval; 
%     elseif length(spring_scale) == 1 && input == 1 && length(penalty_array) == 1
%             optimal(r) = cvx_optval;           
    elseif input == 1 && length(penalty_array) > 1
        optimal(r) = (L_1*cos(theta_init1)*Tdot1_final*(x(N,2)+x(N,4))+...
             L_2*cos(theta_init1+theta_unshift_init2+pi)*Tdot1_final*(x(N,2)+x(N,4))...
             +L_2*cos(theta_init1+theta_unshift_init2+pi)*Tdot2_final*(x(N,6)+x(N,8)));
    end 
    % Period and natural frequency for the continuous system. Used to
    % compare against A_1 and A_lin
    [period_continuous(r), frequency_continuous(r)] = Nonlinear_Dynamics(x, N,m_L1,m_L2);   

    % Store our trajectories 
    if utilize_compliance == 1
        xcvx = x;         
    elseif utilize_compliance == 0 && input == 1
        zcvx = z;
        if save_z_data == 1
            filename = 'nominal_rigid_trajectory.mat'; 
            save(filename,'zcvx')             
        end

    end 

    % Save x data. Typically the data is used for analyzing convergence 
    if save_x_data == 1
        x_compare_convergence = x;
        filename = 'nominal_compliant_trajectory.mat';
        save(filename,'x_compare_convergence')   
    end


    % Compute the upward velocities in order to compare rigid to compliant
    % behavior
    if save_y_velocity_data == 1
        Tdot1_array = zeros(N,1);
        Tdot2_array = zeros(N,1);   
        if utilize_compliance == 1
            z_1_evaluate_y_velocity = x(:,1)+x(:,3);
            z_dot_1_evaluate_y_velocity = x(:,2)+x(:,4); 
            z_dot_2_evaluate_y_velocity = x(:,6)+x(:,8);
            filename = 'yVelocityCompliant.mat';
        elseif utilize_compliance == 0
            z_1_evaluate_y_velocity = z(:,1);
            z_dot_1_evaluate_y_velocity = z(:,2);
            z_dot_2_evaluate_y_velocity = z(:,4);
            filename = 'yVelocityRigid.mat';
        end 
        for w = 1:N
            Tdot1_array(w) = get_ankle_jacobian(z_1_evaluate_y_velocity(w));
            Tdot2_array(w) = solve_knee_jacobians2(theta2_unshift(w)) ;       
        end 
        if utilize_compliance == 1
            y_velocity_compliant = get_y_velocity(theta1(:,:)',theta2_unshift(:,:)',...
                z_dot_1_evaluate_y_velocity(:,:),z_dot_2_evaluate_y_velocity(:,:),...
                L_1,L_2,Tdot1_array(:,:),Tdot2_array(:,:));
            save(filename,'y_velocity_compliant')
        elseif utilize_compliance == 0
             y_velocity_rigid = get_y_velocity(theta1(:,:)',theta2_unshift(:,:)',...
                z_dot_1_evaluate_y_velocity(:,:),z_dot_2_evaluate_y_velocity(:,:),...
                L_1,L_2,Tdot1_array(:,:),Tdot2_array(:,:));
            save(filename,'y_velocity_rigid');     
        end 

    end 
    % Store our input trajectory
    ucvx = u;
end
legend('show')
if input == 1
    title('Trajectories for Jumping Behavior')
else
    title('Trajectories for Zero Input Behavior')
    axis([0.19,.3,.14,.26])
end 
xlabel('z_1 (m)')
ylabel('z_2 (m)')
set(gca,'Fontsize',18)
set(gca,'fontname','times')
grid on


% Plot eigenvalues for varying ficticious mass values 
if utilize_compliance == 1 && Euler_method == 0 && length(mass) > 1
    figure()
    plot(mass, frequency_continuous, 'DisplayName', 'Continuous System','Linewidth',1.5)
    hold on 
    plot(mass, frequency_a_1, 'DisplayName', 'A_1','Linewidth',1.5)
    plot(mass, frequency_a_lin,'k', 'DisplayName', 'A_{lin}','Linewidth',1.5)
    legend('show')
    xlabel('M_p (kg)')
    ylabel('Natural Frequency (Hz)')
    title('Maximum Eigenvalue for Varying M_p')
    set(gca,'Fontsize',18)
    set(gca,'fontname','times')
    grid on

    figure()
    plot(optimal,frequency_continuous)
    xlabel('frequency of continuous system')
    ylabel('optimal velocity (m/s)')
end 

% Plot error vs. number of iterations
if analyze_convergence == 1
    figure()
    semilogy(1:iters,error,'Linewidth',1.5)
    xlabel('Iteration')
    ylabel('Error')
    title('Rate of Convergence')
    set(gca,'Fontsize',18)
    set(gca,'fontname','times')
    grid on
end 


% Spring Stiffness Comparison
if utilize_compliance == 1 && Euler_method == 0 && length(spring_scale) > 1
    figure()
    plot(k_values, -optimal)
    xlabel('Spring Constant (N/m)')
    ylabel('optimal velocity (m/s)')    
end 

% Penalty comparison
if length(penalty_array) > 1
    figure()
    semilogx(penalty_array, optimal,'Linewidth',1.5)
    xlabel('Penalty on u_1 and u_2')
    ylabel('Optimal Velocity (m/s)') 
    title('Optimal Velocity Based on Input Penalization')
    set(gca,'Fontsize',18)
    set(gca,'fontname','times')
    grid on
end 


% figure()
% title('comparison with mL')
% plot(H11_graph, 'DisplayName', 'M11')
% hold on
% plot(H22_graph, 'DisplayName', 'M22')
% hold off
% legend('show')

% Plot dispacements and velocities over time 
figure()
plot([0:N-1]*time_step, real(x(:,1))+real(x(:,3)),'DisplayName','x1')
hold on 
plot([0:N-1]*time_step, real(x(:,2))+real(x(:,4)),'DisplayName','x1dot')
plot([0:N-1]*time_step, real(x(:,5))+real(x(:,7)),'DisplayName','x2')
plot([0:N-1]*time_step, real(x(:,6))+real(x(:,8)),'DisplayName','x2dot')
legend('show')
xlabel('time (s)')
ylabel('trajectory (m)')
hold off 

% Plot Spring Deflection over time 
figure()
plot([0:N-1]*time_step,real(x(:,1)*1000),'Linewidth',1.5,'DisplayName','\delta_1')
hold on
plot([0:N-1]*time_step,real(x(:,5)*1000),'Linewidth',1.5,'DisplayName','\delta_2')
legend('show')
xlabel('Time (s)')
ylabel('Spring Deflection (mm)')
title('Spring Deflection Trajectories')
set(gca,'Fontsize',18)
set(gca,'fontname','times')
grid on
hold off
axis([0,.5,-10,11])

% Plot u's 
figure()
plot([0:N-2]*time_step,ucvx(:,1),'Linewidth',1.5,'DisplayName','u_1')
hold on
plot([0:N-2]*time_step,ucvx(:,2),'Linewidth',1.5,'DisplayName','u_2')
legend('show')
axis([0,.5,-16,16])
xlabel('Time (s)','Interpreter','latex')
ylabel('Input Current (A)')
title('Input Trajectories')
set(gca,'Fontsize',18)
set(gca,'fontname','times')
grid on
hold off 

%% Checks: Kinetic and Potential Energy, and Torque
if utilize_compliance == 1
    T = zeros(N-1,1);
    U = zeros(N-1,1);
    diss = zeros(N-1,1);
    Torque1 = zeros(N,1);
    Torque2 = zeros(N,1);
    for j = 1:N-1

        % kinetic energy and potential energy
         [T(j),U(j),diss(j),Torque1(j),Torque2(j)] = checks(x(j:j+1,:),m_L1,m_L2);
    end 

    % Total Mechanical Energy
    figure()
    plot([1:N-1]'*time_step,real(T),'DisplayName','kinetic energy');
    hold on
    plot([1:N-1]'*time_step,real(U),'DisplayName','potential energy');
    plot([1:N-1]'*time_step,real(T)+real(U)+real(diss),'DisplayName','total energy');
    legend('show')
    xlabel('time (s)')
    ylabel('Energy (J)')
    hold off

    % Torques
    figure()
    plot([1:N]'*time_step, Torque1,'DisplayName','torque1')
    hold on
    title('torque 1 verification')
    xlabel('time (s)')
    ylabel('torque1')
    legend('show')
    hold off

    figure()
    plot([1:N]'*time_step, Torque2,'DisplayName','torque2')
    hold on
    title('torque 2 verification')
    xlabel('time (s)')
    ylabel('torque2')
    legend('show')
    hold off
end
%% Animation
figure()
clf()
v = VideoWriter('new');
% v.framerate=60;
v.FrameRate = 1/time_step; 
v.Quality=98;
disp(v.VideoCompressionMethod)
open(v)
for i=1:N
    % visualization in a plot
    clf()
    axs=plot([0,0],[0,0]);
    if utilize_compliance == 1 
    %delta = [1,0,0,0,0,0,0,0]*xcvx(i,:)';
        % link 1
        x1= [1,0,1,0,0,0,0,0]*xcvx(i,:)';
        % link 2
        x2= [0,0,0,0,1,0,1,0]*xcvx(i,:)';    
    elseif utilize_compliance == 0
        % link 1
        x1= [1,0,0,0]*zcvx(i,:)';
        % link 2
        x2= [0,0,1,0]*zcvx(i,:)'; 
    end
    theta1=get_theta1(x1); 
    theta2_unshift= get_theta_2(x2); 
    theta2_shift = theta2_unshift+pi; 


    l1 = get_l1(theta1);
    l2 = get_l2(theta1,theta2_shift);
    m1 = get_mass1(theta1);
    m2 = get_mass2(theta1,theta2_shift);

    plot([0, l1(1), l2(1)],[0,l1(2),l2(2)],'o-');


    hold on
    plot([-foot_length,0],[0,0])
    hold on
    plot([m1(1)],[m1(2)],'o','MarkerSize',20)
    plot([m2(1)],[m2(2)],'o','MarkerSize',20)
    
    if utilize_compliance == 1 && input == 1
    % plot reaction forces 
%         if i < N
%             reaction = c_constraint(i,:);
% 
%             point = [-foot_length;0];
%             R1 = .005*(reaction(1,1)*[mu;1]+reaction(1,2)*[-mu;1]);
%             plot([point(1),point(1)+R1(1)],[point(2),point(2)+R1(2)])
%             R2 = .005*(reaction(1,3)*[mu;1]+reaction(1,4)*[-mu;1]);
%             plotv(1*R2,'-')
% 
%         end     
    end
    
    if utilize_compliance == 1
        delta1 = [1,0,0,0,0,0,0,0]*xcvx(i,:)';
        delta2 = [0,0,0,0,1,0,0,0]*xcvx(i,:)';
        springscale = 50; 
        points1 = get_A1(theta1,delta1*springscale);
        points2 = get_A2(theta1,theta2_shift,delta1*springscale);
        plot(points1(1,:), points1(2,:),'o','MarkerSize',10);
        plot(points2(1,:), points2(2,:),'o','MarkerSize',10);
    end 
    axis('equal')

    axis([-.5,.4,-.1,1.3])
    drawnow
%     F(i)=getframe;
    writeVideo(v,getframe)
    %pause(0.001)
end
close(v)