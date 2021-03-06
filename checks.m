% Compute torques, kinetic energy, potential energy, and energy dissipated
% by dampers
function [T, U, diss, Torque1, Torque2] = checks(x)
    global m_1 m_2 L_1 L_c1 L_c2 time_step k_S m_M m_S b_S b_M work  
    
    % angular position, velocity, and acceleration of the two links
    [theta_1, theta_dot1,theta_doubledot1, ~,theta_2_shift, theta_dot2,...
        theta_double_dot2, ~,~,~,~,~,~] = get_pos_vel_acc(x);

    % rigid body dynamics
    [H11, H22, H12, H21, h, G1, G2] = get_dynamic_components(theta_1,theta_2_shift);
    
    % Actuator torques 
    [Torque1, Torque2] = get_torques(H11, H22, H12, H21, h, G1, G2,...
    theta_doubledot1, theta_double_dot2, theta_dot2, theta_dot1);
  
    % Accumulation of work over time
    work = work + (b_S*(x(2,2))^2+ b_M*(x(2,4))^2 + b_S*(x(2,6))^2 + b_M*(x(2,8))^2)*time_step;
    % Kinetic energy of actuator from motors and springs
    load_components = .5 * m_M * x(1,4).^2 + 0.5 * m_S * x(1,2).^2 + .5 * m_M * x(1,8).^2 + 0.5 * m_S * x(1,6).^2;
    % Kinetic energy of robot + KE of actuator
    T = (1/2)*H11*theta_dot1^2+H12*theta_dot1*theta_dot2+(1/2)*H22*theta_dot2^2+load_components;
    % Potential energy of robot and spring
    U =  m_1*9.81*L_c1*sin(theta_1)+m_2*9.81*(L_1*sin(theta_1)+L_c2*sin(theta_1+theta_2_shift))+.5*k_S*(x(2,1))^2;     
    diss = work;

end 

    

