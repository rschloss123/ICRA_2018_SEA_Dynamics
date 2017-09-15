% Compute torques, kinetic energy, potential energy, and energy dissipated
% by dampers
function [T, U, diss, Torque1, Torque2] = checks(x,m_L1,m_L2)
   global m_1 m_2 L_1 L_c1 L_c2 time_step k_S m_M m_S b_S b_M work
    % calculate velocity and acceleration of actuator outputs
    % actuator 1
    x_1 = x(1,1)+x(1,3);
    x_dot1 = (x(2,1)+x(2,3)-x(1,1)-x(1,3))/time_step;
    xdouble_dot1 = (x(2,2)+x(2,4)-x(1,2)-x(1,4))/time_step; 
    % derivative of transmission = dtheta/dx
    Tdot1 = get_ankle_jacobian(x_1);
    % (d^2)theta/d(x^2)
    Tddot1 = 2500/(21*(1 - ((1250*x_1^2)/21 - 457/168)^2)^(1/2)) + (6250000*x_1^2*((1250*x_1^2)/21 - 457/168))/(441*(1 - ((1250*x_1^2)/21 - 457/168)^2)^(3/2));
    % get theta
    theta_1 = get_theta1(x_1);
    % compute theta_dot and theta_ddot using transmission component
    theta_dot1 = Tdot1*x_dot1;
    theta_doubledot1 = Tddot1*x_dot1^2+Tdot1*xdouble_dot1; 

    % actuator 2
    x_2 = x(1,5)+x(1,7);
    x_dot2 = (x(2,5)+x(2,7)-x(1,5)-x(1,7))/time_step;
    xdouble_dot2 = (x(2,6)+x(2,8)-x(1,6)-x(1,8))/time_step; 
    theta_2_unshift = get_theta_2(x_2);
    [Tdot2, Tdouble_dot2] = solve_knee_jacobians2(theta_2_unshift);

    
    theta_2_shift = theta_2_unshift+pi;
    % compute theta_dot and theta_double_dot using transmission component
    theta_dot2 = Tdot2*x_dot2;
    theta_double_dot2 = Tdouble_dot2*x_dot2^2+Tdot2*xdouble_dot2; 
    % rigid body dynamics
    [H11, H22, H12, H21, h, G1, G2] = get_dynamic_components(theta_1,theta_2_shift);

    % Torque produced by actuators     
    T1 = H11*theta_doubledot1+H12*theta_double_dot2-h*theta_dot2^2-2*h*theta_dot1*theta_dot2+G1;
    T2 = H22*theta_double_dot2+H21*theta_doubledot1+h*theta_dot1^2+G2;
    
    % Accumulation of work over time
    work = work + (b_S*(x(2,2))^2+ b_M*(x(2,4))^2 + b_S*(x(2,6))^2 + b_M*(x(2,8))^2)*time_step;
    % Kinetic energy of actuator from motors and springs
    load_components = .5 * m_M * x(1,4).^2 + 0.5 * m_S * x(1,2).^2 + .5 * m_M * x(1,8).^2 + 0.5 * m_S * x(1,6).^2;
    % Kinetic energy of robot + KE of actuator
    T = (1/2)*H11*theta_dot1^2+H12*theta_dot1*theta_dot2+(1/2)*H22*theta_dot2^2+load_components;
    % Potential energy of robot and spring
    U =  m_1*9.81*L_c1*sin(theta_1)+m_2*9.81*(L_1*sin(theta_1)+L_c2*sin(theta_1+theta_2_shift))+.5*k_S*(x(2,1))^2;     
    diss = work;
    Torque1 = T1;
    Torque2 = T2;

end 

    

