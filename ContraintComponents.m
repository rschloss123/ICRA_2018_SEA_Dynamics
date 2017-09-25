% Used to compute reaction forces on foot, which is bolted to the ground
function [wrench_nominal,Constraint_Jacobian,checkY,checkT] = ContraintComponents(x,ddx1_est,ddx2_est, A_lin,bias_lin)
    global m_1 m_2 L_1 L_c1 L_c2 time_step

    % angular position, velocity, and acceleration of the two links
    [theta_1, theta_dot1,theta_doubledot1, ~,theta_2_shift, theta_dot2,...
        theta_double_dot2, Am1,Am2,Tdot1,~,Tdot2,~] =...
        get_pos_vel_acc(x);
    
    % rigid body dynamics
    [H11, H22, H12, H21, h, G1, G2] = get_dynamic_components(theta_1,theta_2_shift);

    % Actuator torques 
    [T1,~] = get_torques(H11, H22, H12, H21, h, G1, G2,...
    theta_doubledot1, theta_double_dot2, theta_dot2, theta_dot1);
  
    % check: y force balance
    checkY = m_1*9.81+m_2*9.81;
    % check: torque balance 
    checkT = - m_1*9.81*L_c1*cos(theta_1)-m_2*9.81*(L_1*cos(theta_1)+L_c2*cos(theta_1+theta_2_shift))-G1;

    
    % derivative of F_x with respect to theta_double-dot1
    dFx_dalpha1 = -m_1*L_c1*sin(theta_1)-m_2*L_1*sin(theta_1)-m_2*L_c2*sin(theta_1+theta_2_shift);
    % derivative of F_x with respect to theta_double-dot2
    dFx_dalpha2 = -m_2*L_c2*sin(theta_1+theta_2_shift); 
    % derivative of F_y with respect to theta_double-dot1
    dFy_dalpha1 = m_1*L_c1*cos(theta_1)+m_2*L_1*cos(theta_1)+m_2*L_c2*cos(theta_1+theta_2_shift);
    % derivative of F_y with respect to theta_double-dot2
    dFy_dalpha2 = m_2*L_c2*cos(theta_1+theta_2_shift);
    % derivative of torque with respect to theta_double-dot1
    dT_dalpha1 = -H11;
    % derivative of torque with respect to theta_double-dot2
    dT_dalpha2 = -H12;
    
    Constraint_Jacobian = ...
              [dFx_dalpha1, dFx_dalpha2;...
               dFy_dalpha1, dFy_dalpha2;...
               dT_dalpha1, dT_dalpha2]...
              *[Tdot1 0; 0 Tdot2];

    acceleration_adjust = [ddx1_est; ddx2_est];
    x_nominal = m_1*Am1(1)+m_2*Am2(1);
    y_nominal = m_1*9.81+m_2*9.81+m_1*Am1(2)+m_2*Am2(2);
    torque_nominal = -m_1*9.81*L_c1*cos(theta_1)-m_2*9.81*(L_1*cos(theta_1)+L_c2*cos(theta_1+theta_2_shift))-T1;
    adjust_nominal_acceleration = ...
        Constraint_Jacobian*(acceleration_adjust-[0,1,0,1,0,0,0,0;0,0,0,0,0,1,0,1]*(bias_lin+(A_lin-eye(8,8))*x(1,:)')/time_step);
    wrench_nominal = [x_nominal; y_nominal; torque_nominal]-adjust_nominal_acceleration;     

end

       