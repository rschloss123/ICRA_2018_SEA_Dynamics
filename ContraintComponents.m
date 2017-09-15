% Used to compute reaction forces on foot,
function [wrench_nominal,Constraint_Jacobian,checkY,checkT] = ContraintComponents(x,ddx1_est,ddx2_est, A_lin,bias_lin)
    global m_1 m_2 L_1 L_c1 L_c2 time_step
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
    if Tdouble_dot2 > 1000
        stop = 0;
        
    end 
    
    theta_2_shift = theta_2_unshift+pi;
    % compute theta_dot and theta_double_dot using transmission component
    theta_dot2 = Tdot2*x_dot2;
    theta_double_dot2 = Tdouble_dot2*x_dot2^2+Tdot2*xdouble_dot2; 
    % rigid body dynamics
    [H11, H22, H12, H21, h, G1, G2] = get_dynamic_components(theta_1,theta_2_shift);

    % Link 1 position, velocity, acceleration
    Pm1 = [L_c1*cos(theta_1);...
           L_c1*sin(theta_1)];
    % Vm1 = [-L_c1*sin(theta_1)*theta_dot1;...
    %        L_c1*theta_dot1*cos(theta_1)];
    Am1 = [-L_c1*cos(theta_1)*theta_dot1^2-L_c1*sin(theta_1)*theta_doubledot1;
           -L_c1*(theta_dot1^2)*sin(theta_1)+L_c1*theta_doubledot1*cos(theta_1)];
    
    % Link 2 position, velocity, acceleration

    Pm2 = [L_1*cos(theta_1)+L_c2*cos(theta_1+theta_2_shift);
           L_1*sin(theta_1)+L_c2*sin(theta_1+theta_2_shift)];

    % Vm2 = [-L_1*sin(theta_1)*theta_dot1-L_c2*sin(theta_1+theta_2)*theta_dot1-L_c2*sin(theta_1+theta_2)*theta_dot2;
    %        L_1*cos(theta_1)*theta_dot1+L_c2*cos(theta_1+theta_2)*theta_dot1+L_c2*cos(theta_1+theta_2)*theta_dot2];

    Am2 = [-L_1*cos(theta_1)*theta_dot1^2-L_1*sin(theta_1)*theta_doubledot1-L_c2*cos(theta_1+theta_2_shift)*theta_dot1^2-L_c2*sin(theta_1+theta_2_shift)*theta_doubledot1-2*L_c2*cos(theta_1+theta_2_shift)*theta_dot1*theta_dot2-L_c2*cos(theta_1+theta_2_shift)*theta_dot2^2-L_c2*sin(theta_1+theta_2_shift)*theta_double_dot2;...
           -L_1*sin(theta_1)*theta_dot1^2+L_1*cos(theta_1)*theta_doubledot1-L_c2*sin(theta_1+theta_2_shift)*theta_dot1^2+L_c2*cos(theta_1+theta_2_shift)*theta_doubledot1-2*L_c2*sin(theta_1+theta_2_shift)*theta_dot1*theta_dot2-L_c2*sin(theta_1+theta_2_shift)*theta_dot2^2+L_c2*cos(theta_1+theta_2_shift)*theta_double_dot2];

    % Torque produced by actuators     
    T1 = H11*theta_doubledot1+H12*theta_double_dot2-h*theta_dot2^2-2*h*theta_dot1*theta_dot2+G1;
    %T2 = H22*theta_double_dot2+H21*theta_doubledot1+h*theta_dot1^2+G2;

    % y force balance
    checkY = m_1*9.81+m_2*9.81;
    % torque balance 
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
    
    % derivative of theta_double_dot1 with respect to x_double_dot1
    Trans1 = Tdot1;
    % derivative of theta_double_dot1 with respect to x_double_dot1
    Trans2 = Tdot2; 
    
    Constraint_Jacobian = ...
              [dFx_dalpha1, dFx_dalpha2;...
               dFy_dalpha1, dFy_dalpha2;...
               dT_dalpha1, dT_dalpha2]...
              *[Trans1 0; 0 Trans2];

    acceleration_adjust = [ddx1_est; ddx2_est];
    x_nominal = m_1*Am1(1)+m_2*Am2(1);
    y_nominal = m_1*9.81+m_2*9.81+m_1*Am1(2)+m_2*Am2(2);
    torque_nominal = -m_1*9.81*L_c1*cos(theta_1)-m_2*9.81*(L_1*cos(theta_1)+L_c2*cos(theta_1+theta_2_shift))-T1;
    adjust_nominal_acceleration = ...
        Constraint_Jacobian*(acceleration_adjust-[0,1,0,1,0,0,0,0;0,0,0,0,0,1,0,1]*(bias_lin+(A_lin-eye(8,8))*x(1,:)')/time_step);
    wrench_nominal = [x_nominal; y_nominal; torque_nominal]-adjust_nominal_acceleration;%+Constraint_Jacobian*[0,1,0,1,0,0,0,0;0,0,0,0,0,1,0,1]*((A_lin-eye(8,8))*x(2,:)'+bias_lin)*1/time_step;      

end

       