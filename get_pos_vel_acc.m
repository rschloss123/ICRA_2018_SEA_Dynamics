% Position, velocity, acceleration, and transmission components for links 1
% and 2
function [theta_1, theta_dot1,theta_doubledot1, theta_2_unshift, theta_2_shift,theta_dot2,...
    theta_double_dot2, Am1, Am2, Tdot1, Tddot1, Tdot2, Tdouble_dot2] =...
    get_pos_vel_acc(x)

    global time_step L_1 L_c1 L_c2
    % calculate velocity and acceleration of actuator outputs
    % actuator 1
    x_1 = x(1,1)+x(1,3);
    x_dot1 = (x(2,1)+x(2,3)-x(1,1)-x(1,3))/time_step;
    xdouble_dot1 = (x(2,2)+x(2,4)-x(1,2)-x(1,4))/time_step; 
    % derivative of transmission = dtheta/dx
    % Tddot1 = % (d^2)theta/d(x^2)
    [Tdot1, Tddot1] = get_ankle_jacobian(x_1);
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
    % absolute angle for theta_2
    theta_2_shift = theta_2_unshift+pi;

    % compute theta_dot and theta_double_dot using transmission component
    theta_dot2 = Tdot2*x_dot2;
    theta_double_dot2 = Tdouble_dot2*x_dot2^2+Tdot2*xdouble_dot2; 
    
    % Link 1 position, velocity, acceleration at center of mass
    %Pm1 = [L_c1*cos(theta_1);...
    %       L_c1*sin(theta_1)];
    % Vm1 = [-L_c1*sin(theta_1)*theta_dot1;...
    %        L_c1*theta_dot1*cos(theta_1)];
    Am1 = [-L_c1*cos(theta_1)*theta_dot1^2-L_c1*sin(theta_1)*theta_doubledot1;
           -L_c1*(theta_dot1^2)*sin(theta_1)+L_c1*theta_doubledot1*cos(theta_1)];
    
    % Link 2 position, velocity, acceleration at center of mass

    % Pm2 = [L_1*cos(theta_1)+L_c2*cos(theta_1+theta_2_shift);
    %        L_1*sin(theta_1)+L_c2*sin(theta_1+theta_2_shift)];

    % Vm2 = [-L_1*sin(theta_1)*theta_dot1-L_c2*sin(theta_1+theta_2)*theta_dot1-L_c2*sin(theta_1+theta_2)*theta_dot2;
    %        L_1*cos(theta_1)*theta_dot1+L_c2*cos(theta_1+theta_2)*theta_dot1+L_c2*cos(theta_1+theta_2)*theta_dot2];

    Am2 = [-L_1*cos(theta_1)*theta_dot1^2-L_1*sin(theta_1)*theta_doubledot1-L_c2*cos(theta_1+theta_2_shift)*theta_dot1^2-L_c2*sin(theta_1+theta_2_shift)*theta_doubledot1-2*L_c2*cos(theta_1+theta_2_shift)*theta_dot1*theta_dot2-L_c2*cos(theta_1+theta_2_shift)*theta_dot2^2-L_c2*sin(theta_1+theta_2_shift)*theta_double_dot2;...
           -L_1*sin(theta_1)*theta_dot1^2+L_1*cos(theta_1)*theta_doubledot1-L_c2*sin(theta_1+theta_2_shift)*theta_dot1^2+L_c2*cos(theta_1+theta_2_shift)*theta_doubledot1-2*L_c2*sin(theta_1+theta_2_shift)*theta_dot1*theta_dot2-L_c2*sin(theta_1+theta_2_shift)*theta_dot2^2+L_c2*cos(theta_1+theta_2_shift)*theta_double_dot2];

end 
