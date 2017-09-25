function [l1, l2, m1,m2] = link_animation_components(theta1, theta2_shift)
    global L_1 L_2 L_c1 L_c2 time_step 
    
    % location of end of link 1
    l1 = [cos(theta1)*L_1; sin(theta1)*L_1];
    % location of end of link 2
    l2 = [cos(theta1)*L_1+cos(theta1+theta2_shift)*L_2; sin(theta1)*L_1+sin(theta1+theta2_shift)*L_2];
 

   % COM locations for animation
    m1 = [cos(theta1)*L_c1; sin(theta1)*L_c1];
    m2 = [cos(theta1)*L_1+cos(theta1+theta2_shift)*L_c2; sin(theta1)*L_1+sin(theta1+theta2_shift)*L_c2];
    
end 