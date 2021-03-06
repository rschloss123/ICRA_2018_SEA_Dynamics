% Create the discrete linear time varying update equation for
% x_n+1 = A_lin,n*x_n+B_n*u_n+bias_n
function [A_lin,B_lin,bias_lin] = linearize_matrices(A,B_12,B_34,F_x_components,...
    F_u_components,F_b_components)

    A_lin = A+B_34*F_x_components;
    B_lin= B_12+B_34*F_u_components;
    bias_lin = B_34*F_b_components; 
        
end 