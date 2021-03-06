% Create the discrete matrices A and B from the continuous A_1 and B_1

function[A,B] = get_discrete_matrices(A_1,B_1,actuator_dynamics_considered)
global time_step
% do approximate integration to get a discrete time model
num_iters = 100;
A = expm(A_1*time_step);

if actuator_dynamics_considered == 1 % compliant system
    B = [0,0,0,0;0,0,0,0;0,0,0,0;0,0,0,0;0,0,0,0;0,0,0,0;0,0,0,0;0,0,0,0];
else % rigid system
    B = [0,0;0,0;0,0;0,0];    
end
for i=1:num_iters
    B = B + time_step/num_iters * expm(A_1*(time_step-time_step/num_iters))*B_1;
end