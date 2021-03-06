% We are considering a two link system which contains two actuators 
% states: [delta1; delta_dot1; y1; y_dot1; delta2; delta_dot2; y2; y_dot2]
% input: u = [motor1 current in amps; external force1; motor2 current in amps; external force2]
% where delta is spring deflection, and y is ball screw deflection
% total actuator deflection, x, satisfies x = y + delta for each individual
% actuator

% For actuator dynamics
% we put the results above into a matrix form:
% E_0 chi_dot = A_0 chi + B_0 u
function [E_0,A_0,B_0,A_1,B_1]= get_continuous_matrices(m_L1,m_L2)

global m_M b_M k_M Km_M k_S m_S b_S k_L b_L Km_L

E_0 = [1,       0,       0,       0, 0,       0,       0,       0;
       0, m_S+m_L1,       0,     m_L1, 0,       0,       0,       0;
       0,       0,       1,       0, 0,       0,       0,       0;
       0,     m_L1,       0, m_M+m_L1, 0,       0,       0,       0;
       0,       0,       0,       0, 1,       0,       0,       0;
       0,       0,       0,       0, 0, m_S+m_L2,       0,     m_L2;
       0,       0,       0,       0, 0,       0,       1,       0;
       0,       0,       0,       0, 0,     m_L2,       0, m_M+m_L2];       
   
A_0 = [         0,          1,          0,          0,         0,          0,          0,          0;
       -k_S - k_L, -b_S - b_L,       -k_L,       -b_L,         0,          0,          0,          0; 
                0,          0,          0,          1,         0,          0,          0,          0;
             -k_L,       -b_L, -k_M - k_L, -b_M - b_L,         0,          0,          0,          0;
                0,          0,          0,          0,         0,          1,          0,          0;
                0,          0,          0,          0,-k_S - k_L, -b_S - b_L,       -k_L,       -b_L; 
                0,          0,          0,          0,         0,          0,          0,          1;
                0,          0,          0,          0,      -k_L,       -b_L, -k_M - k_L, -b_M - b_L];

% Below, columns 1 and 3 refer to the inputs to actuator 1
% and columns 2 and 4 refer to the inputs to actuator 2
B_0 = [    0,     0,     0,     0;
           0,     0, -Km_L,     0;
           0,     0,     0,     0;
       -Km_M,     0, -Km_L,     0;
           0,     0,     0,     0; 
           0,     0,     0, -Km_L; 
           0,     0,     0,     0; 
           0, -Km_M,     0, -Km_L];

% do inverse to get standard model 
A_1 = E_0\A_0;
B_1 = E_0\B_0;

end 