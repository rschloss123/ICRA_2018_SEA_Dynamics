% Parameters of the actuator and robot systems
function get_parameters(scale)
global m_M b_M k_M Km_M m_S b_S work k_L b_L Km_L...
    J1 J2 m_1 m_2 L_1 L_2 L_c1 L_c2 
% Motor
m_M = 293; % Kg
b_M = 1680; % N s/m
k_M = 0; % motor subsystem has no spring effect
Km_M = -157; % N/amp, negative due to the way the motor is connected 

% Spring
m_S = 1.7; % Kg
b_S = 0; % N s/m

% Load
b_L = 0; % N s/m
k_L = 0; % N/m
Km_L = 1; % N/m (relation to setpoint)
work = 0; % used in ConstraintComponents.m to calculated accumulation of energy dissipated by the dampers 

% link masses (kg)
m_1 = 3.77;
m_2 = 30.608;
% link lengths (m)
L_1 = .5; 
L_2 = .5; 
% center of mass position in x direction, m 
L_c1 = .29557; 
L_c2 = .4361;  
% inertia of links 1 and 2
% Inertia is about the center of mass of the arm, along axis of rotation. 
% kg-m^2
J1 = .077;
J2 = .1502;
end 