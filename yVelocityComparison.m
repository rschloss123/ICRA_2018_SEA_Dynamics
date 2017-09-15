% Plot upward velocities for rigid and compliant cases
function yVelocityComparison
    N = 50;
    time_step = .0095;
    load('yVelocityCompliant')
    load('yVelocityRigid')
    figure
    plot((0:N-1)*time_step,y_velocity_compliant,'Linewidth',1.5)
    hold on 
    plot((0:N-1)*time_step,y_velocity_rigid','Linewidth',1.5)
    set(gca,'Fontsize',18)
    legend('Compliant', 'Rigid')
    xlabel('Time (s)')
    ylabel('y-Velocity (m/s)')
    title('Comparison of Upward Velocities')
    set(gca,'fontname','times')
    grid on

end 