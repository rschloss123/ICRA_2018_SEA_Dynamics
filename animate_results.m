% Animate the motion of the 2-link Draco leg (bolted to the ground) 
function animate_results(utilize_compliance, trajectory, N, foot_length, input)
    global L_1 L_2 L_c1 L_c2 time_step 
    
    % Animation
    figure()
    clf()
   
    v = VideoWriter('new1');
    % v.framerate=60;
    v.FrameRate = 1/time_step; 
    v.Quality=98;
    disp(v.VideoCompressionMethod)
    open(v)
    for i=1:N
        % visualization in a plot
        clf()
        axs=plot([0,0],[0,0]);
        if utilize_compliance == 1 
            % link 1
            x1= [1,0,1,0,0,0,0,0]*trajectory(i,:)';
            % link 2
            x2= [0,0,0,0,1,0,1,0]*trajectory(i,:)';    
        elseif utilize_compliance == 0
            % link 1
            x1= [1,0,0,0]*trajectory(i,:)';
            % link 2
            x2= [0,0,1,0]*trajectory(i,:)'; 
        end
        theta1=get_theta1(x1); 
        theta2_unshift= get_theta_2(x2); 
        theta2_shift = theta2_unshift+pi; 

        [l1, l2, m1, m2] = link_animation_components(theta1, theta2_shift);

        plot([0, l1(1)],[0,l1(2)],'o-','Linewidth',1.3, 'Color', 'k');
        hold on
        plot([l1(1), l2(1)],[l1(2),l2(2)],'Linewidth',1.3,'Color', 'k')


     
        plot([-foot_length,0],[0,0],'Linewidth',1.3,'Color', 'k')
        hold on
        plot([m1(1)],[m1(2)],'o','MarkerSize',20, 'Color','[.3 .2 .1]')
        plot([m2(1)],[m2(2)],'o','MarkerSize',20, 'Color','[.3 .2 .1]')

        if utilize_compliance == 1 && input == 1
        % plot reaction forces 
    %         if i < N
    %             reaction = c_constraint(i,:);
    % 
    %             point = [-foot_length;0];
    %             R1 = .005*(reaction(1,1)*[mu;1]+reaction(1,2)*[-mu;1]);
    %             plot([point(1),point(1)+R1(1)],[point(2),point(2)+R1(2)])
    %             R2 = .005*(reaction(1,3)*[mu;1]+reaction(1,4)*[-mu;1]);
    %             plotv(1*R2,'-')
    % 
    %         end     
        end

        if utilize_compliance == 1
            delta1 = [1,0,0,0,0,0,0,0]*trajectory(i,:)';
            delta2 = [0,0,0,0,1,0,0,0]*trajectory(i,:)';
            springscale = 25; 
            [points1, points2] = get_spring_locations(theta1, theta2_shift, delta1*springscale, delta2*springscale);
            plot(points1(1,:), points1(2,:),'o','MarkerSize',10, 'Color', 'm');
            plot(points2(1,:), points2(2,:),'o','MarkerSize',10,'Color', 'm');
        end 
        plot([-.5 .4],[-.005 -.005],'-','Color',[.5 .5 .5])
        axis('equal')

        axis([-.5,.4,-.1,1.3])
        set(gca, 'YTickLabel',[]);
        set(gca, 'YTick',[]);
        set(gca, 'XTickLabel',[]);
        set(gca, 'XTick',[]);  
        drawnow
    %     F(i)=getframe;
        writeVideo(v,getframe)
        %pause(0.001)
    end
    close(v)
    
end 
