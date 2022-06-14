%plot and save the time histories of the EQ record vs. the detected pulse
%
%by Yuchuan Tang @SEU, 1/13/2022
%--------------------------------------------------------------------------

function fn_PlotDetdPulse2(folder,RecdNameStr,grd_acc,grd_vel,grd_disp,dt,...
                          pulse_type,TH_istart,TH_istop,a_pulse,v_pulse,d_pulse)

         gravity = 981;     %gravity acceleration (unit: cm/s^2)
             
         recd_time = 0 : dt : (length(grd_vel)-1)*dt;   %time points of the ground motion time series

         L_detd = TH_istop - TH_istart;     %length of the detected pulse, which may be truncated, in the original EQ record 
                                            %(the original EQ record may miss the early part of the full pulse)

         figure('Name', RecdNameStr,'Units','normalized','OuterPosition',[0 0 1 1]);

         subplot(3,1,1);
         plot(recd_time, grd_acc./gravity, 'b-', 'LineWidth', 0.5);
         hold on;
         plot(recd_time(TH_istart:TH_istop), a_pulse(end-L_detd:end)./gravity, 'r--', 'LineWidth', 1);         
         legend('EQ record', 'Pulse');
         ylabel('Accel. (g)');

         subplot(3,1,2)
         plot(recd_time, grd_vel, 'b-', 'LineWidth', 0.5);
         hold on;
         plot(recd_time(TH_istart:TH_istop), v_pulse(end-L_detd:end), 'r--', 'LineWidth', 1);
         title(strcat("Pulse ", pulse_type), 'Interpreter','none');
         ylabel('Vel. (cm/s)');

         subplot(3,1,3)
         plot(recd_time, grd_disp, 'b-', 'LineWidth', 0.5);
         hold on;
         plot(recd_time(TH_istart:TH_istop), d_pulse(end-L_detd:end), 'r--', 'LineWidth', 1);
         xlabel('Time (s)');
         ylabel('Disp. (cm)');

         temp = strsplit(folder, {'/','\'});
         temp = strcat(temp{end-2},'\',temp{end-1},'\',temp{end});
         sgtitle(temp,'FontSize',10,'Interpreter','none');

%          pause
         
         saveas(gcf, strcat(folder,'\DetdPulse_',pulse_type,'.JPEG'), 'jpeg');

         close(gcf);                   

end

