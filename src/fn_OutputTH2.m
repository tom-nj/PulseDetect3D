% Write the time history data of the detected pulse and the corresponding 
% ground motion into *.AT2 and *.VT2 files
%
% by Yuchuan Tang @ SEU, 7/14/2020
%--------------------------------------------------------------------------

function fn_OutputTH2(path_output,RecdNameStr,subfolder,Rot_string,dt,grd_acc,...
                      grd_vel,pulse_type,TH_istart,TH_istop,a_pulse,v_pulse,flag_outputGRM)

        gravity = 981;      %gravity acceleration in cm/s2
        
        temp = strsplit(RecdNameStr, {'\', '/'});
        recd_name = temp{end};
        np = length(grd_vel);
        grd_acc = grd_acc./gravity;         %convert unit to g

        pulse_vth = zeros(size(grd_vel));     
        pulse_ath = zeros(size(grd_acc));
        L_detd = TH_istop - TH_istart;      %length of the detected pulse, which may be truncated, in the original EQ record 
                                            %(the original EQ record may miss the early part of the full pulse waveform)
        pulse_ath(TH_istart:TH_istop) = a_pulse(end-L_detd:end)./gravity;   %pulse (possibly truncated) matching grd_acc time history
        pulse_vth(TH_istart:TH_istop) = v_pulse(end-L_detd:end);            %pulse (possibly truncated) matching grd_vel time history
                
        %--output time histories of the rotated ground motion
        if (flag_outputGRM < 1)
        folder = strcat(path_output, '\', RecdNameStr, '\', Rot_string);
        
        filename1 = strcat(recd_name,'_',Rot_string);
        
        fid = fopen(strcat(folder, '\', filename1, '.VT2'), 'w');
        fprintf(fid, 'PEER Format PulseDetect3D_2022 GRM\n');        
        fprintf(fid, '%s\n', RecdNameStr);
        fprintf(fid, 'VELOCITY TIME SERIES IN UNITS OF CM/S\n');
        fprintf(fid, 'NPTS=  %d, DT=   %.3f SEC\n', np, dt);
        fprintf(fid,'%10.5E %10.5E %10.5E %10.5E %10.5E\n',grd_vel);    
        fclose(fid);
        
        fid = fopen(strcat(folder, '\', filename1, '.AT2'), 'w');
        fprintf(fid, 'PEER Format PulseDetect3D_2022 GRM\n');        
        fprintf(fid, '%s\n', RecdNameStr);
        fprintf(fid, 'ACCELERATION TIME SERIES IN UNITS OF G\n');
        fprintf(fid, 'NPTS=  %d, DT=   %.3f SEC\n', np, dt);
        fprintf(fid,'%10.5E %10.5E %10.5E %10.5E %10.5E\n',grd_acc);    
        fclose(fid);
        
        %--record the names of the generated time history files to a text list       
        fid = fopen(strcat(path_output, '\List_GeneratedTHfiles.txt'), 'a');
        temp = char(subfolder);
        fprintf(fid, '%s\n', strcat(temp(length(path_output)+2:end), '\', filename1));
        fclose(fid);        
        end
        
        %--output time histories of the detected pulse        
        filename2 = strcat(recd_name,'_',Rot_string, '_', pulse_type);
        
        fid = fopen(strcat(subfolder, '\', filename2, '.VT2'), 'w');       
        fprintf(fid, 'PEER Format PulseDetect3D_2022 Output Pulse\n');
        fprintf(fid, '%s\n', RecdNameStr);
        fprintf(fid, 'VELOCITY TIME SERIES IN UNITS OF CM/S\n');
        fprintf(fid, 'NPTS=  %d, DT=   %.3f SEC\n', np, dt);
        fprintf(fid,'%10.5E %10.5E %10.5E %10.5E %10.5E\n',pulse_vth);    
        fclose(fid);

        fid = fopen(strcat(subfolder, '\', filename2, '.AT2'), 'w');
        fprintf(fid, 'PEER Format PulseDetect3D_2022 Output Pulse\n');
        fprintf(fid, '%s\n', RecdNameStr);
        fprintf(fid, 'ACCELERATION TIME SERIES IN UNITS OF G\n');
        fprintf(fid, 'NPTS=  %d, DT=   %.3f SEC\n', np, dt);
        fprintf(fid,'%10.5E %10.5E %10.5E %10.5E %10.5E\n',pulse_ath);    
        fclose(fid);
        
        %--record the names of the generated time history files to a text list       
        fid = fopen(strcat(path_output, '\List_GeneratedTHfiles.txt'), 'a');
        temp = char(subfolder);
        fprintf(fid, '%s\n', strcat(temp(length(path_output)+2:end), '\', filename2));
        fclose(fid);
        

%{
        %--output time histories of the residual motion        
        resid_ath = grd_acc - pulse_ath;    %residual acceleration time history in g
        resid_vth = grd_vel - pulse_vth;    %residual velocity time history in cm/s
        
        filename3 = strcat(filename2, '_resid');
        
        fid = fopen(strcat(subfolder, '\', filename3, '.VT2'), 'w');       
        fprintf(fid, 'PEER Format PulseDetect3D_2022 Output Residual Motion\n');
        fprintf(fid, '%s\n', RecdNameStr);
        fprintf(fid, 'VELOCITY TIME SERIES IN UNITS OF CM/S\n');
        fprintf(fid, 'NPTS=  %d, DT=   %.3f SEC\n', np, dt);
        fprintf(fid,'%10.5E %10.5E %10.5E %10.5E %10.5E\n',resid_vth);    
        fclose(fid);

        fid = fopen(strcat(subfolder, '\', filename3, '.AT2'), 'w');
        fprintf(fid, 'PEER Format PulseDetect3D_2022 Output Residual Motion\n');
        fprintf(fid, '%s\n', RecdNameStr);
        fprintf(fid, 'ACCELERATION TIME SERIES IN UNITS OF G\n');
        fprintf(fid, 'NPTS=  %d, DT=   %.3f SEC\n', np, dt);
        fprintf(fid,'%10.5E %10.5E %10.5E %10.5E %10.5E\n',resid_ath);    
        fclose(fid);
        
        %--record the names of the generated time history files to a text list       
        fid = fopen(strcat(path_output, '\List_GeneratedTHfiles.txt'), 'a');
        temp = char(subfolder);
        fprintf(fid, '%s\n', strcat(temp(length(path_output)+2:end), '\', filename3));
        fclose(fid);
%}
        
end


