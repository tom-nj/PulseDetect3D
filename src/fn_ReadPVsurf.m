% read pseudo-velocity spectrum surface of ground motion from file
%
% by Yuchuan Tang @ SEU, 6/11/2022
%--------------------------------------------------------------------

function PV_surf = fn_ReadPVsurf(filename,damping_desired)

       fid=fopen(filename, 'r');
       tline = fgetl(fid);
       temp = textscan(tline,'%f','Delimiter',' ');
       logTs = temp{1,1}';   %samping points of log10(Ts) for PV surface
       
       tline = fgetl(fid);
       temp = textscan(tline,'%f','Delimiter',' ');
       damping_read = temp{1,1}';   %sampling points of damping ratio for PV surface
       
    [~,index] = ismember(round(damping_desired,4),round(damping_read,4));    %round-off error affects ismember()
    if (min(index) < 1)
        fclose(fid);
        errordlg('The damping ratios in the PV surface data file do not align with those desired!');
        return;
    else
        temp = ftell(fid);
        fseek(fid, temp , 'bof');
        temp = fscanf(fid,'%e',[length(logTs),inf]);   %read log10(PV) data
        PV_surf = temp';
        PV_surf = PV_surf(index,:);     %corresponding to damping ratios in damping_desired       
        fclose(fid);            
     end
end

