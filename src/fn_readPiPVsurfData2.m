% read from file the Pi_PV surface data (the velocity-amplified region only) 
% corresponding to a pulse waveform
%
% by Yuchuan Tang @ SEU, 6/11/2022
%--------------------------------------------------------------------

function [logPi_T,Pi_xi,Pipv_VAR] = fn_readPiPVsurfData2(path_inputPiSpec,flg1,para1,para2,damping_array)

FileLocation=strcat(path_inputPiSpec,'PiPVsurf_',flg1 ,'pulse_',para1,'_',para2,'.dat');
fid=fopen(FileLocation,'r');
tline = fgetl(fid);
temp = textscan(tline,'%f','Delimiter',' ');
logPi_T = temp{1,1}';   %samping points of log10(Pi_T) for Pi_PV surface

tline = fgetl(fid);
temp = textscan(tline,'%f','Delimiter',' ');
Pi_xi = temp{1,1}';     %sampling points of damping ratio for Pi_PV surface

[~,index] = ismember(round(damping_array,4),round(Pi_xi,4));    %round-off error affects ismember()
if (min(index) < 1)
    fclose(fid);
    errordlg('The damping ratios for PV surface do not align with those for Pi_PV surface!');
    return;
else
    Pi_xi = Pi_xi(index);

    temp = ftell(fid);
    fseek(fid, temp , 'bof');
    temp=fscanf(fid,'%e',[length(logPi_T),inf]);  %log10(Pi_PV) data at the sampling grid
    Pipv_VAR=temp';
    Pipv_VAR=Pipv_VAR(index,:);     %corresponding to damping ratios in Pi_xi
    fclose(fid);
end

end

