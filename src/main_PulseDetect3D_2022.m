%**************************************************************************
% Velocity pulse detection and characterization in the ground motion 
% components rotated to different horizontal orientations
% according to the method presented in the paper,
% Y. Tang, C. Wu, J. Wang. A hybrid characterization framework
% for structure-significant pulse-like features in ground motions,
% Soil Dynamics and Earthquake Engineering, 160 (2022) 107325.
% ( https://doi.org/10.1016/j.soildyn.2022.107325 )
%
% developed by Yuchuan Tang at Southeast University, Nanjing, China
% released on June 12, 2022 with copyright reserved.
%**************************************************************************

clear
close all

%% Basic parameters 
%--the list of earthquake records to analyze
EQrecdList_input = '..\Example\Input_GRM\List_GRM.txt';

%--path to locate the AT2 files of the EQ records
path_inputAT2 = '..\Example\Input_GRM';

%--path to locate the standard Pi spectra of the candidate pulse types
path_inputPiSpec = '..\PiPVsurf_database\';

%--path to output the results (figures and list of detected pulses)
path_output = '..\Output';
mkdir(path_output);

%--output file containing a list of detected pulse parameters
pulselist_output = strcat(path_output,'\Pulse_Parameters.txt');

%--rotation angle interval for rotating the ground motion from the 1st input
%--component toward the 2nd input component, which are orthogonal horizontal components 
RotIntv = 45;       %unit: deg; must be a factor of 90 and 0<RotIntv<=180; assign 180 to investigate the 1st component ONLY
num_Rot = floor(180/RotIntv);   %number of horizontal orientations to rotate to
rot_array = (RotIntv/180*pi).*[0:num_Rot-1];    %angles in radian to rotate from the 1st input component toward the 2nd component

%--PGV threshold to skip ground motions with low PGV
PGV_threshold = 20;     %unit: cm/s

%--default parameters for generating response spectrum surfaces
damping = 0.05 : 0.01 : 0.2;        %structure damping ratio (the grid must be aligned with that for Pi_PV surface)
logTs = -1.0 : 0.01 : 2.0;          %must have identical interval as logPi_T of the Pi_PV surface
Ts = 10.^logTs;                     %structure period (unit: sec)
ws = 2*pi./Ts;                      %circular frequency (unit: rad/s), must be a one-row vector
flag_CorR = 1;          %1--to compute PV surface; 0--to read PV surface from a previously generated file
flag_savePV = 1;        %1--to save the computed PV surface to file; 0--not to save to file

%--velocity-amplified region threshold (ratio to PGV) considered 
VelAmp_tol = 0.6;

%--threshold of corr2 that indicates subregion-to-subregion similarity between PV and Pi_PV surfaces
R2_threshold = 0.9;

%--pulse type candidates (Note: the corresponding admissible wavelets must be incorporated into MATLAB beforehand)
PsType_Array = [...
                "MP_1_90";...  %M&P pulse with gamma=1, nu=90
                "MP_1.5_0";... %M&P pulse with gamma=1.5, nu=0 (not admissible wavelet)
                "MP_2_0";...   %M&P pulse with gamma=2, nu=0
                "MP_2_90";...  %M&P pulse with gamma=2, nu=90
                "MP_3_0";...   %M&P pulse with gamma=3, nu=0
                "MP_3_90";...  %M&P pulse with gamma=3, nu=90
                "MP_4_0";...   %M&P pulse with gamma=4, nu=0
                "MP_4_90"...   %M&P pulse with gamma=4, nu=90
                ];

numPsType = length(PsType_Array);       %total number of pulse type candidates

%% Read Pi_PV spectrum surface data associated with candidate pulse types
logPi_T_cell = cell(numPsType, 1);      %dimensionless structural period Pi_T
Pi_xi_cell=cell(numPsType, 1);          %structure damping ratio Pi_xi
Pipv_VAR_cell = cell(numPsType, 1);     %dimensionless pseudo spectral velocity Pi_PV

for i_type = 1 : numPsType
    %---input the standard Pi_PV surface of current pulse type
    temp = char(PsType_Array(i_type));      %convert string to character array
    temp = split(temp, '_');
    
    %---specify the input Pi_PV-surface file name CAREFULLY!!
    [logPi_T_cell{i_type},Pi_xi_cell{i_type},Pipv_VAR_cell{i_type}] = fn_readPiPVsurfData2(path_inputPiSpec,temp{1},temp{2},temp{3},damping);
end

%% Read the list of earthquake records to analyze
fid1=fopen(EQrecdList_input, 'r');
namesRec=textscan(fid1, '%s');      %returns a cell array of strings
fclose(fid1);

%% Loop through three-component EQ record sets in the list
% for indexRec = 3 : 3 : size(namesRec{1},1)
for indexRec = 3 : 3 : 6
    disp('***************************'); 
    indexRec

    RecdNameStr1 = namesRec{1}{indexRec-2};
    if(strcmp(RecdNameStr1(end-3:end),'.AT2'))
        RecdNameStr1 = RecdNameStr1(1:end-4);
    end
    disp(RecdNameStr1);
    RecdNameStr2 = namesRec{1}{indexRec-1};
    if(strcmp(RecdNameStr2(end-3:end),'.AT2'))
        RecdNameStr2 = RecdNameStr2(1:end-4);
    end
    disp(RecdNameStr2);

    %--read accel. data of two orthogonal horizontal GRM components from two AT2 files
    [RecdNameStr,grd_acc1,azim1,grd_acc2,azim2,dt,flag_RotDir]=fn_ReadOrthGRM2(path_inputAT2,RecdNameStr1,RecdNameStr2);
                
    folder = strcat(path_output, '\', RecdNameStr);      %folder to save output files for current EQ record
    mkdir(folder);
    
    %--calculate response spectrum surfaces of the GRM at various rotation angles
    if (flag_CorR)  
        [SD_matrix, ~, ~] = fn_RotRespSpecSurfLS(grd_acc1,grd_acc2,dt,1,Ts,damping,rot_array);
    end
    
  %--loop through the rotation angles to detect pulse in ground velocity in each rotated orientation 
  for i_Rot = 1 : length(rot_array)
    disp('======================');
    Rot_string = strcat('Rot',num2str(RotIntv*(i_Rot-1)))    %current rotation angle (unit: degree) as a string
    
    %--azimuth of current direction rotated from the 1st input GRM component
    if (flag_RotDir > 0)    %rotating from azim1 toward azim2 is clockwise
        azim = azim1 + RotIntv*(i_Rot-1);
        if(azim >= 360)
            azim = azim - 360;
        end
    else     %rotating from azim1 toward azim2 is anti-clockwise
        azim = azim1 - RotIntv*(i_Rot-1);
        if(azim < 0)
            azim = azim + 360;
        end
    end
    
    %--rotated ground motion in current orientation
    grd_acc_rot=grd_acc1.*cos(rot_array(i_Rot))+grd_acc2.*sin(rot_array(i_Rot));      %ground accel. (unit: cm/s2)
    grd_vel_rot=cumtrapz(dt, grd_acc_rot);      %ground velocity (unit: cm/s)
    PGV = max(abs(grd_vel_rot));                %peak ground velocity (unit: cm/s)     
    if(PGV < PGV_threshold) 
        fid2 = fopen(pulselist_output, 'a');
        fprintf(fid2, '%s %s %.1f LowPGV\n', RecdNameStr, Rot_string, PGV);
        fclose(fid2);
        
        disp("Skip ground motion for low PGV");        
        continue;       %skip pulse detection for ground motion with low PGV
    end
    grd_disp_rot=cumtrapz(dt, grd_vel_rot);      %ground displacement (unit: cm)   
    
    if (flag_CorR)
        %--compute PV spectrum surface of ground motion on the fly
        PV_surf = ws.*SD_matrix(:,:,i_Rot);     %PV spectrum surface (unit: cm/s)
        PV_surf = log10(PV_surf);
        
        %--save PV spectrum surface data to file
        if (flag_savePV)
            temp=strfind(RecdNameStr,'\');
            filename=strcat(folder,'\',RecdNameStr(temp+1:end),'_',Rot_string,'_PVsurf','.dat');   
            fn_SavePVsurf(filename,logTs,damping,PV_surf);
        end
    else    
       %--read PV spectrum surface of ground motion from file
       temp = strfind(RecdNameStr,'\');
       filename=strcat(folder,'\',RecdNameStr(temp+1:end),'_',Rot_string,'_PVsurf','.dat');   
       PV_surf = fn_ReadPVsurf(filename,damping);   %log10(PV) data
    end       

      fn_FixOrientation_PulseDetect(path_output,pulselist_output,folder,...
                                    VelAmp_tol,R2_threshold,PsType_Array,...
                                    logPi_T_cell,Pipv_VAR_cell,...
                                    RecdNameStr,Rot_string,...
                                    grd_acc_rot,grd_vel_rot,grd_disp_rot,PGV,dt,...
                                    logTs,PV_surf);

  end
  
end
          