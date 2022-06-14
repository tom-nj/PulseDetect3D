%this file contains the following subroutines:
%   fn_ReadOrthGRM2()
%   fn_ReadAT2File2A()
%   fn_ConvAzimuth()
%
% by Yuchuan Tang@SEU, 6/27/2020 
%--------------------------------------------------------------------------

function [RecdNameStr,grd_acc1,azim1,grd_acc2,azim2,dt,flag_RotDir]= ...
                  fn_ReadOrthGRM2(path_inputAT2,RecdNameStr1,RecdNameStr2)
%
%read the acceleration data of two orthogonal horizontal components
%of the ground motion as recorded from *.AT2 files in the PEER format

    %--check whether the two input files are a pair of orthogonal components   
    RecdNameL = length(RecdNameStr1);
    if( RecdNameL ~=  length(RecdNameStr2) )
        msgbox(strcat('Different length between ', RecdNameStr1, RecdNameStr2));
        return;
    else
        for n = RecdNameL-1 : -1 : RecdNameL-4
            flag_common = strncmp(RecdNameStr1, RecdNameStr2, n);
            if(flag_common)
                RecdNameStr = RecdNameStr1(1:n);    %the common part of the two AT2 file names
                break;
            end
        end
        if(~flag_common)
            msgbox(strcat('NOT orthogonal components: ', RecdNameStr1, RecdNameStr2));
            return;
        end
    end
    
    %--read the accel. data of the EQ record from an AT2 file at FileLocation
    FileLocation = strcat(path_inputAT2, '\', RecdNameStr1, '.AT2');
    [grd_acc1, dt1, NPTS1, azim1] = fn_ReadAT2File2A(FileLocation);

    FileLocation = strcat(path_inputAT2, '\', RecdNameStr2, '.AT2');
    [grd_acc2, dt2, NPTS2, azim2] = fn_ReadAT2File2A(FileLocation);

    %--skip the EQ record when the time steps of the two components differ
    if(dt1 ~= dt2)
        msgbox('The input GRM components have different time steps!');
        return;
    else
        dt = dt1;
    end
    
    %--skip the EQ record when the difference between number of data points 
    %--recorded in the two horizontal GRM components is greater than 20
    if(abs(NPTS1-NPTS2) > 20)
        msgbox('The input GRM components have quite different lengths!');
        return;
    elseif(NPTS1 < NPTS2)
        grd_acc2 = grd_acc2(1:NPTS1);
    elseif(NPTS1 > NPTS2)
        grd_acc1 = grd_acc1(1:NPTS2);
    end
            
    %--check the directional relativity of the two input GRM components
    azim_diff = azim2 - azim1;
    switch azim_diff
        case {90, -270}
            flag_RotDir = 1;    %rotating from azim1 toward azim2 is clockwise
        case {-90, 270}
            flag_RotDir = -1;    %rotating from azim1 toward azim2 is anti-clockwise
        otherwise
            msgbox('The input GRM components are NOT orthogonal!');
            return;
    end
    
end


%**************************************************************************
function [grd_acc, dt, NPTS, azim] = fn_ReadAT2File2A(FileLocation)

  filename=cellstr(FileLocation);         
  fid=fopen(filename{1}, 'r');
        
  header1 = fgetl(fid);
  if( ~strcmp(header1(1:4), 'PEER') && ~strcmp(header1(1:5), 'Baker') )
      msgbox('The earthquake record is not in PEER format!');
      return;
  end
        
  header2 = fgetl(fid);
  temp = strsplit(header2, ',');
  temp = strtrim( temp{end} );
  azim = fn_ConvAzimuth(temp);       %azimuth angle
  
  fgetl(fid);
  temp1 = ftell(fid);
  header4 = fgetl(fid);
  fseek(fid, temp1, 'bof');

  if( strcmp(header4(end-1:end), 'DT') )
     %--PEER format "4111  0.0050  NPTS, DT"            
     %--get number of points and time step in the *.AT2 file        
        temp2=textscan(fid, '%f %f %*s %*s',1);
  elseif( strcmp(header4(1:4), 'NPTS') )
     %--PEER format "NPTS= 4111 , DT= .00500 SEC,"
     %--get number of points and time step in the *.AT2 file               
        temp2=textscan(fid, '%*s %f %*s %*s %f %*s', 1);   
  end       
        
  dt=temp2{1, 2};          % time step, sec

  %--add the next two commands to deal with the 4th row in AT2 file that contains extra tailing characters
  fseek(fid, temp1, 'bof'); 
  fgetl(fid);
        
  %--read the acceleration data
  data = fscanf(fid, '%f');
  fclose(fid);    
    
  gravity=981;                   %gravity acceleration, cm/sec^2       
  grd_acc=gravity.*data;         %unit: cm/sec^2

  NPTS=size(grd_acc,1);          %number of data points
                
end

%**************************************************************************
function [azim] = fn_ConvAzimuth(strtemp)

%convert the input string that contains azimuth information 
%to numerical value of the azimuth within [0, 360) deg.
%------------------------------------------------------------------------

  switch strtemp               
      case {'NS','N','L'}
          azim=0;
      case {'EW','E','T'}
          azim=90;
      case {'SN','S'}
          azim=180;
      case {'WE','W'}
          azim=270;
      case 'UP'
          msgbox('Vertical GRM component is input by mistake!');
          return;          
      otherwise
          numtemp = str2double(strtemp);
          if(~isnan(numtemp))
             if(numtemp>=0 && numtemp<=360)
                 azim = numtemp;
             else
                 msgbox('The azimuth info in *.AT2 is unexpected!');
                 return;
             end
          elseif(length(strtemp)<3)
                 msgbox('The azimuth info in *.AT2 is unexpected!');
                 return;      
          else
              temp = str2double( strtemp(2:end-1) );
              if(isnan(temp))
                 msgbox('The azimuth info in *.AT2 is unexpected!');
                 return;                  
              else
                  switch strcat(strtemp(1),strtemp(end))
                      case 'NE'
                          azim = temp;
                      case 'NW'
                          azim = 360 - temp;
                      case 'SE'
                          azim = 180 - temp;
                      case 'SW'
                          azim = 180 + temp;
                      otherwise
                        msgbox('The azimuth info in *.AT2 is unexpected!');
                        return;                          
                  end
              end
          end
   end

end


