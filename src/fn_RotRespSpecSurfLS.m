%This file contains 3 subroutines:
%  fn_RotRespSpecSurfLS()
%  fn_SincInterp()
%  fn_RotLSFixed()

%**************************************************************************
function [us_maxL,vs_maxL,as_maxL] = fn_RotRespSpecSurfLS(grd_acc1,grd_acc2,dt,flag_Sinc,Ts_array,damping_array,rot_angles)
% Response spectrum surfaces of a linear SDOF system 
% subjected to input motion rotated to various horizontal directions
%
%Input: 
%       grd_acc1, grd_acc2 -- acceleration time histories of the two orthogonal horizontal components of GRM (unit: cm/s2)
%       dt -- time step of grd_acc1 and grd_acc2
%       flag_Sinc -- control whether Sinc interpolation is applied to grd_acc1 and grd_acc2
%       Ts_array -- array of structural periods for generating resp. spectra
%       damping_array -- array of structural damping ratios for generating resp. spectra
%       rot_angles -- angles (unit: radian) rotating from the 1st original GRM component toward the 2nd component
%
%Output: 
%        us_maxL(i,j,k) -- spectral displacement (unit: cm) at damping(i) and Ts(j) corresp. to rot_angles(k)
%        vs_maxL(i,j,k) -- spectral velocity (unit: cm/s) at damping(i) and Ts(j) corresp. to rot_angles(k)
%        as_maxL(i,j,k) -- spectral acceleration (unit: cm/s^2) at damping(i) and Ts(j) corresp. to rot_angles(k)
%
% by Yuchuan Tang @ SEU, 1/12/2022
%--------------------------------------------------------------------------

        nstep=size(grd_acc1,1);
        if(nstep ~= size(grd_acc2,1))
            warndlg('The two orthogonal horizontal components have different lengths!');
            return;
        end
        
        num_Rot = length(rot_angles);

%%---- Sinc interpolation of the original acceleration record if needed
        if(flag_Sinc)
            [acc1_interp, dt_interp] = fn_SincInterp(grd_acc1, dt);
            [acc2_interp, ~] = fn_SincInterp(grd_acc2, dt);
        else
            acc1_interp = grd_acc1;
            acc2_interp = grd_acc2;
            dt_interp = dt;
        end

%%----- Response spectra calculation       
        num_Period = length(Ts_array);        %number of structure period points
        num_Damping = length(damping_array);  %number of structure damping points        
        us_maxL=zeros(num_Damping,num_Period,num_Rot);        % peak relative displacement array
        vs_maxL=zeros(num_Damping,num_Period,num_Rot);        % peak relative velocity array
        as_maxL=zeros(num_Damping,num_Period,num_Rot);        % peak total acceleration array
        
        for j=1:num_Period
          for i=1:num_Damping
              
            ws=2*pi/Ts_array(j);         % structure frequency
            damping = damping_array(i);     %structure damping ratio
            
            %--peak and ending responses of linear SDOF structure during the grd. accel. excitation in various rotated directions
            [u_peak1, v_peak1, a_peak1, u_end, v_end]=fn_RotLSFixed(ws,damping,dt_interp,acc1_interp,acc2_interp,rot_angles);         
                        
            wd=ws*sqrt(1-damping^2);     %damped frequency

            for k = 1 : num_Rot
                temp_after = atan(sqrt(1-damping^2)*v_end(k)/(ws*u_end(k)+damping*v_end(k)));
                if (temp_after < 0 )
                    temp_after = pi + temp_after;
                end
                t_0v_after = temp_after/wd;      %the time at which after-EQ velocity is zero
                u_peak2 = u_end(k)*cos(temp_after)+(v_end(k)+u_end(k)*damping*ws)/wd*sin(temp_after);
                u_peak2 = u_peak2*exp(-damping*ws*t_0v_after);
                u_peak2 = abs(u_peak2);       %the peak displacement after EQ             
            
                temp_after = atan((ws*u_end(k)+2*damping*v_end(k))*wd/(damping*ws*u_end(k)-(1-2*damping^2)*v_end(k))/ws);
                if (temp_after < 0 )
                    temp_after = pi + temp_after;
                end
                t_0a_after = temp_after/wd;     %the time at which after-EQ accel. is zero
                v_peak2 = v_end(k)*cos(temp_after)-(u_end(k)*ws+v_end(k)*damping)*ws/wd*sin(temp_after);
                v_peak2 = v_peak2*exp(-damping*ws*t_0a_after);
                v_peak2 = abs(v_peak2);        %the peak velocity after EQ
            
                temp_after = atan((2*damping*ws*u_end(k)-(1-4*damping^2)*v_end(k))*wd...
                                  /((2*damping^2-1)*ws*u_end(k)-(3-4*damping^2)*damping*v_end(k))/ws);
                if (temp_after < 0 )
                    temp_after = pi + temp_after;
                end
                t_0j_after = temp_after/wd;            %the time at which after-EQ jerk is zero
                a_peak2 = -cos(temp_after)*ws*(2*damping*v_end(k)+ws*u_end(k));
                a_peak2 = a_peak2+sin(temp_after)*ws^2*(damping*ws*u_end(k)-(1-2*damping^2)*v_end(k))/wd;
                a_peak2 = a_peak2*exp(-damping*ws*t_0j_after);
                a_peak2 = abs(a_peak2);        %the peak accel. after EQ 
           
                %--peak transient response is taken as the maximum during and after EQ
                us_maxL(i,j,k)=max(u_peak1(k), u_peak2);      %relative disp. (unit: cm)
                vs_maxL(i,j,k)=max(v_peak1(k), v_peak2);      %true relative vel. (unit: cm/sec)
                as_maxL(i,j,k)=max(a_peak1(k), a_peak2);      %true total accel. (unit: cm/sec^2)
            end
          end
        end
end
  
%**************************************************************************
function [data_interp, dt_interp]=fn_SincInterp(data, time_dt)
% Sinc interpolation of the original acceleration time series
% following RCTC by Wang&Stewart, et al. (2017)
%
%Input variables:
%       data: original acceleration time series (one-column vector)
%       time_dt: time step of the original data (scalar)
%
% by Yuchuan Tang @ SEU, May 16, 2020
%--------------------------------------------------------------------------

dt_desired = 0.005;         % desired time step for find peak response accurately

if (time_dt <= dt_desired)    % no need to interpolate
   data_interp = data;
   dt_interp = time_dt;

else
    
 % calculate interpolation factor (referenced to GM_RotD.R by Wang&Stewart, et al., 2017)
 interp_fac = power(2, min(ceil(log(time_dt/dt_desired)/log(2)), 3) );     % should be a power of 2

 npts = length(data);
 npw2_exp = floor(log(npts)/log(2));
 npw2_temp = power(2, npw2_exp);
  
 if (npw2_temp == npts)
    npw2 = npts;
 else
    npw2 = 2^(npw2_exp+1);
    data(npts+1 : npw2) = zeros(npw2-npts,1);   % ensure length(data) is a power of 2
 end
 
 nyqst = npw2/2 + 1;
 m = length(data);
 a = fft(data);
 b = [a(1:nyqst); zeros((m*interp_fac-m),1); a((nyqst+1):m)];
 b(nyqst) = b(nyqst)/2;
 b(nyqst+m*interp_fac-m) = b(nyqst);
 y = ifft(b);
 y = y*interp_fac;
 y = real(y);
 data_interp = y(1:(npts*interp_fac));      %one-column vector

 dt_interp = time_dt/interp_fac;        %scalar
 
end
end


%**************************************************************************
function [disp_max, vel_max, accel_max, u_end, v_end]=fn_RotLSFixed(wn,damping,dt,acc_input1,acc_input2,rot_angles)
% response of linear fixed-base SDOF structure
% to ground acceleration input rotated to various nonredundant directions
%
% solved with analytical expressions referenced to Appendix to A.K. Gupta's
% Book "Response Spectrum Method" (1990) and originally Nigam&Jennings (1969)
%
%%----input parameters
%wn: natural circular frequency of structure
%damping: damping ratio of structure
%dt: time step of the input acceleration time series
%acc_input1, acc_input2: acceleration time series of the two orthogonal horizontal components of input ground motion
%rot_angles: angles in radian to rotate the ground motion from the input acc_input1 toward acc_input2
%
% by Yuchuan Tang @ SEU, May 16, 2020
%--------------------------------------------------------------------------

if(length(acc_input1) ~= length(acc_input2))
     msgbox('The two input acceleration time series have different lengths!');
     return;
end

%%----calculate the individual responses to acc_input1 and acc_input2 respectively
temp = size(acc_input1);
disp1 = zeros(temp);      %relative structure displacement due to acc_input1
vel1 = zeros(temp);       %relative structure velocity due to acc_input1
accel1 = zeros(temp);     %absolute structure acceleration due to acc_input1
disp2 = zeros(temp);      %relative structure displacement due to acc_input2
vel2 = zeros(temp);       %relative structure velocity due to acc_input2
accel2 = zeros(temp);     %absolute structure acceleration due to acc_input2

wD = wn*sqrt(1-damping^2);
k = power(wn,2);
k1 = exp(-damping*wn*dt);

A11 = k1*(damping*wn/wD*sin(wD*dt) + cos(wD*dt));
A12 = k1*sin(wD*dt)/wD;
B12 = -1/k*(1-2*damping/(wn*dt)+k1*(((2*damping^2-1)/wD/dt)*sin(wD*dt)+2*damping/(wn*dt)*cos(wD*dt)));
B11 = -B12+(A11-1)/k;

A21 = -k*A12;
A22 = A11 - 2*damping*wn*A12;
B21 = -(A11-1)/(wn^2*dt) - A12;
B22 = -B21 - A12;

%--response at the first time step starting from quiet initial condition
disp1(1) = A11*0 + A12*0 + B11*0 + B12*acc_input1(1);  %relative displacement
vel1(1) = A21*0 + A22*0 + B21*0 + B22*acc_input1(1);   %relative velocity
accel1(1) = -2*damping*wn*vel1(1) - power(wn,2)*disp1(1);    %total acceleration

disp2(1) = A11*0 + A12*0 + B11*0 + B12*acc_input2(1);  %relative displacement
vel2(1) = A21*0 + A22*0 + B21*0 + B22*acc_input2(1);   %relative velocity
accel2(1) = -2*damping*wn*vel2(1) - power(wn,2)*disp2(1);    %total acceleration

%--response at the following time steps
for j = 2 : length(acc_input1)
   disp1(j) = A11*disp1(j-1) + A12*vel1(j-1) + B11*acc_input1(j-1) + B12*acc_input1(j);      %relative displacement
   vel1(j) = A21*disp1(j-1) + A22*vel1(j-1) + B21*acc_input1(j-1) + B22*acc_input1(j);   %relative velocity
   accel1(j) = -2*damping*wn*vel1(j) - power(wn,2)*disp1(j);             %absolute acceleration

   disp2(j) = A11*disp2(j-1) + A12*vel2(j-1) + B11*acc_input2(j-1) + B12*acc_input2(j);      %relative displacement
   vel2(j) = A21*disp2(j-1) + A22*vel2(j-1) + B21*acc_input2(j-1) + B22*acc_input2(j);   %relative velocity
   accel2(j) = -2*damping*wn*vel2(j) - power(wn,2)*disp2(j);             %absolute acceleration
end

%%----calculate the responses to ground acceleration rotated to various horizontal directions
num_Rot = length(rot_angles);
disp_max = zeros(num_Rot, 1);
vel_max = zeros(num_Rot, 1);
accel_max = zeros(num_Rot, 1);
u_end = zeros(num_Rot, 1);
v_end = zeros(num_Rot, 1);
a_end = zeros(num_Rot, 1);

time = 0 : dt : (length(acc_input1)-1)*dt;

for i = 1 : num_Rot
    n1 = cos(rot_angles(i));
    n2 = sin(rot_angles(i));
    disp = disp1.*n1 + disp2.*n2;
    vel = vel1.*n1 + vel2.*n2;
    accel = accel1.*n1 + accel2.*n2;
    
    disp_max(i) = max(abs(disp));   %peak relative displacement during GRM
    vel_max(i) = max(abs(vel));        %peak relative velocity during GRM
    accel_max(i) = max(abs(accel));    %peak total acceleration during GRM

    u_end(i) = disp(end);     %relative displacement at the end of GRM duration
    v_end(i) = vel(end);      %relative velocity at the end of GRM duration
    a_end(i) = accel(end);    %absolute acceleration at the end of GRM duration
    
end

end



