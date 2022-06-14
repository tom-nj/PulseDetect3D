% this files contains the following subroutines:
% fn_PulseIndicator
% fn_MilestoneTime
%
% by Yuchuan Tang @ SEU, 6/12/2022

%*************************************************************************
function [PI_Shahi,Arr_Shahi,Ep_Chang,Arr_Chang,PI_Kard] = ...
               fn_PulseIndicator(grd_vel,dt,PGV,TH_istart,TH_istop,v_pulse)
               
%---pulse indicator and late arriving check following Shahi&Baker (2014)
        velTH_pulse = zeros(size(grd_vel));
        L_detd = TH_istop - TH_istart;       
        velTH_pulse(TH_istart:TH_istop) = v_pulse(end-L_detd:end);
        velTH_resid = grd_vel - velTH_pulse;
        PGV_ratio = max(abs(velTH_resid))/PGV;
        CSV_orig = trapz(dt, power(grd_vel, 2));
        CSV_resid = trapz(dt, power(velTH_resid, 2));
        Energy_ratio = CSV_resid/CSV_orig;
        PC_Shahi = 0.63*PGV_ratio + 0.777*Energy_ratio;
        PI_Shahi = 9.384*(0.76-PC_Shahi-0.0616*PGV)*(PC_Shahi+6.914E-4*PGV-1.072)-6.179;    %PGV in cm/s
        
        [t17_orig] = fn_MilestoneTime(grd_vel, dt, 0.17*CSV_orig);
        csv_pulse_t17orig = trapz(dt, power(velTH_pulse(1:round(t17_orig/dt+1)), 2));
        CSV_pulse = trapz(dt, power(velTH_pulse, 2));
        if(csv_pulse_t17orig > 0.05*CSV_pulse)
            Arr_Shahi = "Early";
        else
            Arr_Shahi = "Late";
        end
       
%---relative energy contained within the entire velocity pulse following Chang et al.(2016) & Zhai et al. (2013)
        CSV_orig_PulseWindow = trapz(dt, power(grd_vel(TH_istart:TH_istop), 2));
        Ep_Chang = CSV_orig_PulseWindow / CSV_orig;        
        ts_Chang = dt*(TH_istart-1);
        [t5_orig] = fn_MilestoneTime(grd_vel, dt, 0.05*CSV_orig);
        [t95_orig] = fn_MilestoneTime(grd_vel, dt, 0.95*CSV_orig);
        Omega_mid_Zhai = (t5_orig + t95_orig)/2;
        
        if(ts_Chang < Omega_mid_Zhai)
            Arr_Chang = "Early";
        else
            Arr_Chang = "Late";
        end
        
%---pulse indicator following Kardoutsou et al. (2017)
        f1 = grd_vel - mean(grd_vel);
        f2 = velTH_pulse - mean(velTH_pulse);
        PI_Kard = dot(f1,f2)/sqrt(dot(f1,f1)*dot(f2,f2));

end

%*************************************************************************
function [t_CSVx] = fn_MilestoneTime(time_series, dt, CSVx)

%--input parameters
%  time_series: velocity values of the velocity time series
%  dt: time step of the velocity time series
%  CSVx: desired milestone CSV value of time_series
%--output parameters
%  t_CSVx: time point at which CSV of time_series reaches CSVx
%-------------------------------------------------------------

        csv_working = 0;
        for i = 2 : length(time_series)
            csv_working = csv_working + dt*(time_series(i-1)^2+time_series(i)^2)/2;
            if(csv_working >= CSVx)
                t_CSVx = (i-1)*dt;     %i=1 corresponds to time equal to zero
                break;
            end
        end
            
end

