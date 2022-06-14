%This file contains the subroutines:
%  fn_THDetec3C()
%  fn_THsimilar2E()
%  fun_InterpVelTH1()
%
%**************************************************************************
function [J_star,VelPulse_star_itype]=fn_THDetec3C(pulse_type,Dx_uniq,Dy_uniq,...
                                            dist_uniq,grd_vel,dt,flag_optm)

% to detect i_type pulse in the velocity time history of the EQ record
%
% by Yuchuan Tang@SEU, 6/11/2022
%------------------------------------------------------------------------

    J_star = [];
    VelPulse_star_itype = [];
    
    nPnt_grm = length(grd_vel);
    [temp1, temp2] = max(grd_vel);
    [temp3, temp4] = min(grd_vel);
    PGV = max(temp1, abs(temp3));
    index_core = sort([temp2, temp4]);
    
    for j_seed = length(Dx_uniq) : -1 : 1      %check the possible matches from right to left for velocity pulse
        
        Tp = power(10.0, Dx_uniq(j_seed));     %preliminary pulse period in s  (10.^Dx_uniq(j_seed))
        Ap = power(10.0, Dy_uniq(j_seed));   %preliminary envelope amplitude of pulse velocity in cm/s, sign undetermined yet
    
        %%---generate the time history of the corresponding pulse
        temp = char(pulse_type);
        switch temp(1:2)
            case 'MP'
                [~, ~, v_pulse, ~, ~] = fn_MPPulseTH1(pulse_type, Ap, Tp, dt);
            otherwise
                errordlg('Unrecognized pulse waveform!','Error');
                return;
        end
        
        if(length(grd_vel) < length(v_pulse))
            continue
        end
        
        %---take a ground velocity window for pulse detection
        %---the window extends length(v_pulse) before index_core(1) as well as after index_core(2)
        grd_vel_window = [];
        nPnt_pulse = length(v_pulse);
        temp = nPnt_pulse - (index_core(1)-1);
        if(temp > 0)    %need to pad zeros before grd_vel(1) to fullfill length(v_pulse) before index_core(1)
            nPnt_pad0Bef = min(temp, floor(nPnt_pulse/3));      %number of zeros to pad before grd_vel(1), capped
            grd_vel_window = [zeros(nPnt_pad0Bef,1); grd_vel(1:index_core(2))];
            index_windowL = 1 - nPnt_pad0Bef;   %pseudo-index of grd_vel corresponding to grd_vel_window(1)
        else
            index_windowL = index_core(1)-nPnt_pulse;   %index of grd_vel corresponding to grd_vel_window(1)
            grd_vel_window = grd_vel(index_windowL : index_core(2));
        end
        temp = nPnt_pulse - (nPnt_grm-index_core(2));
        if(temp>0)
            grd_vel_window = [grd_vel_window; grd_vel(index_core(2)+1:end)];  
        else
            grd_vel_window = [grd_vel_window; grd_vel(index_core(2)+1:index_core(2)+nPnt_pulse)];  
        end
        
        %---detect v_pulse in grd_vel_window
        [TH_status,TH_istart,TH_istop,TH_updn,TH_dist,Ap_optm,Tp_optm,t_epoch_optm] = ...
                             fn_THsimilar2E(grd_vel_window,v_pulse,dt,pulse_type,Ap,Tp,-999,flag_optm);
        
        if( TH_status > 0 )     %potential pulse detected in the EQ velocity time history 
            J_star = [J_star, j_seed];
            TH_istart = max(index_windowL+TH_istart-1, 1);         %convert index of grd_vel_window to index of grd_vel
            TH_istop = index_windowL+TH_istop-1;        %convert index of grd_vel_window to index of grd_vel
            VelPulse_star_itype = [VelPulse_star_itype; ...
                                   Ap_optm,Tp_optm,t_epoch_optm,dist_uniq(j_seed),TH_istart,TH_istop,TH_updn,TH_dist];           
            
        end
    end

end

%**************************************************************************

function [status,TH_istart,TH_istop,TH_updn,TH_dist,Ap_optm,Tp_optm,t_epoch_optm] = ...
                 fn_THsimilar2E(grm,pulse,dt,pulse_type,Ap,Tp,t_epoch,flag_optm)

%to detect similarity in time history between pulse and grm (ground motion)
%Output: 
%   status (1--sucess; -1,-2--failures for different causes);
%   Ap_optm,Tp_optm,t_epoch_optm are final pulse parameters w.r.t. velocity TH match
%when the input flag_optm is zero (optimization NOT triggered in this subroutine):
%   TH_updn*pulse(end-(TH_istop-TH_istart):end) is similar to grm(TH_istart:TH_istop)
%           where the pulse may be truncated at beginning
%   TH_dist is the Euclidean distance, normalized by the pulse length, between 
%           TH_updn*pulse(end-(TH_istop-TH_istart):end) and grm(TH_istart:TH_istop).
%when the input flag_optm is nonzero (optimization triggered in this subroutine):
%   v_pulse(end-(TH_istop-TH_istart):end) is similar to grm(TH_istart:TH_istop)
%           where the pulse may be truncated at beginning
%   TH_dist is the Euclidean distance, normalized by the pulse length, between 
%           v_pulse(end-(TH_istop-TH_istart):end) and grm(TH_istart:TH_istop).
%
%NOTE:
%1)loop through various metric options in findsignal()
%2)ignore "pulse" that does not cover the PGV-time of the ground motion
%
% by Yuchuan Tang @SEU,11/12/2021
%--------------------------------------------------------------------------

if(length(grm) < length(pulse))
    status = -1;       %the duration of the grm is shorter than the pulse (Unreasonable!)
    TH_istart = 0;      %assign fake values as follows
    TH_istop = 0;
    TH_dist = 1.0E9;
    TH_updn = 0;
    Ap_optm = -999;
    Tp_optm = -999;
    t_epoch_optm = -999; 

else
    
    [grm_max, i_max] = max(grm);    
    [grm_min, i_min] = min(grm);
    if(grm_max > abs(grm_min))
        PGV = grm_max;
        i_PGV = i_max;
    else
        PGV = abs(grm_min);
        i_PGV = i_min;        
    end
    
    istart = zeros(2,1);
    istop = zeros(2,1);
    dist = zeros(2,1);

    key = ["squared", "euclidean", "absolute"];        %metric options for findsignal function
        
  for j_key = 1 : length(key)               
    [istart(1),istop(1),dist(1)] = findsignal(grm, -1.0.*pulse, 'Metric', key(j_key));    %pulse shape upside down
    [istart(2),istop(2),dist(2)] = findsignal(grm, pulse, 'Metric', key(j_key));          %normal pulse shape   
    
    %--check whether the matched pulse covers the global peak or valley of grm
        if( istart(1)<i_PGV && i_PGV<istop(1) )
            flag1 = true;
        else
            flag1 = false;
        end

        if( istart(2)<i_PGV && i_PGV<istop(2) )
            flag2 = true;
        else
            flag2 = false;
        end
    
        if(flag1 && ~flag2)
            status = 1;
            i_updn = 1;
            break;
        elseif(~flag1 && flag2)
            status = 1;
            i_updn = 2;
            break;
        elseif(flag1 && flag2)
            status = 1;
            [~, i_updn] = min(dist);            
            break;
        else
            status = -2;       %neither +/-pulse covers the global peak/valley
        end

  end
  
  if(status < 0)        %no pulse detected in the EQ velocity time history
    TH_istart = 0;      %assign fake values as follows
    TH_istop = 0;
    TH_dist = 1.0E9;
    TH_updn = 0;
    Ap_optm = -999;
    Tp_optm = -999;
    t_epoch_optm = -999; 
    
  else         %potential pulse detected in the EQ velocity time history
    TH_istart = istart(i_updn);
    TH_istop = istop(i_updn);
    TH_updn = (-1)^i_updn;       %+1: normal pulse shape; -1: upside down
    Ap = TH_updn*Ap;
    [~, temp] = max(abs(pulse));   %locate max absolute value of velocity
    vp = pulse(temp);        %peak velocity with sign
    
    if( ~flag_optm )   %--No optimization on pulse parameters w.r.t. velocity TH match

      temp = TH_updn.*pulse(end-(TH_istop-TH_istart):end) - grm(TH_istart:TH_istop);
      TH_dist = median(abs(temp));      %time history distance measure
      
      Ap_optm = Ap;       %already updated the sign of Ap
      Tp_optm = Tp;       %keep the input value 
      t_epoch_optm = t_epoch;  %keep the input value
      
    else       %--Optimization on pulse parameters w.r.t. velocity TH match
      nstep = length(grm);
      grm_time = dt.*[0:1:nstep-1]';
     
      fun_Vel=@(paramt,grm_time)fun_InterpVelTH1(paramt, grm_time, grm, pulse_type);
      paramt0=[Ap, Tp, grm_time(TH_istop)];       %starting from the already detected pulse
      
      if(vp>0)
          Ap_lb = min(0.8*Ap, 0.9*(PGV/vp)*Ap);     %lower bound of Ap for optimization
          Ap_ub = min(1.2*Ap, (PGV/vp)*Ap);     %upper bound of Ap for optimization
      elseif(vp<0)
          Ap_lb = max(1.2*Ap, (-PGV/vp)*Ap);     %lower bound of Ap for optimization
          Ap_ub = max(0.8*Ap, (-0.9*PGV/vp)*Ap);     %upper bound of Ap for optimization
      else
          msgbox("vp value is zero!!");
      end
      lb = [Ap_lb, 0.8*Tp, max(grm_time(TH_istop)-min(Tp,5),dt*(length(pulse)-1))];       %lower bound for optimization
      ub = [Ap_ub, 1.2*Tp, min(grm_time(TH_istop)+min(Tp,5),grm_time(end))];       %upper bound for optimization
      options = optimoptions('lsqcurvefit', 'Display','off');      
      [paramt,~] = lsqcurvefit(fun_Vel, paramt0, grm_time, grm, lb, ub, options);

      Ap_optm = paramt(1);    %optimized envelope amplitude of velocity pulse w.r.t. velocity TH match
      Tp_optm = paramt(2);    %optimized pulse period w.r.t. velocity TH match
      time_pulseend = paramt(3);    %optimized ending time of pulse in the coordinate of grm_time
      TH_istop = find(grm_time<=time_pulseend, 1, 'last');     %time step in grm_time corresponding to time_pulseend
      
      temp = char(pulse_type);
      switch temp(1:2)
        case 'MP'
            [~, ~, v_pulse, ~, ~] = fn_MPPulseTH1(pulse_type, Ap_optm, Tp_optm, dt);
        otherwise
            errordlg('Unrecognized pulse waveform!','Error');
            return;
      end     
      TH_istart = TH_istop-length(v_pulse)+1;
      if(TH_istart<=0)       %grm does not contain the entire pulse waveform
          TH_istart = 1;
      end
      
      temp = v_pulse(end-(TH_istop-TH_istart) : end) - grm(TH_istart:TH_istop);
      TH_dist = median(abs(temp));      %time history distance measure

    end

  end

end

end

%**************************************************************************
function [vel_new]=fun_InterpVelTH1(paramt, grm_time, grm, pulse_type)

%--this subroutine is called in the subroutine fn_THDetec3C()
%  for refining the pulse time history match to the ground velocity TH
%
% by Yuchuan Tang@SEU,11/12/2021
%-------------------------------------------------------------------------

Ap=paramt(1);    % envelope amplitude of velocity pulse (cm/s)
Tp=paramt(2);    % velocity pulse period (s)
time_pulseend=paramt(3);     % pulse ending time in the coordinate of grm_time

dt_pulse = min(grm_time(2)-grm_time(1), 0.005);       %time step of pulse

temp = char(pulse_type);
switch temp(1:2)
   case 'MP'
        [t_pulse, ~, v_pulse, ~, ~] = fn_MPPulseTH1(pulse_type, Ap, Tp, dt_pulse);
   otherwise
        errordlg('Unrecognized pulse waveform!','Error');
        return;
end
time_pulse = t_pulse - t_pulse(end) + time_pulseend;    %pulse time range in the coordinate of grm_time

vel_new=grm;

index1 =find(grm_time >= time_pulse(1), 1); 
index2 =find(grm_time <= time_pulse(end), 1, 'last'); 
if(isempty(index1) || isempty(index2))
   msgbox('Pulse does not overlap grm in THDetec Optimization!'); 
   return;
end

vel_new(index1:index2)=interp1(time_pulse, v_pulse, grm_time(index1:index2), 'linear');
    
end
