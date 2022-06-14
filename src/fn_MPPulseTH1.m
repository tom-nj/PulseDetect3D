% generate the time history of Mavroeidis&Papageorgiou pulse
%
% by Yuchuan Tang@SEU, 6/11/2022
%--------------------------------------------------------------------------

function [t_pulse, a_pulse, v_pulse, d_pulse, t_start] = fn_MPPulseTH1(pulse_type, A, Tp, dt)

%---input parameters:
% pulse type
% A: M&P velocity pulse envelope amplitude including +/- sign (unit: cm/sec)
% Tp: pulse period (unit: sec)
% dt: time step of the time history (unit: sec)

%---output parameters in units of cm, sec

temp = char(pulse_type);
temp = split(temp, '_');
if (strcmp(temp{1}, 'MP'))
    gamma = str2double(temp{2});
    nu = str2double(temp{3})*(pi/180);     % unit: rad
else
    errordlg(strcat("Pulse type ",pulse_type, " is not recognized for generating TH!"), 'Error');
    return;
end

%%---- generate pulse      
if(gamma >= 1)
    fp = 1.0/Tp;
    tend = ceil(gamma/fp/dt)*dt;    % unit: sec
    
    t_pulse=[0 : dt : tend]';       % time steps (local time coordinate of the pulse)
    t_epoch = gamma/(2*fp);         % unit: sec  
    
    nstep=size(t_pulse,1);
    a_pulse=zeros(nstep,1);
    v_pulse=zeros(nstep,1);
    d_pulse=zeros(nstep,1);    
             
        fp = 1/Tp;            % unit: Hz
        
        nstep=length(t_pulse);
        for i=1:nstep
        
            t=t_pulse(i);
            
            if ( t < t_epoch-gamma/(2*fp) )
                if (gamma~=1)
                    d_pulse(i)=A/(4*pi*fp*(1-gamma^2))*sin(nu-pi*gamma); 
                else
                    d_pulse(i)=0;
                end
            elseif ( t > t_epoch+gamma/(2*fp) )
                if (gamma~=1) 
                    d_pulse(i)=A/(4*pi*fp*(1-gamma^2))*sin(nu+pi*gamma); 
                else
                    d_pulse(i)=A/(4*fp)*cos(nu);
                end
            else
                    temp=sin( 2*pi*fp/gamma*(t-t_epoch) )*cos(2*pi*fp*(t-t_epoch)+nu)...
                             +gamma*sin(2*pi*fp*(t-t_epoch)+nu)*( 1+cos(2*pi*fp/gamma*(t-t_epoch)) );
                    a_pulse(i) = -A*pi*fp/gamma*temp;
                    
                    v_pulse(i) = A/2*( 1+cos(2*pi*fp/gamma*(t-t_epoch)) )*cos(2*pi*fp*(t-t_epoch)+nu);
                    
                 if (gamma~=1)   
                    temp=sin( 2*pi*fp*(t-t_epoch)+nu )+gamma/(2*(gamma-1))*sin(2*pi*fp*(gamma-1)/gamma*(t-t_epoch)+nu)...
                             +gamma/(2*(gamma+1))*sin(2*pi*fp*(gamma+1)/gamma*(t-t_epoch)+nu);
                 else
                    temp=sin(2*pi*fp*(t-t_epoch)+nu)+sin(4*pi*fp*(t-t_epoch)+nu)/4+3*sin(nu)/4+(2*pi*fp*(t-t_epoch)+pi)*cos(nu)/2;
                 end
                    d_pulse(i) = A/(4*pi*fp)*temp;                    
            end

        end
      
      d_pulse = d_pulse - d_pulse(1);   %enforce zero disp. at pulse beginning
      t_start = t_pulse(1);     %start of the pulse in local temporal coordinate t_pulse
        
end

end