function [CWT_coefabs,CWT_vp_abs]=fn_THDetecCWT(wname,Dx,vp0,grd_vel,dt)

%compute CWT coefficients for the pulse candidates associated with waveform 
%of wname based on the pulse periods identified from spectrum similarity
%
%Include the tolerance parameter:
%tol_vpdiff (relative difference in |vp| between the CWT result and the preliminary 
%            estimation based on the vertical shift Dz of response spectrum)
%
% by YTang@SEU, 12/13/2021
%--------------------------------------------------------------------------
    
    tol_vpdiff = 0.3;
    
    CWT_coefabs = zeros(size(Dx));   %initialization
    CWT_vp_abs = zeros(size(Dx));    %initialization
    
    for j_seed = length(Dx) : -1 : 1        %loop through multiple candidates based on sub-to-sub similarity in spectrum
        
        Tp = power(10.0, Dx(j_seed));         %preliminary pulse period in s
        vp_abs0 = vp0(j_seed);    %preliminary estimate of peak velocity |vp| in cm/s (sign undetermined yet)

        iter = 8;       %the number of iterations
        freq = 1/(Tp/dt);   %pseudo-frequency corresp. to sampling interval equal to 1 (refer to scal2frq function in MATLAB Help)
        pulse_scale = centfrq(wname, iter)/freq;
        cwt_coefs = cwt(grd_vel, pulse_scale, wname);        % peform the Continuous Wavelet Transform
        coef_absmax = max(abs(cwt_coefs));        % find the absolute max coefficient
        
        wtype = wavemngr('type',wname);
        switch wtype
            case 1 , [~,psi,~] = wavefun(wname,iter);
            case 2 , [~,psi,~,~,~] = wavefun(wname,iter);
            case 3 , [~,psi,~] = wavefun(wname,iter);
            case 4 , [psi,~] = wavefun(wname,iter);
            case 5 , [psi,~] = wavefun(wname,iter);
        end
        vp_abs = coef_absmax * max(abs(psi))/sqrt(pulse_scale);      %velocity pulse amplitude (absolute value)
        
        diff = abs(vp_abs/vp_abs0 - 1);     %relative difference from the prelimiary estimate of |vp|
        if(diff < tol_vpdiff)
            CWT_coefabs(j_seed) = coef_absmax;
            CWT_vp_abs(j_seed) = vp_abs;
        else
%             msgbox('Big difference in |vp| between the spectrum matching and the CWT!');
            continue;
        end
    end
        
end
