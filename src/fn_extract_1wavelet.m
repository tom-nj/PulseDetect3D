function [pulse_th,vp,t_vp,t_start,coef_peak] = fn_extract_1wavelet(signal,dt,wname,Tp)
% extract one wavelet in time domain from the CWT results
%
% by Yuchuan Tang@SEU, 1/17/2022
% modified to work with waveforms whose negative side has larger amplitude, 5/21/2023
%--------------------------------------------------------------------------

    iter = 8;       %the number of iterations
    freq = 1/(Tp/dt);   %pseudo-frequency corresp. to sampling interval equal to 1 (refer to scal2frq function in MATLAB Help)
    pulse_scale = centfrq(wname, iter)/freq;

    cwt_coefs = cwt(signal, pulse_scale, wname);        % peform the Continuous Wavelet Transform
    [coef_peak, col] = max(abs(cwt_coefs));        % find peak coefficient and the corresponding column index
    coef = cwt_coefs(col);      %peak coefficient with its actual sign
    
    wtype = wavemngr('type',wname);
    switch wtype
        case 1 , [~,psi,x_vals] = wavefun(wname,iter);
        case 2 , [~,psi,~,~,x_vals] = wavefun(wname,iter);
        case 3 , [~,psi,x_vals] = wavefun(wname,iter);
        case 4 , [psi,x_vals] = wavefun(wname,iter);
        case 5 , [psi,x_vals] = wavefun(wname,iter);
    end

    basis = x_vals.*pulse_scale;
    basis = basis + (col - median(basis));
    basis = basis.*dt;
    y_vals = coef * psi/sqrt(pulse_scale);
    [temp1, temp2] = max(abs(y_vals));
    vp = y_vals(temp2);    %peak velocity (+/- sign as is) 
    t_vp = basis(temp2);    %occuring time of vp
    delta = basis(2) - basis(1);
    
    time = 0:dt:(length(signal)-1)*dt;     %must be identical to the vector, recd_time, in the main program 
    prelude = 0:delta:(basis(1)-0.00001*delta);
    num_pads = ceil((time(end)-basis(end))/delta);
    tail = (basis(end)+delta):delta:(basis(end)+delta*num_pads);
    final_basis = [prelude basis tail];
    final_yvals = [zeros(size(prelude)) y_vals zeros(1,num_pads)];
    pulse_th = interp1(final_basis,final_yvals,time);       %pulse velocity values    
    pulse_th(isnan(pulse_th)) = 0;

    t_start = basis(1);   %start of the generated wavelet in the temporal coordinate represented by the vector, time
    
end