% to detect velocity pulse in the ground motion component rotated to
% a specific orientation
%--------------------------------------------------------------------------

function fn_FixOrientation_PulseDetect(path_output,pulselist_output,folder,...
                                    VelAmp_tol,R2_threshold,PsType_Array,...
                                    logPi_T_cell,Pipv_VAR_cell,...
                                    RecdNameStr,Rot_string,...
                                    grd_acc,grd_vel,grd_disp,PGV,dt,...
                                    logTs,PV_surf)
%
%The key steps in the current pulse characterization algorithm are
%(1) loop through multiple pulse types for optimal result;
%(2) subregion-to-subregion similarity between the PV spectrum surface and 
%    the Pi_PV spectrum surface is identified through traversal computation 
%    of corr2 between all allowable overlapping sub-regions of them;
%(3) if the waveform of the pulse type belongs to admissible wavelets, 
%    then pulse detection in time domain is conducted through wavelet analysis; 
%(4) if the waveform of the pulse type is not admissible wavelet, then pulse
%    detection in time domain is conducted using findsignal() in MATLAB;
%(5) for each pulse type, two potential pulses (possibly identical) are identified:
%    the one with |vp| closest to PGV (or optionally the largest |vp|) and
%    the one with the largest |CWT coefficient| if pulse waveform is admissible wavelet. 
%
%Note:
%(1) the response spectrum surface of the EQ record and the Pi-spectrum surface  
%    must have identical discretization interval along the logTs(logPi_T) axis;
%(2) the response spectrum surface of the EQ record and the Pi-spectrum surface 
%    must have identical discretization interval along the xi(Pi_xi) axis;
%(3) peak velocity vp differs from envelope amplitude Ap for M&P pulse with nu>0.
%
% by Yuchuan Tang @ Southeast University, 6/11/2022
%--------------------------------------------------------------------------

        nstep=size(grd_vel,1);      %number of time steps in the EQ record
        
        numPsType = length(PsType_Array);       %total number of pulse-type candidates

       %--get the velocity-amplified region of the PV surface of current rotated GRM
        [~,colidx]=find(PV_surf>log10(VelAmp_tol.*PGV));
        PV_VAR=PV_surf(:,colidx(1):colidx(end));
        logTs=logTs(colidx(1):colidx(end));
       
       %--Initialization for storage across multiple pulse types
        R2_star = zeros(numPsType, 2);
        BandWidth_prmt = zeros(numPsType, 2);
        vp0_star = zeros(numPsType, 1);
        CWTcoef_prmt = zeros(numPsType, 1);
        tral_X_star = zeros(numPsType, 2);
        tral_Z0_star = zeros(numPsType, 2);
        colsPV_match_star = cell(numPsType, 2);
        jtype_supstar = zeros(2,1);
        tral_X_supstar = zeros(2,1);
        tral_Z0_supstar = zeros(2,1);
        R2_supstar = zeros(2,1);
        colsPV_match_supstar = zeros(2,2);
        
       %% --select potential pulses of each pulse type
     for j_type = 1 : numPsType
           %--load the existing standard Pi-spectrum of current pulse type 
            pulse_type = PsType_Array(j_type);         %pulse type

       %--traverse all possible sub-to-sub overlapping between PV_VAR and Pipv_VAR
       clear R2_peaks ishift tral_Z0 colsPV_match;
       [~,R2_peaks,ishift,tral_Z0,~,colsPV_match] = fn_corr2Trav(PV_VAR, Pipv_VAR_cell{j_type}, R2_threshold);  %multiple candidates may exist
       if(isnan(ishift))
           disp('No pulse is detected!'); 
           continue;
       end
       
           tral_X = logTs(ishift) - logPi_T_cell{j_type}(1);

       %--potential j_type pulse corresp. to the desired |vp|  
           vp0 = zeros(size(tral_Z0));            
       
           temp = char(pulse_type);
           pulseType_split = split(temp, '_');
           if(strcmp(pulseType_split{1}, 'MP'))
                MP_gamma = str2double(pulseType_split{2});
                MP_nu = str2double(pulseType_split{3});
           else
                errordlg('Unrecognized pulse waveform!','Error');
                return;
           end
           
           if (strcmp(pulseType_split{3}, '0'))   %M&P pulse with nu=0 whose Ap=vp
                vp0 = 10.^tral_Z0;      %preliminary |vp| based on spectrum matching
           else     %M&P pulse with nu>0 whose Ap different from vp
               for k_jtype = 1:length(tral_Z0)     %loop through the candidates of j_type pulse
                    [~, ~, v_pulse_temp, ~, ~] = fn_MPPulseTH1(pulse_type,10^tral_Z0(k_jtype),10^tral_X(k_jtype),10^tral_X(k_jtype)/100);
                    vp0(k_jtype) = max(abs(v_pulse_temp));      %|vp|<|Ap| for M&P pulse with 0<nu<pi
                end
           end
             
            if (MP_gamma>1 && mod(MP_gamma,1)==0)   %M&P waveform with gamma being a natural number >1
                flag_wavelet = 1;
                wname = strcat('mp',pulseType_split{2},'-',pulseType_split{3});                                              
            elseif (mod(MP_nu,90)==0 && mod(MP_nu,180)~=0)   %M&P waveform with cos(nu)=0
                flag_wavelet = 1;
                wname = strcat('mp',pulseType_split{2},'-',pulseType_split{3});                               
            else
                flag_wavelet = 0;   %current M&P waveform is not admissible wavelet             
            end
            
            if flag_wavelet   
                [CWT_coef,CWT_vp] = fn_THDetecCWT(wname,tral_X,vp0,grd_vel,dt);
                vp0 = CWT_vp;       %supersede the preliminary values based on spectrum matching
            end
            
           [~, star] = min(abs(vp0 - PGV));   %select j_type pulse whose |vp| closest to PGV 
           tral_X_star(j_type,1) = tral_X(star);        %select 1st potential j_type pulse
           tral_Z0_star(j_type,1) = tral_Z0(star);
           vp0_star(j_type) = vp0(star);    %preliminary |vp| of 1st potential j_type pulse 
           R2_star(j_type,1) = R2_peaks(star);       
           colsPV_match_star{j_type,1} = colsPV_match(star,:);

       %--potential j_type pulse corresp. to the largest |CWT coef| if j_type is admissible wavelet       
            if flag_wavelet   
                [CWTcoef_prmt(j_type), star] = max(CWT_coef);     %largest |CWT coefficient|
                tral_X_star(j_type,2) = tral_X(star);         %select 2nd potential j_type pulse
                tral_Z0_star(j_type,2) = tral_Z0(star);  
                R2_star(j_type,2) = R2_peaks(star);                 
                colsPV_match_star{j_type,2} = colsPV_match(star,:);          
            end
                        
     end
     
     %% --select superstar candidates among all pulse types    
        [~, jtype_supstar(1)] = min(abs(vp0_star-PGV));     %1st candidate that has |vp| closest to PGV
        tral_X_supstar(1) = tral_X_star(jtype_supstar(1),1);
        tral_Z0_supstar(1) = tral_Z0_star(jtype_supstar(1),1);
        R2_supstar(1) = R2_star(jtype_supstar(1),1);  
        colsPV_match_supstar(1,:) = colsPV_match_star{jtype_supstar(1),1};
     
     if (sum(CWTcoef_prmt)>0)     %some pulse types are admissible wavelets (flag_wavelet=1)
        [~, jtype_supstar(2)] = max(CWTcoef_prmt);      %2nd candidate that has the largest |CWT coef|
        tral_X_supstar(2) = tral_X_star(jtype_supstar(2),2);   
        tral_Z0_supstar(2) = tral_Z0_star(jtype_supstar(2),2);
        R2_supstar(2) = R2_star(jtype_supstar(2),2);
        colsPV_match_supstar(2,:) = colsPV_match_star{jtype_supstar(2),2};    
     end
         
     temp = [jtype_supstar, tral_X_supstar, tral_Z0_supstar];
     [~, supstarUniq, temp1]  = unique(temp, 'rows', 'stable');   %remove repeated superstar candidate(s)
     supstarUniq_mark = zeros(size(supstarUniq));
     for k = 1 : length(temp1)
         supstarUniq_mark(temp1(k)) = supstarUniq_mark(temp1(k)) + 10^(k-1);    %indicate superstar in terms of the criterion (|vp| closest to PGV or the largest |CWT coef|)
     end
     
     %% --loop through the superstar candidates to detect them in velocity time history
     flag_outputGRM = 0;
     for j_unique = 1 : length(supstarUniq)
            i_supstar = supstarUniq(j_unique);
            if(jtype_supstar(i_supstar) == 0)     %no candidate
               continue; 
            end

            disp('-------------');            
            i_type = jtype_supstar(i_supstar);            
            pulse_type = PsType_Array(i_type);
            Dx_supstar = tral_X_supstar(i_supstar);
            Dz0_supstar = tral_Z0_supstar(i_supstar);      %based on similarity match in resp. spectrum

            %--seek pulse in time domain
            temp = char(pulse_type);
            pulseType_split = split(temp, '_');
            gamma_Mav = str2double(pulseType_split{2});
            nu_Mav = str2double(pulseType_split{3});  %in degree       
            if (gamma_Mav>1 && mod(gamma_Mav,1)==0)   %M&P waveform with gamma being a natural number >1
                flag_wavelet = 1;
                wname = strcat('mp',pulseType_split{2},'-',pulseType_split{3});                                              
            elseif (mod(nu_Mav,90)==0 && mod(nu_Mav,180)~=0)   %M&P waveform with cos(nu)=0
                flag_wavelet = 1;
                wname = strcat('mp',pulseType_split{2},'-',pulseType_split{3});                               
            else
                flag_wavelet = 0;   %current M&P waveform is not admissible wavelet             
            end

            if (flag_wavelet)        %current potential pulse is admissible wavelet
                Tp = power(10.0, Dx_supstar)     %pulse period identified
                [~,vp,t_vp,t_start,~] = fn_extract_1wavelet(grd_vel, dt, wname, Tp);   %detection through CWT
                                        %vp is peak velocity (minus sign indicates turning original waveform upside down)
                if (strcmp(pulseType_split{1}, 'MP') && strcmp(pulseType_split{3}, '0'))     %M&P pulse with nu=0 whose envelope amplitude Ap=vp
                    Ap = vp;        %unit: cm/sec
                elseif (strcmp(pulseType_split{1}, 'MP'))   %M&P pulse with nu>0 whose envelope amplitude Ap different from vp
                    temp1 = 2*pi/Tp*(t_vp-t_start);
                    Ap = 2*vp/((1-cos(temp1/gamma_Mav))*cos(temp1-gamma_Mav*pi+nu_Mav/180*pi));  %per the formula for velocity at t_vp (unit: cm/sec)
                end
                Dz_supstar = log10(abs(Ap));    %based on the result from CWT in time domain
            else        %current potential pulse is not admissible wavelet
                [J_star,VelPulse_star]=fn_THDetec3C(pulse_type,Dx_supstar,Dz0_supstar,-999,grd_vel,dt,0);   %detection thru findsignal
                
                if(isempty(J_star))     %fail to detect the potential pulse in time domain
                    continue;
                else
                 Ap = VelPulse_star(1,1);          %detected envelope amplitude of velocity pulse (unit: cm/sec)
                 Dz_supstar = log10(abs(Ap));     %may have been updated after optimization on TH match
                 Tp = VelPulse_star(1,2)          %detected velocity pulse period (unit: sec)
                 Dx_supstar = log10(Tp);          %may have been updated after optimization on TH match
                 TH_istart = VelPulse_star(1,5);     %starting index of the pulse in the EQ time history
                 TH_istop = VelPulse_star(1,6);      %stopping index of the pulse in the EQ time history
                 TH_updn = VelPulse_star(1,7);       %whether the pulse shape is upside down in the EQ time history 
                end
            end

             %--generate the time history of the detected pulse
             temp = char(pulse_type);
             switch temp(1:2)
               case 'MP'
                    [~, a_pulse, v_pulse, d_pulse, t1] = fn_MPPulseTH1(pulse_type, Ap, Tp, dt);
                    if (~flag_wavelet)        %current pulse waveform is not admissible wavelet
                        [~, temp1] = max(abs(v_pulse));
                        vp = v_pulse(temp1);
                    end
               otherwise
                    errordlg('Unrecognized pulse waveform!','Error');
                    return;
             end

             if flag_wavelet       %current pulse is admissible wavelet        
                temp = round((t_start-t1)/dt);    %t1 in local temporal coordinate should align with t_start in recd_time (possibly truncated pulse starts at 0 in local temporal coordinate) 
                TH_istart = max(1, temp);
                TH_istop = min(nstep, temp+length(v_pulse)-1);
             end
             
             d_pulse = d_pulse + (grd_disp(TH_istart) - d_pulse(1));    %shift pulse disp. to meet grd_disp(TH_istart)
             
             %--pulse indicator and early arriving check
             [PI_Shahi,Arr_Shahi,Ep_Chang,Arr_Chang,PI_Kard] = fn_PulseIndicator(grd_vel,dt,PGV,TH_istart,TH_istop,v_pulse);
               if (PI_Kard < 0.55 || Ep_Chang < 0.3)     % criterion per Tang et al. (2022)
                   continue;        %not pulse-like
               end
             
             %--open output file to record the detected pulse parameters
             fid2 = fopen(pulselist_output, 'a');
             fprintf(fid2, ['%s %s %.2f %d %s',repmat(' %.2f',1,3),repmat(' %.3f %s',1,2),' %.3f\n'], ...
                           RecdNameStr,Rot_string,PGV,supstarUniq_mark(j_unique),pulse_type,Tp,Ap,vp,...
                           PI_Shahi,Arr_Shahi,Ep_Chang,Arr_Chang,PI_Kard);
             fclose(fid2);
          
             %--create a subfolder to save output figures
             subfolder = strcat(folder,'\',Rot_string,'\',num2str(supstarUniq_mark(j_unique)));
             mkdir(subfolder);
                    
        %--plot the time histories of the EQ record vs. the detected pulse                 
        fn_PlotDetdPulse2(subfolder,RecdNameStr,grd_acc,grd_vel,grd_disp,dt,...
                          pulse_type,TH_istart,TH_istop,a_pulse,v_pulse,d_pulse);

        %---record the detected pulse and the corresponding GRM into *.AT2 and *.VT2 files
        fn_OutputTH2(path_output,RecdNameStr,subfolder,Rot_string,dt,grd_acc,grd_vel,...
                     pulse_type,TH_istart,TH_istop,a_pulse,v_pulse,flag_outputGRM);     
        flag_outputGRM = 10;
                       
     end
     
end

          