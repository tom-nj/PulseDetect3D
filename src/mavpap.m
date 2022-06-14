function [PSI,X] = mavpap(LB,UB,N,IN4,IN5) 
%MAVPAP Mavroeidis-Papageorgiou wavelet (Mavroeidis and Papageorgiou, 2003)
%   [PSI,X] = MAVPAP(LB,UB,N) returns wavelet values on an N point regular grid 
%   in the interval [LB,UB] associated with the Mav-Pap wavelet specified
%   by the character vector W.
%   N = 2^iter where iter is the number of iterations specified beforehand. 
%   W = 'mpGamma_Nu' where possible values for Gamma and Nu are:
%       Gamma is integer,  Nu is in degree
%   Output arguments are the wavelet function PSI
%   computed on the grid X, and the grid X.
%
%   This wavelet has [-pi pi] as effective support.
%
%   Implemented by Yuchuan Tang@Southeast Univ., Nov. 7, 2021
%
%   Referenced to 
%              Tang et al. A hybrid characterization framework for
%              structure-significant pulse-like features in ground motions,
%              Soil Dynamics and Earthquake Engineering, 160 (2022) 107325.
%----------------------------------------------------------------------------------------------

% Check arguments.
%-----------------
if nargin > 3
    IN4 = convertStringsToChars(IN4);
end

Gamma = 2;      %default value
Nu = 0;         %default value
nbIn = nargin;
switch nbIn
    case {0,1,2,3}
        error(message('Wavelet:FunctionInput:NotEnough_ArgNum'));

    case 5 , Gamma = IN4; Nu = IN5;
    case 4
        if ischar(IN4)
            label = deblank(IN4);
            ind   = strncmpi('mp',label,2);
            if isequal(ind,1)
                label(1:2) = [];
                len = length(label);
                if len>0
                    ind = strfind(label,'-');
                    if isempty(ind)
                        Gamma = []; % error 
                    else
                        Gamma = str2num(label(1:ind-1));
                        label(1:ind) = [];
                        Nu = str2num(label);    
                    end
                else
                    Nu = []; % error     
                end
            else
                Nu = []; % error 
            end
        else
            Nu = []; % error 
        end
end

err = isempty(Nu) || isempty(Gamma);
if ~err 
    err = ~isnumeric(Nu) || ~isnumeric(Gamma) || (Nu<0) || (Gamma<1);
end
if err
    error(message('Wavelet:WaveletFamily:Invalid_WavNum'))
end

% Compute values of the Mav.-Pap. wavelet.
X = linspace(LB,UB,N);        % wavelet support.
PSI = (1 + cos(X)) .* cos(Gamma * X + Nu/180*pi);   %Nu in degree

if (abs(Gamma)==1)
    L2norm_PSI = pi*(1.5+0.25*cos(Nu/90*pi));
elseif (abs(Gamma)==0.5)
    L2norm_PSI = pi*1.5;
else
    L2norm_PSI = 1.5*pi+3*sin(2*Gamma*pi)*cos(Nu/90*pi)/(4*Gamma*(1-Gamma^2)*(1-4*Gamma^2));
end

PSI = PSI./sqrt(L2norm_PSI);      %make L2 norm of psi(x) equal to 1

end
