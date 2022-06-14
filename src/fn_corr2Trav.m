%%--this file contains the following subroutines
% fn_corr2Trav()
% fn_corr2peaks()

%*************************************************************************
function [R2,R2_peaks,ishift_peaks,tral_Z0,Z_diff,colsA_match] = fn_corr2Trav(A, B, R2_threshold)
%Traversal algorithm: to compute the corr2 value for all possible global 
%and local matches between B and any submatrix Ai of A as long as size(Ai)=size(B).
%A and B have the same number of rows.
%
%by Yuchuan Tang@SEU, 6/16/2021
%--------------------------------------------------------------------------

       minBandRatio = 0.25;  %ratio of narrowest band width considered for similarity match to size(B,2)
       
       mAB = size(B,1);
       if mAB ~= size(A,1)
           errordlg('Inconsistent dimensions between matrices for similarity match!');
       end
       
       nA = size(A,2);
       nB = size(B,2);
       if nB > nA
           errordlg('The size of signal matrix exceeds the size of data matrix!');
       end
       
       A_colsum = sum(A,1);     %sum of the elements in individual columns
       A_colsqsum = sum(A.^2,1);    %square sum of the elements in individual columns
       B_colsum = sum(B,1);     %sum of the elements in individual columns
       B_colsqsum = sum(B.^2,1);    %square sum of the elements in individual columns
       
       num_shift = nA-nB+1;         %number of possible shift steps
       minBand = floor(minBandRatio*nB);     %narrowest band width considered for computating corr2
       numBandwidth = nB-minBand+1;     %number of different band widths considered
       
       R2 = cell(1,num_shift);      %cell array to store corr2 values
       for jstart = 1 : num_shift      %iteration step of jstart must be 1
           %%jstart is index of starting column of submatrix of width nB within A
           jstop = jstart+nB-1;     %index of stopping column of submatrix of width nB within A
           Ai = A(:, jstart:jstop);
           Ai_colsum = A_colsum(jstart:jstop);
           Ai_colsqsum = A_colsqsum(jstart:jstop);
           Ai_dot_B = dot(Ai,B,1);
           R2{1,jstart} = cell(1, numBandwidth);
           
           %--first case: narrowest band width considered
           iBandwidth = 1;   %index of current band width in the array [minBand:nB]
           numBand = nB-minBand+1;    %number of different bands of width minBand within B
           m_nBand = mAB*minBand;     %pre-calculated value of mAB times current band width
           
           AiBand_colsum = zeros(numBand,1);        %initialization
           AiBand_colsqsum = zeros(numBand,1);
           BBand_colsum = zeros(numBand,1);
           BBand_colsqsum = zeros(numBand,1);
           AiBand_dot_BBand = zeros(numBand,1);
           R2{1,jstart}{1,iBandwidth}=zeros(numBand,1);

           Band_head = 1;   %local column index in Ai where current band of width minBand starts
           Band_end = minBand;  %local column index in Ai where current band of width minBand ends
           AiBand_colsum(1) = sum(Ai_colsum(Band_head:Band_end));
           AiBand_colsqsum(1) = sum(Ai_colsqsum(Band_head:Band_end));
           BBand_colsum(1) = sum(B_colsum(Band_head:Band_end));
           BBand_colsqsum(1) = sum(B_colsqsum(Band_head:Band_end));               
           AiBand_dot_BBand(1) = sum(Ai_dot_B(Band_head:Band_end));                       
           numerator = AiBand_dot_BBand(1)*m_nBand - AiBand_colsum(1)*BBand_colsum(1);
           denominator = (AiBand_colsqsum(1)*m_nBand - AiBand_colsum(1)^2)*(BBand_colsqsum(1)*m_nBand - BBand_colsum(1)^2);
           R2{1,jstart}{1,iBandwidth}(1) = numerator/sqrt(denominator);      %corr2(Ai(:,1:minBand), B(:,1:minBand))
               
           for iBand = 2:numBand
               temp1 = iBand-1;
               temp2 = Band_head;
               Band_head = Band_head+1;    %local column index in Ai where current band of width minBand starts
               Band_end = Band_end+1;      %local column index in Ai where current band of width minBand ends                       
               AiBand_colsum(iBand) = AiBand_colsum(temp1) - Ai_colsum(temp2) + Ai_colsum(Band_end);
               AiBand_colsqsum(iBand) = AiBand_colsqsum(temp1) - Ai_colsqsum(temp2) + Ai_colsqsum(Band_end); 
               BBand_colsum(iBand) = BBand_colsum(temp1) - B_colsum(temp2) + B_colsum(Band_end);
               BBand_colsqsum(iBand) = BBand_colsqsum(temp1) - B_colsqsum(temp2) + B_colsqsum(Band_end);                       
               AiBand_dot_BBand(iBand) = AiBand_dot_BBand(temp1) - Ai_dot_B(temp2) + Ai_dot_B(Band_end);
               numerator = AiBand_dot_BBand(iBand)*m_nBand - AiBand_colsum(iBand)*BBand_colsum(iBand);
               denominator = (AiBand_colsqsum(iBand)*m_nBand-AiBand_colsum(iBand)^2)*(BBand_colsqsum(iBand)*m_nBand-BBand_colsum(iBand)^2);
               R2{1,jstart}{1,iBandwidth}(iBand) = numerator/sqrt(denominator);      %corr2(Ai(:,Band_head:Band_end), B(:,Band_head:Band_end))                 
           end               
           
           %--subsequent cases: increasing band width           
           for nBand = minBand+1 : 1 : nB     %nBand is band width
               iBandwidth = nBand-minBand+1;   %index of current band width nBand in the array [minBand:nB]
               numBand = nB-nBand+1;    %number of possible different bands of width nBand within B
               m_nBand = mAB*nBand;     %pre-calculated value of mAB times current band width
               
               R2{1,jstart}{1,iBandwidth}=zeros(numBand,1);     %initialization
               
               temp2 = nBand-1;
               for iBand = 1 : numBand
                   %--current band of width nBand overlays B(:,iBand:iBand+nBand-1) on Ai(:,iBand:iBand+nBand-1)
                   temp2 = temp2+1;
                   AiBand_colsum(iBand) = AiBand_colsum(iBand) + Ai_colsum(temp2);
                   AiBand_colsqsum(iBand) = AiBand_colsqsum(iBand) + Ai_colsqsum(temp2);
                   BBand_colsum(iBand) = BBand_colsum(iBand) + B_colsum(temp2);
                   BBand_colsqsum(iBand) = BBand_colsqsum(iBand) + B_colsqsum(temp2);                       
                   AiBand_dot_BBand(iBand) = AiBand_dot_BBand(iBand) + Ai_dot_B(temp2);
                   numerator = AiBand_dot_BBand(iBand)*m_nBand - AiBand_colsum(iBand)*BBand_colsum(iBand);
                   denominator = (AiBand_colsqsum(iBand)*m_nBand-AiBand_colsum(iBand)^2)*(BBand_colsqsum(iBand)*m_nBand-BBand_colsum(iBand)^2);
                   R2{1,jstart}{1,iBandwidth}(iBand) = numerator/sqrt(denominator);       %corr2(Ai(:,Band_head:Band_end), B(:,Band_head:Band_end))
               end                                                            
           end           
       end
       
       %--obtain the prominent peaks among the local max corr2 per shift of B matrix on A matrix
       flag_centercov = true;       %1-only consider the bands that cover the center of B matrix
       [R2_shiftmax,R2_shiftmax_base] = fn_corr2peaks(R2,nB,num_shift,minBand,flag_centercov);    %local max of corr2 per shift of B on A
       R2_globalmax = max(R2_shiftmax);
       
       if (R2_globalmax < R2_threshold)
           R2_peaks = NaN;
           ishift_peaks = NaN;
           tral_Z0 = NAN;
           Z_diff = NAN;
           colsA_match = NAN;
       else
           temp = [R2_shiftmax; -Inf];
           temp = max(temp, R2_threshold);   %cut off R2_shiftmax at the threshold
           [~, ishift_peaks] = findpeaks(temp);
          
          R2_peaks = R2_shiftmax(ishift_peaks);     %similarity match candidates for subsequent pulse detection
          colsA_match = zeros(length(ishift_peaks),2);      %initialization
          tral_Z0 = zeros(size(ishift_peaks));      %initialization
          Z_diff = cell(size(ishift_peaks));       %initialization
          for j = 1 : length(ishift_peaks)
              colsA_match(j,:)=(ishift_peaks(j)-1) + R2_shiftmax_base(ishift_peaks(j),:);     %start and end column indices in A corresp. to current similarity match
              Z_diff{j}=A(:,colsA_match(j,1):colsA_match(j,2))-B(:,R2_shiftmax_base(ishift_peaks(j),1):R2_shiftmax_base(ishift_peaks(j),2));    %vertical difference between the matched bands of A and horizontally shifted B
              tral_Z0(j)=mean(Z_diff{j}, 'all');      %preliminary vertical translation between the similarity match
          end
          
       end
       
end


%*************************************************************************
function [R2_shiftmax,R2_shiftmax_base] = fn_corr2peaks(R2,nB,num_shift,minBand,flag_centercov)
%Obtain the max corr2 for every shift of B matrix sliding through the columns of A matrix.
%flag_centercov=true leads to consider only the bands that cover the center of B.
%
%by Yuchuan Tang@SEU, 6/1/2021
%--------------------------------------------------------------------------

     R2_shiftmax = -9.*ones(num_shift,1);      %initialization
     R2_shiftmax_base = zeros(num_shift,2);

     for jstart = 1 : num_shift
          for nBand = minBand : nB     %nBand is band width
              iBandwidth = nBand-minBand+1;   %index of current band width nBand in the array [minBand:nB]
              numBand = nB-nBand+1;    %number of possible different bands of width nBand within B matrix
                        
              center = 0.5*nB;
              if (flag_centercov && nBand < center)                        
                   iBand_lb = ceil(center+0.5)-nBand; %lowest index among bands that cover the center of B matrix
                   iBand_ub = floor(center+0.5);    %largest index among bands that cover the center of B
              else
                   iBand_lb = 1;
                   iBand_ub = numBand;
              end
                        
              [nBand_R2max, imax] = max(R2{1,jstart}{1,iBandwidth}(iBand_lb : iBand_ub));     %local max among the considered bands of width nBand
              if(nBand_R2max > R2_shiftmax(jstart))
                    R2_shiftmax(jstart) = nBand_R2max;     %corr2 value
                    R2_shiftmax_base(jstart,1) = iBand_lb-1+imax;        %start column index in B matrix corresp. to nBand_R2max
                    R2_shiftmax_base(jstart,2) = R2_shiftmax_base(jstart,1)+nBand-1;     %end column index in B matrix corresp. to nBand_R2max
              end
              
          end                    
     end
end

