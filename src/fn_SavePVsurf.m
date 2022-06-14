function fn_SavePVsurf(filename,logTs,damping,PV_surf)
%--save pseudo-velocity response surface to file

        fid=fopen(filename, 'w');
        fprintf(fid, '%e ', logTs);
        fprintf(fid, '\n');
        fprintf(fid, '%e ', damping);
        fprintf(fid, '\n');
        for j = 1 : length(damping)
            fprintf(fid, '%e ',PV_surf(j,:) );
            fprintf(fid, '\n'); 
        end
        fclose(fid);
end

