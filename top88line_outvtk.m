% Function to output data from top88 as vtk file
% Note that the structure will be rotated anti-clockwise by 90deg
function top88line_outvtk(name,h,nely,nelx,density,sens,disp)
 of = fopen(name,'w'); % Open file for writing
 %% Print header and grid info 
 fprintf(of,'# vtk DataFile Version 3.0\nPara0\nASCII\nDATASET RECTILINEAR_GRID\n');
 fprintf(of,'DIMENSIONS %i %i 1\n',nelx+1,nely+1);
 fprintf(of,'X_COORDINATES %i float\n',nelx+1);
 for i=0:nelx
     fprintf(of,'%d ',h*i);
 end
 fprintf(of,'\nY_COORDINATES %i float\n',nely+1);
 for i=0:nely
     fprintf(of,'%d ',h*i);
 end
 fprintf(of,'\nZ_COORDINATES 1 int\n0\n');
 %% Print density values
 fprintf(of,'\nCELL_DATA %i\nSCALARS density float 1\nLOOKUP_TABLE default\n',nelx*nely);
 fprintf(of,'%d\n',density);
 %% Print objective sensitivity values
 if(length(sens)>1)
    fprintf(of,'\nSCALARS sens double 1\nLOOKUP_TABLE default\n',nelx*nely);
    fprintf(of,'%12.4e\n',sens);
 end
 %% Print displacements
 if(length(disp)>1)
     nnode = (nelx+1)*(nely+1);
     fprintf(of,'\n\nPOINT_DATA %i\nVECTORS disp double\n',nnode);
     for i=1:nnode
         j = i*2;
         fprintf(of,'%12.4e %12.4e 0.0\n',-disp(j),disp(j-1));
     end
 end
 fclose(of); % Close file
end