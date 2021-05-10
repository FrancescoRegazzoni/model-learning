function pgfplots_write3Ddatafile(XX,YY,ZZ,filename)

    fileID = fopen(filename,'w');    
    nx = size(XX,1);
    ny = size(YY,2);
    for iy = 1:ny
        for ix = 1:nx
            fprintf(fileID,'%0.16e %0.16e %0.16e\n',XX(ix,iy), YY(ix,iy), ZZ(ix,iy));            
        end
        if iy < ny
            fprintf(fileID,'\n');
        end
    end
    fclose(fileID);

end