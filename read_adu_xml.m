function read_adu_xml
% read head file of the ADU07 MT system, from which coordinates are picked out. 
% Then, a .txt file will be generated in the pwd.
nxml = dir([pwd,'\*.xml']);

for isite = 1:length(nxml)
    siten = nxml(isite).name;
    fid = fopen(siten);
    sline = fgetl(fid);
    while 1
        sline = fgetl(fid);        
        if ~ischar(sline); break; end
        sline = strtrim(sline);
        if strncmp(sline,'<Height>',8)
            dind = isstrprop(sline,'digit');
            elev = str2double(sline(dind));
            elev = elev/100; % unit:m
        end
        
        if strncmp(sline,'<Latitude>',10)
            dind = isstrprop(sline,'digit');
            lat = str2double(sline(dind));
            lat = lat/3600/1000; % unit:degree
        end
        
        if strncmp(sline,'<Longitude>',11)
            dind = isstrprop(sline,'digit');
            lon = str2double(sline(dind));
            lon = lon/3600/1000; % unit: degree
        end
    end
    if isite == 1
        fids = fopen('coord_aduxml.txt','wt');
        fprintf(fids,'%s\n','#HEAD ID  LON LAT ELEV');
    else
        fids = fopen('coord_aduxml.txt','a+t');
    end
    fprintf(fids,'%s  ',siten(1:end-4));
    fprintf(fids,'%12.8f  %12.8f  %10.3f\n',[lon, lat, elev]);
    
    if isite == length(siten)
        fclose(fids);
    end
end

end
