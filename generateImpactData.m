% Nicholas Masso
% Meteorite Impact Dataset Generation Program
% Written for integration with DEEDS
% Created 9/14/2020
% Edit 11/16/2020: added export to .mat file 

function out = generateImpactData(xLength, yLength, latitude, longitude, duration, startDate)
    % GENERATEIMPACTDATA  take parameters for a case and generate the XTMV
    % dataset for use with DEEDS. Thsi tool is made to be used on the
    % command line and only 
    %   out = generateImpactData(xLength, yLength, lat, lon, duration, startDate, filename)
    %       creates a csv file in the current directory with the outputs of the sim.
    %       returns 1 on proper exit. The filename is auto-generated.
    %       
    %       also generates a .mat file with a timeseries of the same data.
    %
    %   INPUTS
    %       xLength = area bounding box x side length (m)
    %       yLength = area bounding box y side length (m)
    %       latitude = location of base (degrees)
    %       longitude = location of base (degrees)
    %       duration = sample duration in years
    %       startDate = ISO String of the start date of the case

    
    
    % convert the string into a datetime object
    d = datetime(startDate);

    
    % run a simulation for these conditions
    xtmv = getImpact([0,xLength],[0,yLength],latitude,longitude,[0,duration*365*24*3600],d);
    
    % prepare the data to be written to a file
    labels = ["xloc (m)","yloc (m)", "t (seconds)", "mass (grams)", "vX (Km/s)","vY (Km/s)","vZ (Km/s)"];
    outdata = [xtmv];
    filename = "xtmv_loc_" + latitude + "_" + longitude + "_year_" + year(d) + "_duration_" + duration;
    
    % write to .mat file
    ts = timeseries(xtmv, xtmv(:,3));
    save(filename, 'ts', '-v7.3');
    
    % Write the data to csv file
    fid = fopen(filename + ".csv","w");
    fprintf(fid, "%s,", labels);
    fprintf(fid, "\n");
    for i = 1:length(outdata(:,1))
        for j = 1:length(outdata(1,:))
            fprintf(fid, "%.12f,", outdata(i,j));
        end
        fprintf(fid,"\n");
    end
    fclose(fid);
    
    out = 1;
end
