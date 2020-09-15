% Nicholas Masso
% Meteorite Impact Example Script
% Created 6/7/2020
% updated 9/8/2020 final version of function
% updated 9/14/2020 added write to csv

clc
clear


sideLength = 100; % Square sidelength, m
years = 10; % number of years to run the sim
lat = 0; % latitude of base
lon = 0; % longitude of base
filename = "out.csv";


xtmv = getImpact([0,sideLength],[0,sideLength],lat,lon,[0,years*365*24*3600],datetime("14-Dec-2024 00:00:00"));
d = getCrater(xtmv,1.83,1.6, "Granular");

labels = ["xloc (m)","yloc (m)", "t (seconds)", "mass (grams)", "vX (Km/s)","vY (Km/s)","vZ (Km/s)", "crater diameter (m)"];
%fid = fopen(filename,"w");
%fprintf(fid,"%s,",labels);
%csvwrite(filename, [labels; xtmv, d])
%fclose(fid);


fid = fopen(filename,"w");
fprintf(fid, "%s,", labels);
fprintf(fid, "\n");
out = [xtmv, d];
for i = 1:length(out(:,1))
    for j = 1:length(out(1,:))
        fprintf(fid, "%.12f,", out(i,j));
    end
    fprintf(fid,"\n");
end
fclose(fid);


figure(1)
scatter(xtmv(:,1),xtmv(:,2),d * 1000)
title('Impact Craters - Size and Location')
xlabel('x Location (m)')
ylabel('y Location (m)')
legend('Crater size (m * 1000)')

figure(2)
quiver3(xtmv(:,1),xtmv(:,2),zeros(length(xtmv(:,1)),1),xtmv(:,5),xtmv(:,6),xtmv(:,7));
title('Impactor Velocities (km/s)');
xlabel('x Location (m)')
ylabel('y Location (m)')

figure(3)
histogram(log10(xtmv(:,4)),20)
title('Impactor Mass')
xlabel('log_{10}(Mass) (g)')
ylabel('# of Occurences')