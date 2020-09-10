clc
clear
% updated 9/8/2020

sideLength = 100; % Square sidelength, m
years = 10; % number of years to run the sim
lat = 0; % latitude of base
lon = 0; % longitude of base

xtmv = getImpact([0,sideLength],[0,sideLength],lat,lon,[0,years*365*24*3600],datetime("14-Dec-2024 00:00:00"));
d = getCrater(xtmv,1.83,1.6, "Granular");

labels = ["xloc (m)","yloc (m)", "t (seconds)", "mass (grams)", "vX (Km/s)","vY (Km/s)","vZ (Km/s)", "crater diameter (m)"];
csvwrite( "case1.csv", [xtmv, d])

%{
fid = fopen("case1.csv","w");
fprintf(fid,"%s,",labels)
out = [xtmv, d];
for i = length(out)
    for j = length(out(1,:))
        fprintf(fid, "%f,", out(i,j))
    end
    fprintf(fid,"\n")
end
%}

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