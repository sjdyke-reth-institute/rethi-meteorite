% this is how to format a GenerateImpactData line
% this generates the required .mat files, with appropriate names
% also generates a .csv file with the same data

xlength = 200000; % meters
zlength = 200000; % meters
latitude = 0; % degrees
longitude = 0; % degrees
duration = 1/(365*24); % years
startdate = "1-Dec-2020"; % datestring

generateImpactData(zlength, xlength, latitude, longitude, duration, startdate);

clear