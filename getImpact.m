% Nicholas Masso
% Meteorite Impact Modeling Program
% Adapted to MATLAB from Dr. Ilias Bilionis' python script
% with several additions
% Created 6/7/2020
% Edited 8/4/2020
% edit 7/17/2020: added shower modeling code
% edit 8/8/2020: added documentation

function xtmv = getImpact(areaX, areaY, latitude, longitude, timeHorizon, startDate)
    % GETIMPACT  sample for all meteor impacts in a given area in a time window
    %   xtmv = getImpact(areaX, areaY, latitude, longitude, timeHorizon)
    %       If startdate is omitted, the current date and time is used
    %
    %   xtmv = getImpact(areaX, areaY, latitude, longitude, timeHorizon, startDate)
    %       Start date is the zero on the timehorizon - if you are using a 
    %       later section of the timeHorizon, keep the startdate the same.
    %
    %   INPUTS
    %       areaX = [left, right] (meters)
    %       areaY = [top, bottom] (meters)
    %       latitude = scalar (degrees)
    %       longitude = scalar (degrees)
    %       timeHorizon = [start, end] (seconds)
    %       startDate = datetime object of start of entire test period (zero on
    %           timehorizon)
    % 	OUTPUTS
    %       returns a set of xtmv tuples (#events x 7 matrix)
    %           [xloc,yloc (meters), t (seconds), m (grams), vX,vY,vZ (Km/s)]
    % 
    %   the events can then be fed into the crater size function when its
    %   determined which part of the structure is hit
    %   the startDate is the day of "zero" time in the timeHorizon
    %   this will determine if meteor showers will be active
    
    if nargin < 6
        % if the start date is not known, it can be omitted by not
        % including it. However, the program still needs some start date,
        % and will use the computer's current time.
        startDate = datetime("today");
    end
    
    % Determine area from the given bounding box
    width = areaX(2) - areaX(1);
    height = areaY(2) - areaY(1);
    area = width * height;
    
    % Determine the total time this section of the simulation is running
    period = timeHorizon(2) - timeHorizon(1);
    
    % This section determines the distribution of objects in space.
    % The integral does not converge, so we have to pick reasonable bounds.
    % Any impactor below 10^-8 grams will leave a crater less than a
    % millimeter in diameter, even at highest velocities. 
    % 10^25 is the mass of the moon in grams (reference)
    toalMass = integral(@getFlux, 10^(-8), 10^18, 'ArrayValued', true, 'Waypoints', 10 .^ [-9:0.1:0]);
    
    % The Poisson function is XTM, and we will do the velocity part later.
    rate = area * period * toalMass;
    numEvents = poissrnd(rate);
    
    % Create an empty array to put values in
    xtmv = zeros(numEvents,7);
    
    % This is a distribution for the possible masses to consider. If this
    % is weighed towards higher values, then the program will take longer,
    % but the average mass will increase.
    pdLog = makedist('Exponential', 'mu', 2);
    
    % handle all sporadic impacts
    for i = 1:numEvents
        while true
            % Pick a mass
            m = 10^((random(pdLog) - 2.9) * 2.8);
            
            % Determine how likely it is
            ratio = getFlux(m) / toalMass;
            
            % Check if that likeliness is enough here
            u = rand;
            if u < ratio
                break
            end
        end
        
        % Pick a random point in space and time to place this impact
        x = [areaX(1) + width * rand, areaY(1) + height * rand];
        t = timeHorizon(1) +  period * rand;
        
        % Velocity and direction are dependent on location
        v = getVelocity(latitude,1);
        
        % Save the value
        xtmv(i,:) = [x,t,m,v];
        % fprintf("%d %d %d %d\n", x, t, m ,v);
    end
    
    % Now to handle all impacts that are caused by meteor showers
    
    try
        % If a CSV that denotes all the shower data is present, read it.
        % This is where more showers could be added, or values changed
        showers = importdata("showers_noSpo.csv");
    catch E
        % If not, just use the data here. 
        fprintf("\nNOTE: Shower constants file not found - using default values for the 6 most intense showers\n");
        showers.data = [41,2.46,99.4,124.9,150.4,339,-17,25,2.33;
                61,1.22,127,139.5,152,46,58,110,0.5;
                67,0.95,192.9,207.9,222.9,95,16,20,0.31;
                71,0.41,230.4,234.4,238.4,153,22,15,0.24;
                36,2.96,234.4,261.4,288.4,112,32,140,4.56;
                41.7,1.92,253,283,313,231.5,48.5,110,2.33];
        showers.textdata = {'Code','v','rho','Lstart','Lmax','Lend','RA','DEC','ZHR','m0';
            'Southern Delta Aquariids','','','','','','','','','';
            'Perseids','','','','','','','','','';
            'Orionids','','','','','','','','','';
            'Leonids','','','','','','','','','';
            'Geminids','','','','','','','','','';
            'Quadrantids','','','','','','','','',''};
    end
    
    % The best way I could figure out how to do Poisson dists for periodic
    % occurances is this:
    
    %   ^
    %   |
    % r |         ----
    % a |        |    |
    % t |      --      --
    % e | ____|          |______ 
    %   |-------------------------> time
  
    % how much to multiply the ZHR for the low and high probability zones
    lowMultiplier = 0.3;
    highMultiplier = 0.8;
    % how far from the center the cutoff for the low/hi border is
    loHiRadius = 0.5; 

    % Get the starting orbit right ascension from the date
    [startRA, ~] = getOrbitPos(startDate + seconds(timeHorizon(1)));

    e = numEvents + 1;
    
    for s = 1:length(showers.data(:,1))
        % for each shower...
        % calculate cutoffs
        sstart = deg2rad(showers.data(s,3));
        smid = deg2rad(showers.data(s,4));
        sstop = deg2rad(showers.data(s,5));
        sq1 = smid - (smid - sstart) * loHiRadius;
        sq2 = smid + (sstop - smid) * loHiRadius;
        points = [sstart, sq1, sq2, sstop];

        % check start position
        if startRA < sstart
            % outside shower
            i = 4;
        elseif startRA < sq1
            % starts in low part
            i = 1;
        elseif startRA < sq2
            % starts in high part
            i = 2;
        elseif startRA < sstop
            % starts in low part
            i = 3;
        else
            % outside shower
            i = 4;
        end

        % now to "traverse" through the total duration
        timeLo = 0;
        timeHi = 0;
        ctime = timeHorizon(1);
        loArray = [];
        lPos = 1;
        hiArray = [];
        hPos = 1;

        while ctime < timeHorizon(2)
            if i == 4
                % then there wont be any time added to the totals
                % also we are going to wraparound the year
                i = 1;
                if ctime == timeHorizon(1)
                    ctime = ctime + getTimeFromOrbitPos(startRA,sstart);
                else
                    ctime = ctime + getTimeFromOrbitPos(sstop,sstart);
                end
            else
                % then each duration will add some to the time sections
                addTime = getTimeFromOrbitPos(points(i), points(i+1));
                if ctime + addTime > timeHorizon(2)
                    addTime = timeHorizon(2) - ctime;
                end

                if i == 2
                    timeHi = timeHi + addTime;
                    hiArray(hPos,1) = ctime;
                    hiArray(hPos,2) = ctime + addTime;
                else
                    timeLo = timeLo + addTime;
                    loArray(lPos,1) = ctime;
                    loArray(lPos,2) = ctime + addTime;
                end

                ctime = ctime + addTime;
                i = i + 1;
            end
        end

        % fprintf("%s %d %d\n", string(showers.textdata(s+1,1)), timeLo, timeHi);

        % so now that we know the total time, we need to generate some random
        % events to fit in these windows. Might want to create the windows also

        scaleLo = lowMultiplier * showers.data(s,8) / 10;
        scaleHi = highMultiplier * showers.data(s,8) / 10;

        scaleLo = scaleLo * 100;
        scaleHi = scaleHi * 100;

        uBound = 10^18;

        % for lower intensity
        loMassFlux = integral(@(x)getFlux(x,scaleLo), 10^(-8), ...
            uBound, 'Waypoints', 10 .^ [-8:0.01:0]);
        numLowEvents = poissrnd(loMassFlux * timeLo * area);

        % higher intensity
        HiMassFlux = integral(@(x)getFlux(x,scaleHi), 10^(-8), ...
            uBound, 'Waypoints', 10 .^ [-8:0.01:0]);
        numHighEvents = poissrnd(HiMassFlux * timeHi * area);

        % fprintf("Events: %d %d\n", numLowEvents, numHighEvents);

        % Now to distribute them along the different time periods
        for i = 1:numLowEvents
            while true
                m = 10^((random(pdLog) - 2.9) * 2.8);
                ratio = getFlux(m) / loMassFlux;
                u = rand;
                if u < ratio
                    break
                end
            end
            x = [areaX(1) + width * rand, areaY(1) + height * rand];
            t = timeLo * rand;
            % place it at the right point in the total timehorizon
            j = 1;
            while 1==1
                treal = times(j,1) + t;
                if treal > times(j,2)
                    t = t - (times(j,2) - times(j,1));
                    j = j+1;
                else
                    break
                end
            end

            t = treal;
            v = showers.data(s,1) * getMeteorVector(startDate + seconds(t), latitude, longitude, showers.data(s,6), showers.data(s,7));
            if v(3) < 0
                xtmv(e,:) = [x,t,m,v'];
                e = e + 1;
            end
        end

        % Unfortunately, we have to do this again for the high portion. 
        % We could make an array for numHighEvents and HiMassFlux etc.
        % But for now this works
        for i = 1:numHighEvents
            while true
                m = 10^((random(pdLog) - 2.9) * 2.8);
                ratio = getFlux(m) / HiMassFlux;
                u = rand;
                if u < ratio
                    break
                end
            end
            x = [areaX(1) + width * rand, areaY(1) + height * rand];
            t = timeHi * rand;
            % place it at the right point in the timehorizon
            j = 1;
            while 1==1
                treal = times(j,1) + t;
                if treal > times(j,2)
                    t = t - (times(j,2) - times(j,1));
                    j = j+1;
                else
                    break
                end
            end

            t = treal;
            v = showers.data(s,1) * getMeteorVector(startDate + seconds(t), latitude, longitude, showers.data(s,6), showers.data(s,7));
            if v(3) < 0
                xtmv(e,:) = [x,t,m,v'];
                e = e + 1;
            end
        end
    end
end

function [results] = getVelocity(latitude, numSamples)
    % GETVELOCITY  Given a latitude, get semi-random velocity vectors
    %   [r] = GETVELOCITY(latitude, numSamples)
    %       returns a list of velocity vectors that has numSamples many
    %       entries
    %
    %       Vectors are in Km/s
    % 
    %       Higher latitudes (absolute value) have lower normal impact velocities
    %       to reflect the different average angle of incidence
    
    % create an empty array to store results
    results = zeros(numSamples,3);
    
    % Create a distribution to pull "random" values from. We are centering
    % on the mean velocity in near-moon space
    pd = makedist('Exponential','mu',17); % 17 km/s
    pdLat = makedist('Normal', 'mu', 1, 'sigma', 0.39);
    
    % We can scale this based on the latitude
    probLoc = (90 - abs(latitude) + 1) * 0.5 / 90;
    latMultiplier = icdf(pdLat,probLoc);

    for i = 1:numSamples
        z = rand(1, 3);
        z(1) = 2 * z(1) - 1;
        z(2) = 2 * z(2) - 1;
        z(3) = -1 * z(3);
        u = z / norm(z);
        v = random(pd) * u;
        v(3) = v(3) * latMultiplier;
        results(i,:) = v;
    end
end

function N = getFlux(m, scale)
    % GETFLUX  given a mass (grams) return the flux per square meter per second.
    % 	N = GETFLUX(m) returns the flux of that mass (grams)
    %
    %   N = GETFLUX(m, scale) scales the flux of that mass based on the
    %       scale given. By default the scale is 1.
    % 
    % Also works with arrays of masses. Can be integrated. Integral does
    % not converge.
    
    % Necessary to have the optional scale parameter
    if nargin < 2
        scale = 1;
    end
    
    N = m;
    
    for j = 1:length(m)
        % These are the values from the paper, and my own interpolation of
        % the data given. This is also a piecewise function.
        N_param = [-4.74951621401954e-06,-0.000166361975378819,-0.00136220536929310,0.00452350312025261,0.0246648712783819,-1.34658478023318,-14.7105254223861];
        if m(j) <= 1000
            s = polyval(N_param,log10(m(j)));
        else
            s = -15.42 - 0.8 * log10(m(j));
        end
        N(j) = 10 ^ s;
        N(j) = N(j) * scale;
    end
end


%%% METEOR SHOWER CODE %%%

function relativeMeteor = getMeteorVector(time, lat, lon, ra, dec)
    % GETMETEORVECTOR  finds the vector from a shower origin to the surface
    %   vectors = GETMETEORVECTOR(time, lat, lon, ra, dec)
    %
    %   INPUTS:
    %       time: a datetime object for the time you want to find the
    %       vector for
    %       lat: latitude of the location on the moon surface (degrees)
    %       lon: longitude of the location on the moon surface(degrees)
    %       ra: Right Ascension of the shower origin point on the celestial
    %           sphere (degrees) 
    %       dec: declination of the shower origin point on the celestial
    %           sphere (degrees)
    %
    %   OUTPUTS:
    %       relativeMeteor: a vector for a meteor impacting the surface.
    %       The z axis points upwards, the y axis points towards the north
    %       pole. The x axis points directly east.
    
    [meteor(1),meteor(2),meteor(3)] = sph2cart(deg2rad(ra), deg2rad(dec), 1);
    meteor = -1 .* meteor;
    normV = getSurfNorm(time,lat,lon);
    % now to compute the euler angles
    % roll is just how far off the z axis is from the celestial north
    % this is just the latitude, if the moon has no tilt.
    roll = deg2rad(90 - lat);
    
    pitch = 0; % we are not rotating around the y axis
    
    % now to make the y' axis point towards the north, we have to
    % point the z' 180 degrees around the z axis (make x and y both negative)
    % project it on the xy plane
    normV2 = [-normV(1), -normV(2), 0];

    % find angle between our modified z' and y axis
    oY = [0,1,0];
    yaw = atan2(norm(cross(oY,normV2)),dot(oY,normV2));
    if normV(1) < 0
        % angle between two vectors will be the smaller distance
        % around the unit circle...if we know we are > pi, we have to 
        % account for that
        yaw = 2 * pi - yaw;
    end
    
    dcm = e2m([yaw,pitch,roll],[3,2,1]);
    relativeMeteor = dcm * meteor';
end

function normV = getSurfNorm(time, lat, lon)
    % GETSURFNORM  gets the unit normal vector from the surface of the 
    % moon to the celestial sphere. Is time dependent.
    %   normV = GETSURFNORM(time, lat, lon)
    %
    %   INPUTS: 
    %       time: datetime object for the time you want to get the vector
    %       lat: latitude of the location you want the vector for (degrees)
    %       lon: longitude of the location you want the vector (degrees)
    %   
    %   OUTPUTS: 
    %       normV: vector eminating from the surface of the moon. 
    %           z points outward from the ground
    %           y points towards the north pole
    %           x points east

    zero = datetime("1970-3-21"); % Jan 1 1970 is the unix epoch, but we want the equinox
    secs = posixtime(time) - posixtime(zero);

    % siderial month is 27.321661 days long. this returns the moon to the same
    % point in the stars, aka relative to the celestial sphere.
    sidMonth = 27.321661 * 24 * 3600;

    % at the eqionox in 1970, the RA of the moon was 10h 36m 29s, and thats relative
    % to earth so we have to subtract 180 so that the normal vector points
    % outward.
    RAzero = 2 * pi * (10 + (36 + (29/60))/60)/24 - pi;

    modSecs = 2 * pi * mod(secs,sidMonth) / sidMonth;

    RAmoonSurface = mod(modSecs + RAzero, 2 * pi);
    % where are you on the moon? lets change the RA and DEC of our moon vector
    % to match your location. 

    RAyou = RAmoonSurface + deg2rad(lon);
    DECyou = deg2rad(lat);
    [x,y,z] = sph2cart(RAyou, DECyou, 1);
    normV = [x,y,z];
end

function [RA, DEC] = getOrbitPos(time)
    % GETORBITPOS  find the location of the earth around the sun.
    %   [RA, DEC] = GETORBITPOS(time)
    %       time: datetime object for the time you want to know the earths
    %       position for.
    %
    %       RA, DEC: right ascension and declination of the earth from
    %       around the sun, in degrees. declination is always zero.

    year = 365.24 * 24 * 3600;
    zero = datetime("1970-3-21"); % Jan 1 1970 is the unix epoch, but we want the equinox
    secs = posixtime(time) - posixtime(zero);
    RA = 2 * pi * mod(secs,year) / year;
    DEC = 0;
end

function seconds = getTimeFromOrbitPos(startRA, endRA)
    % GETTIMEFROMORBITPOS  from a difference in RA around the sun,
    %       determine how much time has passed. Helper function for other
    %       stuff.
    %   seconds = GETTIMEFROMORBITPOS(startRA, endRA)
    %       start and end RA are in RADIANS!!!
    %
    % Function does not handle durations over a year!!!
    
    year = 365.24 * 24 * 3600;
    if endRA < startRA
        endRA = endRA + 2 * pi;
    end
    arc = endRA - startRA;
    seconds = arc * year / (2 * pi);
end

function M = e2m(E, A)
    %==========================================================================
    %   Description:    E2M
    %
    %       Converts 3-Vector Of Body Rotation Euler Angles, E, To a 3x3
    %       Direction Cosine Matrix, M.  The Rotation Sequence Is
    %       E(1) Degrees About A(1) Axis, Then E(2) Degrees About A'(2) Axis,
    %       And Finally E(3) Degrees About A''(3) Axis.
    %
    %       Note:  If The Rotation Sequence E(1) About A(1), E(2) About A'(2),
    %       And E(3) About A''(3) Rotates Frame A Into Frame B, Then M Maps
    %       Frame A To Frame B.  M And E Represent The Orientation Of 
    %       Frame B With Respect To Frame A.
    %
    %       Examples:
    %           a)  1 -> 2 -> 3 Euler Rotation Sequence:    E = [e1 e2 e3], A = [1 2 3]
    %                   e1 Degrees About 1-Axis, Then e2 Degrees About 2'-Axis,
    %                   And Finally e3 Degrees About 3''-Axis
    %           b)  2 -> 1 -> 3 Euler Rotation Sequence:    E = [e1 e2 e3], A = [2 1 3]
    %                   e1 Degrees About 2-Axis, Then e2 Degrees About 1'-Axis,
    %                   And Finally e3 Degrees About 3''-Axis
    %           c)  3 -> 2 -> 1 Euler Rotation Sequence:    E = [e1 e2 e3], A = [3 2 1]
    %                   e1 Degrees About 3-Axis, Then e2 Degrees About 2'-Axis,
    %                   And Finally e3 Degrees About 1''-Axis
    %
    %   Inputs:
    %
    %       E           3-Vector Of Euler Body Rotation Angles (In Radians) Denoting
    %                   The Mapping Of Frame A Into Frame B Through Two Intermediate
    %                   Frames A' and A''.  Note:  Canonically, -pi/2 <= e(2) <= pi/2,
    %                   And In Any Case, There Is A Singularity At e(2) = +/- pi/2.
    %       A           3-Vector Denoting The Sequential Rotation Axes.
    %
    %   Outputs:
    %
    %       M           3x3 Direction Cosine Matrix Mapping Frame A To Frame B.
    %
    %   Reference:  Wertz, "Spacecraft Attitude Determination and Control"
    %                   Pages 763-764.
    %
    %
    %  Created By: P. R. Shah
    %
    %==========================================================================

    % --- Compute Direction Cosine Matrix Corresponding To The First Rotation

    if (A(1) == 1)
       M1 = [1     0         0
             0  cos(E(1)) sin(E(1))
             0 -sin(E(1)) cos(E(1))];
    elseif (A(1) == 2)
       M1 = [cos(E(1)) 0 -sin(E(1))
                0      1     0
             sin(E(1)) 0  cos(E(1))];
    elseif (A(1) == 3)
       M1 = [ cos(E(1)) sin(E(1)) 0 
             -sin(E(1)) cos(E(1)) 0
                 0         0      1];
    else
       M1 = eye(3);
    end

    % --- Compute Direction Cosine Matrix Corresponding To The Second Rotation

    if (A(2) == 1)
       M2 = [1     0         0
             0  cos(E(2)) sin(E(2))
             0 -sin(E(2)) cos(E(2))];
    elseif (A(2) == 2)
       M2 = [cos(E(2)) 0 -sin(E(2))
                0      1     0
             sin(E(2)) 0  cos(E(2))];
    elseif (A(2) == 3)
       M2 = [ cos(E(2)) sin(E(2)) 0 
             -sin(E(2)) cos(E(2)) 0
                 0         0      1];
    else
       M2 = eye(3);
    end

    % --- Compute Direction Cosine Matrix Corresponding To The Third Rotation

    if (A(3) == 1)
       M3 = [1     0         0
             0  cos(E(3)) sin(E(3))
             0 -sin(E(3)) cos(E(3))];
    elseif (A(3) == 2)
       M3 = [cos(E(3)) 0 -sin(E(3))
                0      1     0
             sin(E(3)) 0  cos(E(3))];
    elseif (A(3) == 3)
       M3 = [ cos(E(3)) sin(E(3)) 0 
             -sin(E(3)) cos(E(3)) 0
                 0         0      1];
    else
       M3 = eye(3);
    end

    % --- Compute Direction Cosine Matrix Corresponding To The Euler Rotation Sequence
    M = M3 * M2 * M1;
end
