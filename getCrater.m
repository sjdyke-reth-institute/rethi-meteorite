% Nicholas Masso
% Crater Sizing Equation (basic)
% Created 6/10/2020
% Edited 8/8/2020

function d = getCrater(xtmv, impactorDensity, surfaceDensity, surfaceType)
    % GETCRATER  determine crater diameter based on several parameters.
    %   d = GETCRATER(xtmv, impactorDensity, surfaceDensity, surfaceType)
    %
    %   INPUTS:
    %       xtmv: the xtmv vector retrieved from GETIMPACT()
    %       impactorDensity: density of impacting body in g/cm^3
    %       surfaceDensity: density of surface in g/cm^3
    %       surfaceType: can be "Solid" or "Granular"
    %
    %   OUTPUTS:
    %       d: diameter of crater in meters.
    %           returns NaN on incorrect parameters.
    % 
    % Details on depth are unknown, but radius could probably be used.

    d = zeros(length(xtmv(:,1)),1);
    
    if surfaceType ~= "Solid" && surfaceType ~= "Granular"
        d = NaN;
        return
    end
    
    b = [0,0,-1];
    
    for i = 1:length(xtmv(:,1))
        a = xtmv(i,5:7);
        v = norm(a);
        theta = 90 - atan2(norm(cross(a,b)), dot(a,b));
        if theta < 0 || theta > 90
            d(i) = NaN;
        else
            KE = 0.5 * xtmv(i,4) / 1000 * v ^2;
            if surfaceType == "Granular"
                d(i) = 0.25 * impactorDensity^(1/6) * surfaceDensity ^ (-1/2) * 1.67^ (-0.165) * ...
                    KE ^ 0.29 * sin(theta) ^ (1/3);
            else
                d(i) = 0.015 * impactorDensity^(1/6) * surfaceDensity ^ (-1/2) * 1.67^ (-0.165) * ...
                    KE ^ 0.37 * sin(theta) ^ (2/3);
            end
        end
    end
end
