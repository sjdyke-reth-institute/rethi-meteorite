# RETHi Meteorite Impact Modeling
A set of MATLAB tools for retrieving impacts on the Lunar surface.

## Usage

### For use in Scripts and SIMULINK Models
The tools were originally created as pure MATLAB functions that could be used in a larger script or SIMULINK model. These can be simply downloaded and added in the SIMULINK model, taking time inputs.

See the `testcase.m` file for example usage.

### For use in Dataset Generation
As RETHi transitions to using the DEEDS database system, a script is used to retrieve various parameters from the function and export a dataset directly into the database.

This is currently a work in progress.

## Documentation
### getImpact
GETIMPACT  sample for all meteor impacts in a given area in a time window
  `xtmv = getImpact(areaX, areaY, latitude, longitude, timeHorizon)`
      If startdate is omitted, the current date and time is used

  `xtmv = getImpact(areaX, areaY, latitude, longitude, timeHorizon, startDate)`
      Start date is the zero on the timehorizon - if you are using a 
      later section of the timeHorizon, keep the startdate the same.

  INPUTS
      areaX = [left, right] (meters)
      areaY = [top, bottom] (meters)
      latitude = scalar (degrees)
      longitude = scalar (degrees)
      timeHorizon = [start, end] (seconds)
      startDate = datetime object of start of entire test period (zero on
          timehorizon)
	OUTPUTS
      returns a set of xtmv tuples (#events x 7 matrix)
          [xloc,yloc (meters), t (seconds), m (grams), vX,vY,vZ (Km/s)]

  the events can then be fed into the crater size function when its
  determined which part of the structure is hit
  the startDate is the day of "zero" time in the timeHorizon
  this will determine if meteor showers will be active
    
### getVelocity
GETVELOCITY  Given a latitude, get semi-random velocity vectors
  `[r] = GETVELOCITY(latitude, numSamples)`
      returns a list of velocity vectors that has numSamples many
      entries

      Vectors are in Km/s

      Higher latitudes (absolute value) have lower normal impact velocities
      to reflect the different average angle of incidence
    

GETFLUX  given a mass (grams) return the flux per square meter per second.
	`N = GETFLUX(m)` returns the flux of that mass (grams)

    `N = GETFLUX(m, scale)` scales the flux of that mass based on the
      scale given. By default the scale is 1.

Also works with arrays of masses. Can be integrated. Integral does
not converge.
    
### getMeteorVector
GETMETEORVECTOR  finds the vector from a shower origin to the surface
  `vectors = GETMETEORVECTOR(time, lat, lon, ra, dec)`

  INPUTS:
      time: a datetime object for the time you want to find the
      vector for
      lat: latitude of the location on the moon surface (degrees)
      lon: longitude of the location on the moon surface(degrees)
      ra: Right Ascension of the shower origin point on the celestial
          sphere (degrees) 
      dec: declination of the shower origin point on the celestial
          sphere (degrees)

  OUTPUTS:
      relativeMeteor: a vector for a meteor impacting the surface.
      The z axis points upwards, the y axis points towards the north
      pole. The x axis points directly east.

### getSurfNorm
GETSURFNORM  gets the unit normal vector from the surface of the 
moon to the celestial sphere. Is time dependent.
  `normV = GETSURFNORM(time, lat, lon)`

  INPUTS: 
      time: datetime object for the time you want to get the vector
      lat: latitude of the location you want the vector for (degrees)
      lon: longitude of the location you want the vector (degrees)
  
  OUTPUTS: 
      normV: vector eminating from the surface of the moon. 
          z points outward from the ground
          y points towards the north pole
          x points east

### getOrbitPos
GETORBITPOS  find the location of the earth around the sun.
  `[RA, DEC] = GETORBITPOS(time)`
      time: datetime object for the time you want to know the earths
      position for.

      RA, DEC: right ascension and declination of the earth from
      around the sun, in degrees. declination is always zero.

### getTimeFromOrbitPos
GETTIMEFROMORBITPOS  from a difference in RA around the sun,
      determine how much time has passed. Helper function for other
      stuff.
  `seconds = GETTIMEFROMORBITPOS(startRA, endRA)`
      start and end RA are in RADIANS!!!

Function does not handle durations over a year!!!


