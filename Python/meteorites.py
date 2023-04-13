"""Stochastic meteorite-impact event model

Author:
    Ilias Bilionis
    R Murali Krishnan

Date:
    05.15.2021
    03.01.2023

"""

import numpy as np
import operator as op
import numpy.typing as npt
import scipy.stats as st
import scipy.integrate as integrate

from typing import Union, List, Callable, Dict, NamedTuple
from pprint import pprint

class ImpactEvent(NamedTuple):
    t: float
    x: npt.NDArray
    m: float
    v: npt.NDArray

    def __str__(self):
        _str_ = f"""
Impact at time: {self.t: 0.2f},
location: {self.x},
of mass: {self.m:0.2e},
with velocity: {self.v} [{np.linalg.norm(self.v):.2f}]
"""
        return _str_


Number = Union[int, float]

def sample_velocity(speed_average: Number, size: Number=10) -> npt.NDArray:
    """
    Sample velocity vector

    :param speed_average:   The observed average speed
    :param size:            Number of samples
    """
    z = np.random.rand(size, 3)
    z[:, :2] = 2. * z[:, :2] - 1.
    z[:, 2] = -2. * z[:, 2]
    u = z / np.linalg.norm(z, axis=1)[:,None]
    v_s = st.expon(scale=speed_average).rvs(size)[:,None]
    v = v_s * u
    return v


def fn_mass_flux(m: Number) -> Number:
    """According to data, the rate of meteorite impact is the mass flux
    at the Moon's surface

    Calculates the flux of meteorite with a particular mass in hr^-1 * m^-2
    
    :param m:   Mass of the meteorite in grams (gm)
    """

    def model_logN_small_mass(m) -> Number:
        """Flux model for small meteorites

        Reference:
            Grun et al (1985), https://doi.org/10.1016/0019-1035(85)90121-6
        """
        theta = [
            -4.74951621401954e-06,
            -0.000166361975378819,
            -0.00136220536929310,
            0.00452350312025261,
            0.0246648712783819,
            -1.34658478023318,
            -14.7105254223861
        ]
        return np.polyval(theta, np.log10(m))

    def model_logN_large_mass(m):
        """Flux model for large meteorites
        
        Reference:
            Dycus, Robert (1969), https://www.jstor.org/stable/40674750
        """
        return -15.42 - 0.8 * np.log10(m)

    logN = model_logN_large_mass if m > 1000 else model_logN_small_mass
    return 10 ** logN(m) * 3600 # s^-1 * m^-2 * s * hr^-1


def sample_xtmv(
    area_box: npt.NDArray, 
    time_horizon: List[float], 
    speed_average: float, 
    mass_flux_multiplier: float=1.0,
    verbose: bool=False,
    print_interval: float=1000) -> List[ImpactEvent]:
    """Sample impact events given parameters about the location
    
    :param area_box:            Area considered for the habitat
    :param time_horizon:        Horizon of time considered for the analysis
    :param speed_average:       Average speed of the meteorite impact
    """
    area = np.prod(area_box[:, 1] - area_box[:, 0])
    time_period = time_horizon[1] - time_horizon[0]
    rate_per_area_time = integrate.quad(fn_mass_flux, 0.0, np.inf)
    rate = rate_per_area_time[0] * area * time_period * mass_flux_multiplier
    num_events = st.poisson(rate).rvs()

    xtmv = []
    if verbose: print(f"Area experiences `{num_events}` meteorite strikes")
    # using acceptance-rejection sampling to pick `x,t,m`
    for i in range(num_events):
        while True:
            # Propose an `m`
            # Sampling from an exponential
            m = st.expon(scale=1e-3).rvs()
            ratio = fn_mass_flux(m) / rate_per_area_time[0]
            u = np.random.rand()
            if u <= ratio:
                break
        x = area_box[:,0] + (area_box[:,1] - area_box[:,0]) * np.random.rand(2)
        t = time_horizon[0] + time_period * np.random.rand()
        v = sample_velocity(speed_average, size=1)
        xtmv.append(ImpactEvent(t, x, m, v))
        if verbose and ((i > 0) and (i % print_interval) == 0): print(f" sampled **{i + 1}** events")
    sorted_xtmv = sorted(xtmv, key=op.attrgetter('t')) 
    return sorted_xtmv

if __name__ == "__main__":

    # Simulation time
    operation_time = 10 * 365 * 24. # hours
    dt = 0.5 # hours
    max_size = int(operation_time // dt)
    #########################################################################
    ## Probability of finding any meteorite within a given area and a given
    ## time horizon
    #########################################################################
    area_box = np.array([[-200., 200.],
                         [-200., 200.]]) # m
    time_horizon = [0, operation_time]
    speed_average = 20
    num_sims = 100
    impact_sim = 0
    mass = []
    v_s = []
    num_impacts = []
    impacts = []
    for _ in range(num_sims):
        impact_events = sample_xtmv(area_box, time_horizon, speed_average)
        if impact_events:
            n_impacts = len(impact_events)
            num_impacts.append(n_impacts)
            print(f"@sim={_}, num_impacts = {n_impacts}")
            mass.append(sum([impact.m for impact in impact_events]) / n_impacts)
            v_s.append(np.average([np.linalg.norm(impact.v) for impact in impact_events]))
            impact_sim += 1
            impacts.append(impact_events)
    
    print(f"{impact_sim}/{num_sims} have impacts!!")
    print(f"Ave Impacts: {np.average(num_impacts)}, Average mass: {np.average(mass)}, Average velocity: {np.average(v_s)}")
    pprint(impacts[-1])
