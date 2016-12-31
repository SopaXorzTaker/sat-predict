import math
import time
from datetime import datetime, timedelta

HOUR = 60 * 60
DAY = 24 * HOUR
EARTH_RADIUS = 6378137.0
EARTH_ECCENTRICITY = 8.1819190842622e-2
EARTH_MU = 3.986004418e14
EARTH_ROTATION = 7.2921158553e-5
EARTH_J2 = 1.0826269e-3
DEG_TO_RAD = math.pi / 180
RAD_TO_DEG = 180 / math.pi
MIN_PRECISION = 1e-3


def tle_to_params(line1, line2):
    epoch_year, epoch_day, ballistic_coefficient = int(line1[18:20]), float(line1[20:32]), float(line1[33:43])
    inclination, ascension, eccentricity, argument_of_perigee, mean_anomaly, mean_motion =\
        float(line2[8:16]), float(line2[17:25]), float("0." + line2[26:33]), float(line2[34:42]), float(line2[43:51]),\
        float(line2[52:63])
    return epoch_year, epoch_day, inclination, ascension, eccentricity, argument_of_perigee,\
        mean_anomaly, mean_motion, ballistic_coefficient


def julian_date(year, day):
    a = (year - 1) // 100
    b = 2 - a + a//4

    julian = int((year - 1) * 365.25) + int(30.6001 * 14) + 1720994.5 + b + day
    return julian


def jd_theta(jd):
    ut = (jd + 0.5) % 1
    tu = (jd - (ut + 2451545.0)) / 36525
    gmst = 24110.54841 + tu * (8640184.812866 + tu * (0.093104 - tu * 6.2e-6))
    gmst = (gmst + DAY * 1.00273790934 * ut) % DAY
    theta = gmst / DAY * math.pi * 2
    return theta


def ecef_to_coords(x, y, z):
    asq, esq = EARTH_RADIUS**2, EARTH_ECCENTRICITY**2
    b = math.sqrt(asq * (1-esq))
    bsq = b**2
    ep = math.sqrt((asq-bsq)/bsq)
    p = math.sqrt(x**2 + y**2)
    th = math.atan2(EARTH_RADIUS*z, b*p)
    lon = math.atan2(y, x)
    lat = math.atan2(z + (ep**2)*b*(math.sin(th)**3), p - esq*EARTH_RADIUS*(math.cos(th)**3))
    # n = EARTH_RADIUS/math.sqrt(1-esq*(math.sin(lat)**2))
    # alt = p / math.cos(lat) - n

    lat *= RAD_TO_DEG
    lon *= RAD_TO_DEG

    return lat, lon


def find_eccentric_anomaly(mean_anomaly, eccentricity):
    eccentric_anomaly = mean_anomaly

    while abs(mean_anomaly - (eccentric_anomaly - eccentricity * math.sin(eccentric_anomaly))) > MIN_PRECISION:
        eccentric_anomaly -= (eccentric_anomaly - (mean_anomaly + eccentricity * math.sin(eccentric_anomaly))) / \
                             (1 - eccentricity * math.cos(eccentric_anomaly))

    return eccentric_anomaly


class Satellite(object):
    def __init__(self, params):
        self.epoch_year, self.epoch_day, self.inclination, self.ascension, self.eccentricity,\
            self.argument_of_perigee, self.mean_anomaly, self.mean_motion, self.ballistic_coefficient = params

        self.epoch_timestamp = (datetime(2000 + self.epoch_year, 1, 1) +
                                timedelta(self.epoch_day - 1) - datetime(1970, 1, 1)).total_seconds()

        self.epoch_julian = julian_date(2000 + self.epoch_year, self.epoch_day - 1.0)

    def predict(self, timestamp, last_timestamp=0.0, last_mean_motion=0.0, last_mean_anomaly=0.0, last_ascension=0.0, last_argument=0.0):
        if last_timestamp == 0.0:
            last_mean_motion, last_mean_anomaly, last_ascension, last_argument = self.mean_motion, self.mean_anomaly, self.ascension, self.argument_of_perigee

        delta_t = timestamp - last_timestamp
        current_mean_motion = last_mean_motion + (delta_t / DAY * self.ballistic_coefficient * 2)
        orbital_period = DAY / current_mean_motion

        semi_major_axis = (((orbital_period / (math.pi * 2)) ** 2) * EARTH_MU) ** (1 / 3)

        current_mean_anomaly = ((last_mean_anomaly * DEG_TO_RAD +
                                 delta_t * math.sqrt(EARTH_MU / semi_major_axis ** 3)) * RAD_TO_DEG) % 360

        eccentric_anomaly = find_eccentric_anomaly(current_mean_anomaly * DEG_TO_RAD, self.eccentricity) * RAD_TO_DEG

        true_anomaly = 2 * math.atan2(math.sqrt(1 + self.eccentricity) * math.sin(eccentric_anomaly * DEG_TO_RAD / 2),
                                      math.sqrt(1 - self.eccentricity) * math.cos(
                                          eccentric_anomaly * DEG_TO_RAD / 2)) * RAD_TO_DEG

        radius = (semi_major_axis * (1 - self.eccentricity**2)) /\
                 (1 + self.eccentricity * math.cos(true_anomaly * DEG_TO_RAD))

        # from https://smallsats.org/2013/01/20/j2-propagator/
        ascension_change = -(1.5*EARTH_MU**0.5*EARTH_J2*EARTH_RADIUS**2/((1-self.eccentricity**2)*semi_major_axis**3.5))*math.cos(self.inclination*DEG_TO_RAD) * RAD_TO_DEG
        argument_change = ascension_change*(2.5*math.sin(self.inclination*DEG_TO_RAD)**2 - 2) / math.cos(self.inclination*DEG_TO_RAD)

        current_ascension = last_ascension + ascension_change * delta_t
        current_argument = last_argument + argument_change * delta_t

        cos_arg_tru, sin_arg_tru = math.cos(current_argument * DEG_TO_RAD + true_anomaly * DEG_TO_RAD), \
            math.sin(current_argument * DEG_TO_RAD + true_anomaly * DEG_TO_RAD)

        cos_inc, sin_inc = math.cos(self.inclination * DEG_TO_RAD), math.sin(self.inclination * DEG_TO_RAD)
        cos_asc, sin_asc = math.cos(current_ascension * DEG_TO_RAD), math.sin(current_ascension * DEG_TO_RAD)

        # from http://ccar.colorado.edu/asen5070/handouts/kep2cart_2002.doc
        position = (
            radius * (cos_asc * cos_arg_tru - sin_asc * sin_arg_tru * cos_inc),
            radius * (sin_asc * cos_arg_tru + cos_asc * sin_arg_tru * cos_inc),
            radius * (sin_inc * sin_arg_tru),
        )

        theta = jd_theta(self.epoch_julian + timestamp/DAY)

        cos_theta, sin_theta = math.cos(theta), math.sin(theta)
        x = position[0]*cos_theta + position[1]*sin_theta
        y = position[0]*-sin_theta + position[1]*cos_theta
        z = position[2]

        """
                ap_plus_pe = semi_major_axis * 2
                ap_minus_pe = self.eccentricity * ap_plus_pe
                apogee = (ap_plus_pe + ap_minus_pe) / 2
                perigee = apogee - ap_minus_pe

                apogee -= EARTH_RADIUS
                perigee -= EARTH_RADIUS
        """

        latitude, longitude = ecef_to_coords(x, y, z)

        return current_mean_motion, current_mean_anomaly, current_ascension, current_argument, latitude, longitude


#sat_params = tle_to_params("1 27460U 02035A   16353.76973984 -.00000051  00000-0  00000+0 0  9996",
#                           "2 27460   0.0395 354.2710 0006616 279.5515  86.0990  1.00272668 52981")

sat_params = tle_to_params("1 25544U 98067A   16358.18798762  .00001342  00000-0  27737-4 0  9991",
                           "2 25544  51.6428 195.6612 0006387   5.3135 104.2861 15.53921298 34384")

#sat_params = tle_to_params("1 07276U 74026A   16361.60125193  .00000194  00000-0  56065-3 0  9994",
#                           "2 07276  62.7229 202.5803 6911586 286.6968  12.4802  2.45094949200390")

#sat_params = tle_to_params("1 37158U 10045A   16359.31422074 -.00000086  00000-0  00000+0 0  9995",
#                           "2 37158  40.7711 161.8701 0746301 269.9385 273.1385  1.00272864 23007")

#sat_params = tle_to_params("1 33591U 09005A   16364.50642449  .00000064  00000-0  59538-4 0  9993",
#                           "2 33591  99.0688 323.1466 0014603  32.6768 327.5304 14.12139775406705")

sat = Satellite(sat_params)

print(sat_params)
print(sat.predict(0))

lats, lons = [], []

m, ma, asc, arg = 0.0, 0.0, 0.0, 0.0

for i in range(DAY):
    m, ma, asc, arg, lat, lon = sat.predict(i+1, i, m, ma, asc, arg)
    lats.append(lat)
    lons.append(lon)


import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

map = Basemap(projection='cyl', lon_0=0, resolution='l')
map.drawcoastlines(linewidth=0.25)
map.drawcountries(linewidth=0.25)
map.drawmapboundary(fill_color='aqua')
map.fillcontinents(color='coral', lake_color='aqua', alpha=0.5)
points = map(lons, lats)
map.scatter(*points)
plt.show()
