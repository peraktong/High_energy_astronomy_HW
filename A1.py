from astropy.table import Table
from astropy.io import fits
from astropy import units as u
from astropy.coordinates import SkyCoord
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import math

from astropy.time import Time
from astropy import units as u
from astropy.coordinates import SkyCoord
from math import sin, asin, cos, acos, pi

def celestial_to_galactic(ra,dec):
    c_icrs = SkyCoord(ra=ra * u.degree, dec=dec * u.degree, frame='icrs')
    c_galactic = c_icrs.galactic

    # [l,b]

    return c_galactic.to_string('decimal')[0],c_galactic.to_string('decimal')[1]

celestial_to_galactic = np.vectorize(celestial_to_galactic)

def ra_dec_to_zenith(time,ra,dec):


    lon = -105.820417
    lat = 32.780361

    c = SkyCoord(ra=ra*u.degree, dec=dec*u.degree, frame='icrs')

    # in degree
    ra = c.ra.hour/24*360

    ha = time.sidereal_time("mean").hour/24*360 - ra;
    print(ra,time.sidereal_time("mean").hour/24*360,ha)

    if (ha < 0):
        ha = ha + 360
        #    print "ha:", ha
        # convert degrees to radians
    ha = ha * pi / 180
    dec = dec * pi / 180
    lat = lat * pi / 180

    # compute altitude in radians
    sin_alt = sin(dec) * sin(lat) + cos(dec) * cos(lat) * cos(ha);
    alt = asin(sin_alt)

    # compute azimuth in radians
    # divide by zero error at poles or if alt = 90 deg
    cos_az = (sin(dec) - sin(alt) * sin(lat)) / (cos(alt) * cos(lat))
    az = acos(cos_az)

    # convert radians to degrees
    hrz_altitude = alt * 180 / pi;
    hrz_azimuth = az * 180 / pi;

    zenithAngle = 90.0 - alt

    return zenithAngle


class plot():

    def plot_skymap_P1(self):

        table = Table.read("vanvelzen12-v1.0.fits")

        l_array = []
        b_array = []
        D_array = []
        L_array = []

        for i in range(0, len(table)):
            ra_i = table[i]["ra"]
            dec_i = table[i]["dec"]

            D_array.append(table[i]["D"])
            L_array.append(table[i]["Lsyn"])

            l_i, b_i = celestial_to_galactic(ra_i, dec_i)

            l_array.append(l_i)
            b_array.append(b_i)

            print("Doing %d of %d" % (i, len(table)))
            # print(table[i]["Lsyn"])

        L_array = np.array(L_array)
        l_array = np.array(l_array)
        b_array = np.array(b_array)
        D_array = np.array(D_array)

        mask = (D_array < 100) & (L_array > 2.6 * 10 ** 40)

        # print(mask)


        # plot:


        font = {'weight': 'bold', 'size': 20}
        matplotlib.rc('font', **font)

        f, ax = plt.subplots()


        #ax6


        pl = ax.scatter(l_array[mask], b_array[mask], marker='o', c=L_array[mask],
                    vmin=np.min(L_array[mask]), vmax=np.max(L_array[mask]), alpha=1)


        ax.set_xlabel('l', fontsize=20)
        ax.set_ylabel('b', fontsize=20)

        # cbar_ax = f.add_axes([1.05, 0.15, 0.02, 0.7])
        cb = f.colorbar(pl)

        cb.set_label("vLv(1.1Ghz) $erg/s$", fontsize=20)
        f.suptitle('Sky map for selected data from vanvelzen12', fontsize=30)

        # save it

        fig = matplotlib.pyplot.gcf()

        fig.set_size_inches(17.5, 14.5)

        save_path = "Sky_map" + ".png"
        fig.savefig(save_path, dpi=500)

        plt.close()

    def A2_zenith_angle(self):
        t_jd = Time(2455197.5, format='jd', scale='utc', location=(-105.820417, 32.780361))


model = plot()

model.plot_skymap_P1()





















