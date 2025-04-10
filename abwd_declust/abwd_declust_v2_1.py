import numpy as np
import pandas as pd
from shapely.prepared import prep
from shapely import MultiPoint, Polygon, geometry
from itertools import compress

def abwd_crjb(catalogue_pd, rupture_area_poly, rjb_cutoff, window_method, fs_time_prop):
    from datetime import datetime as dt
    import time
    def decimalyear(catalogue_pd):
        def toYearFraction(date):
            def sinceEpoch(date):  # returns seconds since epoch
                return time.mktime(date.timetuple())

            s = sinceEpoch

            year = date.year
            startOfThisYear = dt(year=year, month=1, day=1)
            startOfNextYear = dt(year=year + 1, month=1, day=1)

            yearElapsed = s(date) - s(startOfThisYear)
            yearDuration = s(startOfNextYear) - s(startOfThisYear)
            fraction = yearElapsed / yearDuration

            return date.year + fraction

        time_pd = pd.to_datetime(catalogue_pd.datetime, format='ISO8601')
        time_pd.reset_index(drop=True, inplace=True)
        dtime_orig = []

        for i in range(len(time_pd)):
            temp_orig = toYearFraction(time_pd[i])
            dtime_orig.append(temp_orig)
        return dtime_orig

    def haversine_oq(lon1, lat1, lon2, lat2, radians=False, earth_rad=6371.227):
        """
        Allows to calculate geographical distance
        using the haversine formula.

        :param lon1: longitude of the first set of locations
        :type lon1: numpy.ndarray
        :param lat1: latitude of the frist set of locations
        :type lat1: numpy.ndarray
        :param lon2: longitude of the second set of locations
        :type lon2: numpy.float64
        :param lat2: latitude of the second set of locations
        :type lat2: numpy.float64
        :keyword radians: states if locations are given in terms of radians
        :type radians: bool
        :keyword earth_rad: radius of the earth in km
        :type earth_rad: float
        :returns: geographical distance in km
        :rtype: numpy.ndarray
        """
        if not radians:
            cfact = np.pi / 180.
            lon1 = cfact * lon1
            lat1 = cfact * lat1
            lon2 = cfact * lon2
            lat2 = cfact * lat2

        # Number of locations in each set of points
        if not np.shape(lon1):
            nlocs1 = 1
            lon1 = np.array([lon1])
            lat1 = np.array([lat1])
        else:
            nlocs1 = np.max(np.shape(lon1))
        if not np.shape(lon2):
            nlocs2 = 1
            lon2 = np.array([lon2])
            lat2 = np.array([lat2])
        else:
            nlocs2 = np.max(np.shape(lon2))
        # Pre-allocate array
        distance = np.zeros((nlocs1, nlocs2))
        i = 0
        while i < nlocs2:
            # Perform distance calculation
            dlat = lat1 - lat2[i]
            dlon = lon1 - lon2[i]
            aval = (np.sin(dlat / 2.) ** 2.) + (np.cos(lat1) * np.cos(lat2[i]) *
                                                (np.sin(dlon / 2.) ** 2.))
            distance[:, i] = (2. * earth_rad * np.arctan2(np.sqrt(aval),
                                                          np.sqrt(1 - aval))).T
            i += 1
        return distance

    def polyrupcentroid(rupture_area_poly):
        rupture_area_centroid = []
        for i in range(len(rupture_area_poly)):
            rupture_area_centroid_temp = [rupture_area_poly[i].centroid.x, rupture_area_poly[i].centroid.y]
            rupture_area_centroid.append(rupture_area_centroid_temp)
        return rupture_area_centroid

    def polyres1km(rupture_area_poly):
        rupture_poly_multipoint = []

        for i in range(len(rupture_area_poly)):
            rupture_area_poly_ith = rupture_area_poly[i]
            x, y = rupture_area_poly_ith.exterior.xy
            xy = np.transpose(np.array([np.array(x), np.array(y)]))
            n = int(rupture_area_poly_ith.length * 111.2) + 8
            d = np.cumsum(np.r_[0, np.sqrt((np.diff(xy, axis=0) ** 2).sum(axis=1))])
            d_sampled = np.linspace(0, d.max(), n)
            xy_interp = np.c_[np.interp(d_sampled, d, xy[:, 0]), np.interp(d_sampled, d, xy[:, 1]),]
            xy_resampled_temp = MultiPoint(list(zip(xy_interp[:, 0], xy_interp[:, 1])))

            rupture_poly_multipoint.append(xy_resampled_temp)
        return rupture_poly_multipoint

    def crjb(rupture_area_poly,rup1km_multipoint, rupture_area_centroid):
        crjb_final = []
        for i in range(len(rupture_area_centroid)):
            rup1km_centroid = rupture_area_centroid[i]
            crjb_distance_initial = []
            for j in range(len(rup1km_multipoint.geoms) - 1):
                crjb0_filter = rupture_area_poly.contains(geometry.Point(rup1km_centroid))
                if crjb0_filter:
                    distance_temp = [[0.0]]
                else:
                    distance_temp = haversine_oq(rup1km_centroid[0], rup1km_centroid[1], rup1km_multipoint.geoms[j].x,
                                             rup1km_multipoint.geoms[j].y)
                crjb_distance_initial.append(distance_temp[0][0])
            crjb_final.append(min(crjb_distance_initial))
        return np.array(crjb_final)

    decimal_year = decimalyear(catalogue_pd)
    decimal_yearnp = np.array(decimal_year)
    rupture_area_centroid = polyrupcentroid(rupture_area_poly)
    rup1km_multipoint = polyres1km(rupture_area_poly)

    neq = len(catalogue_pd.mag)

    DAYS = 364.75
    if window_method == 'GardnerKnopoff':
        sw_time = np.power(10.0, 0.032 * catalogue_pd.mag + 2.7389) / DAYS
        sw_time[catalogue_pd.mag < 6.5] = np.power(10.0, 0.5409 * catalogue_pd.mag[catalogue_pd.mag < 6.5] - 0.547) / DAYS
        sw_space = np.power(10.0, 0.1238 * catalogue_pd.mag + 0.983)
        sw_space[sw_space>=0] = rjb_cutoff
    elif window_method == 'Gruenthal':
        sw_time = np.abs((np.exp(-3.95 + np.sqrt(0.62 + 17.32 * catalogue_pd.mag))) / 364.75)
        sw_time[catalogue_pd.mag >= 6.5] = np.power(10, 2.8 + 0.024 * catalogue_pd.mag[catalogue_pd.mag >= 6.5]) / 364.75
        sw_space = np.exp(1.77 + np.sqrt(0.037 + 1.02 * catalogue_pd.mag))
        sw_space[sw_space>=0] = rjb_cutoff
    elif window_method == 'Urhammer':
        sw_time = np.exp(-2.87 + 1.235 * catalogue_pd.mag) / 364.75
        sw_space = np.exp(-1.024 + 0.804 * catalogue_pd.mag)
        sw_space[sw_space>=0] = rjb_cutoff
    else:
        sw_time = np.power(10.0, 0.032 * catalogue_pd.mag + 2.7389) / DAYS
        sw_time[catalogue_pd.mag < 6.5] = np.power(10.0, 0.5409 * catalogue_pd.mag[catalogue_pd.mag < 6.5] - 0.547) / DAYS
        sw_space = np.power(10.0, 0.1238 * catalogue_pd.mag + 0.983)
        sw_space[sw_space>=0] = rjb_cutoff

    eqid = np.arange(0, neq, 1)
    vcl = np.zeros(neq, dtype=int)
    id0 = np.flipud(np.argsort(catalogue_pd.mag, kind="heapsort"))

    sw_time = sw_time[id0]
    sw_space = sw_space[id0]
    year_dec = decimal_yearnp[id0]
    eqid = eqid[id0]

    rupture_area_poly_sort = [rupture_area_poly[i] for i in id0]
    rupture_area_centroid_sort = [rupture_area_centroid[i] for i in id0]
    rup1km_multipoint_sort = [rup1km_multipoint[i] for i in id0]

    flagvector = np.zeros(neq, dtype=int)

    clust_index = 0

    for i in range(0, neq - 1):
        if vcl[i] == 0:
            # Find Events inside both fore- and aftershock time windows
            dt = year_dec - year_dec[i]
            vsel = np.logical_and(vcl == 0, np.logical_and(dt >= (-sw_time[i] * fs_time_prop), dt <= sw_time[i]))
            aftershock_centroids = list(compress(rupture_area_centroid_sort, vsel))
            crjb_dist = crjb(rupture_area_poly_sort[i],rup1km_multipoint_sort[i], aftershock_centroids)
            vsel1 = crjb_dist <= sw_space[i]
            vsel[vsel] = vsel1
            temp_vsel = np.copy(vsel)
            temp_vsel[i] = False

            if any(temp_vsel):
                # Allocate a cluster number
                vcl[vsel] = clust_index + 1
                flagvector[vsel] = 1
                # For those events in the cluster before the main event,
                # flagvector is equal to -1
                temp_vsel[dt >= 0.0] = False
                flagvector[temp_vsel] = -1
                flagvector[i] = 0
                clust_index += 1

    id1 = np.argsort(eqid, kind="heapsort")
    eqid = eqid[id1]
    vcl = vcl[id1]
    flagvector = flagvector[id1]

    return flagvector, vcl

def abwd_rclosestp2h(catalogue_pd, rupture_area_poly, rjb_cutoff, window_method, fs_time_prop):
    from datetime import datetime as dt
    import time
    def decimalyear(catalogue_pd):
        def toYearFraction(date):
            def sinceEpoch(date):  # returns seconds since epoch
                return time.mktime(date.timetuple())

            s = sinceEpoch

            year = date.year
            startOfThisYear = dt(year=year, month=1, day=1)
            startOfNextYear = dt(year=year + 1, month=1, day=1)

            yearElapsed = s(date) - s(startOfThisYear)
            yearDuration = s(startOfNextYear) - s(startOfThisYear)
            fraction = yearElapsed / yearDuration

            return date.year + fraction

        time_pd = pd.to_datetime(catalogue_pd.datetime, format='ISO8601')
        time_pd.reset_index(drop=True, inplace=True)
        dtime_orig = []

        for i in range(len(time_pd)):
            temp_orig = toYearFraction(time_pd[i])
            dtime_orig.append(temp_orig)
        return dtime_orig

    def haversine_oq(lon1, lat1, lon2, lat2, radians=False, earth_rad=6371.227):
        """
        Allows to calculate geographical distance
        using the haversine formula.

        :param lon1: longitude of the first set of locations
        :type lon1: numpy.ndarray
        :param lat1: latitude of the frist set of locations
        :type lat1: numpy.ndarray
        :param lon2: longitude of the second set of locations
        :type lon2: numpy.float64
        :param lat2: latitude of the second set of locations
        :type lat2: numpy.float64
        :keyword radians: states if locations are given in terms of radians
        :type radians: bool
        :keyword earth_rad: radius of the earth in km
        :type earth_rad: float
        :returns: geographical distance in km
        :rtype: numpy.ndarray
        """
        if not radians:
            cfact = np.pi / 180.
            lon1 = cfact * lon1
            lat1 = cfact * lat1
            lon2 = cfact * lon2
            lat2 = cfact * lat2

        # Number of locations in each set of points
        if not np.shape(lon1):
            nlocs1 = 1
            lon1 = np.array([lon1])
            lat1 = np.array([lat1])
        else:
            nlocs1 = np.max(np.shape(lon1))
        if not np.shape(lon2):
            nlocs2 = 1
            lon2 = np.array([lon2])
            lat2 = np.array([lat2])
        else:
            nlocs2 = np.max(np.shape(lon2))
        # Pre-allocate array
        distance = np.zeros((nlocs1, nlocs2))
        i = 0
        while i < nlocs2:
            # Perform distance calculation
            dlat = lat1 - lat2[i]
            dlon = lon1 - lon2[i]
            aval = (np.sin(dlat / 2.) ** 2.) + (np.cos(lat1) * np.cos(lat2[i]) *
                                                (np.sin(dlon / 2.) ** 2.))
            distance[:, i] = (2. * earth_rad * np.arctan2(np.sqrt(aval),
                                                          np.sqrt(1 - aval))).T
            i += 1
        return distance

    def polysubrupcentroid(rupture_area_poly, delta):
        sqkm_deg = (111.1949 ** 2)

        def grid_bounds(geom, delta):
            minx, miny, maxx, maxy = geom.bounds
            nx = int((maxx - minx) / delta)
            ny = int(abs((abs(maxy) - abs(miny)) / delta))
            if ny <= 2:
                ny = 2
            else:
                ny = int(abs((abs(maxy) - abs(miny)) / delta))
            if nx <= 2:
                nx = 2
            else:
                nx = int((maxx - minx) / delta)
            gx, gy = np.linspace(minx, maxx, nx), np.linspace(miny, maxy, ny)
            grid = []
            for i in range(len(gx) - 1):
                for j in range(len(gy) - 1):
                    poly_ij = Polygon([[gx[i], gy[j]], [gx[i], gy[j + 1]], [gx[i + 1], gy[j + 1]], [gx[i + 1], gy[j]]])
                    grid.append(poly_ij)
            return grid

        def partition(geom, delta):
            grid_centroid = []
            prepared_geom = prep(geom)
            grid = list(filter(prepared_geom.intersects, grid_bounds(geom, delta)))
            for i in range(len(grid)):
                grid_centroid_temp = [grid[i].centroid.x, grid[i].centroid.y]
                grid_centroid_filter = geom.contains(geometry.Point(grid_centroid_temp))
                if grid_centroid_filter:
                    grid_centroid.append(grid_centroid_temp)
                else:
                    continue
            return grid, grid_centroid

        # Define Empty Lists
        grid = []
        grid_centroid = []
        grid_centroid_length = []
        grid_centroid_multi = []

        for i in range(len(rupture_area_poly)):
            if rupture_area_poly[i].area * sqkm_deg <= 3:
                grid_temp = rupture_area_poly[i]
                grid_centroid_temp = [[grid_temp.centroid.x, grid_temp.centroid.y]]
            else:
                grid_temp, grid_centroid_temp = partition(rupture_area_poly[i], delta)

            grid.append(grid_temp)
            grid_centroid.append(grid_centroid_temp)
            grid_centroid_length.append(len(grid_centroid_temp))

        # Convert Centroid List to Shapely Multipoints for ease of reading
        for i in range(len(grid_centroid)):
            temp_multi = MultiPoint(points=list(grid_centroid[i]))
            grid_centroid_multi.append(temp_multi)
        return grid_centroid_multi

    def rclosest_p2h(rupture_area_poly, grid_centroid_multi, epicentre_vsel):
        rclosest_p2h_final = []
        for i in range(len(epicentre_vsel)):
            epicentre_ref = epicentre_vsel[i]
            rclosest_p2h_distance_initial = []
            for j in range(len(grid_centroid_multi.geoms)):
                rjb0_filter = rupture_area_poly.contains(geometry.Point(epicentre_ref))
                if rjb0_filter:
                    distance_temp = [[0.0]]
                else:
                    distance_temp = haversine_oq(epicentre_ref[0], epicentre_ref[1], grid_centroid_multi.geoms[j].x,
                                             grid_centroid_multi.geoms[j].y)
                rclosest_p2h_distance_initial.append(distance_temp[0][0])
            rclosest_p2h_final.append(min(rclosest_p2h_distance_initial))
        return np.array(rclosest_p2h_final)

    decimal_year = decimalyear(catalogue_pd)
    decimal_yearnp = np.array(decimal_year)
    grid_centroid_multi = polysubrupcentroid(rupture_area_poly, 0.008)

    neq = len(catalogue_pd.mag)

    DAYS = 364.75
    if window_method == 'GardnerKnopoff':
        sw_time = np.power(10.0, 0.032 * catalogue_pd.mag + 2.7389) / DAYS
        sw_time[catalogue_pd.mag < 6.5] = np.power(10.0, 0.5409 * catalogue_pd.mag[catalogue_pd.mag < 6.5] - 0.547) / DAYS
        sw_space = np.power(10.0, 0.1238 * catalogue_pd.mag + 0.983)
        sw_space[sw_space>=0] = rjb_cutoff
    elif window_method == 'Gruenthal':
        sw_time = np.abs((np.exp(-3.95 + np.sqrt(0.62 + 17.32 * catalogue_pd.mag))) / 364.75)
        sw_time[catalogue_pd.mag >= 6.5] = np.power(10, 2.8 + 0.024 * catalogue_pd.mag[catalogue_pd.mag >= 6.5]) / 364.75
        sw_space = np.exp(1.77 + np.sqrt(0.037 + 1.02 * catalogue_pd.mag))
        sw_space[sw_space>=0] = rjb_cutoff
    elif window_method == 'Urhammer':
        sw_time = np.exp(-2.87 + 1.235 * catalogue_pd.mag) / 364.75
        sw_space = np.exp(-1.024 + 0.804 * catalogue_pd.mag)
        sw_space[sw_space>=0] = rjb_cutoff
    else:
        sw_time = np.power(10.0, 0.032 * catalogue_pd.mag + 2.7389) / DAYS
        sw_time[catalogue_pd.mag < 6.5] = np.power(10.0, 0.5409 * catalogue_pd.mag[catalogue_pd.mag < 6.5] - 0.547) / DAYS
        sw_space = np.power(10.0, 0.1238 * catalogue_pd.mag + 0.983)
        sw_space[sw_space>=0] = rjb_cutoff

    eqid = np.arange(0, neq, 1)
    vcl = np.zeros(neq, dtype=int)
    id0 = np.flipud(np.argsort(catalogue_pd.mag, kind="heapsort"))

    epicentre = np.transpose(np.array([catalogue_pd.hyp_lon, catalogue_pd.hyp_lat]))
    epicentre_sort = epicentre[id0]
    rupture_area_poly_sort = [rupture_area_poly[i] for i in id0]
    grid_centroid_multi_sort = [grid_centroid_multi[i] for i in id0]
    sw_time = sw_time[id0]
    sw_space = sw_space[id0]
    year_dec = decimal_yearnp[id0]
    eqid = eqid[id0]
    flagvector = np.zeros(neq, dtype=int)

    clust_index = 0

    for i in range(0, neq - 1):
        if vcl[i] == 0:
            # Find Events inside both fore- and aftershock time windows
            dt = year_dec - year_dec[i]
            vsel = np.logical_and(vcl == 0, np.logical_and(dt >= (-sw_time[i] * fs_time_prop), dt <= sw_time[i]))
            epicentre_vsel = (epicentre_sort[vsel])
            rclosest_p2h_dist = rclosest_p2h(rupture_area_poly_sort[i],grid_centroid_multi_sort[i], epicentre_vsel)
            vsel1 = rclosest_p2h_dist <= sw_space[i]
            # print(vsel1.shape)
            vsel[vsel] = vsel1
            temp_vsel = np.copy(vsel)
            temp_vsel[i] = False

            if any(temp_vsel):
                # Allocate a cluster number
                vcl[vsel] = clust_index + 1
                flagvector[vsel] = 1
                # For those events in the cluster before the main event,
                # flagvector is equal to -1
                temp_vsel[dt >= 0.0] = False
                flagvector[temp_vsel] = -1
                flagvector[i] = 0
                clust_index += 1

    id1 = np.argsort(eqid, kind="heapsort")
    eqid = eqid[id1]
    vcl = vcl[id1]
    flagvector = flagvector[id1]

    return flagvector, vcl

def abwd_rclosestp2p(catalogue_pd, rupture_area_poly, rjb_cutoff, window_method, fs_time_prop):
    from datetime import datetime as dt
    import time
    def decimalyear(catalogue_pd):
        def toYearFraction(date):
            def sinceEpoch(date):  # returns seconds since epoch
                return time.mktime(date.timetuple())

            s = sinceEpoch

            year = date.year
            startOfThisYear = dt(year=year, month=1, day=1)
            startOfNextYear = dt(year=year + 1, month=1, day=1)

            yearElapsed = s(date) - s(startOfThisYear)
            yearDuration = s(startOfNextYear) - s(startOfThisYear)
            fraction = yearElapsed / yearDuration

            return date.year + fraction

        time_pd = pd.to_datetime(catalogue_pd.datetime, format='ISO8601')
        time_pd.reset_index(drop=True, inplace=True)
        dtime_orig = []

        for i in range(len(time_pd)):
            temp_orig = toYearFraction(time_pd[i])
            dtime_orig.append(temp_orig)
        return dtime_orig

    def haversine_oq(lon1, lat1, lon2, lat2, radians=False, earth_rad=6371.227):
        """
        Allows to calculate geographical distance
        using the haversine formula.

        :param lon1: longitude of the first set of locations
        :type lon1: numpy.ndarray
        :param lat1: latitude of the frist set of locations
        :type lat1: numpy.ndarray
        :param lon2: longitude of the second set of locations
        :type lon2: numpy.float64
        :param lat2: latitude of the second set of locations
        :type lat2: numpy.float64
        :keyword radians: states if locations are given in terms of radians
        :type radians: bool
        :keyword earth_rad: radius of the earth in km
        :type earth_rad: float
        :returns: geographical distance in km
        :rtype: numpy.ndarray
        """
        if not radians:
            cfact = np.pi / 180.
            lon1 = cfact * lon1
            lat1 = cfact * lat1
            lon2 = cfact * lon2
            lat2 = cfact * lat2

        # Number of locations in each set of points
        if not np.shape(lon1):
            nlocs1 = 1
            lon1 = np.array([lon1])
            lat1 = np.array([lat1])
        else:
            nlocs1 = np.max(np.shape(lon1))
        if not np.shape(lon2):
            nlocs2 = 1
            lon2 = np.array([lon2])
            lat2 = np.array([lat2])
        else:
            nlocs2 = np.max(np.shape(lon2))
        # Pre-allocate array
        distance = np.zeros((nlocs1, nlocs2))
        i = 0
        while i < nlocs2:
            # Perform distance calculation
            dlat = lat1 - lat2[i]
            dlon = lon1 - lon2[i]
            aval = (np.sin(dlat / 2.) ** 2.) + (np.cos(lat1) * np.cos(lat2[i]) *
                                                (np.sin(dlon / 2.) ** 2.))
            distance[:, i] = (2. * earth_rad * np.arctan2(np.sqrt(aval),
                                                          np.sqrt(1 - aval))).T
            i += 1
        return distance

    def polyres1km(rupture_area_poly):
        rupture_poly_multipoint = []

        for i in range(len(rupture_area_poly)):
            rupture_area_poly_ith = rupture_area_poly[i]
            x, y = rupture_area_poly_ith.exterior.xy
            xy = np.transpose(np.array([np.array(x), np.array(y)]))
            n = int(rupture_area_poly_ith.length * 111.2) + 12
            d = np.cumsum(np.r_[0, np.sqrt((np.diff(xy, axis=0) ** 2).sum(axis=1))])
            d_sampled = np.linspace(0, d.max(), n)
            xy_interp = np.c_[np.interp(d_sampled, d, xy[:, 0]), np.interp(d_sampled, d, xy[:, 1]),]
            xy_resampled_temp = MultiPoint(list(zip(xy_interp[:, 0], xy_interp[:, 1])))

            rupture_poly_multipoint.append(xy_resampled_temp)
        return rupture_poly_multipoint

    def rclosest_p2p(rupture_area_poly, rup1km_multipoint, rup1km_multipoint_ref):
        rclosest_p2p_final = []
        for k in range(len(rup1km_multipoint_ref)):
            rclosest_p2p_initial2 = []
            for i in range(len(rup1km_multipoint_ref[k].geoms)):
                rclosest_p2p_initial = []
                for j in range(len(rup1km_multipoint.geoms)):
                    rjb0_filter = rupture_area_poly.contains(
                        geometry.Point([rup1km_multipoint_ref[k].geoms[i].x, rup1km_multipoint_ref[k].geoms[i].y]))
                    if rjb0_filter:
                        distance_temp = [[0.0]]
                    else:
                        distance_temp = haversine_oq(rup1km_multipoint_ref[k].geoms[i].x,
                                                 rup1km_multipoint_ref[k].geoms[i].y, rup1km_multipoint.geoms[j].x,
                                                 rup1km_multipoint.geoms[j].y)
                    rclosest_p2p_initial.append(distance_temp[0][0])
                rclosest_p2p_initial2.append(min(rclosest_p2p_initial))
            rclosest_p2p_final.append(min(rclosest_p2p_initial2))
        return rclosest_p2p_final

    decimal_year = decimalyear(catalogue_pd)
    decimal_yearnp = np.array(decimal_year)
    rup1km_multipoint = polyres1km(rupture_area_poly)

    neq = len(catalogue_pd.mag)

    DAYS = 364.75
    if window_method == 'GardnerKnopoff':
        sw_time = np.power(10.0, 0.032 * catalogue_pd.mag + 2.7389) / DAYS
        sw_time[catalogue_pd.mag < 6.5] = np.power(10.0, 0.5409 * catalogue_pd.mag[catalogue_pd.mag < 6.5] - 0.547) / DAYS
        sw_space = np.power(10.0, 0.1238 * catalogue_pd.mag + 0.983)
        sw_space[sw_space >= 0] = rjb_cutoff
    elif window_method == 'Gruenthal':
        sw_time = np.abs((np.exp(-3.95 + np.sqrt(0.62 + 17.32 * catalogue_pd.mag))) / 364.75)
        sw_time[catalogue_pd.mag >= 6.5] = np.power(10,
                                                    2.8 + 0.024 * catalogue_pd.mag[catalogue_pd.mag >= 6.5]) / 364.75
        sw_space = np.exp(1.77 + np.sqrt(0.037 + 1.02 * catalogue_pd.mag))
        sw_space[sw_space >= 0] = rjb_cutoff
    elif window_method == 'Urhammer':
        sw_time = np.exp(-2.87 + 1.235 * catalogue_pd.mag) / 364.75
        sw_space = np.exp(-1.024 + 0.804 * catalogue_pd.mag)
        sw_space[sw_space >= 0] = rjb_cutoff
    else:
        sw_time = np.power(10.0, 0.032 * catalogue_pd.mag + 2.7389) / DAYS
        sw_time[catalogue_pd.mag < 6.5] = np.power(10.0, 0.5409 * catalogue_pd.mag[catalogue_pd.mag < 6.5] - 0.547) / DAYS
        sw_space = np.power(10.0, 0.1238 * catalogue_pd.mag + 0.983)
        sw_space[sw_space >= 0] = rjb_cutoff

    eqid = np.arange(0, neq, 1)
    vcl = np.zeros(neq, dtype=int)
    id0 = np.flipud(np.argsort(catalogue_pd.mag, kind="heapsort"))

    rupture_area_poly_sort = [rupture_area_poly[i] for i in id0]
    rup1km_multipoint_sort = [rup1km_multipoint[i] for i in id0]
    sw_time = sw_time[id0]
    sw_space = sw_space[id0]
    year_dec = decimal_yearnp[id0]
    eqid = eqid[id0]
    flagvector = np.zeros(neq, dtype=int)

    clust_index = 0

    for i in range(0, neq - 1):
        if vcl[i] == 0:
            # Find Events inside both fore- and aftershock time windows
            dt = year_dec - year_dec[i]
            vsel = np.logical_and(vcl == 0, np.logical_and(dt >= (-sw_time[i] * fs_time_prop), dt <= sw_time[i]))
            rup1km_multipoint_ref = list(compress(rup1km_multipoint_sort, vsel))
            rclosest_p2p_dist = rclosest_p2p(rupture_area_poly_sort[i], rup1km_multipoint_sort[i], rup1km_multipoint_ref)
            vsel1 = rclosest_p2p_dist <= sw_space[i]
            vsel[vsel] = vsel1
            temp_vsel = np.copy(vsel)
            temp_vsel[i] = False

            if any(temp_vsel):
                # Allocate a cluster number
                vcl[vsel] = clust_index + 1
                flagvector[vsel] = 1
                # For those events in the cluster before the main event,
                # flagvector is equal to -1
                temp_vsel[dt >= 0.0] = False
                flagvector[temp_vsel] = -1
                flagvector[i] = 0
                clust_index += 1

    id1 = np.argsort(eqid, kind="heapsort")
    eqid = eqid[id1]
    vcl = vcl[id1]
    flagvector = flagvector[id1]

    return flagvector, vcl
