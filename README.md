#GKLatLonTransformer
Transformation of latitude and longitude to plane Gauss-Kruger and vise versa.
```
usage: gaussKrugerTransform.py [-h] -t {GeodeticToPlane,PlaneToGeodetic}
                               [-e {wgs84,CGCS2000}]
                               latitude longitude

Tool to transform geodedic to plane and reverse in Gauss-Kruger projection

positional arguments:
  latitude              Latitude(in GG MM SS.S format) or X(in meters) of
                        point to transform
  longitude             Longitude(in GG MM SS.S format) or Y(in meters) of
                        point to convert in GG MM SS.S format

optional arguments:
  -h, --help            show this help message and exit
  -t {GeodeticToPlane,PlaneToGeodetic}, --transformationType {GeodeticToPlane,PlaneToGeodetic}
                        Type of transformation
  -e {wgs84,CGCS2000}, --ellipsoid {wgs84,CGCS2000}
                        Ellipsoid
```

**Examples**

1. To transform lat lon to XY:
```
C:\Python35\python.exe gaussKrugerTransform.py -t GeodeticToPlane "22 15 58.98294" "111 28 52.15387" -e wgs84

Geodetic to plane transformation:
x : 2463376.6501757316
y : 19549592.072119053
```
2. To transform XY to lat lon
```
C:\Python35\python.exe gaussKrugerTransform.py -t PlaneToGeodetic 2562038.2708 19512837.2851 -e CGCS2000

Plane to geodetic transformation:
lat(dms) : 23.092871681917128
lon(dms) : 111.07312980400093
```

23.092871681917128 means 23°09'28''.71681917128