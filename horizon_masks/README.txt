Horizon masks of the ground stations
====================================

One file per station, named exactly as the station in ground_stations.txt:
horizon_masks\<station>.txt (or .csv), e.g. horizon_masks\Svalbard.txt.

Two columns per line: azimuth [deg from North, clockwise] and the minimum
elevation [deg] at that azimuth, separated by a comma, tabs, spaces or a
semicolon (then a decimal comma is accepted). Header lines and # comments
are skipped. The mask is linear between the points and wraps around North.

A satellite is in contact only above both the mask and the Min Elevation of
the Orbit & Passes tab. Example:

# azimuth, min elevation
0, 3.0
90, 6.5
180, 2.0
270, 4.0
