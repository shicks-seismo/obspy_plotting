#!/usr/bin/env python
"""
Script to download and plot RaspberryShake station data
Also computes and plots theoretical phase arrival times and raypaths

Stephen Hicks
Imperial College London
Feb 2020
"""

from obspy.clients.fdsn import Client
from obspy import UTCDateTime
from obspy.taup import TauPyModel
from obspy.geodetics.base import locations2degrees
import matplotlib.pyplot as plt

# Start of parameters to define
NETWORK = "AM"   # AM = RaspberryShake network
STATION = "RAEBE"  # Station code of station to get data for
CHANNEL = "EHZ"  # channel to grab data for (e.g. EHZ, SHZ, EHE, EHN)
EQ_TIME = "2020-01-07T08:24:26"  # origin time of earthquake
T_START = 0 # Length in seconds of data to plot before origin time
T_END = 1250  # Length in seconds of data to plot after origin time
PHASES = ["P", "S"]   # list of phases to compute theoretical times for
EVT_LAT = 17.916  # Latitude of event
EVT_LON = -66.813  # Longitude of event 
EVT_Z = 10  # Depth of event
STA_LAT = 51.33  # Latitude of station  
STA_LON = -0.49  # Longitude of station
F1 = 0.3  # High-pass filter corner
F2 = 0.7  # Low-pass filter corner 
LABEL = "M 6.4 Puerto Rico"  # Title to plot on figure
MODEL = 'iasp91'  # Velocity model to predict travel-times through
# End of parameters to define

# Define fdsn client to get data from
client = Client('http://fdsnws.raspberryshakedata.com')

# Define start and end time
orig_time = UTCDateTime(EQ_TIME)
t1 = orig_time - T_START
t2 = orig_time + T_END
# Download and filfter data
st = client.get_waveforms(NETWORK, STATION, "00", CHANNEL,
                          starttime=t1, endtime=t2, attach_response=True)
st.merge()
st.detrend(type="demean")
st.remove_response()
st.filter("bandpass", freqmin=F1, freqmax=F2, corners=4)
st.trim(t1, t2)

# Set-up figure
fig = plt.figure(figsize=(12, 8))
plt.suptitle(LABEL)
ax = plt.subplot(121)

# Set-up taup travel-time model
dist = locations2degrees(EVT_LAT, EVT_LON, STA_LAT, STA_LON)
model = TauPyModel(model=MODEL)

# Now plot the waveform data
ax.plot(st[0].times(reftime=orig_time), st[0].data*1000, linewidth=0.2,
        color="darkred")
ymin, ymax = ax.get_ylim()

# Now plot the theoretical arrival times
for phase in PHASES:
    phase = [phase]
    tt = model.get_travel_times(source_depth_in_km=EVT_Z,
                                distance_in_degree=dist,
                                phase_list=phase)
    ax.vlines(tt[0].time, ymin, ymax, color="blue",
              linewidth=1.2, zorder=3, linestyle="--", alpha=0.5)
    ax.text(tt[0].time*1.02, ymax, phase[0], fontsize=12,
            horizontalalignment="left", verticalalignment="top")
ax.set_xlabel("Time after earthquake (s)")
ax.set_title("{:}.{:}.{:}.{:}\nBandpass filter: {:}-{:} Hz".format(
    st[0].stats.network, st[0].stats.station, st[0].stats.location,
    st[0].stats.channel, F1, F2))
ax.set_ylabel("Ground velocity (mm/s)")

# Now plot the raypaths through the Earth
ax2 = plt.subplot(122, projection='polar')
arrivals = model.get_ray_paths(
    source_depth_in_km=EVT_Z, distance_in_degree=dist,
    phase_list=PHASES)
ax3 = arrivals.plot_rays(phase, legend=False, ax=ax2, show=False,
                         label_arrivals=True)
ax3.set_title("Epicentral distance: {:3.1f}$^\circ$".format(dist))

# Save and plot the figure
plt.tight_layout(rect=[0, 0.03, 1, 0.95])
plt.savefig("traces.png")
plt.show()
