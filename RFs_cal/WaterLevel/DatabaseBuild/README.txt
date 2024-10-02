This file shows how to request EQ data with JWEED, and compile the waveform in Antelope software.

Data request with JWEED:
(1) request EQ events (MB>=5.4, MHDF and WHDF from NEIC) and save the event file;
(2) request station list;
(3) generate request forms and send requests to IRIS (send_request.m)

Build up Antelope tables:
(1) rdseed2sac / seed2db / sd2db: rdseed2db_RFs.csh, seed2db_RFs.csh
    (NOTE: for older data, need sd2db)
(2)
% change channels 1C* to BH* in *.wfdisc
dbset WALLOWA.wfdisc chan BHZ_01 BHZ
dbset WALLOWA.wfdisc chan 1C2 BHN
dbset WALLOWA.wfdisc chan 1C3 BHE

% change station names from seriel numbers to real names
% dbset WALLOWA.wfdisc sta serial_nuumber sta_name

% fix channel id:
dbfixchanids dbname
(3)
% download event informataion (EHDF files) from website (ehdf USGS)
% cat all the ehdf files into one big file (cat command)
% cat HDF_weekly_9908 HDF_weekly_9909  HDF_weekly_9910 HDF_weekly_9911  HDF_weekly_9912 >> all_HDF_weekly

% create origin table using ehdf filef
pde2origin *.dat WALLOWA

(4) ae

pick up arrivals: dbpick Wallowa
sa on
pal on
ph P (P phase)
rec
cm 500 (change maximum # of traces to display)
sc *:BHZ

se (the orid to any specific event)
fe (first event)
ne (next event)
pe (previous event)
le (last event)

ae
wa

