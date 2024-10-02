#/bin/csh
# extract EQ waveform seed files into Antelope database tables

set wkdir = /Users/congli/Research/RFs_NewEngland
set network = RFs_1995_now_nopick
set db = $network

set dbdir = $wkdir/$db
if ( ! ( -e $dbdir ) ) mkdir $dbdir

cd $dbdir

foreach seedfile ( `ls -f $wkdir/seed/Events_*.seed`)
echo $seedfile

seed2db -v $seedfile $db

end # of each seed file



