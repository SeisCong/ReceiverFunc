#/bin/csh
##################################################################
# Extract EQ waveform seed files into Antelope database tables
#
#  Cong Li
# 
# 11/2/2017
##################################################################

set wkdir = $1 # workdir
set dbName = $2 #databasename
set seeddir = $3 # raw data dir

set dbdir = $wkdir/$dbName
if ( ! ( -e $dbdir ) ) mkdir $dbdir

cd $dbdir
echo  "Read seed files" > $wkdir/log_mkdb.txt
foreach seedfile ( `ls -f $seeddir/Events_*.seed`)
echo $seedfile 
echo  $seedfile >> $wkdir/log_mkdb.txt

seed2db -v $seedfile $dbName

end # of each seed file



