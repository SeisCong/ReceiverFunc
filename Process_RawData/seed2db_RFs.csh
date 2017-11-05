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
if ( ! ( -e $dbdir ) ) then
	mkdir $dbdir
endif
cd $dbdir
if (!(-e $wkdir/log_mkdb.txt)) then
	echo  "Read seed files" > $wkdir/log_mkdb.txt
else
	foreach seedfile ( `ls -f $seeddir/*.seed`)
 	  echo $seedfile 
 	  set nfile = `echo | grep -c $seedfile $wkdir/log_mkdb.txt`
 	  echo "find" $nfile "seed files"
 	  if ($nfile != 0) then
 	    echo "Seed file has been added into the database"
      else
      	seed2db -v $seedfile $dbName
        echo  $seedfile >> $wkdir/log_mkdb.txt
	  endif 
 	end # of each seed file
endif



