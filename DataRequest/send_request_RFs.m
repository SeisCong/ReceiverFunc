%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Send data request to BREQ_FAST
%
%
%%  Cong Li     04/15/16
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clear all;
Porjdir=['/Users/congli/Project/ReceiverFunc_Appalachians/'];
Workdir=['DataRequest/EQs_request/'];    
filelist = dir([Porjdir Workdir 'Events*']);
 for ii=1:length(filelist)
fname=[Porjdir Workdir char(filelist(ii).name)]

%   fname=['./EQs_request/Events_QM15_199']

unix( [' mail breq_fast@iris.washington.edu < ' fname ]); 

display(['request ' fname ' is sent']);

pause(2);   % sleep the request for half an hour

 end


