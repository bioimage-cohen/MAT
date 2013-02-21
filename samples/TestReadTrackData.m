% [objTracks gConnect] = SampleReadTrackData(DatasetDir, DatasetName) - Sample LEVer track data writer.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%    Copyright 2013 Andrew Cohen, Eric Wait and Mark Winter
% 
%    This file is part of LEVer - the tool for stem cell lineaging. See
%    http://bioimage.coe.drexel.edu for details
% 
%    LEVer is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
% 
%    LEVer is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
% 
%    You should have received a copy of the GNU General Public License
%    along with LEVer in file "gnu gpl v3.txt".  If not, see 
%    <http://www.gnu.org/licenses/>.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [objTracks gConnect] = SampleReadTrackData(DatasetDir, DatasetName)
global CellHulls

bValidCell = ~[CellHulls.deleted];
backMap = find(bValidCell);

objTracks = struct('Label',cell(1,length(CellHulls)), 'ccc',cell(1,length(CellHulls)),...
                   'outID',{0}, 'inID',{0});

fname = fullfile(DatasetDir,['Tracked_' DatasetName '.txt']);
fid=fopen(fname,'rt');
bDone=0;
TrackList=[];

while ~bDone
    dd=fscanf(fid,'%d,%d,%d\n',3);
    if -1==dd(1)
       bDone=1;
       break
    end
    TrackList=[TrackList;dd'];
end

bDone=0;
InList=[];
dd=textscan(fid,'%f,%f,%f');
InList=[dd{1},dd{2},dd{3}];

fclose(fid);

for i=1:size(TrackList,1)
    o1 = backMap(TrackList(i,2));
    o2 = backMap(TrackList(i,3));
    
    objTracks(o1).Label=TrackList(i,1);
    objTracks(o2).Label=TrackList(i,1);
    objTracks(o1).outID=o2;
    objTracks(o2).inID=o1;
end

nLabel=max([objTracks.Label])+1;
for n=1:length(objTracks)
    if (objTracks(n).Label>0)
        continue;
    end
    
    objTracks(n).Label = nLabel;
    nLabel=nLabel+1;
end

gConnect=sparse(InList(:,2),InList(:,1),InList(:,3),length(CellHulls),length(CellHulls));
end