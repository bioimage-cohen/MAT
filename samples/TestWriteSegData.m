% SampleWriteSegData(DatasetDir, DatasetName) - Sample LEVer segmentation data writer.

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

function SampleWriteSegData(DatasetDir, DatasetName)
global CONSTANTS CellHulls

sz = CONSTANTS.imageSize;

th = max([CellHulls.time]);
hashedHulls = cell(th,1);

% reset tracking info
for n=1:length(CellHulls)
    hashedHulls{CellHulls(n).time} = [hashedHulls{CellHulls(n).time}; n];
end

fname = fullfile(DatasetDir, ['SegObjs_Test_' DatasetName '.txt']);
fid=fopen(fname,'wt');

fprintf(fid, '%d %d\n', sz(2), sz(1));
fprintf(fid, '%d %d\n\n', th, length(CellHulls));

for t=1:length(hashedHulls)
    fprintf(fid, '%d\n', length(hashedHulls{t}));
    
    for i=1:length(hashedHulls{t})
        [r c]=ind2sub(sz,CellHulls(hashedHulls{t}(i)).indexPixels);
        
        COM = mean([c r], 1);
        fprintf(fid, '%f %f %d:', COM(1), COM(2), length(r));
        
        for k=1:length(r)
            fprintf(fid, ' (%d,%d)', c(k), r(k));
        end
        
        fprintf(fid,'\n');
    end
end

fclose(fid);
end