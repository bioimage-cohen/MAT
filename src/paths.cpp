//***********************************************************************
//
//    Copyright 2013 Andrew Cohen and Mark Winter
// 
//    This file is part of the Multitemporal Association Tracker. See
//    http://bioimage.coe.drexel.edu for details
// 
//    LEVer is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
// 
//    LEVer is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
// 
//    You should have received a copy of the GNU General Public License
//    along with LEVer in file "gnu gpl v3.txt".  If not, see 
//    <http://www.gnu.org/licenses/>.
//
//
//***********************************************************************

#include "tracker.h"

int BuildHistoryPath(CSourcePath* historyPath, CSourcePath* path, int occlLookback)
{
	int pathStartGIdx = path->index[0];

	// Find the assigned in-edge
	int histTrackID = -1;
	int histIdx = gAssignedConnectIn[pathStartGIdx];
	if ( histIdx >=0 )
	{
		CSourcePath* histpath = gConnectIn[pathStartGIdx][histIdx];
		// Single segment agreement requirement
		if ( histpath->index.size() > 2	)
			histTrackID = histpath->trackletID;
	}

	if ( histTrackID >= 0 )
	{
		tPathList::iterator histIter = gAssignedTracklets[histTrackID].begin();
		while ( histIter != gAssignedTracklets[histTrackID].end() )
		{
			CSourcePath* curPath = *histIter;

			historyPath->PushPoint(curPath->index[0]);
			++histIter;
		}
	}

	for ( int i=0; i < path->index.size(); ++i )
		historyPath->PushPoint(path->index[i]);

	return histTrackID;
}

int DepthFirstBestPathSearch(CSourcePath path, int bestGIdx, int t, int tEnd, int occlLookback)
{
	bool bFinishedSearch = true;
	if ( t < tEnd )
	{
		for ( int nextPt=0; nextPt < gFrameHash[t].size(); ++nextPt )
		{
			int nextGIdx = GetGlobalIdx(t, nextPt);
			path.PushPoint(nextGIdx);
			double chkCost = GetCost(path.index, 0, 1);
			path.PopPoint();

			if ( chkCost == dbltype::infinity() )
				continue;

			bFinishedSearch = false;

			path.PushPoint(nextGIdx);
			bestGIdx = DepthFirstBestPathSearch(path, bestGIdx, t+1, tEnd, occlLookback);
			path.PopPoint();
		}
	}

	if ( bFinishedSearch && (path.index.size() > 1) )
	{
		CSourcePath historyPath;

		int startGIdx = path.index[0];
		int nextGIdx = path.index[1];

		int historyTrackID = BuildHistoryPath(&historyPath, &path, occlLookback);
		int srcPathIdx = historyPath.index.size() - path.index.size();

		double newPathCost = GetCost(historyPath.index, srcPathIdx, 0);
		if ( newPathCost == dbltype::infinity() )
			return bestGIdx;

		path.trackletID = historyTrackID;
		path.cost = newPathCost;

		if ( gConnectOut[startGIdx].count(nextGIdx) == 0 )
		{
			CSourcePath* newPath = new CSourcePath(path);
			gConnectOut[startGIdx].insert(std::pair<int,CSourcePath*>(nextGIdx, newPath));
			gConnectIn[nextGIdx].insert(std::pair<int,CSourcePath*>(startGIdx, newPath));
		}
		else if ( newPathCost < gConnectOut[startGIdx][nextGIdx]->cost )
		{
			*(gConnectOut[startGIdx][nextGIdx]) = path;
		}

		if ( bestGIdx < 0 || newPathCost < gConnectOut[startGIdx][bestGIdx]->cost )
		{
			bestGIdx = nextGIdx;
		}
	}

	return bestGIdx;
}

void BuildBestPaths(std::map<int,int>& bestOutEdges, int t, int occlLookback)
{
	if ( t-occlLookback < 0 )
		return;

	int numDetections = gFrameHash[t-occlLookback].size();

	int tEnd = std::min<int>(t+gWindowSize, gNumFrames);

	for ( int srcIdx=0; srcIdx < numDetections; ++srcIdx )
	{
		int startGIdx = GetGlobalIdx(t-occlLookback, srcIdx); 
		if ( occlLookback > 0 && gAssignedConnectOut[startGIdx] > 0 )
			continue;

		CSourcePath srcPath;
		srcPath.PushPoint(startGIdx);
		int bestGIdx = DepthFirstBestPathSearch(srcPath, -1, t+1, tEnd, occlLookback);
		if ( bestGIdx >= 0 )
		{
			bestOutEdges[startGIdx] = bestGIdx;
		}
	}
}

