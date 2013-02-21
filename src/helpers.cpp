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

#include <stdio.h>

double DebugGetCost(double& forwardCost, std::vector<double>& partialCosts, const std::vector<int>& path, int srcFrameIdx)
{
	std::vector<int> tempPath;
	tempPath.reserve(path.size());

	partialCosts.reserve(path.size() - srcFrameIdx - 1);

	for ( int i=srcFrameIdx; i < path.size(); ++i )
		tempPath.push_back(path[i]);

	forwardCost = GetCost(tempPath, 0, 0);

	tempPath.clear();

	for ( int i=0; i <= srcFrameIdx; ++i )
		tempPath.push_back(path[i]);

	for ( int i=srcFrameIdx+1; i < path.size(); ++i )
	{
		tempPath.push_back(path[i]);

		double testCost = GetCost(tempPath, srcFrameIdx, 0);
		partialCosts.push_back(testCost);
	}

	return GetCost(path, srcFrameIdx, 0);
}

int maxDigits(const std::vector<int>& index, int startIdx)
{
	int maxElem = 0;
	for ( int i=startIdx; i < index.size(); ++i )
		maxElem = std::max<int>(maxElem, index[i]);

	const int maxAllowedDigits = 20;
	for ( int i=1; i < maxAllowedDigits; ++i )
	{
		maxElem /= 10;

		if ( maxElem == 0 )
			return i;
	}

	return maxAllowedDigits;
}

void PrintDebugPathInfo(const CSourcePath& path, int srcFrameIdx)
{
	FILE* fDebug = fopen("debugPathInfo.txt", "at");

	if ( fDebug == NULL )
		return;

	double memorylessCost;
	std::vector<double> partialCosts;
	double totalCost = DebugGetCost(memorylessCost, partialCosts, path.index, srcFrameIdx);

	int pathEntrySize = maxDigits(path.index, srcFrameIdx);

	if ( pathEntrySize < 7 )
		pathEntrySize = 7;

	char buildPathPrintStr[20];
	char buildCostPrintStr[20];

	_snprintf(buildPathPrintStr, 20, "%%-%dd ", pathEntrySize);
	_snprintf(buildCostPrintStr, 20, "%%-%dg ", pathEntrySize);

	fprintf(fDebug, "------------------------------\n");
	fprintf(fDebug, "Total Cost: %g\n", totalCost);
	fprintf(fDebug, "Memoryless Cost: %g\n\n", memorylessCost);

	int offset = 0;
	for ( int i=0; i <= srcFrameIdx; ++i )
		offset += fprintf(fDebug, "%-d ", path.index[i]);

	for ( int i=srcFrameIdx+1; i < path.index.size(); ++i )
		fprintf(fDebug, buildPathPrintStr, path.index[i]);

	fprintf(fDebug, "\n");
	for ( int i=0; i < offset; ++i )
		fprintf(fDebug, " ");

	for ( int i=srcFrameIdx+1; i < path.index.size(); ++i )
		fprintf(fDebug, buildCostPrintStr, partialCosts[i-srcFrameIdx-1]);

	fprintf(fDebug, "\n\n\n");

	fclose(fDebug);
}

