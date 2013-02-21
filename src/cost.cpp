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

#undef max
#undef min

const double costEpsilon = 1e-3;

double GetCCDistance(int node0, int node1)
{
	if ( gDetect[node0].connCompDist.count(node1) == 0 )
		return (gCCMax + 1.0);

	return (gDetect[node0].connCompDist[node1]);
} // GetCCDistance

double GetCOMDistance(int node0, int node1)
{
	double distSq = (SQR(gDetect[node0].X-gDetect[node1].X) + SQR(gDetect[node0].Y-gDetect[node1].Y));
				
	return sqrt(distSq);
} // GetCOMDistance

double GetWeightedDistance(int node0, int node1, double vmax, double ccmax)
{
	double hd,sd;
	double nmax,nmin,cd;

	hd = GetCOMDistance(node0, node1);
	if ( hd > vmax )
		return dbltype::infinity();

	nmax = std::max<int>(gDetect[node0].pixelCoords.size(), gDetect[node1].pixelCoords.size());
	nmin = std::min<int>(gDetect[node0].pixelCoords.size(), gDetect[node1].pixelCoords.size());

	cd = GetCCDistance(node0, node1);
	if ( cd > ccmax )
		return dbltype::infinity();
		
	sd = (nmax - nmin) / nmax;

	return (10.0*hd + 100.0*sd + 1000.0*cd);

} // GetWeightedDistances

double GetCost(const std::vector<int>& path, int srcFrameIdx, int bCheck)
{
	int k;
	double dlcd;
	int startIdx;
	double dlocnX,dlocnY;

	double LocalCost = 0.0;
	double LocationCost = 0.0; 
	double TotalCost = 0.0;
	double OcclusionCost=1.0;


	if ( (path.size() - srcFrameIdx) < 2 )
		return dbltype::infinity();
	
	if (bCheck)
		startIdx = path.size() - 2;
	else
	{				
		int tStart;
			
		tStart = GetTime(path[srcFrameIdx]) - gWindowSize+1;
		tStart = std::max<int>(0, tStart);
		startIdx = srcFrameIdx;
		while ( (GetTime(path[startIdx]) > tStart) && (startIdx > 0) )
			startIdx--;
	}

	for ( k=startIdx; k < path.size()-1; ++k )
	{
		dlcd = GetCOMDistance(path[k], path[k+1]);
		if ( dlcd > gVMax )
			return dbltype::infinity();		
	}

	OcclusionCost = 1.0;
	for ( k=startIdx; k < path.size()-1; ++k )
		OcclusionCost += GetTime(path[k+1]) - GetTime(path[k]) - 1;
	
	if (bCheck)
		return 1.0;

	LocalCost = 3 * GetWeightedDistance(path[srcFrameIdx], path[srcFrameIdx+1], gVMax, gCCMax);
		
	if ( LocalCost == dbltype::infinity() )			
		return dbltype::infinity();

	if ( srcFrameIdx > 0 )
		LocalCost += GetWeightedDistance(path[srcFrameIdx-1], path[srcFrameIdx+1], 2*gVMax, 2*gCCMax);
	else
		LocalCost *= 2;
	
	if ( LocalCost == dbltype::infinity() )			
		return dbltype::infinity();
	
	if ( (srcFrameIdx < path.size()-2) )
		LocalCost += GetWeightedDistance(path[srcFrameIdx], path[srcFrameIdx+2], 2*gVMax, 2*gCCMax);
	else
		LocalCost *= 2;
	
	if ( LocalCost == dbltype::infinity() )
			return dbltype::infinity();

	dlocnX = gDetect[path[srcFrameIdx]].X;
	dlocnY = gDetect[path[srcFrameIdx]].Y;
	for (  k=startIdx; k < srcFrameIdx; ++k )	
	{ 
		dlocnX += gDetect[path[k]].X;
		dlocnY += gDetect[path[k]].Y;
	}
	dlocnX /= (srcFrameIdx-startIdx+1);
	dlocnY /= (srcFrameIdx-startIdx+1);

	for ( k=srcFrameIdx; k < path.size(); ++k )
	{
		LocationCost += SQR(gDetect[path[k]].X - dlocnX) + SQR(gDetect[path[k]].Y - dlocnY);
	}
	LocationCost /= (path.size()-srcFrameIdx);
	LocationCost = sqrt(LocationCost);

	TotalCost = LocalCost + LocationCost;
	if ( path.size() < 2*gWindowSize + 1 )
	{
		double LengthPenalty;
		LengthPenalty = (2*gWindowSize + 1) - path.size();
		TotalCost = 2*TotalCost*LengthPenalty; 
	}

	if (OcclusionCost > 1.0)
		OcclusionCost *= 2;

	TotalCost *= OcclusionCost;

	if ( TotalCost < costEpsilon )
		TotalCost = costEpsilon;

	return TotalCost;
}

