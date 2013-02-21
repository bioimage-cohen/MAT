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

int gImageWidth;
int gImageHeight;

int GetGlobalIdx(int t, int idx)
{
	if ( t >= gFrameHash.size() )
		return -1;

	if ( idx >= gFrameHash[t].size() )
		return -1;

	return (gFrameHash[t][idx]);
}

int GetTime(int globalIdx)
{
	if ( globalIdx >= gDetect.size() )
		return -1;

	return (gDetect[globalIdx].time);
}

unsigned int GetLinearCoordIdx(int X, int Y)
{
	return X + gImageWidth * Y;
}

void MakeCoordSet(std::set<int>& coordSet, int detIdx)
{
	coordSet.clear();

	int numCoord = gDetect[detIdx].pixelCoords.size();
	for ( int i=0; i < numCoord; ++i )
	{
		unsigned int uniqueLinIdx = GetLinearCoordIdx(gDetect[detIdx].pixelCoords[i].X, gDetect[detIdx].pixelCoords[i].Y);
		coordSet.insert(uniqueLinIdx);
	}
}

double CalcConnectedDist(const std::set<int>& coordSet, int detIdx, int nextDetIdx)
{
	SDetection& curDet = gDetect[detIdx];
	SDetection& nextDet = gDetect[nextDetIdx];

	// Don't bother to calculate CC distance when past center of mass velocity threshold
	double comDistSq = SQR(curDet.X - nextDet.X) + SQR(curDet.Y - nextDet.Y);
	if ( comDistSq > SQR(gVMax) )
		return dbltype::infinity();

	// Overlap distance |A intersect B| / min(|A|,|B|);
	int overlapCount = 0;
	for ( int i=0; i < nextDet.pixelCoords.size(); ++i )
		overlapCount += coordSet.count(GetLinearCoordIdx(nextDet.pixelCoords[i].X, nextDet.pixelCoords[i].Y));

	if ( overlapCount > 0 )
	{
		int compSize = std::min(curDet.pixelCoords.size(), nextDet.pixelCoords.size());
		return (1.0 - ((double) overlapCount) / compSize);
	}

	// If no overlap then find the minimum pixel distance
	double minDistSq = dbltype::infinity();
	for ( int i=0; i < curDet.pixelCoords.size(); ++i )
	{
		for ( int j=0; j < nextDet.pixelCoords.size(); ++j )
		{
			double chkDistSq = SQR(curDet.pixelCoords[i].X - nextDet.pixelCoords[j].X) + SQR(curDet.pixelCoords[i].Y - nextDet.pixelCoords[j].Y);
			if ( chkDistSq > minDistSq )
				continue;

			minDistSq = chkDistSq;
		}
	}

	return sqrt(minDistSq);
}

void CalcConnectedDistToFrame(const std::set<int>& coordSet, int detIdx, int time)
{
	SDetection& curDet = gDetect[detIdx];

	if ( time >= gFrameHash.size() )
		return;

	for ( int i=0; i < gFrameHash[time].size(); ++i )
	{
		int nextDetIdx = GetGlobalIdx(time, i);

		double connDist = CalcConnectedDist(coordSet, detIdx, nextDetIdx);
		if ( connDist == dbltype::infinity() )
			continue;

		curDet.connCompDist.insert(std::pair<int,double>(nextDetIdx, connDist));
	}
}

// Precalculates connected component distances for each detection
// up to two frames into the future
int CalcConnectedDistances()
{
	std::set<int> coordSet;

	for ( int t=0; t < gNumFrames-1; ++t )
	{
		for ( int i=0; i < gFrameHash[t].size(); ++i )
		{
			int detIdx = GetGlobalIdx(t,i);

			MakeCoordSet(coordSet, detIdx);
			CalcConnectedDistToFrame(coordSet, detIdx, t+1);
			CalcConnectedDistToFrame(coordSet, detIdx, t+2);
		}
	}

	return 0;
}

// Read text segmentation data from a file
int ReadSegmentationData(char* filename)
{
	FILE* fp;

	fp = fopen(filename, "r");
	if ( !fp )
		return -1;

	fscanf(fp, "%d %d\n", &gImageWidth, &gImageHeight);
	fscanf(fp, "%d %d\n\n", &gNumFrames, &gNumDetections);

	gDetect.resize(gNumDetections);
	gFrameHash.resize(gNumFrames);

	int totalRead = 0;

	for ( int t=0; t < gNumFrames; ++t )
	{
		int frameDetections;
		fscanf(fp, "%d\n", &frameDetections);

		gFrameHash[t].resize(frameDetections);

		for ( int ptItr = 0; ptItr < frameDetections; ++ptItr )
		{
			gFrameHash[t][ptItr] = totalRead + ptItr;

			int numPixels;
			SDetection* curPt = &gDetect[totalRead + ptItr];
			fscanf(fp, "%lf %lf %d:", &(curPt->X), &(curPt->Y), &numPixels);

			curPt->time = t;

			// Read in pixel coordinates for each pixel in the detection
			// This is used to calculate connected component distances
			curPt->pixelCoords.resize(numPixels);
			int xPix, yPix;
			for ( int pixItr = 0; pixItr < numPixels; ++pixItr )
			{
				fscanf(fp, " (%d,%d)", &xPix,&yPix);
				curPt->pixelCoords[pixItr].X = xPix;
				curPt->pixelCoords[pixItr].Y = yPix;
			}

			fscanf(fp,"\n");
		}

		totalRead += frameDetections;
	}

	fclose(fp);

	return gNumFrames;
}

int ReadDetectionData(int argc, char* argv[])
{
	int checkResult;

	if ( argc < 1)
		return -1;

	checkResult = ReadSegmentationData(argv[0]);
	if ( checkResult < 0 )
		return -1;

	checkResult = CalcConnectedDistances();
	if ( checkResult < 0 )
		return -1;

	gConnectOut.resize(gDetect.size());
	gConnectIn.resize(gDetect.size());

	gAssignedConnectIn.resize(gDetect.size());
	gAssignedConnectOut.resize(gDetect.size());

	gAssignedTrackID.resize(gDetect.size());

	for ( int i=0; i < gDetect.size(); ++i )
	{
		gAssignedConnectIn[i] = -1;
		gAssignedConnectOut[i] = -1;
		gAssignedTrackID[i] = -1;
	}

	return 1;
}

