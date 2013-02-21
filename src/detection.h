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

// Detection related types and variables
struct SPixelCoord
{
	int X;
	int Y;
};

struct SDetection
{
	// Frame the detection is in
	int time;

	// centroid coordinates
	double X;
	double Y;
	
	// Pixel coordinates
	std::vector<SPixelCoord> pixelCoords;

	// Connected component distances
	std::map<int,double> connCompDist;
};

// Read and initialize detection related data.  All globals listed below must be filled or initialized by the end of this routine.

// Globals:
//  gDetect - filled with detection data for all frames
//  gFrameHash - lists indices (into gDetect) for each frame of movie

//  gConnectOut,gConnectIn - initialized to numPts(total detections) empty std::maps each
//  gAssignedConnectIn - same size as gConnectIn, initialized to -1
//  gAssignedConnectOut - same as gAssignedConnectIn, these are for quick lookup of assigned paths
int ReadDetectionData(int argc, char* argv[]);

int GetGlobalIdx(int t, int idx);
int GetTime(int globalIdx);
int GetLocalIdx(int globalIdx);

