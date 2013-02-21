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

#include <stdio.h>

#include <list>
#include <vector>
#include <set>
#include <map>
#include <limits>

#include "detection.h"
#include "cost.h"
#include "paths.h"

// Convenience defines
#define SQR(x) ((x)*(x))
#define DOT(x1,y1,x2,y2) ((x1)*(x2) + (y1)*(y2))
#define LENGTH(x,y) (sqrt((SQR(x))+(SQR(y))))
#define SIGN(x) (((x) >= 0.0) ? (1.0) : (-1.0) )

typedef char int8;
typedef short int16;
typedef long int32;
//typedef long long int64;
typedef unsigned char uint8;
typedef unsigned short uint16;
typedef unsigned long uint32;
//typedef unsigned long long uint64;

typedef std::list<CSourcePath*> tPathList;

// Detection related global variables
extern int gNumFrames;
extern int gNumDetections;

extern std::vector<SDetection> gDetect;
extern std::vector<std::vector<int>> gFrameHash;

// Path-graph global variables
extern std::vector<std::map<int,CSourcePath*>> gConnectOut;
extern std::vector<std::map<int,CSourcePath*>> gConnectIn;

// For quick edge lookup from point (like inID/outID)
extern std::vector<int> gAssignedConnectOut;
extern std::vector<int> gAssignedConnectIn;
extern std::vector<int> gAssignedTrackID;

// Final tracklet assignments
extern std::vector<tPathList> gAssignedTracklets;

// Global tracking parameters
extern int gWindowSize;
extern double gVMax;
extern double gCCMax;

