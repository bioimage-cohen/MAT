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

#include <vector>
#include <list>
#include <map>

typedef std::numeric_limits<double> dbltype;

class CSourcePath
{
public:
	CSourcePath()
	{
		trackletID = -1;
		cost = dbltype::infinity();
	}

	void PushPoint(int globalIdx)
	{
		index.push_back(globalIdx);
	}

	void PopPoint()
	{
		if ( index.size() <= 1 )
			return;

		index.pop_back();
	}

public:

	int trackletID;
	double cost;

	std::vector<int> index;
};

int GetGlobalIdx(int t, int idx);
int GetTime(int globalIdx);
int GetLocalIdx(int globalIdx);
void BuildBestPaths(std::map<int,int>& bestOutEdges, int t, int occlLookcback = 0);

