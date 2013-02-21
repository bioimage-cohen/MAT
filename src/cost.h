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



// Get cost from a path list.
// path - vector of indices into gDetect, indicates current path being explored.
// 
// srcFrameIdx - the index into path index vector of the source point (start of new path).
//		srcFrameIdx is trivially 0 if there is no history being used.
// 
// bCheck - If bCheck is true then only do a check to verify the path doesn't violate constraints.
double GetCost(const std::vector<int>& path, int srcFrameIdx,int bCheck);

