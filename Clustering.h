#pragma once

#include "SpatialStruct.h"

namespace Altair
{
	Integer scaleCluster2DPoints(std::vector<Point> Points, const double scale, const bool stableCC = false, const bool verbose = false)
	{
		SpatialStruct spatial(Points, scale, verbose);
		return spatial.computeClusters(stableCC);
	}

	std::set<std::set<Integer>> getScaleCluster2DPoints(std::vector<Point> Points, const double scale, const bool stableCC = false, const bool verbose = false)
	{
		SpatialStruct spatial(Points, scale, verbose);
		auto clusters = spatial.computeClusters(stableCC);
		return spatial.getClusters();
	}

	void printScaleCluster2DPoints(std::vector<Point> Points, const double scale, const bool stableCC = false, const bool verbose = false, std::ostream& outStream = std::cout)
	{
		SpatialStruct spatial(Points, scale, verbose);
		auto clusters = spatial.computeClusters(stableCC);
		spatial.printClusters(outStream);
	}
}