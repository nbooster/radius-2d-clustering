#pragma once

#include <set>

#include "Point.h"
#include "Chunk.h"

/*
 * We can increase maxPointValue to DOUBLE_MAX,
 * but then we have to change how we compute the euclidean distance,
 * so that no overflow happens
 * (we divide xd and yd by a factor, compute sqrt(xd * xd + yd + yd),
 * and then multiply the result by that factor).
 *
 * Average distance of 2 points inside a square with length 2R ---> 2R * 0.52... ~ R = maxPointValue
 * AvgDist = (maxPointValue - minPointValue) * (2 + sqrt(2) + 5 * log(sqrt(2) + 1)) / 15;
 */

constexpr double maxPointValue = 1.0e100;
constexpr double minPointValue = -maxPointValue;
constexpr double THRESHOLD = 4.0e9; // <--- change it based on your available RAM

class SpatialStruct
{
public:

	SpatialStruct(std::vector<Point>& data, const double scale, const bool verbose = true);

	SpatialStruct(const std::vector<Point>&& data, const double scale) = delete;

	Integer computeClusters(const bool byXY = true);

	void printClusters(std::ostream& outStream = std::cout) const;

	std::set<std::set<Integer>> getClusters(void) const;

private:

/* Methods */

	void printMessage(const std::string_view message) const;

	template<char parent = 'C'>
	inline Integer getParent(const Integer x) const noexcept;

	void initialize(const std::vector<Point>& data, const double scale);

	bool checkDistance(const Point& A, const Point& B) const noexcept;

	bool compareChunkPoints(const Chunk& A, const Chunk& B) const noexcept;

	template<char Execution = 'S'>
	std::vector<std::pair<Integer, Integer>> neighbourChunks(const Integer index) const;

	template<char Execution = 'S'>
	void visitChunk(const Integer index);

	template<char Execution = 'S'>
	void chunkedSpaceMethod(void);

	template<bool byX = true>
	void axisConnectedComponents(std::mutex& firstMutex, bool& stopThread, bool& firstIsX);

	void connectedComponentsMethod(const bool byXY = false);

/* Fields */

	bool printMessages_{ true };
	bool initialized_{ false };
	bool useChunks_{ false };

	double minX_{ maxPointValue };
	double maxX_{ minPointValue };
	double minY_{ maxPointValue };
	double maxY_{ minPointValue };
	double scale_{ 0.0 };
	double minusScale_{ 0.0 };
	long double scaleSquared_{ 0.0 };

	Integer rows_{ 0 };
	Integer columns_{ 0 };
	Integer columnsMinusOne_{ 0 };
	Integer rowsMinusOne_{ 0 };
	Integer clusters_{ 0 };

	static inline const double threshold_{ THRESHOLD };
	std::vector<Point>& Points_;

	std::vector<Integer> indices_;
	std::vector<Integer> parents_;
	std::vector<Integer> indicesY_;
	std::vector<Integer> parentsY_;

	std::vector<Chunk> chunks_;
	std::vector<Integer> indexList_;
	std::vector<Integer> chunkParents_;

	std::mutex chunkParentMutex_;
};