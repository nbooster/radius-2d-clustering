#define _USE_MATH_DEFINES

#include <unordered_map>
#include <execution>
#include <algorithm>
#include <iomanip>
#include <sstream>
#include <string>
#include <format>

#include "SpatialStruct.h"

static const auto availableThreads = std::jthread::hardware_concurrency() / 2;

template<class T> requires std::is_arithmetic_v<T>
static std::string formatNumber(T value)
{
	std::stringstream ss{};
	ss.imbue(std::locale(""));
	ss << std::fixed << std::showpoint << std::setprecision(3) << value;
	return ss.str();
}

void SpatialStruct::printMessage(const std::string_view message) const
{
	if (printMessages_)
		std::cout << "\n" << message << "\n";
}

template<char parent>
inline Integer SpatialStruct::getParent(const Integer x) const noexcept
{
	Integer localX{ x };
	auto temp = localX;

	while (true)
	{
		if (parent == 'C')
			temp = chunkParents_[localX];
		else if (parent == 'X')
			temp = parents_[localX];
		else
			temp = parentsY_[localX];

		if (temp == localX)
			return localX;

		localX = temp;
	}
}

SpatialStruct::SpatialStruct(std::vector<Point>& data, const double scale, const bool verbose) :
	printMessages_{ verbose }, scale_{ scale }, minusScale_{ -scale }, scaleSquared_{ scale * scale }, Points_{ data }
{
	if (!data.size() || data.size() >= maxValue || scale <= 0.0L || !std::isfinite(scale * scale))
	{
		std::cout << "\nScale lenght is not positive or is too big, or point set is empty or too big. Structure wasn't created...\n";
		return;
	}

	initialize(data, scale);
	initialized_ = true;
}

void SpatialStruct::initialize(const std::vector<Point>& data, const double scale)
{
	const auto& [minXP, maxXP] = std::minmax_element(std::execution::par_unseq, data.cbegin(), data.cend());
	const auto& [minYP, maxYP] = std::minmax_element(std::execution::par_unseq, data.cbegin(), data.cend(), 
									[](const Point& a, const Point& b) { return a.y() < b.y(); });

	minX_ = (*minXP).x();
	minY_ = (*minYP).y();
	maxX_ = (*maxXP).x();
	maxY_ = (*maxYP).y();

	if (!std::isfinite(std::pow(maxX_ - minX_, 2) + std::pow(maxY_ - minY_, 2)))
	{
		std::cout << "\nNumber range is too big. Result may be incorrect. Structure wasn't created...\n";
		initialized_ = false;
		return;
	}

	printMessage(std::format("Minimum X: {}", formatNumber(minX_)));
	printMessage(std::format("Maximum X: {}", formatNumber(maxX_)));
	printMessage(std::format("Minimum Y: {}", formatNumber(minY_)));
	printMessage(std::format("Maximum Y: {}", formatNumber(maxY_)));

	const auto chunkLength{ scale / M_SQRT2 };
	const auto tmpRows{ ceil((maxY_ - minY_) / chunkLength) };
	const auto tmpColumns{ ceil((maxX_ - minX_) / chunkLength) };

	printMessage(std::format("Scale length : {}", formatNumber(scale)));

	if (tmpRows < threshold_ / tmpColumns)
	{
		useChunks_ = true;
		rows_ = static_cast<Integer>(tmpRows);
		columns_ = static_cast<Integer>(tmpColumns);
		rowsMinusOne_ = rows_ - 1;
		columnsMinusOne_ = columns_ - 1;

		printMessage(std::format("Will use the Chunked Space Method.\n\nNumber of chunks: {}", formatNumber(rows_ * columns_)));
		printMessage(std::format("Rows: {}\n\nColumns: {}", formatNumber(rows_), formatNumber(columns_)));
	}
	else
	{
		printMessage("Will use the Connecteed Components Method.");
		printMessage(std::format("Structure was created for {} points.", formatNumber(data.size())));

		return;
	}

	indexList_.resize(data.size(), maxValue);
	Chunk::setIndices(indexList_);

	Integer numberOfChunks{ rows_ * columns_ };

	chunks_.resize(numberOfChunks);

	Integer index{ 0 };

	for (const auto& point : data)
	{
		Integer y { static_cast<decltype(y)>((point.x() - minX_) / chunkLength) };
		Integer x { static_cast<decltype(x)>((point.y() - minY_) / chunkLength) };

		x -= (x == rows_);
		y -= (y == columns_);

		chunks_[static_cast<std::vector<Chunk, std::allocator<Chunk>>::size_type>(x) * columns_ + y].addIndex(index++);
	}

	chunkParents_.resize(numberOfChunks, 0U);
	std::iota(chunkParents_.begin(), chunkParents_.end(), 0U);

	printMessage(std::format("Structure was created for {} points.", formatNumber(data.size())));
}

Integer SpatialStruct::computeClusters(const bool byXY)
{
	if (!initialized_)
		return 0;

	if (clusters_)
		return clusters_;

	if (useChunks_)
	{
		if (const auto chunks = rows_ * columns_; chunks < 1'000'000 || chunks / Points_.size() > 10)
			chunkedSpaceMethod<'P'>();
		else
			chunkedSpaceMethod();
	}
	else
	{
		connectedComponentsMethod(byXY);
	}

	return clusters_;
}

bool SpatialStruct::checkDistance(const Point& A, const Point& B) const noexcept
{
	const long double xd{ A.x() - B.x() };
	if (xd > scale_ || xd < minusScale_)
		return false;
	const long double yd{ A.y() - B.y() };
	return minusScale_ <= yd && yd <= scale_ && xd * xd + yd * yd <= scaleSquared_;
}

bool SpatialStruct::compareChunkPoints(const Chunk& A, const Chunk& B) const noexcept
{
	for (const auto indexA : A)
	{
		const auto& pointA{ Points_[indexA] };

		for (const auto indexB : B)
			if (checkDistance(pointA, Points_[indexB]))
				return true;
	}

	return false;
}

template<char Execution>
std::vector<std::pair<Integer, Integer>> SpatialStruct::neighbourChunks(const Integer index) const
{
	const Integer i{ index / columns_ };
	const Integer j{ index % columns_ };

	if (Execution == 'S')
	{
		if (i == 0)
		{
			if (j == 0)
				return { {0, 1}, {1, 0}, {1, 1}, {0, 2}, {2, 0}, {2, 1}, {1, 2} };
			else if (j == columnsMinusOne_)
				return { {1, j}, {1, j - 1}, {2, j}, {1, j - 2}, {2, j - 1} };
			else [[likely]]
				return { {0, j + 1}, {1, j}, {1, j - 1}, {1, j - 2}, {1, j + 1}, 
						{0, j + 2}, {1, j + 2}, {2, j}, {2, j - 1}, {2, j + 1} };
		}
		else if (i == rowsMinusOne_)
		{
			if (j == columnsMinusOne_)
				return {};
			else [[likely]]
				return { {i, j + 1}, {i, j + 2} };
		}
		else [[likely]]
		{
			if (j == 0)
				return { {i + 1, 0}, {i + 1, 1}, {i, 1}, {i + 2, 0}, {i + 2, 1}, {i, 2}, {i + 1, 2} };
			else if (j == columnsMinusOne_)
				return { {i + 1, j}, {i + 1, j - 1}, {i + 2, j}, {i + 2, j - 1},  {i + 1, j - 2} };
			else [[likely]]
				return { {i, j + 1}, {i + 1, j}, {i + 1, j - 1}, {i + 1, j + 1}, {i + 2, j}, {i, j + 2}, 
						{i + 2, j - 1}, {i + 2, j + 1}, {i + 1, j + 2}, {i + 1, j - 2} };
		}
	}
	else
	{
		if (i == 0)
		{
			if (j == 0)
				return { {0, 1}, {1, 0}, {1, 1}, {0, 2}, {2, 0}, {2, 1}, {1, 2} };
			else if (j == columnsMinusOne_)
				return { {1, j}, {1, j - 1}, {2, j}, {1, j - 2}, {2, j - 1}, {0, j - 1}, {0, j - 2} };
			else [[likely]]
				return { {0, j + 1}, {1, j}, {1, j - 1}, {1, j - 2}, {1, j + 1}, {0, j + 2}, {1, j + 2}, {2, j}, 
						{2, j - 1}, {2, j + 1}, {0, j - 1}, {0, j - 2} };
		}
		else if (i == rowsMinusOne_)
		{
			if (j == 0)
				return { {i - 1, 0}, {i - 2, 0}, {i - 1, 1}, {i - 2, 1}, {i - 1, 2} , {i, j + 1}, {i, j + 2} };
			else if (j == columnsMinusOne_)
				return { {i, j - 1}, {i, j - 2}, {i - 1, j - 1}, {i - 1, j - 2}, {i - 1, j}, {i - 2, j}, {i - 2, j - 1} };
			else [[likely]]
				return { {i - 1, j - 2}, {i - 1, j - 1}, {i - 1, j}, {i - 1, j + 1}, {i - 1, j + 2}, {i - 2, j - 1}, 
						{i - 2, j}, {i - 2, j + 1}, {i, j - 1}, {i, j - 2}, {i, j + 1}, {i, j + 2} };
		}
		else [[likely]]
		{
			if (j == 0)
				return { {i - 1, 0}, {i - 2, 0}, {i - 1, 1}, {i - 2, 1}, {i - 1, 2}, {i + 1, 0}, {i + 1, 1}, {i, 1}, 
						{i + 2, 0}, {i + 2, 1}, {i, 2}, {i + 1, 2} };
			else if (j == columnsMinusOne_)
				return { {i - 1, j}, {i - 2, j}, {i - 1, j - 1}, {i - 2, j - 1}, {i - 1, j - 2}, {i, j - 2}, {i, j - 1}, 
						{i + 1, j}, {i + 1, j - 1}, {i + 2, j}, {i + 2, j - 1},  {i + 1, j - 2} };
			else [[likely]]
				return { {i, j - 1}, {i, j - 2}, {i - 1, j - 2}, {i - 1, j - 1}, {i - 1, j}, 
						{i - 1, j + 1}, {i - 1, j + 2}, {i - 2, j - 1}, {i - 2, j}, {i - 2, j + 1},
						{i, j + 1}, {i + 1, j}, {i + 1, j - 1}, {i + 1, j + 1}, {i + 2, j}, 
						{i, j + 2}, {i + 2, j - 1}, {i + 2, j + 1}, {i + 1, j + 2}, {i + 1, j - 2 } };
		}
	}
}

template<char Execution>
void SpatialStruct::visitChunk(const Integer index)
{
	if (Execution == 'S')
	{
		const auto& chunk{ chunks_[index] };

		const auto chunkParent{ getParent(index) };

		for (const auto& [first, second] : neighbourChunks(index))
		{
			if (std::cmp_greater_equal(first, rows_) || std::cmp_greater_equal(second, columns_)) [[unlikely]]
				continue;

			const Integer temp{ first * columns_ + second };

			const auto& neighbour{ chunks_[temp] };

			if (neighbour.isEmpty())
				continue;

			if (compareChunkPoints(chunk, neighbour))
				chunkParents_[getParent(temp)] = chunkParent;
		}
	}
	else
	{
		const auto& chunk{ chunks_[index] };

		for (const auto& [first, second] : neighbourChunks<'P'>(index))
		{
			if (std::cmp_greater_equal(first, rows_) || std::cmp_greater_equal(second, columns_)) [[unlikely]]
				continue;

			const Integer temp{ first * columns_ + second };

			const auto& neighbour{ chunks_[temp] };

			if (neighbour.isEmpty())
				continue;

			if (compareChunkPoints(chunk, neighbour))
			{
				chunkParentMutex_.lock();
				chunkParents_[getParent(temp)] = getParent(index);
				chunkParentMutex_.unlock();
			}
		}
	}
}

template<char Execution>
void SpatialStruct::chunkedSpaceMethod(void)
{
	const Integer numberOfChunks{ static_cast<Integer>(chunkParents_.size()) };

	auto chunkTraverseLambda = [this](const Integer fromIndex, const Integer toIndex)
	{
		for (Integer index{ fromIndex }; index < toIndex; ++index)
			if (!chunks_[index].isEmpty())
				visitChunk<Execution>(index);
	};

	if (Execution == 'S')
		chunkTraverseLambda(0, numberOfChunks);
	else
	{
		{
			printMessage(std::format("Using {} threads...", formatNumber(availableThreads)));

			std::vector<std::jthread> threadPool;
			threadPool.reserve(availableThreads);

			const auto chunksPerThread = rows_ * columns_ / availableThreads;

			for (Integer index{ availableThreads }; index--;)
				threadPool.emplace_back(chunkTraverseLambda, index * chunksPerThread, (index + 1) * chunksPerThread - 1);
		}
	}

	clusters_ = 0;
	Integer sum = 0;

	for (Integer index{ numberOfChunks }; index--;)
	{
		if (!chunks_[index].isEmpty())
		{
			const auto chunkParent = &chunkParents_[index];

			clusters_ += (*chunkParent == index);

			if (Execution == 'P')
				*chunkParent = getParent(*chunkParent);

			++sum;
		}
	}

	printMessage(std::format("Empty Chunks are: {} % of total.", formatNumber(100.0 * (numberOfChunks - sum) / numberOfChunks)));
}

std::set<std::set<Integer>> SpatialStruct::getClusters(void) const
{
	if (!initialized_ || !clusters_)
		return {};

	std::unordered_map<Integer, std::set<Integer> > clustersMap;

	if (useChunks_)
	{
		const Integer N{ static_cast<Integer>(chunkParents_.size()) };

		for (Integer i{ N }; i--;)
			for (const auto& index : chunks_[i])
				clustersMap[chunkParents_[i]].insert(index);
	}
	else
	{
		const Integer N{ static_cast<Integer>(Points_.size()) };

		for (Integer index{ N }; index--;)
			clustersMap[parents_[index]].insert(indices_[index]);
	}

	std::set<std::set<Integer>> clustersSet;

	for (auto& [first, second] : clustersMap)
		clustersSet.insert(std::move(second));

	return clustersSet;
}

void SpatialStruct::printClusters(std::ostream& outStream) const
{
	const auto& clusters{ getClusters() };

	outStream << "\n";
	for (auto& cluster : clusters)
	{
		for (auto pointIndex : cluster)
			outStream << formatNumber(pointIndex) << " ";
		outStream << "\n\n";
	}
}

template<bool byX>
void SpatialStruct::axisConnectedComponents(std::mutex& firstMutex, bool& stopThread, bool& firstIsX)
{
	const auto N{ static_cast<Integer>(Points_.size()) };
	const auto constScale{ scale_ };
	const auto minusConstScale{ -scale_ };
	const auto constScaleSquared{ scaleSquared_ };

	auto* indicesPtr{ &indices_ };
	auto* parentsPtr{ &parents_ };

	if (!byX)
	{
		indicesPtr = &indicesY_;
		parentsPtr = &parentsY_;
	}

	auto& indicesRef{ *indicesPtr };
	auto& parentsRef{ *parentsPtr };

	indicesRef.resize(N);
	parentsRef.resize(N);

	std::iota(parentsRef.begin(), parentsRef.end(), 0);
	std::iota(indicesRef.begin(), indicesRef.end(), 0);

	if (byX)
		std::sort(std::execution::par_unseq, indicesRef.begin(), indicesRef.end(), [&](Integer a, Integer b) { return Points_[a] < Points_[b]; });
	else
		std::sort(std::execution::par_unseq, indicesRef.begin(), indicesRef.end(), [&](Integer a, Integer b) { return Points_[a].y() < Points_[b].y(); });

	std::mutex parentsMutex;

	auto threadLambda = [&](Integer k)
	{
		for (Integer i{ k }; i < N && !stopThread; i += availableThreads)
		{
			const auto& indexIPoint{ Points_[indicesRef[i]] };

			for (Integer j{ i + 1 }; j < N; ++j)
			{
				const auto& indexJPoint{ Points_[indicesRef[j]] };

				if (byX)
				{
					const long double xd{ indexJPoint.x() - indexIPoint.x() };

					if (xd > constScale)
						break;

					const long double yd{ indexJPoint.y() - indexIPoint.y() };

					if (minusConstScale <= yd && yd <= constScale && xd * xd + yd * yd <= constScaleSquared)
					{
						parentsMutex.lock();
						parentsRef[getParent<'X'>(j)] = getParent<'X'>(i);
						parentsMutex.unlock();
					}
				}
				else
				{
					const long double yd{ indexJPoint.y() - indexIPoint.y() };

					if (yd > constScale)
						break;

					const long double xd{ indexJPoint.x() - indexIPoint.x() };

					if (minusConstScale <= xd && xd <= constScale && xd * xd + yd * yd <= constScaleSquared)
					{
						parentsMutex.lock();
						parentsRef[getParent<'Y'>(j)] = getParent<'Y'>(i);
						parentsMutex.unlock();
					}
				}
			}
		}
	};

	{
		std::vector<std::jthread> threadPool;
		threadPool.reserve(availableThreads);

		for (Integer index{ availableThreads }; index--;)
			threadPool.emplace_back(threadLambda, index);
	}

	{
		std::lock_guard lock(firstMutex);

		if (stopThread)
			return;

		stopThread = true;
	}

	if (!byX)
		firstIsX = false;

	clusters_ = 0;

	for (Integer index{ N }; index--;)
	{
		clusters_ += (parentsRef[index] == index);

		if (byX)
			parentsRef[index] = getParent<'X'>(parentsRef[index]);
		else
			parentsRef[index] = getParent<'Y'>(parentsRef[index]);
	}
}

void SpatialStruct::connectedComponentsMethod(const bool byXY)
{
	bool firstIsX{ true };
	bool stopThread{ false };
	std::mutex firstMutex;

	auto lambdaByX = [&]()
	{
		axisConnectedComponents<true>(firstMutex, stopThread, firstIsX);
	};

	auto lambdaByY = [&]()
	{
		axisConnectedComponents<false>(firstMutex, stopThread, firstIsX);
	};

	if (!byXY)
	{
		std::jthread threadAxisX(lambdaByX);

		printMessage(std::format("Using {} threads...", formatNumber(availableThreads)));

		return;
	}

	{
		std::jthread threadAxisX(lambdaByX);
		std::jthread threadAxisY(lambdaByY);

		printMessage(std::format("Using all {} threads...", formatNumber(std::jthread::hardware_concurrency())));
	}

	if (!firstIsX)
	{
		std::copy(std::execution::par_unseq, indicesY_.begin(), indicesY_.end(), indices_.begin());
		std::copy(std::execution::par_unseq, parentsY_.begin(), parentsY_.end(), parents_.begin());
	}
}