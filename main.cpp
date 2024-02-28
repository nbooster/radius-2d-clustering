#include <chrono>
#include <random>

#include "Clustering.h"

/*
 * For numberOfPoints = 1 Billion and:
 * factor = 1.0e3 --> (20 GB RAM, 12 mins) --> chunkedSpaceMethod (for factor < 1.0e3 --> less time and RAM)
 * factor = 1.0e5 --> (24.6 GB RAM, 26 mins) --> connectedComponentsMethod (for factor > 1.0e5 --> less time, same RAM)
 * factor = 1.0e4 or 1.0e5 --> (less than 32 GB RAM, more time) with either method.
 * [ My CPU (AMD FX 8350 - 32nm - 8 cores) and RAM (32 GB Dual-Channel DDR3 669MHz) are more than 10 years old... ]
 */

constexpr double factor = 1.0e10;
constexpr long double scaleLength = maxPointValue / factor;
constexpr Integer numberOfPoints = 100'000'000;
constexpr std::uint32_t randomSeed = 42;
constexpr bool verbose = false;
constexpr bool stableCC = true; /* <--makes the connected components method more stable on the input,
									by sorting in both axis, and using all available threads,
									but it may slow dows if other things are running at the same time. */

std::vector<Point> randomPoints(const Integer N, const Integer seed, const double minVal, const double maxVal)
{
	std::default_random_engine generator;
	generator.seed(seed);
	std::uniform_real_distribution distribution(minVal, maxVal);

	std::vector<Point> Points;
	Points.resize(N);

	for (Integer index{ N }; index--;)
		Points[index] = { distribution(generator), distribution(generator) };

	return Points;
}

int main(int argc, char* argv[])
{
	std::cout << "\nCreating Points...\n";
	std::vector<Point> Points = randomPoints(numberOfPoints, randomSeed, minPointValue, maxPointValue);
	//{ {0,10}, {2,0}, {2,10}, {13,10}, {14,10} };  // --> similar to the picture example, output is: 2 clusters (with scaleLength 10)

	std::cout << "\nComputing Clusters...\n";
	auto begin = std::chrono::steady_clock::now();
	auto result = RadiusClustering::scaleCluster2DPoints(Points, scaleLength, stableCC, verbose);
	auto end = std::chrono::steady_clock::now();

	std::cout << "\nNumber of clusters: " << result;
	auto duration = std::chrono::duration_cast<std::chrono::seconds>(end - begin).count();

	std::cout << "\n\nComputing time was: = " << duration << "[secs]\n" << '\a';
}
