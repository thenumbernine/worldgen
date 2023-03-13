#include <stdlib.h>
#include <time.h>
#include "Tensor/Tensor.h"
#include "Tensor/Matrix.h"
#include "Tensor/Grid.h"
#include "WorldGen/Noise.h"
#include "Image/Image.h"
#include "Common/Macros.h"

auto rad(auto x) { return x * M_PI / 180.; }
auto deg(auto x) { return x * 180. / M_PI; }
auto round(auto x) { return floor(x + .5); }

template<typename T, int N>
auto getrange(Tensor::Grid<T, N> const & m) {
	T vmin = INFINITY;
	T vmax = -INFINITY;
	for (auto const & x : m) {
		vmin = std::min(vmin, x);
		vmax = std::max(vmax, x);
	}
	return Tensor::vec<T,2>{vmin, vmax};
};


// weighted by surface area
// hs holds the value to be binned
// weights holds how much to contribute to the bin, or none means weight of 1
template<typename real>
auto histogram(
	Tensor::Grid<real, 2> const & hs,
	int numBins = 200,
	std::function<real(Tensor::int2)> weight = [](Tensor::int2) { return 1.; }
) {
	std::vector<real> bins(numBins);
	auto hrange = getrange(hs);
	auto [hmin, hmax] = hrange;
	for (auto i : hs.range()) {
		bins[ (int)((hs(i) - hmin) / (hmax - hmin) * (real)numBins * .999) ] += weight(i);
	}
	return std::make_tuple(bins, hrange);
};

template<typename real>
auto lookupLinear(auto m, real f) {
	real i = f*(m.size()-1);
	int j = (int)i;
	real s = i - j;
	if (j == m.size() - 1) {
		s = 1;
		j = m.size()-2;
	}
	return m[j] * (1 - s) + m[j+1] * s;
};


int main() {
	using namespace WorldGen;
	using namespace Tensor;
	using real = double;
	using real2 = vec<real,2>;
	using real3 = vec<real,3>;

	srand(time(0));

	int n = 1080;
	auto size = int2{2*n, n};
	
	auto lonlat = Grid<real2, 2>(size, [&](int2 i) {
		real lon = ((i.x-.5)/n - 1) * M_PI;// [-pi, pi]
		real lat = ((i.y-.5)/n - .5) * M_PI;// [-pi/2, pi/2]
		return real2{lon, lat};
	});

	auto dxs = Grid<real2, 2>(size, [&](int2 i) -> real2 {
		auto [lon, lat] = lonlat(i);
		real dlon = (2.*M_PI)/(real)size.x;
		real dlat = M_PI/(real)size.y * cos(lat);
		return {dlon, dlat};
	});

	auto areas = Grid<real,2>(size, [&](int2 i) -> real {
		return dxs(i).product();
	});

	auto pts = Grid<real3, 2>(size, [&](int2 is) {
		auto [lon, lat] = lonlat(is);
		return real3{
			cos(lat) * cos(lon),
			cos(lat) * sin(lon),
			sin(lat)
		};
	});

	auto seed = real3([](int i) { return (real)(rand() & 0xffff); });
	auto hs = Grid<real, 2>(size, [&](int2 is) {
		auto pt = pts(is);
		int d = size.x;
		real ampl = 1;
		real h = 0;
		while (d > 1) {
			h += ampl * Noise<real>::noise(pt+seed);
			ampl *= .5;
			pt *= 2;
			d >>= 1;
		}
		return h / 1.5 - .2;
	});

std::cout << "initial h range " << getrange(hs) << std::endl;

	{
		// ok now erosion
		// 1) get grad at point
		// 2) move some height along grad
		// 3) ???
		// 4) profit
		int numPasses = 1;
		//int numPasses = 10;
		//int numPasses = 100;  // too much water
		for (int pass = 0; pass < numPasses; ++pass) {
			auto nhs = hs;

			std::cout << "pass " << pass << "/" << numPasses << std::endl;
			for (auto ij : hs.range()) {
				// TODO error when sampling across the poles
				// needs proper gradient / look up mirror / repeat operations
				int2 ijL = (ij - 1 + size) % size;
				int2 ijR = (ij + 1) % size;
				auto [dlon, dlat] = dxs(ij);
				auto grad = real2{
					(hs(int2{ijR.x, ij.y}) - hs(int2{ijL.x, ij.y})) / (2. * dlon),
					(hs(int2{ij.x, ijR.y}) - hs(int2{ij.x, ijL.y})) / (2. * dlat),
				}.normalize();
				auto h = hs(ij);
				auto some = h * frand();
				nhs(ij) -= some;
				int2 dstij = ((int2)((real2)ij + grad) + size) % size;
				nhs(dstij) += some;
			}
			hs = nhs;
		}
	}

#if 0	// now smooth a bit
	// hmm this is always looking bad ... 
	auto gaussian = [](Grid<real, 2> const & A, int kernelSize, real sigma, auto calcMetric) 
	-> Grid<real, 2>
	{
		auto size = A.size;
		return Grid<real, 2>(size, [&](int2 ij) -> real {
			real sum = 0;
			real total = 0;
			for (int u = -kernelSize; u <= kernelSize; ++u) {
				for (int v = -kernelSize; v <= kernelSize; ++v) {
					auto [i,j] = ij;
					int srci = (int)((real)i+u)%size.x;
					int srcj = (int)((real)j+v)%size.y;
					int2 srcij = {srci,srcj};
					auto g = calcMetric(srcij);
					auto uv = real2{(real)u, (real)v};
					real dsq = uv * g * uv;
					real area = determinant(g);
					real infl = exp(-dsq / (sigma*sigma)) * area;
					sum += A(srcij) * infl;
					total += infl;
				}
			}
			return sum / total;
		});
	};
	hs = gaussian(hs, 10, 1/100, [&](int2 ij) -> mat<real,2,2> {
		auto [dlon, dlat] = dxs(ij);
		return mat<real,2,2>{{dlon, 0},{0,dlat}};
	});
#endif

	// use histogram to determine sealevel at 70% lowest height of all land
#if 0
	{
		int numbins = 200;
		auto [bins, hrange] = histogram<real>(hs, numbins, [&](int2 i) { return areas(i); });
std::cout << "histogram range " << hrange << std::endl;	
//std::cout << "histogram bins " << bins << std::endl;	
		//now add up bins from lowest until we reach 
		real total = std::accumulate(bins.begin(), bins.end(), (real)0, std::plus());
		// total should be 4 pi
std::cout << "histogram total " << total << std::endl;	
		real target = total * .71;	// % water
		real sofar = 0;
		int i = 0;
		for (; i < numbins; ++i) {
			sofar += bins[i];
			if (sofar > target) {
				++i;
				break;
			}
		}
std::cout << "based on bins, land % should be: " << ((1. - sofar / total) * 100.) << "%" << std::endl;
		real hsealevel = (real)i / (real)numbins * (hrange[1] - hrange[0]) + hrange[0];
std::cout << "hsealevel " << hsealevel << std::endl;		
// why is this putting random water dots everywhere?
		//for (auto & h : hs) h -= hsealevel;
		for (auto i : hs.range()) hs(i) -= hsealevel;
	}
#endif 

	{
		real totalArea = 0;
		real landArea = 0;
		for (auto i : hs.range()) {
			real dA = areas(i);
			totalArea += dA;
			if (hs(i) > 0) {
				landArea += dA;
			}
		}
		std::cout << "total area " << totalArea << std::endl;
		std::cout << "total area error " << fabs(1. - totalArea / (4 * M_PI)) << std::endl;
		std::cout << "land area " << landArea << std::endl;
		std::cout << "land area percent " << (landArea / totalArea * 100) << "%" << std::endl;
	}

#if 1
	{
		auto [hmin, hmax] = getrange(hs);
		std::cout << "h range " << hmin << " " << hmax << std::endl;
		for (auto & h : hs) {
			// map above zero to [1,0] and below 0 to [0,-1]
			if (h > 0) {
				h /= hmax;
			} else {
				h = -h / hmin;
			}
		}
	}
#endif

	std::array landGrad = {
		real3(0x27, 0xa5, 0x2a),
		real3(0x34, 0x88, 0x3b),
		real3(0x9a, 0xa7, 0x35),
		real3(0xf2, 0xb3, 0x04),
		real3(0xbf, 0x4a, 0x06),
		real3(0x87, 0x09, 0x00),
		real3(0x73, 0x19, 0x02),
	};
	
	std::array seaGrad = {
		real3(0x00, 0x00, 0xff),
		real3(0x00, 0x00, 0x3f),
	};

	std::array snowGrad = {
		real3(0xcf, 0xcf, 0xcf),
		real3(0xff, 0xff, 0xff),
	};

	// hmm, I've got my Image library, and I've got my matrix lua library, and I've got my matrix ffi library ... hmmmmmmm
	real arcticCircleLat = 80;
	real maxAlt = 9000;
	//real snowAlt = 500;
	real snowAlt = 3500;
	auto img = std::make_shared<Image::Image>(size);
	for (auto i : hs.range()) {
		auto [lon, lat] = lonlat(i);
		real h = hs(i);
		real3 c;
		if (h < 0) {
			c = lookupLinear(seaGrad, -h);
		} else {
			if (fabs(lat) > rad(arcticCircleLat)) {
				c = lookupLinear(snowGrad, h);
			} else if (h * maxAlt > snowAlt) {
				c = lookupLinear(snowGrad, (h - snowAlt / maxAlt) / (1. - snowAlt / maxAlt) );
			} else {
				c = lookupLinear(landGrad, h / (snowAlt / maxAlt));
			}
		}
		for (int k = 0; k < 3; ++k) {
			(*img)(i.x, i.y, k) = (uint8_t)c(k);
		}
	}
	Image::system->write("tmp.png", img);
}
