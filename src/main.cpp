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
auto sqr(auto x) { return x*x; };

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

template<typename real, typename Grid>
auto gaussian(
	Grid const & A,
	int kernelSize,
	real sigma,
	auto calcMetric
) {
	using Type = Grid::Type;
	using int2 = Tensor::int2;
	using real2 = Tensor::vec<real,2>;
	auto size = A.size;
	return Grid(size, [&](int2 ij) -> Type {
		Type sum = {};
		Type total = {};
		for (int u = -kernelSize; u <= kernelSize; ++u) {
			for (int v = -kernelSize; v <= kernelSize; ++v) {
				auto uv = real2{(real)u, (real)v};
				int2 srcij = ((int2)((real2)ij + uv) + size) % size;
				auto g = calcMetric(srcij);
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
	if (!std::isfinite(f)) {
		//throw Common::Exception() << "lookupLinear got a NaN, this is gonna segfault if i try to look it up";
		return Tensor::vec<real,3>(127,127,127);
	}
	real i = f*(m.size()-1);
	int j = (int)i;
	real s = i - j;
	if (j == m.size() - 1) {
		s = 1;
		j = m.size()-2;
	}
	return m[j] * (1 - s) + m[j+1] * s;
};


template<typename real>
auto makeEquirectangular(Tensor::int2 size) {
	return Tensor::Grid<Tensor::vec<real,2>, 2>(size, [&](Tensor::int2 i) {
		real lon = ((i.x - .5)/(real)size.x - .5) * M_PI * 2.;// [-pi, pi]
		real lat = ((i.y - .5)/(real)size.y - .5) * M_PI;		// [-pi/2, pi/2]
		return Tensor::vec<real,2>{lon, lat};
	});
}

auto clamp(auto x, auto a, auto b) {
	return std::max(std::min(x, b), a);
}

template<typename real, bool squareTheCircle = false>
auto makeAzimuthalEquiDist(Tensor::int2 size) {
	using real2 = Tensor::vec<real,2>;
	return Tensor::Grid<real2, 2>(size, [&](Tensor::int2 i) {
		real2 f = (((real2)i + .5)/(real2)size - .5) * 2.;	//[-1,1]^2
		if (squareTheCircle) {	//square the circle	
			real maxp = std::max(fabs(f.x), fabs(f.y));
			real maxf = maxp == 0 ? 0 : (f / maxp).length();
			f /= maxf;
		}
		real r = f.length();
		r = std::min(r, .999);
		real lat = (r - .5) * M_PI;
		real lon = atan2(f.y, f.x);
		return Tensor::vec<real,2>{lon, lat};
	});
}

template<typename real>
auto makeLambertAzimuthalEquiArea(Tensor::int2 size) {
	using real2 = Tensor::vec<real,2>;
	return Tensor::Grid<real2, 2>(size, [&](Tensor::int2 i) {
		real2 f = (((real2)i + .5)/(real2)size - .5) * 2.;	//[-1,1]^2
		real r = f.length();
		r = std::min(r, .999);
		r = asin(r) / (.5 * M_PI);
		real lat = (r - .5) * M_PI;
		real lon = atan2(f.y, f.x);
		return Tensor::vec<real,2>{lon, lat};
	});
}

template<typename real>
auto makePolesAtCorners(Tensor::int2 size) {
	using real2 = Tensor::vec<real,2>;
	return Tensor::Grid<real2, 2>(size, [&](Tensor::int2 i) {
		real2 f = (((real2)i + .5)/(real2)size - .5) * 2.;	//[-1,1]^2
		real lat = f * real2(-.5, .5);	// flip y either here or in image write
		real lon = f * real2(.5, .5);
		real div = 1. - fabs(lat);
		if (div != 0) lon /= div;
		lat = clamp(lat, -1., 1.);
		lon = clamp(lon, -1., 1.);
		lat *= M_PI * .5;
		lon *= M_PI;
		return Tensor::vec<real,2>{lon, lat};
	});
}

int main() {
	using namespace WorldGen;
	using namespace Tensor;
	using real = double;
	using real2 = vec<real,2>;
	using real3 = vec<real,3>;
	using real2x2 = Tensor::mat<real,2,2>;

	srand(time(0));

#if 1
	auto size = int2(2,1) * 1080;	// 2:1 for equirectangular
#endif
#if 0
	auto size = int2(1080);			// 1:1 for everything else
#endif

#if 0	//equirectangular 2x1 (same as nasa bluemarble pics?)
	auto lonlat = makeEquirectangular<real>(size);
#endif
#if 1	//equi-square
	auto lonlat = makeEquirectangular<real>(size);
#endif
#if 0	// azimuthal equi-dist circle
	auto lonlat = makeAzimuthalEquiDist<real, false>(size);
#endif
#if 0
	auto lonlat = makeLambertAzimuthalEquiArea<real>(size);
#endif
#if 0	// azimuthal equi-dist square
	auto lonlat = makeAzimuthalEquiDist<real, true>(size);
#endif
#if 0 // square with poles at corners
	auto lonlat = makePolesAtCorners<real>(size);
#endif

	//sphere pts
	auto pts = Grid<real3, 2>(size, [&](int2 is) {
		auto [lon, lat] = lonlat(is);
		return real3{
			cos(lat) * cos(lon),
			cos(lat) * sin(lon),
			sin(lat)
		};
	});

	auto dxs = Grid<real2, 2>(size, [&](int2 ij) -> real2 {
		//// analytical ...
		//auto [lon, lat] = lonlat(ij);
		//real dlon = (2.*M_PI)/(real)size.x;
		//real dlat = M_PI/(real)size.y * cos(lat);
		//return {dlon, dlat};
		//// numerical ...
		int2 ijL = (ij - 1 + size) % size;
		int2 ijR = (ij + 1) % size;
		return {
			(pts(int2{ijR.x, ij.y}) - pts(int2{ijL.x, ij.y})).length() / 2.,
			(pts(int2{ij.x, ijR.y}) - pts(int2{ij.x, ijL.y})).length() / 2.,
		};
	});

	auto metrics = Grid<real2x2, 2>(size, [&](int2 ij) {
		auto [dlon, dlat] = dxs(ij);
		// TODO proper inner product from the delta vectors
		return real2x2{{dlon*dlon, 0},{0,dlat*dlat}};
	});

	auto areas = Grid<real,2>(size, [&](int2 i) -> real {
		return dxs(i).product();
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

// nah cuz blur undoes the erosion ...
//	hs = gaussian<real>(hs, 10, 100, [&](int2 ij) { return metrics(ij); });

	// use histogram to determine sealevel at 70% lowest height of all land
// why is this putting random water dots everywhere?
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

#if 1	// normalize altitude to [-1,1] while maintaining 0 (since it is at our desired land covering %age)
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

	// 1) calc sealevel temp by latitude
	auto temps = Grid<real, 2>(size, [&](int2 i) -> real {
		auto [lon, lat] = lonlat(i);
		// https://commons.wikimedia.org/wiki/File:Relationship_between_latitude_vs._temperature_and_precipitation.png
		return 25. - 30. * sqr(lat/(M_PI*.5));	// in Celsius
	});

	// ... decrease temp at land (dont at sea)
	for (auto i : temps.range()) {
		if (hs(i) < 0) {
			temps(i) += 10; // land temps about 10 degrees C less than sea temps
		}
	}

	// blur to get rid of hard edges from land/sea boundary
	temps = gaussian<real>(temps, 10, 100, [&](int2 ij) { return metrics(ij); });

	// TODO blur
	// then blur/smooth/diffuse temp
	// then push temp east thx to prevailing winds

	// ... decrease temp by elevation
	real maxAlt = 9000;	// ht of mt everest ish
	for (auto i : temps.range()) {
		real h = hs(i);

		// snow level ...
		// https://www.researchgate.net/publication/232976650_Variation_in_Temperature_With_Altitude_and_Latitude
		// " One degree increase in latitude is roughly equal to a 122 m decrease in elevation, and to a 0.55 o C temperature decrease."
		// http://www.ces.fau.edu/nasa/module-3/why-does-temperature-vary/elevation.php
		// " For every 100-meter increase in elevation, the average temperature decreases by 0.7°C."
		temps(i) -= std::max(0., h) * maxAlt / 100.;
	}

	{
		real totalArea = 0;
		real landArea = 0;
		real landAreaBy2D = 0;
		for (auto i : hs.range()) {
			real dA = areas(i);
			totalArea += dA;
			if (hs(i) > 0) {
				landArea += dA;
				landAreaBy2D++;
			}
		}
		std::cout << "total area by sphere " << totalArea << std::endl;
		std::cout << "land area by sphere " << landArea << std::endl;
		std::cout << "total area error by sphere " << fabs(1. - totalArea / (4 * M_PI)) << std::endl;
		std::cout << "land area percent by sphere " << (landArea / totalArea * 100) << "%" << std::endl;
		std::cout << "land area percent by 2D projection " << (landAreaBy2D / (real)size.product() * 100) << "%" << std::endl;
	}

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

	std::array iceGrad = {
		real3(0x7f, 0x7f, 0xff),
		real3(0xcf, 0xcf, 0xff),
	};

	std::array tempGrad = {
		real3(0, 0, 0xff),
		real3(0xff, 0, 0),
	};

	{
		auto tempRange = getrange(temps);
		auto img = std::make_shared<Image::Image>(size);
		for (auto i : temps.range()) {
			real temp = temps(i);
			real3 c = lookupLinear(tempGrad, (temp - tempRange[0]) / (tempRange[1] - tempRange[0]));
			*(uchar3*)&(*img)(i.x, i.y) = (uchar3)c;
		}
		Image::system->write("temp.png", img);
	}
	
	// hmm, I've got my Image library, and I've got my matrix lua library, and I've got my matrix ffi library ... hmmmmmmm
	auto img = std::make_shared<Image::Image>(size);
	for (auto i : hs.range()) {
		real h = hs(i);
		real temp = temps(i);
		
		real3 c;
		if (h < 0) {
			if (temp < 0) {
				c = lookupLinear(iceGrad, -h);
			} else {
				c = lookupLinear(seaGrad, -h);
			}
		} else {
			if (temp < 0) {
				c = lookupLinear(snowGrad, h);
			} else {
				c = lookupLinear(landGrad, h);
			}
		}
		*(uchar3*)&(*img)(i.x, i.y) = (uchar3)c;
	}
	Image::system->write("out.png", img);
}
