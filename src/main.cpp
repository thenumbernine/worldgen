#include <stdlib.h>
#include <time.h>
#include "Tensor/Tensor.h"
#include "Tensor/Grid.h"
#include "WorldGen/Noise.h"
#include "Image/Image.h"
#include "Common/Macros.h"

auto round(auto x) { return floor(x + .5); }

int main() {
	using namespace WorldGen;
	using namespace Tensor;
	using real = double;
	using real2 = vec<real,2>;
	using real3 = vec<real,3>;

	srand(time(0));

	int n = 432;
	auto size = int2{2*n,n};
	
	auto lonlat = Grid<real2, 2>(size, [&](int2 is) {
		auto [i,j] = is;
		real lon = ((i-.5)/n - 1) * M_PI;// [-pi, pi]
		real lat = ((j-.5)/n - .5) * M_PI;// [-pi/2, pi/2]
		return real2{lon, lat};
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
		return h;
	});

	// ok now erosion
	// 1) get grad at point
	// 2) move some height along grad
	// 3) ???
	// 4) profit
	int numPasses = 0;
	//local numPasses = 10
	//local numPasses = 100  // too much water
	for (int pass = 0; pass < numPasses; ++pass) {
		auto nhs = hs;

		std::cout << "pass " << pass << "/" << numPasses << std::endl;
		for (auto is : hs.range()) {
			auto [i,j] = is;
			auto [lon, lat] = lonlat(is);
			int iR = (i+1)%size.x;
			int iL = (i-1+size.x)%size.x;
			int jR = (j+1)%size.y;
			int jL = (j-1+size.y)%size.y;
			real dlon = fabs((2*M_PI)/(real)size.x);
			real dlat = fabs(M_PI/(real)size.y * cos(lon));
			auto grad = real2{
				(hs(int2{iR,j}) - hs(int2{iL,j})) / dlon,
				(hs(int2{i,jR}) - hs(int2{i,jL})) / dlat,
			};// * 3
			grad = grad.normalize();
			auto h = hs(is);
			auto some = h * frand();
			nhs(is) -= some;
			int dsti = round(i*grad.x);
			dsti = (dsti+size.x)%size.x;
			int dstj = round(j*grad.y);
			dstj = (dstj+size.y)%size.y;
			nhs(int2{dsti,dstj}) += some;
		}
		hs = nhs;
	}

#if 0
// [[ now smooth a bit
// hmm this is always looking bad ... 
local function gaussian(A, kernelSize, sigma, calcMetric)
	local size = A:size()
	kernelSize = kernelSize or 10
	sigma = sigma or 1/100
	return size:lambda(function(i,j)
		local sum = 0
		local total = 0
		for u=-kernelSize,kernelSize do
			for v=-kernelSize,kernelSize do
				local srci = ((i-1+u)%size.x)+1
				local srcj = ((j-1+v)%size.y)+1
				local metric = calcMetric(srci, srcj)
				local fuv = matrix{u, v}
				local dsq = fuv * metric * fuv
				local area = metric:det()
				local infl = exp(-dsq / (sigma*sigma)) * area
				sum = sum + A[srci][srcj] * infl
				total = total + infl
			end
		end
		return sum / total
	end)
end
hs = gaussian(hs, 10, 1/100, function(i,j)
	local lon = m[i][j][2]
	local dlon = fabs((2*M_PI)/size.x)
	local dlat = fabs(M_PI/size.y * cos(lon))
	return matrix{{dlon, 0},{0,dlat}}
end)
//]]


#endif

	auto gethrange = [&]() -> real2 {
		real hmin = INFINITY;
		real hmax = -INFINITY;
		for (auto const h : hs) {
			hmin = std::min(hmin, h);
			hmax = std::max(hmax, h);
		}
		return {hmin, hmax};
	};

	real totalArea = 0;
	real landArea = 0;
	for (auto i : hs.range()) {
		auto h = hs(i);
		auto [lon, lat] = lonlat(i);
		real dlon = (2.*M_PI)/(real)size.x;
		real dlat = M_PI/(real)size.y * cos(lat);
		real dA = dlat * dlon;
		totalArea += dA;
		if (h > 0) {
			landArea += dA;
		}
	}
	auto [hmin, hmax] = gethrange();
	std::cout << "h range " << hmin << " " << hmax << std::endl;
	std::cout << "total area " << totalArea << std::endl;
	std::cout << "total area error " << fabs(1. - totalArea / (4 * M_PI)) << std::endl;
	std::cout << "land area " << landArea << std::endl;
	std::cout << "land area percent " << (landArea / totalArea * 100) << "%" << std::endl;

	for (auto & h : hs) {
		// map above zero to [1,0] and below 0 to [0,-1]
		if (h > 0) {
			h /= hmax;
		} else {
			h = -h / hmin;
		}
	}

	std::array landGrad = {
		real3(0x73, 0x19, 0x02),
		real3(0x87, 0x09, 0x00),
		real3(0xbf, 0x4a, 0x06),
		real3(0xf2, 0xb3, 0x04),
		real3(0x9a, 0xa7, 0x35),
		real3(0x34, 0x88, 0x3b),
		real3(0x27, 0xa5, 0x2a),
	};
	std::array seaGrad = {
		real3(0x00, 0x00, 0xff),
		real3(0x00, 0x00, 0x3f),
	};
	auto lookupLinear = [](auto m, real f) -> real3 {
		real i = f*(m.size()-1);
		int j = (int)i;
		real s = i - j;
		if (j == m.size() - 1) {
			s = 1;
			j = m.size()-2;
		}
		return m[j] * (1 - s) + m[j+1] * s;
	};

	// hmm, I've got my Image library, and I've got my matrix lua library, and I've got my matrix ffi library ... hmmmmmmm
	auto img = std::make_shared<Image::Image>(size);
	for (auto i : hs.range()) {
		real h = hs(i);
		real3 rgb = h > 0
			? lookupLinear(landGrad, h)
			: lookupLinear(seaGrad, -h);
		for (int ch = 0; ch < 3; ++ch) {
			(*img)(i.x, i.y, ch) = (uint8_t)rgb(ch);
		}
	}
	Image::system->write("tmp.png", img);
}
