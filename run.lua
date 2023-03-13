#!/usr/bin/env luajit

local Image = require 'image'
local matrix = require 'matrix'
local math = require 'ext.math'
local range = require 'ext.range'
local noise = require 'simplexnoise.3d'

math.randomseed(os.time())

local n = 432
local size = matrix{2*n,n}
m = size:lambda(function(i,j)
	local c = matrix()
	-- cell-centered calculations
	local lon = ((i-.5)/n - 1) * math.pi		-- [-pi, pi]
	local lat = ((j-.5)/n - .5) * math.pi	-- [-pi/2, pi/2]
	c[1] = lon
	c[2] = lat
	c[3] = math.cos(lat) * math.cos(lon)	-- x
	c[4] = math.cos(lat) * math.sin(lon)	-- y
	c[5] = math.sin(lat)					-- z
	return c
end)


local seed = range(3):mapi(function() return math.random(65536) end)
local function f(x,y,z,seed)
	local h = 0
	local d = size[1]
	local A = 1
	while d > 1 do
		h = h + A * noise(x+seed[1],y+seed[2],z+seed[3])
		d = d / 2
		A = A / 2
		x = x * 2
		y = y * 2
		z = z * 2
	end
	return h
end

local hs = size:lambda(function(i,j)
	local c = m[i][j]
	return f(c[3],c[4],c[5],seed) - .2
end)

-- ok now erosion
-- 1) get grad at point
-- 2) move some height along grad
-- 3) ???
-- 4) profit
local numPasses = 1
--local numPasses = 10
--local numPasses = 100  -- too much water
for pass=1,numPasses do
	local nhs = matrix(hs)
	print('pass '..pass..'/'..numPasses)
	for i=1,size[1] do
		for j=1,size[2] do
			local c = m[i][j]
			local iR = (i%size[1])+1
			local iL = ((i-2)%size[1])+1
			local jR = (j%size[2])+1
			local jL = ((j-2)%size[2])+1
			local dlon = math.abs((2*math.pi)/size[1])
			local dlat = math.abs(math.pi/size[2] * math.cos(c[2]))
			local grad = matrix{
				(hs[iR][j] - hs[iL][j]) / (2 * dlon),
				(hs[i][jR] - hs[i][jL]) / (2 * dlat),
			}-- * 3
			grad = grad:unit()
			local h = hs[i][j]
			local some = h * math.random()
			nhs[i][j] = nhs[i][j] - some
			local dsti = math.round(i+1*grad[1])
			dsti = ((dsti-1)%size[1])+1
			local dstj = math.round(j+1*grad[2])
			dstj = ((dstj-1)%size[2])+1
			nhs[dsti][dstj] = nhs[dsti][dstj] + some
		end
	end
	hs = nhs
end

--[[ now smooth a bit
-- hmm this is always looking bad ... 
local function gaussian(A, kernelSize, sigma, calcMetric)
	local size = A:size()
	kernelSize = kernelSize or 10
	sigma = sigma or 1/100
	return size:lambda(function(i,j)
		local sum = 0
		local total = 0
		for u=-kernelSize,kernelSize do
			for v=-kernelSize,kernelSize do
				local srci = ((i-1+u)%size[1])+1
				local srcj = ((j-1+v)%size[2])+1
				local metric = calcMetric(srci, srcj)
				local fuv = matrix{u, v}
				local dsq = fuv * metric * fuv
				local area = metric:det()
				local infl = math.exp(-dsq / (sigma*sigma)) * area
				sum = sum + A[srci][srcj] * infl
				total = total + infl
			end
		end
		return sum / total
	end)
end
hs = gaussian(hs, 10, 1/100, function(i,j)
	local lon = m[i][j][2]
	local dlon = math.abs((2*math.pi)/size[1])
	local dlat = math.abs(math.pi/size[2] * math.cos(lon))
	return matrix{{dlon, 0},{0,dlat}}
end)
--]]

local function gethrange()
	local hmin = math.huge
	local hmax = -math.huge
	for i=1,size[1] do
		for j=1,size[2] do
			local h = hs[i][j]
			hmin = math.min(hmin, h)
			hmax = math.max(hmax, h)
		end
	end
	return hmin, hmax
end

local function getHistogram(n, hmin, hmax)
	n = n or 200
	local h = matrix{n}
end

local totalArea = 0
local landArea = 0
for i=1,size[1] do
	for j=1,size[2] do
		local c = m[i][j]
		local dlon = (2*math.pi)/size[1]
		local dlat = math.pi/size[2] * math.cos(c[2])
		local dA = dlat * dlon
		totalArea = totalArea + dA
		if hs[i][j] > 0 then
			landArea = landArea + dA
		end
	end
end
local hmin, hmax = gethrange()
print('h range', hmin, hmax)
print('total area', totalArea)
print('total area error', totalArea / (4 * math.pi))
print('land area', landArea)
print('land area percent', (landArea / totalArea * 100)..'%')

-- [[
for i=1,size[1] do
	for j=1,size[2] do
		local h = hs[i][j]
		-- map above zero to [1,0] and below 0 to [0,-1]
		if h > 0 then
			h = h / hmax
		else
			h = -h / hmin
		end
		hs[i][j] = h
	end
end
--]]

local function color(r,g,b) return matrix{r,g,b} / 255 end
local landGrad = matrix{
	color(0x73, 0x19, 0x02),
	color(0x87, 0x09, 0x00),
	color(0xbf, 0x4a, 0x06),
	color(0xf2, 0xb3, 0x04),
	color(0x9a, 0xa7, 0x35),
	color(0x34, 0x88, 0x3b),
	color(0x27, 0xa5, 0x2a),
}
local seaGrad = matrix{
	color(0x00, 0x00, 0xff),
	color(0x00, 0x00, 0x3f),
}
local function lookupLinear(m, f)
	local i = f*(#m-1) + 1
	local j = math.floor(i)
	local s = i - j
	if i == #m then
		s = 1
		j = #m-1
	end
	return m[j] * (1 - s) + m[j+1] * s
end

-- hmm, I've got my Image library, and I've got my matrix lua library, and I've got my matrix ffi library ... hmmmmmmm
local hmin, hmax = gethrange()
local img = Image(size[1], size[2], 3, 'double', function(i,j)
	local h = hs[i+1][j+1]
	if h > 0 then
		return lookupLinear(landGrad, h):unpack()
	else
		return lookupLinear(seaGrad, -h):unpack()
	end
end)
img:save'out.png'
