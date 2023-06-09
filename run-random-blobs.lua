#!/usr/bin/env luajit

local Image = require 'image'
local matrix = require 'matrix'
local math = require 'ext.math'
local range = require 'ext.range'

math.randomseed(os.time())

local function randomunitvec()
	local lon = (math.random() - .5) * 2 * math.pi
	local lat = math.acos(math.random() * 2 - 1) 
	return matrix{
		math.sin(lat) * math.cos(lon),
		math.sin(lat) * math.sin(lon),
		math.cos(lat),
	}
end

local n = 432
m = matrix{2*n,n}:lambda(function(i,j)
	local c = matrix()
	-- cell-centered calculations
	local lon = ((i-.5)/n - 1) * math.pi		-- [-pi, pi]
	local lat = ((j-.5)/n - .5) * math.pi	-- [-pi/2, pi/2]
	c[1] = lon
	c[2] = lat
	c[3] = math.cos(lat) * math.cos(lon)	-- x
	c[4] = math.cos(lat) * math.sin(lon)	-- y
	c[5] = math.sin(lat)					-- z
	c[6] = 0								-- h
	return c
end)

-- ok so
-- earth height is antipodal
-- if one point on earth is above water then the opposite is below water
-- so how about we pick N (# of continents) random planes
-- and nudge the surface in those directions
-- and then with whats left, calc the altitude

local numPlanes = 7
local vs = range(numPlanes):mapi(function()
	local v = randomunitvec()
	v[4] = (math.random() * 3 - 1) * .01	-- how much to push? log?
	v[5] = (math.random() * .5 + .5) * .25 * math.pi	-- radius in radians
	return v
end)
-- [[
-- space them each out a bit?
local vavg = vs:mapi(function(v)
	return matrix{v[1], v[2], v[3], 1} * v[4] * v[5] * v[5]	-- TODO weight by push and size
end):sum()
vavg = vavg / vavg[4]
vavg = matrix{
	vavg[1] * (math.random() * .5 + 1),
	vavg[2] * (math.random() * .5 + 1),
	vavg[3] * (math.random() * .5 + 1),
}
vs = vs:mapi(function(v)
	local is = {1,2,3}
	local nv = (v(is) - vavg(is)):unit()
	v[1] = nv[1]
	v[2] = nv[2]
	v[3] = nv[3]
	return v
end)
--]]
vs:sort(function(a,b)
	local is = {1,2,3}
	return math.deg(math.acos(a(is) * vs[1](is)))
		< math.deg(math.acos(b(is) * vs[1](is)))
end)
print'normal dots:'
print(matrix{#vs,#vs}:lambda(function(i,j)
	local is = {1,2,3}
	return math.deg(math.acos(vs[i](is) * vs[j](is)))
end))
for _,v in ipairs(vs) do
	print('normal', v)
	local eps = v[4]
	local theta = v[5]
	--local costh = math.cos(theta)
	--local dpow = math.random() * 7 + 1
	--local dpow = 1/(math.random() * 7 + 1)
	for i=1,2*n do
		for j=1,n do
			local c = m[i][j]
			local d = c[3]*v[1] + c[4]*v[2] + c[5]*v[3]
			-- d is the dot product with the sealevel surface coordinate, so it's 1 aligned with the normal, -1 opposite the normal
			--d = math.abs(d) ^ dpow * math.sign(d)
			--d = d > .5 and .1 or -.1
			--d = d > .3 and .1 or -.1
			--d = d > 0 and .1 or -.1
			--d = d > costh and eps or (d < -costh and -eps or 0)
			d = math.acos(math.clamp(d, -1, 1))
			d = d < theta and 1 or 0
			--d = 1 + .2*((1./((d*math.deg(theta)-3.)^2+1.) + .5 * math.tanh(d*math.deg(theta)-3.)))
			--d = 1 + .2*math.exp(-(d*math.deg(theta))^2)
			-- power d so it is more focused.  odd power to preserve sign.
			d = d * eps
			-- scale by epsilon
			c[6] = c[6] + d
		end
	end
end

-- [[ now smooth a bit
local s = 10
for i=1,2*n do
	for j=1,n do
		local sum = 0
		local total = 0
		for u=-s,s do
			local f = 2 * u / s
			local csrc = m[ ((i-1+u)%(2*n))+1 ][j]
			local dlon = math.abs((2*math.pi)/(2*n))
			local dlat = math.abs(math.pi/n * math.cos(csrc[2]))
			local infl = math.exp(-f*f) * dlat * dlon
			sum = sum + csrc[6] * infl
			total = total + infl
		end
		m[i][j][6] = sum / total

		local sum = 0
		local total = 0
		for u=-s,s do
			local f = 2 * u / s
			local csrc = m[i][((j-1+u)%n)+1]
			local dlon = math.abs((2*math.pi)/(2*n))
			local dlat = math.abs(math.pi/n * math.cos(csrc[2]))
			local infl = math.exp(-f*f) * dlat * dlon
			sum = sum + csrc[6] * infl
			total = total + infl
		end
		m[i][j][6] = sum / total
	end
end
--]]

-- now sink things a bit... up to numPlanes since thats the max things could overlap and superposition
-- or maybe I should normalize by range here, or histogram, or something

local function gethrange()
	local hmin = math.huge
	local hmax = -math.huge
	for i=1,2*n do
		for j=1,n do
			local c = m[i][j]
			local h = c[6]
			hmin = math.min(hmin, h)
			hmax = math.max(hmax, h)
		end
	end
	return hmin, hmax
end

local waterPercent = .5
local maxAlt = .01
local hmin, hmax = gethrange()
for i=1,2*n do
	for j=1,n do
		local c = m[i][j]
		c[6] = ((c[6] - hmin) / (hmax - hmin) - waterPercent) * maxAlt
	end
end

local totalArea = 0
local landArea = 0
for i=1,2*n do
	for j=1,n do
		local c = m[i][j]
		local dlon = (2*math.pi)/(2*n)
		local dlat = math.pi/n * math.cos(c[2])
		local dA = dlat * dlon
		totalArea = totalArea + dA
		if c[6] > 0 then
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

-- hmm, I've got my Image library, and I've got my matrix lua library, and I've got my matrix ffi library ... hmmmmmmm
local img = Image(2*n, n, 3, 'double', function(i,j)
	local h = m[i+1][j+1][6]
	local l = (h - hmin) / (hmax - hmin)
	if h < 0 then
		return 0,l,1
	else
		return 1,l,0
	end
end)
img:save'out.png'
