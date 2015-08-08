-- Adapted from https://github.com/attractivechaos/plb/blob/master/matmul/matmul_v1.lua
-- dummy cast

local x,err=package.loadlib('/github/ravi-matrix/build/Debug/ravimatrix.dll',
                      'luaopen_ravimatrix')
assert(not err)
local t=x()
local matrix, dim, T = t.matrix, t.dim, t.transpose
local print = print

local function mnew(m: integer, n: integer) 
  local A = matrix(m, n, 0)
  print('using ' .. getmetatable(A).type)
  return A
end

-- this version avoids using slices - we operate on the 
-- one dimensional array
local function mmul(a: number[], b: number[])
  local t1 = os.clock()
  local m: integer, n: integer = dim(a)
  local o: integer, p: integer = dim(b)
  assert(n == p)
  local x: number[] = mnew(m,n)
  local c: number[] = T(b); -- transpose for efficiency
  local xdata: number[] = x
  local adata: number[] = a
  local cdata: number[] = c
  local sum: number
  local cj: integer
  local xi: integer
  local t,s
  for i = 1, m do
    xi = (i-1)*m
    for j = 1, p do
      sum = 0.0;
      cj = (j-1)*p
      for k = 1, n do 
        sum = sum + adata[xi+k] * cdata[cj+k] 
      end
      xdata[xi+j] = sum;
    end
  end
  local t2 = os.clock()
  print("mmul: time ", t2-t1)
  return x;
end


local function gen(n: integer)
  local t1 = os.clock()
  local data = mnew(n, n)
  local tmp: number = 1.0 / n / n;
  local ri: integer 
  for i = 1, n do
    ri = (i-1)*n
    for j = 1, n do
      data[ri+j] = tmp * (i - j) * (i + j - 2) 
    end
  end
  local t2 = os.clock()
  print("gen: time ", t2-t1)
  return data;
end

local function comp(a, b)
  local abs = math.abs
  for i = 1, #a do
    if abs(a[i] - b[i]) > 1e-10 then
      return false
    end
  end
  print 'comp ok'
  return true
end

if ravi.jit() then
  --ravi.dumplua(gen2)
  assert(ravi.compile(gen))
  assert(ravi.compile(mmul, {omitArrayGetRangeCheck=1}))
  --ravi.dumpir(gen2)
end

local n = arg[1] or 1000;
local use_ravi_matrix = arg[2] or false

if use_ravi_matrix then
  matrix=t.matrixR
  T=t.transposeR
  dim=t.dimR
end

n = math.floor(n/2) * 2;

local t1 = os.clock()
local a = gen(n) * gen(n);
local t2 = os.clock()
print("total time taken ", t2-t1)

if use_ravi_matrix then
  local t11 = os.clock()
  local a1 = mmul(gen(n), gen(n));
  local t21 = os.clock()
  print("total time taken ", t21-t11)
  assert(comp(a,a1))
end
