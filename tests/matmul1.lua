-- Adapted from https://github.com/attractivechaos/plb/blob/master/matmul/matmul_v1.lua

-- There are two implementations of matrices
-- 1. Using Ravi number arrays
-- 2. Using userdata 

-- The matrix library offers functions for both implementations
-- The functions that use Ravi number arrays are suffixed with letter R
-- E.g. matrix() and matrixR()

t=require 'ravimatrix'

-- We first set the locals to point to Lua compatible
-- userdata implementation
local matrix, dim, T = t.matrix, t.dim, t.transpose
local print = print

-- Create a new matrix
local function mnew(m: integer, n: integer) 
  local A = matrix(m, n, 0)
  print('using ' .. getmetatable(A).type)
  return A
end

-- Matrix multiplication using Ravi arrays
-- This version avoids using slices - we operate on the 
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

-- Matrix multiplication using userdata (Lua compatible) matrices
-- this version uses userdata implementation
local function mmul_userdata(a, b)
  local t1 = os.clock()
  local m: integer, n: integer = dim(a)
  local o: integer, p: integer = dim(b)
  assert(n == p)
  local x = mnew(m,n)
  local c = T(b); -- transpose for efficiency
  local xdata = x
  local adata = a
  local cdata = c
  local sum
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

-- Generates a matrix 
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

-- Compares two matrices
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
  assert(ravi.compile(gen))
  assert(ravi.compile(mmul, {omitArrayGetRangeCheck=1}))
  assert(ravi.compile(mmul_userdata))
end

local n = arg[1] or 1000;
n = math.floor(n/2) * 2;

-- First we run a benchmark using userdata implementation
local t111 = os.clock()
local a = mmul_userdata(gen(n), gen(n))
local t211 = os.clock()
local control = gen(n) * gen(n)
assert(comp(a,control))
print("Matrix multiplcation using Lua Userdata Matrices: total time taken ", t211-t111)

-- Switch to Ravi number implementation
-- As these locals are captured as upvalues in the above
-- function, we can switch the implementation just by
-- changing the variables!
matrix = t.matrixR
T = t.transposeR
dim = t.dimR

local t1 = os.clock()
local a = gen(n) * gen(n);
local t2 = os.clock()
print("Matrix multiplication using BLAS: total time taken ", t2-t1)

local t11 = os.clock()
local a1 = mmul(gen(n), gen(n));
local t21 = os.clock()
assert(comp(a,a1))
print("Matrix multiplication using Ravi Number arrays: total time taken ", t21-t11)

print 'Test OK'
