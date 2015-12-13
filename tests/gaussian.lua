-- Gaussian elimination

local ravimatrix = require 'ravimatrix'
local assert = assert
local vector, matrix, dim = ravimatrix.vectorR, ravimatrix.matrixR, ravimatrix.dimR
local slice, numarray = table.slice, table.numarray
local write = io.write

local function cast(a)
	return a
end

local function comp(a: number[], b: number[])
  local abs = math.abs
  for i = 1, #a do
    if abs(a[i] - b[i]) > 1e-10 then
      return false
    end
  end
  return true
end

local function gaussian_elimination(A: number[], b: number[])

  local m: integer, n: integer = dim(A)
  assert(m == n)
  assert(#b == m)
  local columns: table = {}

  -- first get the column slices 
  for k = 1,n do
    columns[k] = slice(A, (k-1)*m+1, m)
  end

  for k = 1,n-1 do -- k is the column
    for i = k+1,m do -- i is the row
      -- obtain the column k
      local column: number[] = @number[] (columns[k]) 
      local multiplier: number = column[i]/column[k]
      write('Performing R(' .. i .. ') = R(' .. i .. ') - m(' .. i .. ',' .. k .. ') * R(' .. k .. '); ', 
            'm(' .. i .. ',' .. k .. ') = ', multiplier, "\n")
      -- For the row i, we need to 
      -- do row(i) = row(i) - multipler * row(k)
      for q = k,n do
      	local col: number[] = @number[] (columns[q])
        col[i] = col[i] - multiplier*col[k]
      end
      b[i] = b[i] - multiplier*b[k]
  	end
  end
end

local A: number[] = matrix { {4,8,-4}, {2,5,1}, {1,5,10} }
local b: number[] = vector { 3,8,5 }

-- ravi.dumplua(gaussian_elimination)
gaussian_elimination(A, b)

local expectedA: number[] = matrix { {4,0,0}, {2,1,0}, {1,3,2} }
local expectedb: number[] = vector { 3,2,2 }

assert(comp(A, expectedA))
assert(comp(b, expectedb))