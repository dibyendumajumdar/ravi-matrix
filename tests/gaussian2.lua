-- Gaussian elimination

local ravimatrix = require 'ravimatrix'
local assert = assert
local vector, matrix, dim = ravimatrix.vectorR, ravimatrix.matrixR, ravimatrix.dimR
local slice, numarray, intarray = table.slice, table.numarray, table.intarray
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

-- i = column
local function partial_pivot(columns: table, nrow: integer[], i: integer, n: integer)
  local p: integer = i
  local max: number = 0.0
  local a: number[] = @number[]( columns[i] )

  -- find the row from i to n that has
  -- max absolute value in column[i]
  for row=i, n do
    local value: number = a[row]
    if value < 0.0 then value = -value end
    if value > max then
      p = i
      max = value
    end
  end
  if a[p] == 0.0 then 
    error("no unique solution exists")
  end
  if nrow[i] ~= nrow[p] then
    local temp: integer = nrow[i]
    nrow[i] = nrow[p]
    nrow[p] = temp
  end
end

local function gaussian_elimination(A: number[], b: number[])

  local m: integer, n: integer = dim(A)
  assert(m == n)
  assert(#b == m)

  -- nrow will hold the order of the rows allowing
  -- easy interchange of rows
  local nrow: integer[] = intarray(n)

  -- As ravi matrices are column major we 
  -- create slices for each column for easy access
  -- the vector b can also be treated as an additional
  -- column thereby creating the augmented matrix 
  local columns: table = {}

  -- we use i as the row and j a the column

  -- first get the column slices 
  for j = 1,n do
    columns[j] = slice(A, (j-1)*m+1, m)
  end
  columns[n+1] = b

  -- initialize the nrow vector
  for i = 1,n do
    nrow[i] = i
  end

  for j = 1,n-1 do -- j is the column
    partial_pivot(columns, nrow, j, m)

    -- debug
    for t=j,#nrow do 
      write('row order for column ' .. j .. ' is nrow[' .. nrow[t] .. ']', "\n")
    end

    for i = j+1,m do -- i is the row
      -- obtain the column j
      local column: number[] = @number[] (columns[j]) 
      local multiplier: number = column[nrow[i]]/column[nrow[j]]
      write('Performing R(' .. i .. ') = R(' .. i .. ') - m(' .. nrow[i] .. ',' .. j .. ') * R(' .. nrow[j] .. '); ', 
            'm(' .. nrow[i] .. ',' .. j .. ') = ', multiplier, "\n")
      -- For the row i, we need to 
      -- do row(i) = row(i) - multipler * row(j)
      for q = j,n+1 do
      	local col: number[] = @number[] (columns[q])
        col[nrow[i]] = col[nrow[i]] - multiplier*col[nrow[j]]
      end
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