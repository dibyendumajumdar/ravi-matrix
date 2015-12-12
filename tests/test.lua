t=require 'ravimatrix'

v=t.vector(0)
assert(getmetatable(v).type == "Lua matrix")
assert(#v == 0)

v=t.vectorR(0)
assert(getmetatable(v).type == "Ravi matrix")
assert(#v == 0)

v=t.vector(2, 2.0)
assert(getmetatable(v).type == "Lua matrix")
assert(#v == 2)
assert(v[1] == 2.0)
assert(v[2] == 2.0)

v=t.vectorR(2, 2.1)
assert(getmetatable(v).type == "Ravi matrix")
assert(#v == 2)
assert(v[1] == 2.1)
assert(v[2] == 2.1)

a=t.matrix(1,2,1.0)
a[2] = 2.0
b=t.matrix(2,1,1.0)
b[1] = 4.0
b[2] = 3.0
c=a*b
assert(#c == 1)
assert(c[1] == 10.0)

a=t.matrixR(1,2,1.0)
a[2] = 2.0
b=t.matrixR(2,1,1.0)
b[1] = 4.0
b[2] = 3.0
c=a*b
assert(#c == 1)
assert(c[1] == 10.0)

-- 3x3 matrix A
-- 76 25 11
-- 27 89 51
-- 18 60 32
local m: number[] = { 76,27,18; 25,89,60; 11,51,32 }
A = t.matrixR(3,3,m)
assert(getmetatable(m).type == "Ravi matrix")

-- 3*1 matrix bx
-- 10
-- 7
-- 43
local m: number[] = { 10, 7, 43 }
bx = t.vectorR(m)
assert(getmetatable(m).type == "Ravi matrix")

C = A*bx;
assert(#C == 3)
assert(C[1] == 1408.0)
assert(C[2] == 3086.0)
assert(C[3] == 1976.0)

-- 3x3 matrix A
-- 76 25 11
-- 27 89 51
-- 18 60 32

A = t.matrix { {76,27,18}, {25,89,60}, {11,51,32} }

-- 3*1 matrix bx
-- 10
-- 7
-- 43
bx = t.vector { 10, 7, 43 }

C = A*bx;
assert(#C == 3)
assert(C[1] == 1408.0)
assert(C[2] == 3086.0)
assert(C[3] == 1976.0)

A = t.matrixR { {1,2,3}, {4,5,6} }
B = t.matrixR { {-1,-2,-3}, {-4,-5,-6} }
C = A+B
assert(#C == 6)
for i = 1,#C do assert(C[i] == 0.0) end
C = A-B
assert(#C == 6)
for i = 1,#C do assert(C[i] == 2*A[i]) end

function comp(a: number[], b: number[])
  local abs = math.abs
  for i = 1, #a do
    if abs(a[i] - b[i]) > 1e-10 then
      return false
    end
  end
  return true
end
--ravi.dumplua(comp)
if ravi.jit() then
  assert(ravi.compile(comp))
end

Z=t.copyR(B)
assert(comp(B, Z))

--ravi.dumpllvmasm(comp)

-- Solve square matrix (linear equation)
A=t.matrix { {1,4}, {2,5} }
b=t.vector { 3, 6 }
x=t.solve(A, b)
assert(x[1] == -1.0 and x[2] == 2.0)

A=t.matrixR { {1,4}, {2,5} }
b=t.vectorR { 3, 6 }
x=t.solveR(A, b)
assert(x[1] == -1.0 and x[2] == 2.0)

x1=t.solveR(A, b, 'Q')
assert(comp(x, x1))

x1=t.solveR(A, b, 'S')
assert(comp(x, x1))

A=t.matrixR { {1,4}, {2,5} }
B=t.transposeR(A)
C=t.matrixR { {1,2}, {4,5} }

assert(comp(B, C))

print 'test OK'
