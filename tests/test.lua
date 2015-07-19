x,err=package.loadlib('/github/ravi-matrix/build/Debug/ravimatrix.dll',
                      'luaopen_ravimatrix')
assert(not err)
t=x()

v=t.vector(0)
assert(getmetatable(v).type == "Lua matrix")
assert(#v == 0)

v=t.vectorx(0)
assert(getmetatable(v).type == "Ravi matrix")
assert(#v == 0)

v=t.vector(2, 2.0)
assert(getmetatable(v).type == "Lua matrix")
assert(#v == 2)
assert(v[1] == 2.0)
assert(v[2] == 2.0)

v=t.vectorx(2, 2.1)
assert(getmetatable(v).type == "Ravi matrix")
assert(#v == 2)
assert(v[1] == 2.1)
assert(v[2] == 2.1)

a=t.matrixx(1,2,1.0)
a[2] = 2.0
b=t.matrixx(2,1,1.0)
b[1] = 4.0
b[2] = 3.0
c=a*b
assert(#c == 1)
assert(c[1] == 10.0)

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

print 'test OK'
