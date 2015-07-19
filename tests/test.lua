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


print 'test OK'
