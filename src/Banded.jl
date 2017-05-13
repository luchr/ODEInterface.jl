# Functions and types for banded matrices

import Base: setindex!, getindex
import Base: dump, hash, size, ==, isequal, fill!, convert

"""macro, for importing Banded Matrix types."""
macro import_bandedmatrix()
  :(
    using ODEInterface: BandedMatrix
  )
end

"""macro, for importing Banded Matrix types and functions."""
macro import_bandedmatrixfuncs()
  :(
    using ODEInterface: BandedMatrix, full, setdiagonal!, setdiagonals!,
                        createBandedMatrix
  )
end

"""
  Type for storing and changing a banded matrix.
  
  ## Introduction (connections to proposal #8240)
  
  This is a simple type for storing banded matrices in the band storage
  format used by LAPACK:
  
       http://www.netlib.org/lapack/lug/node124.html
  
  The emphasis here is on *storing* and *changing* the matrix elements,
  because this is needed in the ODE-context.
  
  This is *not* an attempt to fully implement efficient arithmetic
  for banded matrices. This problem is discussed in proposal #8240:
  
       https://github.com/JuliaLang/julia/issues/8240
  
  ## What are banded matrices?
  
  Take a look at the diaognals of a matrix, e.g.
  
       ⎛1 4 2    ⎞
       ⎜5 2 3 1  ⎟
       ⎜  4 3 0 0⎟
       ⎜    3 4 8⎟
       ⎝      2 5⎠
  
  Then we give the following names to the diagonals:
  
                 0  1 2 3  upper diagonals
                   ↘ ↘ ↘ ↘
                 1  1╲4╲2╲0╲0
         lower     ↘╲ ╲ ╲ ╲ ╲          here: l = 1
       diagonals 2  5╲2╲3╲1╲0                u = 2
                   ↘╲ ╲ ╲ ╲ ╲
                 3  0╲4╲3╲0╲0
                   ↘╲ ╲ ╲ ╲ ╲
                    0╲0╲3╲4╲8
                    ╲ ╲ ╲ ╲ ╲
                    0╲0╲0╲2╲5
  
  The 0 diagonal is the main diagonal. The 1st diagonal above this main
  diagonal is called the 1st upper diagonal, etc. 
  The 1st digonal below the main diagonal is called the 1st lower diagonal.
  
  A matrix has lower bandwidth l, if all diagonals below the l lower diagonal
  have only zeros. In the above example the matrix has lower bandwidth l=1
  because the 2nd, 3rd and 4th lower diagonals are all zero.
  
  A matrix has upper bandwidth u, if all diagonals above the u upper diagonal
  have only zeros. In the above example the matrix has upper bandwidth u=2.
  
  A m×n matrix (with m,n≥2) is called banded, if it has a 
  lower bandwidth l<m-2 and/or a upper bandwidth u<n-2.

  ## What is this type for?
  
  This types stores banded matrices. There are functions for banded matrices
  to make it easy to query and change the elements of a banded matrix.
  
  ## Why is this type immutable? Can the entries be changed?
  
  A banded matrix is immutable, i.e. the structure (the upper and lower
  bandwidths cannot be changed). The entries of a banded matrix are *not*
  immutable. 
  
  For an explanation of the storage format, see the help of
  `BandedMatrix_storage`.
  
  There is a constructur, where you can give the entries array for
  the non-zero elements (diagonals) as input argument.
  """
immutable BandedMatrix{T}    <: AbstractMatrix{T}
  m            :: Integer               # number of rows of banded matrix
  n            :: Integer               # number of colums
                                        #    currently: n == size(entries,2)
  l            :: Integer               # lower bandwidth 
  u            :: Integer               # upper bandwidth
  entries      :: Matrix{T}             # array with entries
  
    function (::Type{BandedMatrix{T}}){T}(m::Integer, n::Integer, l::Integer, u::Integer,
                                          entries::Matrix{T})
        m≥1 || throw(ArgumentErrorODE("requirement m≥1", :m))
        n≥1 || throw(ArgumentErrorODE("requirement n≥1", :n))
        l≥0 || throw(ArgumentErrorODE("requirement l≥0", :l))
        u≥0 || throw(ArgumentErrorODE("requirement u≥0", :u))
        l+1 > m &&
            throw(ArgumentErrorODE("requirement: l=$l ≤ m-1 = $m -1", :l))
        u+1 > n &&
            throw(ArgumentErrorODE("requirement: u=$u ≤ n-1 = $n - 1", :u))
        es = size(entries)
        es ≠ (1+l+u, n) &&
            throw(ArgumentErrorODE("requirement: $es=size(entries) == (1+l+u,n); l=$l, u=$u, n=$n",
                                   :entries))
        fill!(entries, zero(T))
        new{T}(m, n, l, u, entries)
    end
end

"""
  ## How are the non-zero entries stored (in this type)?
  
  If M is a m×n matrix with upper bandwidth u and lower bandwith l, then
  the non-zero entries are saved in a (1+l+u)×n matrix N with
  
       N(i-j+u+1,j) = M(i,j) with
  
        max(1,j-u) ≤ i ≤ min(m,l+j)  ⎫   ⎧ -u ≤ i-j ≤ l
                                     ⎬ ⇔ ⎨  1 ≤ i   ≤ m
        max(1,i-l) ≤ j ≤ max(n,u+i)  ⎭   ⎩  1 ≤ j   ≤ n
  
  This is the storage format described in the follwing URL.
  
       http://www.netlib.org/lapack/lug/node124.html
  
  This is not the place to discuss, if a (1+l+u)×min(m,n) matrix is a better
  format for storing the diagonals.

       M(i,j) = N(i-j+u+1,j)  
       N(r,s) = M(r+s-1-u,s)
  
  There can be entries in the matrix N that do *not* correspond to entries
  in M. In this type they are also zero. Great care is taken not to change
  this zeros in order not to violate the `==` and `hash` "contract": Because
  this uncorresponding entries cannot be seen in M, they should have no
  effect on the hash and on the `==` relation.
  
  Keep in mind: There are two different ways of numbering the diagonals.
  The (internal) form ranging from 1 to 1+u+l and the user-friendly form
  ranging von -l ≤ d ≤ u, see `BandedMatrix`.
  """
const BandedMatrix_storage = nothing

function BandedMatrix(m::Integer,n::Integer,l::Integer,u::Integer, initval::Number)
  me = max(0,1+l+u)
  ne = max(0,n)
  bm =  BandedMatrix{typeof(initval)}(m,n,l,u,Matrix{typeof(initval)}(me,ne))
  fill!(bm,initval)
  return bm
end

function hash(bm::BandedMatrix, h::UInt)
  return hash(bm.l,hash(bm.u,hash(bm.m,hash(bm.n,hash(bm.entries,h)))))
end

function ==(bm1::BandedMatrix, bm2::BandedMatrix)
  return (bm1.m==bm2.m) && (bm1.n==bm2.n) && 
         (bm1.l==bm2.l) && (bm1.u==bm2.u) && (bm1.entries==bm2.entries)
end

function isequal(bm1::BandedMatrix, bm2::BandedMatrix)
  return (bm1.m==bm2.m) && (bm1.n==bm2.n) && 
         (bm1.l==bm2.l) && (bm1.u==bm2.u) && 
         isequal(bm1.entries,bm2.entries)
end

function size(bm::BandedMatrix)
  return (bm.m,bm.n)
end

function size(bm::BandedMatrix,d::Integer)
  if d<1
    throw(ArgumentErrorODE("requirement: d≥1, d=$d",:d))
  end
  return  1==d ? bm.m  : 2==d ? bm.n  : 1
end

"""
  The columns in `bm.entries` for diagonal 1 ≤ d ≤ 1+l+u.
  """
function rangeofdiag(bm::BandedMatrix,d::Integer)
  return  max(1,2+bm.u-d):min(bm.n,bm.m+bm.u+1-d) 
end

"""
  tests if `(i,j)` is in the matrix and in the diagonal bands.
  """
function isvalidinband(bm,i,j)
  return (1 ≤ i ≤ bm.m) && (1 ≤ j ≤ bm.n) && (-bm.u ≤ i-j ≤ bm.l)
end

"""
  tests if `(i,j)` is in the diagonal bands.
  
  Caution: This methods does *not* test, if `i` and `j` are
  inside the matrix bounds. Only the diagonal bounds are checked.
  
  This method should only be used, if the matrix bounds/dimensions
  have been checked before.
  """
function isinband(bm,i,j)
  return (-bm.u ≤ i-j ≤ bm.l)
end

function fill!(bm::BandedMatrix,x)
  for d = 1:1+bm.u+bm.l
    bm.entries[d, rangeofdiag(bm,d) ]=x
  end
  return bm
end

"""
  sets the diagonal with number `d` (0 is the main diagonal,see 
  `BandedMatrix`) to the given values.
  """
function setdiagonal!(bm::BandedMatrix,d::Integer,value)
  (-bm.l ≤ d ≤ bm.u) || throw(ArgumentErrorODE(string(
    "requirement: -l ≤ d ≤ u; d=$d, u=",bm.u," l=",bm.l),:d))
  dIndex = 1+bm.u-d
  ra = rangeofdiag(bm,dIndex)
  if isa(value,AbstractArray)
    length(ra) == length(value) || throw(ArgumentErrorODE(string(
      "for d=$d expected length ",length(ra),"; but found ",length(value)),
      :value))
  end
  bm.entries[dIndex,ra] = value
  return value
end

"""
  sets all diagonals at once. `values` must be a 1+u+l cell-array with
  the diagonals starting with the upper right one.
  """
function setdiagonals!(bm::BandedMatrix,values::AbstractArray)
  if length(values) ≠ 1+bm.u+bm.l
    throw(ArgumentErrorODE(string("expected a 1+u+l=",
      1+bm.u+bm.l," array, got length ",length(values)),:values))
  end
  d = bm.u
  for value in values
    setdiagonal!(bm,d,value)
    d -= 1
  end
  return values
end

"""
  sets all diagonals at once copy the diagonals from an other BandedMatrix.
  The other Bandedmatrix must have the same size and same upper and 
  lower bandwidth.
  """
function setdiagonals!{T}(bm::BandedMatrix{T},other::BandedMatrix{T})
  if bm.m ≠ other.m || bm.n ≠ other.n || 
     bm.u ≠ other.u || bm.l ≠ other.l
    throw(ArgumentErrorODE(string("bm has size ",size(bm),
      " and l=",bm.l,", u=",bm.u,"; but other has size ",size(other),
      " and l=",other.l,", u=",other.u),:other))
  end
  bm.entries[:] = other.entries
  return other
end

"""
  For banded matrices: checks if inds are in bands.
  """
function setindex!(bm::BandedMatrix,value,i::Integer,j::Integer)
  isvalidinband(bm,i,j) || throw(BoundsError(bm,(i,j)))
  bm.entries[1+bm.u+i-j,j]=value
end

function setindex!(bm::BandedMatrix,value,ra::UnitRange,j::Integer)
  (isvalidinband(bm,ra.start,j) && isvalidinband(bm,ra.stop,j)) || 
                                   throw(BoundsError(bm,(ra,j))) 
  bm.entries[(ra.start+bm.u+1-j):(ra.stop+bm.u+1-j),j] = value
end

function setindex!(bm::BandedMatrix,value,i::Integer,ra::UnitRange)
  (isvalidinband(bm,i,ra.start) && isvalidinband(bm,i,ra.stop)) ||
                                   throw(BoundsError(bm,(i,ra))) 
  ni = i+bm.u+1-ra.start; nj = ra.start
  for k = 1:length(ra)
    bm.entries[ni,nj] = isa(value,Number)?value:value[k]
    ni -= 1; nj += 1
  end
end

function setindex!(bm::BandedMatrix,value,zra::UnitRange,sra::UnitRange)
  (isvalidinband(bm,zra.start,sra.start) && 
   isvalidinband(bm,zra.stop,sra.stop)) || throw(BoundsError(bm,(zra,sra))) 
  s = 1
  for nj in sra
    ni = zra.start+bm.u+1-nj; 
    for z = 1:length(zra)
      bm.entries[ni,nj] = isa(value,Number)?value:value[z,s]
      ni += 1;
    end
    s += 1
  end
end

function getindex{T}(bm::BandedMatrix{T},i::Integer,j::Integer)
  ((1 ≤ i ≤ bm.m) && (1 ≤ j ≤ bm.n)) || throw(BoundsError(bm,(i,j)))
  return isinband(bm,i,j)?bm.entries[i-j+bm.u+1,j]:zero(T)
end

function getindex{T}(bm::BandedMatrix{T},ind::Integer)
  (j,i) = divrem(ind-1,bm.m)
  return getindex(bm,i+1,j+1)
end

"""
       function isdiagonalempty(A::AbstractArray{T,2},d::Integer)

  method to check if all entries in a diagonal are zero.
  """
function isdiagonalempty{T}(A::AbstractMatrix{T},d::Integer)
  (m,n) = size(A)
  i = (d ≥ 0) ? 1 : 1-d
  j = (d ≥ 0) ? 1+d : 1
  z = zero(T)
  while true
    (A[i,j] ≠ z) && return false
    i+=1; j+=1
    (i>m || j>n) && break
  end
  return true
end

"""
  convert full matrix to BandedMatrix.
  """
function createBandedMatrix{T}(A::AbstractMatrix{T})
  (m,n) = size(A)
  # find last upper diagonal with nonzero entry
  d = n-1
  while d>0
    !isdiagonalempty(A,d) && break
    d-=1
  end
  u = d

  # find last lower diagonal with nonzero entry
  d = -(m-1)
  while d<0
    !isdiagonalempty(A,d) && break
    d+=1
  end
  l = -d
  
  bm = BandedMatrix(m,n,l,u,zero(T))
  for j in 1:bm.n
    for i in  max(1,j-bm.u):min(bm.l+j,bm.m)
      bm.entries[i+bm.u+1-j,j] = A[i,j]
    end
  end

  return bm
end


"""
  Convert banded matrix `bm` to full matrix. Save full matrix values in `f`.
  
  `f` must have the right size.
  """
function fullToArray{T}(bm::BandedMatrix{T},f::AbstractArray{T})
  if size(f) ≠ (bm.m,bm.n)
    throw(ArgumentErrorODE(string("f has wrong size ",size(f),
          "extexted ",(bm.m,bm.n)),:f))
  end
  fill!(f,0)
  for j in 1:bm.n
    for i in  max(1,j-bm.u):min(bm.l+j,bm.m)
      f[i,j] = bm.entries[i+bm.u+1-j,j]
    end
  end
  return f
end

"""
  For banded matrices: generate and return full/dense matrix.
  """
function full{T}(bm::BandedMatrix{T})
  return fullToArray(bm, Array{T}((bm.m,size(bm.entries,2),)))
end

function dump(io::IO, bm::BandedMatrix, n::Integer, indent)
  println(io,typeof(bm)," ",bm.m,"×",bm.n,"; l=",bm.l," and u=",bm.u)
  if n>0
    print(io,indent,"  entries in diagonals ")
    dump(io,bm.entries,n-1,string(indent,"  "))
  end
end


# vim:syn=julia:cc=79:fdm=indent:
