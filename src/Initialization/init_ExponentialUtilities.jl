eval(quote
         using .ExponentialUtilities: expv!, phiv!, arnoldi!, KrylovSubspace

         # Compute out <- exp(A * NSTEPS * dt) * b
         function _expv(A, b, NSTEPS, dt; hermitian=false, m=min(30, size(A, 1)), tol=1e-7)
             # initialization of the krylov subspace
             TA, Tb = eltype(A), eltype(b)
             T = promote_type(TA, Tb)
             Ks = KrylovSubspace{T,real(T)}(length(b), m)
             arnoldi!(Ks, A, b; m=m, ishermitian=hermitian, tol=tol)

             out = similar(b)
             expv!(out, NSTEPS * dt, Ks)

             return out
         end

         function _phiv(A, b, NSTEPS, dt; hermitian=false, m=min(30, size(A, 1)), tol=1e-7)
             # initialization of the krylov subspace
             TA, Tb = eltype(A), eltype(b)
             T = promote_type(TA, Tb)
             Ks = KrylovSubspace{T,real(T)}(length(b), m)
             arnoldi!(Ks, A, b; m=m, ishermitian=hermitian, tol=tol)

             out = Matrix{Float64}(undef, size(A, 1), 3)
             phiv!(out, NSTEPS * dt, Ks, 2)

             return view(out, :, 3) .* dt^2
         end
     end)

eval(load_krylov_LGG09_homog())
eval(load_krylov_LGG09_inhomog())
