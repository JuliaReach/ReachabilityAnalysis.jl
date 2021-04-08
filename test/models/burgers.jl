# 1D Burgers equation: Eq. (6.3) from [1]
# Space discretization: Eq. (6.4) from [1]

# Reference
# [1] Liu, J. P., Kolden, H. Ø., Krovi, H. K., Loureiro, N. F., Trivisa, K., &
# Childs, A. M. (2020). Efficient quantum algorithm for dissipative nonlinear
# differential equations. arXiv preprint arXiv:2011.03185.


# Burgers equation:
@taylorize function burgers!(du, u, p, t)

    np = 4 # discretization points
    nu = 0.05 # kinematic viscosity
    Δx = 1.0/(np-1)
    c1 = nu/Δx^2 # c1 = ν/Δx²
    c2 = 0.25/Δx # c2 = 1/(4Δx)

    du[1] = zero(u[1])
    du[np] = zero(u[np])
    for i in 2:np-1
        du[i] = c1*(u[i-1]-2*u[i]+u[i+1]) - c2*(u[i+1]^2 - u[i-1]^2)
    end

end
