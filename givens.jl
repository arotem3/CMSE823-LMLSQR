using LinearAlgebra

m = 3
n = 2
τ = 0.1

# J = randn(m, n)
# b = randn(m)
J = [0.5 0.5; -sqrt(2) sqrt(2); 0.5 0.5]
b = [1;1;1]

xe = (J'*J + τ^2 * I) \ (J'*b)

Q,R = qr(J)
Q = Matrix(Q)
b1 = Q' * b
b2 = zeros(n)
G = Matrix(τ * I(n)) # we will perform the ellimination of tau*I on this matrix

# h = sqrt(G[n,n]^2 + R[n,n]^2)
# c = R[n,n] / h
# s = -G[n,n] / h
# rows = [R[n,:]'; G[n,:]']
# rows = [c -s; s c] * rows
# R[n, :] = rows[1,:]
# G[n, :] = rows[2,:]
# G

# h = sqrt(G[n-1,n-1]^2 + R[n-1,n-1]^2)
# c = R[n-1,n-1] / h
# s = -G[n-1,n-1] / h
# rows = [R[n-1,:]'; G[n-1,:]']
# rows = [c -s; s c] * rows
# R[n-1,:] = rows[1,:]
# G[n-1,:] = rows[2,:]
# G

# h = sqrt(G[n-1,n]^2 + R[n,n]^2)
# c = R[n,n] / h
# s = -G[n-1,n] / h
# rows = [c -s; s c] * [R[n, :]'; G[n-1, :]']
# R[n,:] = rows[1,:]
# G[n-1,:] = rows[2,:]
# G

for i in n:-1:1
    for j in i:n
        h = sqrt(G[i, j]^2 + R[j, j]^2)
        c = R[j, j] / h
        s = -G[i, j] / h
        rows = [c -s; s c] * [R[j,:]'; G[i,:]']
        R[j, :] = rows[1,:]
        G[i, :] = rows[2,:]

        # apply givens rotation to rhs
        rows = [c -s; s c] * [b1[j]; b2[i]]
        b1[j] = rows[1]
        b2[i] = rows[2]
    end
end

x = R \ b1

norm(x - xe) / norm(xe)