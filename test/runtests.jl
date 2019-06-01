using SRCWA,Test

@testset "reflection" begin

nx,ny,dnx,dny=grid_n(0,0)
lambda=1000
theta=1E-30
phi=0
k0,Kx,Ky,kin=grid_k(nx,ny,theta,phi,lambda,100,100)
upper=halfspace(Kx,Ky,4)
lower=halfspace(Kx,Ky,1)
V0,W0,Kz0=modes_freespace(Kx,Ky)
Su=matrix_ref(upper,V0,W0)
Sl=matrix_tra(lower,V0,W0)
S=concatenate(Su,Sl)
a0te,a0tm=prepare_source(kin,upper.W,4,0,0)
cref=S.S11*a0tm#+a0tm
ctra=S.S21*a0tm
R=a2p(cref,upper,Kx,Ky,kin[3])
T=a2p(ctra,lower,Kx,Ky,kin[3]) 
    
    @test 1/9>R
    @test 8/9>T
end