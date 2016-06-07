def bump(n):
    if n == 0:
        return 1;
    else:
        return n;

def general_discriminant(self,p=0):
    
    T = PolynomialRing(QQ,self.degree(),names='t')
    U.<x> = T[]
    U.<x> = U.quotient_ring(U.ideal(self.polynomial()))
    int_basis = [U(omega.lift()) for omega in self.integral_basis()]
    theta = sum( map(mul,zip(int_basis,T.gens())))  #dot product
    
    m_theta = matrix([(theta*x^i).lift().coefficients() for i in range(self.degree())])
    traces = [(m_theta^i).trace() for i in range(2*self.degree()-1)]
    trace_matrix = matrix([traces[i:i+self.degree()] for i in range(self.degree())])
    disc = trace_matrix.det()
    rel_disc = disc/self.discriminant()
    
    return rel_disc.mod([t^(bump(p))-t for t in T.gens()]).map_coefficients(lambda c: Mod(c,p))