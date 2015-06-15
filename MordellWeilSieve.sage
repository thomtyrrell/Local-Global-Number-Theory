# This code provides some tools for studying the Brauer-Manin obstruction.  In the preprint 
# "The Brauer-Manin obstruction for curves", V. Scharaschkin outlines a method/algorithm for computing the 
# Brauer-Manin obstruction set for a curve by intersecting it with the Mordell-Weil group of its Jacobian 
# modulo primes.  We implement a simple case of this algorithm - computing intersections one prime at a 
# time - for curves of genus 2.

# HOW TO USE THIS CODE

# Choose your favorite elliptic curve, two parameters a,d (a!=0), and call MordellWeilSieve(E,a,d).  This 
# will create a genus 2 curve and compute a number of pieces of data associated to its Jacobian.  To begin 
# collecting the intersection data associated to this curve, simply call generateData().  If there are 
# 0 intersections at a prime, this provides a proof that the curve has no rational points.  

# UPDATES

# 5.28.15:  Rewrites need testing:  eliminated loop counters, utilized list comprehension to replace 
# CartesianProduct.  Remove "range(2)" junk and replace with dictionaries key=elliptic curve value=MWgroup
# Further speed improvements are possible with the use of yield statements/generators

# Still to come:  a=0, removal of "range(2)", generalization of is_quartic_residue(r,p)

# EXAMPLES

# 1. Bremner's curve y^2 = qx^6 - p can be created with MordellWeilSieve(EllipticCurve([0,-p*q^2]),0,1/q).
# 2. A genus 5 curve y^2 = x^12+... can be made by calling MordellWeilSieve with the option genus=5.

class MordellWeilSieve:		
	# Given an elliptic curve E: y^2 = f(x) = x^3 + bx + c and rational numbers a,d, we
	# create a class to study the Jacobian of the genus 2 curve X:y^2 = f(x^2/d + a).  If
	# genus is set to 5, we instead consider X:y^2 = f(x^4/d + a).
	def __init__(self,E,a,d,genus=2):
		
		#Error Checking.  a=0 currently not implemented.
		if d==0:
			print "Please choose a value for d not equal to 0";
			return None;
		elif a!=0:
			print "Not Implemented";
			return None;
		if genus!=2 or genus!=5:
			print "Chosen value of genus not implemented"
			return None;
		
		self.g = genus;
		genusHack = lambda g: 1 if g==5 else 2;
		
		#Here we record the polynomial f, its coefficients, and the parameters a,d.
		x,y = PolynomialRing(QQ,2,'xy').gens();
		self.f = (1/4)*E.torsion_polynomial(2,x);
		self.a = a;
		self.b = self.f.coefficient(x);
		self.c = self.f.subs(x=0)
		self.d = d;
		
		self.badprimes = [p for (p,k) in list(factor(6*E.discriminant()*d))];
		
		#Below we store data related to the Jacobian of the curve X	
		from sage.schemes.elliptic_curves.jacobian import Jacobian;
		self.factors = [E,Jacobian(y^2 - d*(x-a)*self.f)];

		self.ranks = [self.factors[k].rank(only_use_mwrank = false) for k in range(2)];
		self.torsion = [self.factors[k].torsion_subgroup().points() for k in range(2)];
		self.gens = [self.factors[k].gens() + self.torsion[k] for k in range(2)];
		
		self.positiveRanks = ranks[0]>0 and ranks[1]>0;
		
		#See below for a definition.  We run this for initialization purposes.
		self._clear();

		#We use the collections package to count points on the Jacobian modulo p.
		import collections
		self.data = collections.Counter();
		self.density = dict();

	def _clear(self):
		self._p = None;	  #a prime for use in computation
		self._K = None;
		self._factors1 = [None,None];
		self._ident1 = [None,None];

	def _initialize(self,p):
		self._p = p;
		self._K = GF(p);
		self.cp = self._K(self.c);
		self.dp = self._K(self.d);
		self._factors1 = [self.factors[k].reduction(self._p) for k in range(2)];
		self._ident1 = [self._factors1[k](0,1,0) for k in range(2)];

	# reduceMWGroup(p):
	# Given a prime p, this procedure will take the elements of self.gens, reduce them
	# modulo p, and compute the (finite!) abelian group they generate.
	def reduceMWGroup(self,p):
		self._initialize(p);
		mw = [[self._ident1[k]] for k in range(2)];
		for k in range(2):
			for R in self.gens[k]:
				R1 = self._factors1[k](R);
				genR1 = [R1];
				while genR1[-1]!=self._ident1[k]:
					genR1.append(genR1[-1]+R1);
				mw[k] = [R+S for R in mw[k] for S in genR1];
		return mw;

	# is_quartic_residue(a,p):
	# Given a quadratic residue a and a prime p, determines whether a is a biquadratic 
	# residue.  If p = 3 (mod 4) this is a given (if sqrt(a) is a nonsquare, -sqrt(a) 
	# is).  If p = 1 (mod 4) then we compute a^((p-1)/4).
	#
	# This could be easily generalized to handle an arbitrary map X -> E.
	def _is_quartic_residue(self,r):
		R = r^genusHack(self.g)
		if Mod(self._p,4)==3:
			return R.is_square();
		else:
			return R^((self._p-1)/4)==1

	# _isPointOnX(R,S)
	# Checks whether (R,S) lies on the image of X mod p.
	# As the method requires a prime p to be specified, it should only be called from numPointsOnX(p)
	def _isPointOnX(self,R,S):
		if R==self._ident1[0]:
			return S[0]==0 and S[2]!=0;
		elif S==self._ident1[1]:
			return R[0]==0 and R[2]!=0;
		else:
			return  R.xy()[0]*S.xy()[0] == self.dp*self.cp and self._is_quartic_residue(self.dp*R.xy()[0]);
		
    # numPointsOnX(p):
	# Building off of the previous method numPointsOnX(...) counts for a given 
	# prime how many rational points mod p are in the image of X mod p.
	#
	# With zeroCount=true the loop will terminate as soon as a single point is found.  This option is
	# useful when searching for primes with 0 intersection (thus the name).
	def numPointsOnX(self,p,verbose=true,zeroCount=false):
		mwE, mwF = self.reduceMWGroup(p);
		self.data[p]=0;
		for R in mwE:
			for S in mwF:
				if self._isPointOnX(R,S):
					self.data[p]+=1
					if zeroCount:
						return None;
		if not zeroCount:
			self.density[p]=RR(self.data[p])/RR(len(mwE)*len(mwF));
			if verbose:
				if self.data[p]==1:
					print "There is 1 point on X mod %d with density %f" % (p, self.density[p]);
				else:
					print "There are %d points on X mod %d with density %f" % (self.data[p],p,RR(self.density[p]));

	# generateData(): 
	# Will run a sequence of computations, counting points at each prime p between L and U.
	def generateData(self,L=0,U=500,verbose=true,zeroCount=false):	
		for p in Primes():
			if p >= L and  p <= U and not(p in self.badprimes):
				self.numPointsOnX(p,verbose,zeroCount);
			elif p > U:
				self._clear();
				break;
				
	# printData()
	# displays data/densities for inclusion in thesis
	def printData(self,density=false):
		print "Printing data...";
		for p in self.data:
			print p, self.data[p];
			
		if density:
			print "Printing densities..."
			for p in self.density:
				print p, self.density[p];