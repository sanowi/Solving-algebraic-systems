//%
//%##########################################################################
//%
//%     Copyright (C) 2011 - Iwona Skalna
//%     AGH University of Science and Technology
//%     Department of Applied Computer Science
//%
//%     Module: Neumaier-Pownuk method for solving large scale problems
//%             of the form (K + B*D(p)*A)x= a + F*b
//%
//%##########################################################################
//%

#include "../vector_matrix/vector_matrix_base.h"
#include "../utils/inverse.h"
#include "../utils/gausselim.h"
#include "../utils/svd.h"

const real NEUM_EPS = 1.0e-10;

bool
//%----------------------------------------------------------------------------
//% Iterative method for solving systems of the form
//% K + B*D(p)*A = a + F*b
//% from paper "Linear systems with large uncertainties..."
//% of Neumaier and Pownuk
//%----------------------------------------------------------------------------
neum(const dmatrix& K, 
	const dmatrix& A, 
	const dmatrix& B, 
	const imatrix& D,
	const dmatrix& F, 
	const dvector& a, 
	const ivector& b, 
	ivector& w, 
	int& niter)
{
	int	i, n = A.num_rows(); //% ne - number of edges, dof = ne + 1
	real err;
	dmatrix	C(n, n), D0(n, n), ACB(n, n), AC(n, n);
	imatrix	DD(n, n);
	ivector	u(n), v(n), v1(n), v2(n), d(n), d1(n), d2(n), z(n), cc(n), g(n);

	mid(D, D0);
	C = (K + B * D0 * A);
	inv(C); //% C = (K + B * D0 * A)^{-1}
	AC = A * C;
	ACB = (AC * B);
	DD = D0 - D;
	cc = (AC * F) * b;
	g = AC * a + cc; //% v = A*C*a + A*C*F*b

	//% computing initial enclosures [v] and [d] (Neumaier, Pownuk, pp. 12)
	dvector w0(n), wp(n), wpp(n), mg(n);
	dmatrix mDD(n, n), mACB(n, n);
	
	mag(DD, mDD); //% |D0-D|
	mag(ACB, mACB); //% |A*C*B|
	mag(g, mg); //% |ACa+ACFb|

	w0.fill_in(1.0);
	wp = w0 - mDD * mACB * w0; //% w - |D0-D||ACB|w
	wpp = mDD * mg; //% |D0-D||ACa+ACFb|

	if (wp <= 0.0) return false;
	real alpha = wpp[0] / wp[0];
	for (i = 1; i < n; i++) {
		alpha = ISLIB_MAX(alpha, wpp[i] / wp[i]);
	}
	d.fill_in(interval(-alpha, alpha)); //% initial d
	v = g + ACB * d; //% initial v
	d = DD * v; //% improved d
	//% iteration
	niter = 0;
	do {		
		v1 = (g + ACB * d); //% v = ACa + ACFb + ACBd
		for (i = 0; i < n; i++) {
			v1[i] = v[i] & v1[i];
		}
		d1 = DD * v1; // d1 = ( D0 - D ) * v;
		for (i = 0; i < n; i++) {
			d1[i] = d[i] & d1[i];
		}
		err = 0.0;
		for (i = 0; i < n; i++) {
			//% err is used in a stopping criterion
			err = err + abs(v1[i].inf() - v[i].inf()) +
				abs(v1[i].sup() - v[i].sup()) +
				abs(d1[i].inf() - d[i].inf()) +
				abs(d1[i].sup() - d[i].sup());
		}
		v = v1;
		d = d1;
		niter++;
	} while (err > NEUM_EPS);
	w = C * a + (C * F) * b + (C * B) * d;
	return true;
} //% >>> neum <<<

bool
//%----------------------------------------------------------------------------
//% Iterative method for solving systems of the form
//% A^T*D(p)*A = F*b(p)
//% from paper "Linear systems with large uncertainties..."
//% of Neumaier and Pownuk
//%----------------------------------------------------------------------------
neum(const dmatrix& A, 
	const imatrix& D, 
	const dmatrix& F, 
	const ivector& b, 
	ivector& w)
{
	int	i, j, l = 0, n = A.num_rows();; // ne - number of edges, dof = ne + 1
	real s, err;
	dmatrix	C(n, n), D0(n, n), Atr(n, n), ACAtr(n, n), AC(n, n);
	imatrix	DD(n, n);
	ivector	u(n), v(n), v1(n), v2(n), d(n), d1(n), d2(n), c(n), z(n), cc;
	
	mid(D, D0);
	Atr = A.transpose();
	C = (Atr * D0 * A);
	inv(C);
	AC = A * C;
	ACAtr = (AC * Atr);
	DD = D0 - D;
	cc = (AC * F) * b;
	c = D0 * cc;
	//% initial enclosures [v] and [d] (Neumaier, Pownuk, pp. 12)
	for (i = 0; i < n; i++) {
		s = 0;
		for (j = 0; j < n; j++) {
			s += pow(c[j].mag(), 2) / D(j, j).inf();
		}
		z[i] = 0.5 * (c[i] + interval(-1, 1) *
			sqrt(pow(c[i].mag(), 2) * (1 - D(i, i).sup() / D(i, i).inf()) + D(i, i).sup() * s));
		v[i] = z[i] / D(i, i);
		d[i] = (D0(i, i) / D(i, i) - 1.0) * z[i];
	}
	v = (cc + ACAtr * d); //% initial v
	d = DD * v; //% improved d
	//% iteration
	do {
		v1 = (cc + ACAtr * d);
		for (i = 0; i < n; i++) {
			v1[i] = v[i] & v1[i];
		}
		d1 = DD * v1; //% d1 = ( D0 - D ) * v1;
		for (i = 0; i < n; i++) {
			d1[i] = d[i] & d1[i];
		}
		err = 0.0;
		for (i = 0; i < n; i++) {
			err = err + abs(v1[i].inf() - v[i].inf()) +
				abs(v1[i].sup() - v[i].sup()) +
				abs(d1[i].inf() - d[i].inf()) +
				abs(d1[i].sup() - d[i].sup());
		}
		v = v1;
		d = d1;
		l++;
	} while (err > NEUM_EPS);
	w = (C * F) * b + (C * Atr) * d;
	return true;
} //% >>> neum <<<

