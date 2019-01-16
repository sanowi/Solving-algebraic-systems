void
kolev(const omvector& M1, const ovector& b1, const ivector& p1, ivector& w,
		dvector& t1, dvector& t3, dmatrix& t2, dvector& x0)
// Original Kolev's method; calculates only outer solution
{
	int n = b1[0].size();
	int niter = 0;
	int K = p1.size();
	real eps = 1.0e-6;

	dvector b0(n);
	dmatrix A0(n, n), R(n, n);
	ivector p(K);
	omvector M(K + 1);
	ovector b(K + 1);

	ktransform(M1, b1, p1, M, b, p); // the system is transformed so that all [p]=[-1,1]

	A0 = M[0]; // midpoint matrix of the new system
	b0 = b[0]; // midpoint vector of the new system
	R = A0;
	x0 = dvector(n);
	if (inv(R)) { // if the mipoint matrix is nonsingular, the computation continues
		ivector v(n), v1(n);
		dvector t3i(n);
		dmatrix C0(n, K), C(n, K);
		omvector B(K), D(K), D0(K);
		ovector d(K), G(K);

		gausse(A0, b0, x0); // computing x0 - the midpoint solution
		
		//x0 = R * b0; 
		for (int k = 0; k < K; k++) {
			d[k] = R * (b[k + 1] - M[k + 1] * x0); // d^k = R*(b^k-A^k*x0)
			B[k] = R * M[k + 1]; // B^k=R*A^k, tez chyba OK
			for (int j = 0; j < n; j++) {
				C0(j, k) = d[k][j];
			}
		}
		// v^(0)=0
		v = C0 * p; // v^(1)
		for (int k = 0; k < K; k++) {
			D[k] = B[k] * C0;
		}
		// v^(2)=t1+t2*p+t3	
		t1 = dvector(n);
		t3 = dvector(n);
		t2 = dmatrix(n, K);

		t1.fill_in(0.0);
		t2.fill_in(0.0);
		t3.fill_in(0.0);

		for (int i = 0; i < n; i++) { // n is the size of v
			real s = 0.0;
			for (int k = 0; k < K; k++) {
				for (int j = k + 1; j < K; j++) {
					// perhaps this can be estimated better
					// strony 243, 244 z artykulu Koleva
					s += fabs(D[k](i, j) + D[j](i, k));
				}
				t1[i] = t1[i] - D[k](i, k) / 2.0;
				t2(i, k) = C0(i, k); // this is the L operator, which is multiplied by p
				t3[i] = t3[i] + fabs(D[k](i, k) / 2.0);
			}
			t3[i] = t3[i] + s;
		}
		////////////////////////////////////////////////////////////////////////
		// ITERATION
		// computing v^(k), k>=3
		//
		reduce(t1, t2, t3, v);
		do {
			niter++;
			v1 = v;
			for (int k = 0; k < K; k++) {
				G[k] = B[k] * t1;
				for (int j = 0; j < n; j++) {
					C(j, k) = C0(j, k) - G[k][j];
				}
				D0[k] = B[k] * t2;
			}
			t1.fill_in(0.0);
			t3i = t3;
			t3.fill_in(0.0);
			for (int i = 0; i < n; i++) { // n is the size of v
				real s = 0.0;
				real s1 = 0.0;
				for (int k = 0; k < K; k++) {
					for (int j = 0; j < n; j++) {
						s1 += fabs(B[k](i, j)) * t3i[j];
					}
					for (int j = k + 1; j < K; j++) {
						s += fabs(D0[k](i, j) + D0[j](i, k));
					}
					t1[i] = t1[i] - D0[k](i, k) / 2.0;
					t2(i, k) = C(i, k); // this is the L operator, which is multiplied by p
					t3[i] = t3[i] + fabs(D0[k](i, k) / 2.0);
				}
				t3[i] = t3[i] + s + s1;
				// the p-solution is in the form: t1 + t2*p + [t3]
			}
			reduce(t1, t2, t3, v);
		} while (stopCrit(v, v1) > eps);

		reduce(t1, t2, t3, v);
		for (int i = 0; i < n; i++) { // final result for outer solution
			w[i] = v[i] + interval(x0[i]);
		}
	}
	std::cout << "Iterations: " << niter << std::endl;
}


void
kolevorg(const omvector& M1, const ovector& b1, const ivector& p1, ivector& w,
dvector& t1, dvector& t3, dmatrix& t2)
// Kolev's method; calculates both outer and inner solution
{
	int n = b1[0].size();
	int niter = 0;
	int K = p1.size();
	real eps = 1.0e-6;

	dvector b0(n);
	dmatrix A0(n, n), R(n, n);
	ivector p(K);
	omvector M(K + 1);
	ovector b(K + 1);

	ktransform(M1, b1, p1, M, b, p); // the system is transformed so that all [p]=[-1,1]

	A0 = M[0]; // midpoint matrix of the new system
	b0 = b[0]; // midpoint vector of the new system
	R = A0;
	if (inv(R)) { // if the mipoint matrix is nonsingular, the computation continues
		ivector v(n), v1(n);
		dvector x0(n), t3i(n);
		dmatrix C0(n, K), C(n, K);
		omvector B(K), D(K), D0(K);
		ovector d(K), G(K);

		gausse(A0, b0, x0); // computing x0 - the midpoint solution
		for (int k = 0; k < K; k++) {
			d[k] = R * (b[k + 1] - M[k + 1] * x0); // d^k = R*(b^k-A^k*x0)
			B[k] = R * M[k + 1]; // B^k=R*A^k, tez chyba OK
			for (int j = 0; j < n; j++) {
				C0(j, k) = d[k][j];
			}
		}
		// v^(0)=0
		v = C0 * p; // v^(1)
		for (int k = 0; k < K; k++) {
			D[k] = B[k] * C0; // wyglada na to, ze D tez jest dobrze obliczone
		}

		// v^(2)=t1+t2*p+t3	
		t1 = dvector(n);
		t3 = dvector(n);
		t2 = dmatrix(n, K);

		t1.fill_in(0.0);
		t2.fill_in(0.0);
		t3.fill_in(0.0);

		for (int i = 0; i < n; i++) { // n is the size of v
			real s = 0.0;
			for (int k = 0; k < K; k++) {
				for (int j = k + 1; j < K; j++) {
					// perhaps this can be estimated better
					// strony 243, 244 z artykulu Koleva
					s += fabs(D[k](i, j) + D[j](i, k));
				}
				t1[i] = t1[i] - D[k](i, k) / 2.0;
				t2(i, k) = C0(i, k); // this is the L operator, which is multiplied by p
				t3[i] = t3[i] + fabs(D[k](i, k) / 2.0);
			}
			t3[i] = t3[i] + s;
		}

		////////////////////////////////////////////////////////////////////////
		// ITERATION
		// computing v^(k), k>=3
		//
		reduce(t1, t2, t3, v);
		do {
			niter++;
			v1 = v;
			for (int k = 0; k < K; k++) {
				G[k] = B[k] * t1;
				for (int j = 0; j < n; j++) {
					C(j, k) = C0(j, k) - G[k][j];
				}
				D0[k] = B[k] * t2;
			}
			t1.fill_in(0.0);
			t3i = t3;
			t3.fill_in(0.0);
			for (int i = 0; i < n; i++) { // n is the size of v
				real s = 0.0;
				real s1 = 0.0;
				for (int k = 0; k < K; k++) {
					for (int j = 0; j < n; j++) {
						s1 += fabs(B[k](i, j)) * t3i[j];
					}
					for (int j = k + 1; j < K; j++) {
						s += fabs(D0[k](i, j) + D0[j](i, k));
					}
					t1[i] = t1[i] - D0[k](i, k) / 2.0;
					t2(i, k) = C(i, k); // this is the L operator, which is multiplied by p
					t3[i] = t3[i] + fabs(D0[k](i, k) / 2.0);
				}
				t3[i] = t3[i] + s + s1;
				// the p-solution is in the form: t1 + t2*p + [t3]
				// wiec wystarczy zrobic optymalizacje t2*p
			}
			reduce(t1, t2, t3, v);
		} while (stopCrit(v, v1) > eps);

#ifdef INNER_KOL
		///// ************************************
		///// COMPUTING INNER ESTIMATION
		dmatrix mt2(n, K);
		ivector a(n), innsol(n);
		dvector au(n), al(n), pv(K), innu(n), innl(n);
		mag(t2, mt2);
		pv.fill_in(1.0);
		pv = mt2 * pv;
		a = x0 + t1 + t3 * interval(-1.0, 1.0);
		sup(a, au);
		inf(a, al);
		innl = au - pv;
		innu = al + pv;
		// final result for inner solution
		for (int i = 0; i < n; i++) {
			innsol[i] = interval(innl[i], innu[i]);
		}
		std::cout << "Inner:" << std::endl << innsol << std::endl;
		///// ************************************
#endif
		reduce(t1, t2, t3, v);
		for (int i = 0; i < n; i++) { // final result for outer solution
			w[i] = v[i] + interval(x0[i]);
		}
	}
	//std::cout << "Iterations: " << niter << std::endl;
}