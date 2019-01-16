#include <iostream>
#include <fstream>
#include <conio.h>
#include <cmath>
#include <ctime>
#include "config.h"
#include "inverse.h"
#include "gausselim.h"
#include "zVector.h"
#include "mism.h"
#include "TrussStructure.h"
#include "DervSol.h"
#include "SensAnalysis.h"
#include "ceic.h"

//ivector SAnal( const ivector& pv, 
//			   const omatrix& oA, 
//			   const ovector& ob )
//{
//	int			i, 
//				j, 
//				n = oA.num_rows(), 
//				p_n = pv.size();
//	imatrix		A(n,n);
//	ivector		b(n), 
//				x1(n), 
//				x2(n), 
//				Wynik(n), 
//				pvint(p_n);
//	dvector		xsr(n), 
//				xx(n), 
//				bsr(n), 
//				bnew(n);
//	dmatrix		Asr(n,n);
//	dmatrix		dervoA(n, n);
//	ovector		dervsmin(n), 
//				dervsmax(n);
//	dvector		pvreal(p_n);
//
//	pvreal.fill_in(0.0);
//	pvreal[ 0 ] = 1.0;
//	for ( i = 0; i < n; i++ )
//	{
//		dervsmin[ i ] = pvreal;
//		dervsmax[ i ] = pvreal;
//	}
//
//	A = m_from_dep( pv, oA );
//	b = v_from_dep( pv, ob );
//
//	mid( A, Asr );
//	mid( b, bsr );
//	inv(Asr);
//	xsr = Asr * bsr; // midpoint solution
//	for ( i = 0; i < p_n; i++ )
//	{
//		dervoA = PartialDerivative( oA, i ) * (-1.0);
//		bnew = dervoA * xsr;
//		xx = Asr * bnew; // pochodne dx/dp_j
//		for(j = 0; j < n; j++)
//		{
//			if (xx[j] >= 0)
//			{
//				dervsmin[j][i] = pv[i].inf();
//				dervsmax[j][i] = pv[i].sup();
//				continue;
//			}
//			if (xx[j] <= 0)
//			{
//				dervsmin[j][i] = pv[i].sup();
//				dervsmax[j][i] = pv[i].inf();
//			}
//		}
//	}
//
//	// obliczanie ekstremow
//	int licznik = 0;
//	for(i = 0; i < n; i++)
//	{
//		for(j = i + 1; j < n; j++)
//		{
//			if ((dervsmin[ i ] - dervsmax[ j ]) == 0)
//			{
//				licznik++;
//			}
//		}
//	}
//	for(i = 0; i < n; i++)
//	{
//		for(j = i + 1; j < n; j++)
//		{
//			if ((dervsmin[i] - dervsmin[ j ]) == 0)
//			{
//				licznik++;
//			}
//		}
//	}
//	for(i = 0; i < n; i++)
//	{
//		for(j = i + 1; j < n; j++)
//		{
//			if ((dervsmax[ i ] - dervsmax[ j ]) == 0)
//			{
//				licznik++;
//			}
//		}
//	}
//	for(i = 0; i < n; i++)
//	{
//		to_int(dervsmin[ i ], pvint);
//        x1 = ISM(pvint, oA, ob);
//        to_int(dervsmax[ i ], pvint);
//        x2 = ISM(pvint, oA, ob);
//        Wynik[ i ] = interval(x1[ i ].inf(), x2[ i ].sup());
//	}
//	return Wynik;
//}

ivector SAnalF( const ivector& pv, const ivector& qv, const omatrix& oA, const ovector& ob )
{
	int		i, j,
			n = oA.num_rows(),
			np = pv.size(), 
			nq = qv.size();
	imatrix A(n,n);
	ivector b(n), x1(n), 
			x2(n), 
			Wynik(n), 
			pvint(np), 
			qvint(nq);
	dvector xsr(n), 
			xx(n), 
			bsr(n), 
			bnew(n), 
			dervob(n);
	dmatrix Asr(n,n);
	dmatrix dervoA(n, n);
	ovector pdervsmin(n), 
			pdervsmax(n);
	ovector qdervsmin(n), 
			qdervsmax(n);
	dvector pvreal(np), 
			qvreal(nq);

	pvreal.fill_in(0.0);
	pvreal[0]=1.0;
	qvreal.fill_in(0.0);
	qvreal[0]=1.0;
	for(i = 0; i < n; i++)
	{
		pdervsmin[i] = pvreal;
		pdervsmax[i] = pvreal;
		qdervsmin[i] = qvreal;
		qdervsmax[i] = qvreal;
	}

	dvector pmid(np), qmid(nq);
	mid(pv,pmid);
	mid(qv,qmid);
	Asr = m_from_dep(pmid, oA);
	bsr = v_from_dep(qmid, ob);

	inv(Asr); // inverse of a midpoint matrix
	xsr = Asr * bsr;

	// sprawdzanie pochodnych po p
	for(i = 0; i < np; i++)
	{
		dervoA = PartialDerivative(oA, i) * (-1.0);
		bnew = dervoA * xsr;
		xx = Asr * bnew; // pochodne dx/dp_j
		for(j = 0; j < n; j++)
		{
			if (xx[j] > 0)
			{
				pdervsmin[j][i] = pv[i].inf();
				pdervsmax[j][i] = pv[i].sup();
				continue;
			}
			if (xx[j] < 0)
			{
				pdervsmin[j][i] = pv[i].sup();
				pdervsmax[j][i] = pv[i].inf();
				continue;
			}
			if (xx[j] == 0)
			{
				pdervsmin[j][i] = pv[i].mid();
				pdervsmax[j][i] = pv[i].mid();
				continue;
			}
		}
	}
	//std::cout << pdervsmin << std::endl << std::endl;
	//std::cout << pdervsmax << std::endl << std::endl;

	// sprawdzanie pochodnych po q
	for(i = 0; i < nq; i++)
	{
		dervob = PartialDerivative(ob, i);
		bnew = dervob;
		xx = Asr * bnew; // pochodne dx/dq_j
		for(j = 0; j < n; j++)
		{
			if (xx[j] > 0)
			{
				qdervsmin[j][i] = qv[i].inf();
				qdervsmax[j][i] = qv[i].sup();
				continue;
			}
			if (xx[j] < 0)
			{
				qdervsmin[j][i] = qv[i].sup();
				qdervsmax[j][i] = qv[i].inf();
				continue;
			}
			if (xx[j] == 0)
			{
				qdervsmin[j][i] = qv[i].mid();
				qdervsmax[j][i] = qv[i].mid();
				continue;
			}
		}
	}
	for(i = 0; i < n; i++)
	{
		to_int( pdervsmin[i], pvint );
		to_int( qdervsmin[i], qvint );
        ISM(pvint, qvint, oA, ob, x1);
        to_int( pdervsmax[i], pvint );
		to_int( qdervsmax[i], qvint );
        ISM(pvint, qvint, oA, ob, x2);
        Wynik[i] = interval(x1[i].inf(), x2[i].sup());
	}
	return Wynik;
}

void SAnalMonotTest(const ivector& pv, const ivector& qv, const omatrix& oA, const ovector& ob)
{
	int n = oA.num_rows(), np = pv.size();
	dmatrix dervoA(n, n), A(n, n);
	dvector el(n+1), pmin(np), pmax(np), xmin(n), xmax(n);
	ivector x, qnew(1+n), b(n);
	ivector bnew, xw(n), w(n);
	ovector obnew(n);
	oivector xm(np);
	CombinatorialApproach(oA,ob,pv,qv,x);
	qnew[0] = 1.0;
	for(int i = 1; i < n+1; i++)
	{
		qnew[i] = x[i - 1];
	}
	el[0] = 0.0;
	for(int i = 1; i < np; i++)
	{
		dervoA = PartialDerivative(oA, i) * (-1.0);
		for(int j = 0; j < n; j++)
		{
			for(int k = 1; k < n + 1; k++)
			{
				el[k] = dervoA(j, k - 1);
			}
			obnew[j] = el;
		}
		CombinatorialApproach(oA,obnew,pv,qnew,x);
		xm[i] = x;
	}
	std::cout << xm << std::endl;
	std::cout << "End of 1 stage" << std::endl;
	xmin.fill_in(100000000.0);
	xmax.fill_in(-100000000.0);
	pmin[0] = 1.0;
	pmax[0] = 1.0;
	for(int i = 0; i < n; i++)
	{
		for(int j = 1; j < np; j++)
		{
			if (xm[j][i].inf() * xm[j][i].sup() > 0)
			{
				if (xm[j][i].inf() > 0)
				{
					pmin[j] = pv[j].inf();
					pmax[j] = pv[j].sup();
					continue;
				}
				if (xm[j][i].sup() < 0)
				{
					pmin[j] = pv[j].sup();
					pmax[j] = pv[j].inf();
					continue;
				}
			}
			else
			{
				pmin[j] = pv[j].mid();
				pmax[j] = pv[j].mid();
			}
		}
		b = v_from_dep(qv,ob);
		A = m_from_dep(pmin,oA);
		inv(A);
		w = A * b;
		xmin[i] = ISLIB_MIN(xmin[i],w[i].inf());
		A = m_from_dep(pmax,oA);
		inv(A);
		w = A * b;
		xmax[i] = ISLIB_MAX(xmax[i],w[i].sup());
	}
	for(int i = 0; i < n; i++)
	{
		xw[i]=interval(xmin[i],xmax[i]);
	}
	std::cout << xw << std::endl;
}

void readSensResults(dmatrix&);

ivector SAnalFW(const ivector& pv, const ivector& qv, const omatrix& oA, const ovector& ob)
{
	int		i, j,
			n = oA.num_rows(),
			np = pv.size(), 
			nq = qv.size();
	imatrix A(n,n);
	ivector b(n), x1(n), 
			x2(n), 
			Wynik(n), 
			pvint(np), 
			qvint(nq);
	dvector xsr(n), 
			xx(n), 
			bsr(n), 
			bnew(n), 
			dervob(n);
	dmatrix Asr(n,n);
	dmatrix dervoA(n, n);
	ovector pdervsmin(n), 
			pdervsmax(n);
	ovector qdervsmin(n), 
			qdervsmax(n);
	dvector pvreal(np), 
			qvreal(nq);

	std::ofstream fout;
	fout.open("sensitivity.txt",std::ofstream::out);

	pvreal.fill_in(0.0);
	pvreal[0]=1.0;
	qvreal.fill_in(0.0);
	qvreal[0]=1.0;
	for(i = 0; i < n; i++)
	{
		pdervsmin[i] = pvreal;
		pdervsmax[i] = pvreal;
		qdervsmin[i] = qvreal;
		qdervsmax[i] = qvreal;
	}

	dvector pmid(np), qmid(nq);
	mid(pv,pmid);
	mid(qv,qmid);
	Asr = m_from_dep(pmid, oA);
	bsr = v_from_dep(qmid, ob);

	inv(Asr); // inverse of a midpoint matrix
	xsr = Asr * bsr;

	// sprawdzanie pochodnych po p
	for(i = 0; i < np; i++)
	{
		dervoA = PartialDerivative(oA, i) * (-1.0);
		bnew = dervoA * xsr;
		xx = Asr * bnew; // pochodne dx/dp_j
		for(j = 0; j < n; j++)
		{
			if (xx[j] > 0)
			{
				pdervsmin[j][i] = pv[i].inf();
				pdervsmax[j][i] = pv[i].sup();
				continue;
			}
			if (xx[j] < 0)
			{
				pdervsmin[j][i] = pv[i].sup();
				pdervsmax[j][i] = pv[i].inf();
				continue;
			}
			if (xx[j] == 0)
			{
				pdervsmin[j][i] = pv[i].mid();
				pdervsmax[j][i] = pv[i].mid();
				continue;
			}
		}
	}

	// sprawdzanie pochodnych po q
	for(i = 0; i < nq; i++)
	{
		dervob = PartialDerivative(ob, i);
		bnew = dervob;
		xx = Asr * bnew; // pochodne dx/dq_j
		for(j = 0; j < n; j++)
		{
			if (xx[j] > 0)
			{
				qdervsmin[j][i] = qv[i].inf();
				qdervsmax[j][i] = qv[i].sup();
				continue;
			}
			if (xx[j] < 0)
			{
				qdervsmin[j][i] = qv[i].sup();
				qdervsmax[j][i] = qv[i].inf();
				continue;
			}
			if (xx[j] == 0)
			{
				qdervsmin[j][i] = qv[i].mid();
				qdervsmax[j][i] = qv[i].mid();
				continue;
			}
		}
	}

	// obliczanie ekstremow
	/*int licznik = 0;
	for(i = 0; i < n; i++)
	{
		for(j = i + 1; j < n; j++)
		{
			if (pdervsmin[i] - pdervsmax[j] == 0)
			{
				licznik++;
			}
		}
	}*/
	std::cout.precision(12);
	for(i = 0; i < n; i++)
	{
		to_int( pdervsmin[i], pvint );
		to_int( qdervsmin[i], qvint );
        ISM(pvint, qvint, oA, ob, x1);
        to_int( pdervsmax[i], pvint );
		to_int( qdervsmax[i], qvint );
        ISM(pvint, qvint, oA, ob, x2);
        Wynik[i] = interval(x1[i].inf(), x2[i].sup());
	}
	std::cout << pv << std::endl;
	fout << n << " " << np + nq - 2 << std::endl;
	fout.precision(12);
	for(i = 0; i < n; i++)
	{
		for(j = 1; j < np; j++)
		{
			fout << pdervsmin[i][j] << " ";
		}
		for(j = 1; j < nq; j++)
		{
			fout << qdervsmin[i][j] << " ";
		}
		fout << std::endl;
	}
	for(i = 0; i < n; i++)
	{
		for(j = 1; j < np; j++)
		{
			fout << pdervsmax[i][j] << " ";
		}
		for(j = 1; j < nq; j++)
		{
			fout << qdervsmax[i][j] << " ";
		}
		fout << std::endl;
	}
	fout.close();
	dmatrix w;
	readSensResults(w);
	return Wynik;
}

void readSensResults(dmatrix& w)
{
	int n = 0, nn = 0;
	std::ifstream fin;
	fin.open("sensitivity.txt",std::ofstream::in);

	if (!fin.eof())
	{
		fin >> n;
		fin >> nn;
	}
	w = dmatrix(2*n,nn);
	char line[200];
	int i=0, j=0;
	while (!fin.eof())
	{
		fin >> line;
		if (line[0] == '\0') continue;
		w(i,j++)=atof(line);
		if (j==nn) 
		{
			i++;
			j = 0;
		}
	}
	for(i = 0; i < 2*n; i++)
	{
		for(j = 0; j < nn; j++)
		{
			std::cout << w(i,j) << " ";
		}
		std::cout << std::endl;
	}
	fin.close();
}
