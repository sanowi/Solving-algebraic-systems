#include <iostream>
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
#include "SensAnal.h"

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

ivector SAnalF( const ivector& pv, 
			    const ivector& qv, 
				const omatrix& oA, 
				const ovector& ob )
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

	inv(Asr);
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
	for(i = 0; i < n; i++)
	{
		to_int( pdervsmin[ i ], pvint );
		to_int( qdervsmin[ i ], qvint );
        ISM(pvint, qvint, oA, ob, x1);
        to_int( pdervsmax[ i ], pvint );
		to_int( qdervsmax[ i ], qvint );
        ISM(pvint, qvint, oA, ob, x2);
        Wynik[i] = interval(x1[i].inf(), x2[i].sup());
	}
	return Wynik;
}
