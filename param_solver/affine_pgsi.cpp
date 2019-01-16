#include <iostream>
#include <fstream>
#include <conio.h>
#include <math.h>
#include <time.h>
#include "config.h"
#include "inverse.h"
#include "aaf.h"
#include "aafr.h"
#include "PGSI.h"

#define PGSIeps 0.1e-12

bool pgsi( const aafmatrix& A, const aafvector& b, ivector& wek )
/* Preconditioned Gauss-Seidel method - OK */
{
	int			n= b.size();
	AAF			s;
	REAL		ss;
	interval	tmpwek, ttmp;
	dvector		bmid(n);
	dvector		x0(n);
	aafmatrix	Rf(n,n);
	dmatrix		R(n, n);
	aafvector	tmp(n), bb(n);
	aafmatrix	AA(n,n);
    //-------------------------
	std::cout << "aaf version" << std::endl;
	mid( A,R ); // R - inverse of the midpoint matrix = preconditioner
	inv( R );
	for ( int i=0;i<n;i++ )
	{
		for ( int j=0;j<n;j++ )
		{
			Rf(i,j)=AAF( R(i,j) );
		}
	}
	mid( b,bmid );
	AA=Rf*A; // preconditioning
	bb=Rf*b; // preconditioning
	//firstSolution(R, AA, bb, wek); // computing initial solution
	//wek=R*bmid;
	for ( int i=0; i<n; i++ )
	{
		tmp[i]=AAF( wek[i] );
	}
	int iter=0;
	do
	{
		iter++;
		for ( int i=0; i<n; i++ ) 
		{
			s=0.0;
			for ( int j=0; j<n; j++ ) 
			{
				if (i!=j) 
				{
					s=s+tmp[j]*AA(i,j);// tu na pewno jest roznica z metoda Jacobiego
				}
			}
			ttmp=(1.0/AA(i,i)*(bb[i]-s)).reduce();
			ttmp.intersects(tmp[i].reduce(),tmpwek);
			if ( tmpwek.is_empty() )
			{
				return false;
			}
			wek[i]=tmpwek;
			ss=0.0;
			for ( int k=0; k<n; k++ ) 
			{
				ss+=ISLIB_ABS((tmp[k]).reduce().inf()-wek[k].inf())+
					ISLIB_ABS((tmp[k]).reduce().sup()-wek[k].sup());
			}
			if ( ss<PGSIeps ) // approximate equality
			{
				std::cout << "Iter=" << iter << std::endl;
				return true;
			}
			tmp[i]=AAF( wek[i] );
		}   
	}
	while( true );
}

bool pgsi(const aafrmatrix& A, const aafrvector& b, ivector& wek)
/* Preconditioned Gauss-Seidel method - OK */
{
	int			n= b.size();
	AAFR		s;
	REAL		ss;
	interval	tmpwek, ttmp;
	dvector		bmid(n);
	dvector		x0(n);
	aafrmatrix	Rf(n,n);
	dmatrix		R(n, n);
	aafrvector	tmp(n), bb(n);
	aafrmatrix	AA(n,n);
    //-------------------------
	mid(A,R); // R - inverse of the midpoint matrix = preconditioner
	inv(R);
	for(int i=0;i<n;i++)
	{
		for(int j=0;j<n;j++)
		{
			Rf(i,j)=AAFR(R(i,j));
		}
	}
	mid(b,bmid);
	AA=Rf*A; // preconditioning
	bb=Rf*b; // preconditioning
	//firstSolution(R, AA, bb, wek); // computing initial solution
	//wek=R*bmid;
	for(int i=0; i<n; i++)
	{
		tmp[i]=AAFR(wek[i]);
	}
	int iter=0;
	do
	{
		iter++;
		for (int i=0; i<n; i++) 
		{
			s=0.0;
			for (int j=0; j<n; j++) 
			{
				if (i!=j) 
				{
					s=s+tmp[j]*AA(i,j);// tu na pewno jest roznica z metoda Jacobiego
				}
			}
			ttmp=(1.0/AA(i,i)*(bb[i]-s)).reduce();
			ttmp.intersects(tmp[i].reduce(),tmpwek);
			if (tmpwek.is_empty())
			{
				return false;
			}
			wek[i]=tmpwek;
			ss=0.0;
			for(int k=0; k<n; k++) 
			{
				ss+=ISLIB_ABS((tmp[k]).reduce().inf()-wek[k].inf())+
					ISLIB_ABS((tmp[k]).reduce().sup()-wek[k].sup());
			}
			if(ss<PGSIeps) // approximate equality
			{
				std::cout << "Iter=" << iter << std::endl;
				return true;
			}
			tmp[i]=AAFR(wek[i]);
		}   
	}
	while(true);
}
