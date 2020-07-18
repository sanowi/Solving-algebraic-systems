#include <iostream>
#include "op_asm.h"

REAL add(REAL x, REAL y) 
{
	//REAL z;
	short cw;
	short rw;

	__asm {

		fld			x
		fld			y
		fnstcw		cw //store control word 
		mov			ax, cw
		and			ax, 0x03ff //0F3FFh // CRoundNear
		or			ax, 0x0800 //00800h // CRoundUp
		mov			rw, ax
		fldcw		rw
		fadd
		fldcw		cw //load stored control word
	}
	//std::cout << "add: z = " << z << std::endl;
}

REAL sub(REAL x, REAL y) 
{
	short cw;
	short rw;

	__asm {

		fld			x
		fld			y
		fnstcw		cw
		mov			ax, cw
		and			ax, 0x03ff //0F3FFh // CRoundNear
		or			ax, 0x0800 //00800h // CRoundUp
		mov			rw, ax
		fldcw		rw
		fsub
		fldcw		cw
	}
}

REAL mult(REAL x, REAL y) 
{
	short cw;
	short rw;
	REAL z, z2;

	__asm {

		fld			x
		fnstcw		cw
		mov			ax, cw
		and			ax, 0x03ff //0F3FFh // CRoundNear
		or			ax, 0x0800 //00800h // CRoundUp
		mov			rw, ax
		fldcw		rw
		fmul		y
		fstp		z
		fld			x
		fmul		y
		fstp		z2
		fldcw		cw
	}
	printf("aaa %4.20lf bbb %4.20lf\n", z, z2);
	return z;
}

REAL div(REAL x, REAL y) 
{
	short cw;
	short rw;

	__asm {

		fld			x
		fld			y
		fnstcw		cw
		mov			ax, cw
		and			ax, 0x03ff //0F3FFh // CRoundNear
		or			ax, 0x0800 //00800h // CRoundUp
		mov			rw, ax
		fldcw		rw
		fdiv
		fldcw		cw
	}
}