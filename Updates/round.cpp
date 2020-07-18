// Setting rouding modes

#include <iostream>
#include <math.h>
#include <stdio.h>
#include "Round.h"

static unsigned short oldcw, newcw;
static unsigned short CRoundDown	= 0x0400; //0x0500; //% the second values set as well precision bit to extended double
static unsigned short CRoundUp		= 0x0800; //0x0900; //%
static unsigned short CRoundNear	= 0xF3FF; //0xF2FF; //%
static unsigned short CTrunc		= 0x0C00; //0x0D00; //%

#define RTMP
#undef RTMP

//https://www.cs.put.poznan.pl/adanilecki/inline_asm/doda.php#fldcw

//% The link to explanation how to force rounding in VS2013 (IA32 architecture must be set)
//% http://stackoverflow.com/questions/24483382/different-results-when-using-different-cl-exe-compiler-options
//% zmiana - zamiast eax probuje ax, jest to 16 bitowa czesc rejestru eax - niestety bez zmian

//% std::nextafter(x.inf(), -ISLIB_INFINITY)
//% std::nextafter(x.sup(), ISLIB_INFINITY)
//%printf("%.17lf\n", nextafter(3.0 / 7.0, -1.0 / 0.0));
//%printf("%.17lf\n", nextafter(3.0 / 7.0, 1.0 / 0.0));

void RoundDown()
{
#ifdef _MSC_VER // if Visual C++ compiler
    // get the current Control Word to retain all setting bits
    // not related to the rounding control (RC) bits
#ifndef RTMP

	__asm {	fstcw oldcw
    // to insure the storage instruction is completed           
    fwait
    mov ax, oldcw
    // clears only the RC (11 and 10) bits, leaving all other bits unchanged
    // not necessary here because both bits will be set
    and ax, CRoundNear
    // this will set both bits of the RC field to the rounding down mode
    // without affecting any of the other field's bits 
    or ax, CRoundDown
    // use the stack to store the modified Control Word in memory
       //"push %ax\n\t"
    mov newcw, ax //??? oldcw, eax // mov newcw, eax
    //__asm mov eax, _sp
	//__asm push eax, sp
    // load the modified Control Word
	fldcw newcw //oldcw // fldcw newcw
	}
#endif // !RTMP
#endif
} //% >>> RoundDown <<<

void RoundUp()
{
#ifdef _MSC_VER
    // get the current Control Word to retain all setting bits
    // not related to the rounding control (RC) bits
#ifndef RTMP
	__asm{ fstcw oldcw
    // to insure the storage instruction is completed           
    fwait
    mov ax, oldcw
    // clears only the RC bits, leaving all other bits unchanged
    // not necessary here because both bits will be set
    and ax, CRoundNear
    // this will set both bits of the RC field to the truncating mode
    // without affecting any of the other field's bits 
    or ax, CRoundUp
    mov newcw, ax
    // load the modified Control Word
    fldcw newcw	
	}
#endif // !RTMP
#endif
} //% >>> RoundUp <<<

void RoundNear()
{
#ifdef _MSC_VER
    // get the current Control Word to retain all setting bits
    // not related to the rounding control (RC) bits
	__asm{ fstcw oldcw
    // to insure the storage instruction is completed           
    fwait
    mov ax, oldcw
    // clears only the RC bits, leaving all other bits unchanged
    // not necessary here because both bits will be set
    and ax, CRoundNear
    mov newcw, ax
    // load the modified Control Word
    fldcw newcw }
#endif
} //% >>> RoundNear <<<

void Trunc()
{
#ifdef _MSC_VER
    // get the current Control Word to retain all setting bits
    // not related to the rounding control (RC) bits
	__asm{ fstcw oldcw
    // to insure the storage instruction is completed           
    fwait
    mov ax, oldcw
    // clears only the RC bits, leaving all other bits unchanged
    // not necessary here because both bits will be set
    and ax, CRoundNear
    // this will set both bits of the RC field to the truncating mode
    // without affecting any of the other field's bits 
    or ax, CTrunc
    mov newcw, ax
    // load the modified Control Word
	fldcw newcw }
#endif
} //% >>> Trunc <<<

void Trunc2()
{
#ifdef _MSC_VER
    // get the current Control Word to retain all setting bits
    // not related to the rounding control (RC) bits
	__asm{ fstcw oldcw
    // to insure the storage instruction is completed           
    fwait
    mov ax, oldcw
    // clears only the RC bits, leaving all other bits unchanged
    // not necessary here because both bits will be set
    and ax, CRoundNear
    // this will set both bits of the RC field to the truncating mode
    // without affecting any of the other field's bits 
    or ax, CTrunc
    mov newcw, ax
    // load the modified Control Word
    fldcw newcw }
#endif
} //% >>> Trunc2 <<<
