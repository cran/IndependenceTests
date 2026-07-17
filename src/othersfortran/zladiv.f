*> \brief \b ZLADIV performs complex division in real arithmetic, avoiding unnecessary overflow.
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at 
*            http://www.netlib.org/lapack/explore-html/ 
*
*> \htmlonly
*> Download ZLADIV + dependencies 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zladiv.f"> 
*> [TGZ]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zladiv.f"> 
*> [ZIP]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zladiv.f"> 
*> [TXT]</a>
*> \endhtmlonly 
*
*  Definition:
*  ===========
*
*       DOUBLE COMPLEX     FUNCTION ZLADIV( X, Y )
* 
*       .. Scalar Arguments ..
*       DOUBLE COMPLEX         X, Y
*       ..
*  
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> ZLADIV := X / Y, where X and Y are complex.  The computation of X / Y
*> will not overflow on an intermediary step unless the results
*> overflows.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] X
*> \verbatim
*>          X is DOUBLE COMPLEX
*> \endverbatim
*>
*> \param[in] Y
*> \verbatim
*>          Y is DOUBLE COMPLEX
*>          The complex scalars X and Y.
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee 
*> \author Univ. of California Berkeley 
*> \author Univ. of Colorado Denver 
*> \author NAG Ltd. 
*
*> \date September 2012
*
*> \ingroup complex16OTHERauxiliary
*
*  =====================================================================
      COMPLEX(KIND(0D0))     FUNCTION ZLADIV( X, Y )
*
*  -- LAPACK auxiliary routine (version 3.4.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     September 2012
*
*     .. Scalar Arguments ..
      COMPLEX(KIND(0D0))         X, Y
*     ..
*
*  =====================================================================
*
*     .. Local Scalars ..
      DOUBLE PRECISION   ZI, ZR
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLADIV
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE, CMPLX, AIMAG
*     ..
*     .. Executable Statements ..
*
      CALL DLADIV( DBLE( X ), AIMAG( X ), DBLE( Y ), AIMAG( Y ), 
     $             ZR, ZI )
      ZLADIV = CMPLX( ZR, ZI, KIND(0D0))
*
      RETURN
*
*     End of ZLADIV
*
      END
