plug_in, "yfftw";

func rollarr(x, ljoff, rjoff, fcplx=)
/* DOCUMENT rollarr (x, ljoff, rjoff, fcplx=)
         or rollarr, x, ljoff, rjoff
         or rollarr(x)
         or rollarr, x

     "rolls" selected dimensions of the array X.  The roll offsets
     LJOFF and RJOFF (both optional) work in the same fashion as the
     LJDIR and RJDIR arguments to the fft function:

        A scalar LJDIR (and nil RJDIR) rolls all dimensions of X by
        the specified offset.
        Otherwise, the elements of the LJDIR vector [ljoff1, ljoff2, ...]
        are used as the roll offsets for the first, second, etc.
        dimensions of X.
        Similarly, the elements of the RJDIR vector [..., rjoff1, rjoff0]
        are matched to the final dimensions of X, so the next to last
        dimension is rolled by rjoff1 and the last dimension by rjoff0.

        As a special case (mostly for use with the fft function), if
        both LJDIR and RJDIR are nil, every dimension is rolled by
        half of its length.  Thus,
           roll(x)
        it equivalent to
           roll(x, dimsof(x)(2:0)/2)

     The result of the roll function is complex if X is complex, otherwise
     double (i.e.- any other array type is promoted to type double).  If
     roll is invoked as a subroutine, the operation is performed in place.

   SEE ALSO: fft
 */
{
  dims= dimsof(x);
  /* get array to be transformed */
  if (fcplx==1) {
    dims=_(dims(1)-1,dims(3:));
    stds1= 2;
    if (!am_subroutine()) x= x;
  } else if (structof(x)==complex) {
    stds1= 2;
    if (!am_subroutine()) x= x;  /* copy to avoid clobbering */
  } else {
    stds1= 1;
    if (!am_subroutine()) x= x;  /* copy to avoid clobbering */
  }
  ndims= dims(1);
  for (n=1,i=1;i<=ndims;i++) n*= dims(i+1);
  if (ndims<1 || n<2) return x;
  dims= dims(2:0);

  /* get offset arguments */
  if (is_void(ljoff) && is_void(rjoff)) ljoff= dims/2;
  else ljoff= fft_dirs(ndims, ljoff, rjoff);
  ljoff%= dims;
  if (noneof(ljoff)) return x;

  /* get strides */
  stds= array(stds1, ndims);   /* note extra factor of 2 for complex */
  for (i=1 ; i<ndims ; i++) stds(i+1)= stds(i)*dims(i);
  tops= (stds1*n)/(stds*dims);

  /* perform roll operation */
  px= &x;   /* avoids type conversion problems */
  if (structof(x)==complex||structof(x)==double) {
    ws= array(0.0, 2*max(dims(where(ljoff))));
    for (i=1 ; i<=ndims ; i++)
      _rollarr_d, px, ljoff(i), stds(i), dims(i), tops(i), ws;
  } else if (structof(x)==float) {
    ws= array(0.0f, max(dims(where(ljoff))));
    for (i=1 ; i<=ndims ; i++)
      _rollarr_f, px, ljoff(i), stds(i), dims(i), tops(i), ws;
  } else {
    error,"unknown type";
  }
  return x;
}

extern _rollarr_d;
/* PROTOTYPE
   void rollarr_d(pointer a, long ioff, long istd, long n, long n2,
              double array ws)
 */
extern _rollarr_f;
/* PROTOTYPE
   void rollarr_f(pointer a, long ioff, long istd, long n, long n2,
              float array ws)
 */

