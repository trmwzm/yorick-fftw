plug_in, "yfftw";

/* yfftw wrappers */
require, "yfftw.i";

scratch= save(scratch, tmp);
tmp= save(_eval);

func offt (base, void, ljdir=, rjdir=, dims=, usefftw=, inplace=, fcplx=)
/* DOCUMENT obj= offt(ljdir=,rjdir=,dims=,usefftw=,inplace=,fcplx=)
 *
 *   Return a new FFT operator object and initialize workspace
 *   to compute FFT using Yorick's builtin Swarztrauber FFT, or
 *   using the yfftw FFTW plugin -- the "fastest FFT in the West."
 *
 *   Keyword DIMS can be used to pre-specify the dimension list of the
 *   arrays to be transformed.  If left unspecified, the dimension
 *   list is initialized with the workspace the first time the
 *   transform operator is applied.  In all cases, an
 *   operator can only be used on arrays with same dimension lists,
 *   and the same target dimension (same as FFTW, and more restrictive
 *   than builtin fft for which workspace in not direction or array- 
 *   geometry specific.)  
 *   
 *   with USEFFTW==1 and  FCPLX==1 all arrays with a leading dimension of 
 *   length 2 will be taken to be of type complex{float re;float im}
 *
 *   with USEFFTW==1, the workspace has to be manually cleaned up to
 *   free up memory (a shortcoming of yfftw plug_in):  OBJ, reset=1;
 *
 *   Examples:
 *
 *       d= [2,12,32]; x=random_n(d)+1i*random_n(d);
 *       fto= offt(usefftw=1);
 *       write,avg(abs(fft(x)-fto(x))),format="%15.8g\n";
 *       fto, reset=1;
 *
 *       fto= offt(dims=d,ljdir=1,usefftw=1,inplace=1)
 *       fto, x;
 *       fto, reset=1;
 *
 *   NOTE: don't understand why usefftw cannot be in use statement
 *   at line 1 of _eval func.
 * SEE ALSO: fftw, fft. 
 */
{
  ob= base(:);

  save, ob, ljdir=ljdir, rjdir=rjdir, dims=dims,\
              usefftw=usefftw, inplace=inplace, fcplx=fcplx;

  if (fcplx==1 && usefftw!=1)
    error,"fcplx=1 only possible with fftw";

  if (!is_void(dims)) {
    if (usefftw==1) {
      if (fcplx==1 && dims(2)!=2)
        error,"fcplx=1 requires leading length-2 dimension";
      if (inplace==1) {
        local p;
        fftw, dims, ljdir, rjdir, p, fcplx=fcplx;
        save, ob, setup=p;
        p= [];
      } else {
        save, ob, setup= fftw(dims, ljdir, rjdir, fcplx=fcplx);
      }
    } else {
      save, ob, setup= fft_setup(dims, ljdir, rjdir);
    }
  } else {
    save, ob, setup=[];
  }

  return closure(ob, _eval);
}

func _eval (&x, ljdir0, rjdir0, reset=, fcplx=)
{
  use, ljdir, rjdir, dims, setup, inplace;

  if (!is_void(use(fcplx)) && !is_void(fcplx) && fcplx!=use(fcplx))
    error,"fcplx offt(eval) keyword inconsitent with init";
  if (is_void(fcplx)) restore, use, fcplx;
  if (is_void(use(fcplx))) save, use, fcplx;

  if (reset==1) {
    if (use(usefftw)==1)
      fftw_clean, setup, fcplx=fcplx;
    else
      setup= [];
    return;
  }

  dims0= dimsof(x);
  if (!is_void(dims) && anyof(dims0!=dims))
    error,"init with diff. dims";
  if (is_void(dims)) 
    dims= dims0;

  if (!is_void(ljdir0) && !is_void(ljdir) && anyof(ljdir0!=ljdir)) 
    error,"init with diff. ljdir";
  if (is_void(ljdir))
    ljdir= ljdir0;

  if (!is_void(rjdir0) && !is_void(rjdir) && anyof(rjdir0!=rjdir)) 
    error,"init with diff. rjdir";
  if (is_void(rjdir))
    rjdir= rjdir0;

  if (!is_void(inplace)) {
    if (inplace!=am_subroutine()) 
      error,"setup/plan in/out-place, called as out/in place.";
  } else {
    inplace= am_subroutine();
  }

  if (is_void(setup)) {
    if (use(usefftw)==1) {
      if (inplace==1) {
        local p;
        fftw, dims, ljdir, rjdir, p, fcplx=fcplx;      
        setup= p;      
      } else {
        setup= fftw(dims, ljdir, rjdir, fcplx=fcplx);      
      }
    } else {
      setup= fft_setup(dims, ljdir, rjdir);
    }
  }

  if (inplace==1) {
    if (use(usefftw)==1) {
      fftw, x, ljdir, rjdir, setup, fcplx=fcplx, keep=1;
    } else {
      fft_inplace, x, ljdir, rjdir, setup=setup;
    }
  } else {
    if (use(usefftw)==1) {
      return fftw(x, ljdir, rjdir, setup, fcplx=fcplx, keep=1);
    } else {
      return fft(x, ljdir, rjdir, setup=setup);
    }
  }
}
offt= closure(offt, restore(tmp));
restore, scratch;

