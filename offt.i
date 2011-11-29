plug_in, "yfftw";

/* yfftw wrappers */
require, "yfftw.i";

scratch= save(scratch, tmp);
tmp= save(_eval);

func offt (base, ljdir=, rjdir=, dims=, usefftw=, inplace=)
/* DOCUMENT obj= offt(ljdir=,rjdir=,dims=,usefftw=,inplace=)
 *
 *   Return a new FFT operator object and initialize workspace
 *   to compute FFT using Yorick's builtin Swarztrauber FFT, or
 *   using the yfftw FFTW plugin-- the "fastest FFT in the West."
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
 *   with usefftw==1, all arrays with a leading dimension of 
 *   length 2 will be taken to be of type complex{float re;float im}
 *
 *   with usefftw==1, the workspace has to be manually cleaned up to
 *   free up memory (a shortcoming of yfftw plug_in):  OBJ, reset=1;
 *
 *   Examples:
 *
 *       OBJ = offt(ljdir=, rjdir=, usefftw=1)
 *
 *       OBJ(x)
 *
 * SEE ALSO: fftw, fft. 
 */
{
  ob= base(:);

  save, ob, ljdir=ljdir, rjdir=rjdir, dims=dims,
              usefftw=usefftw, inplace=inplace;

  if (!is_void(dims)) {
    if (usefftw) {
      if (inplace) {
        local p;
        fftw, dims, ljdir, rjdir, p;
        save, ob, setup=p;
        p= [];
      } else {
        save, ob, setup= fftw(dims, ljdir, rjdir);
      }
    } else {
      save, ob, setup= fft_setup(dims, ljdir, rjdir);
    }
  } else {
    save, ob, setup=[];
  }

  return closure(ob, _eval);
}

func _eval (&x, ljdir0, rjdir0, reset=)
{
  use, setup;
  restore, use, ljdir, rjdir, dims, inplace, usefftw;
  
  if (reset==1) {
    if (usefftw)
      if (dims(2)==2)
        fftwf_clean, setup;
      else
        fftw_clean, setup;
    else
      setup= [];
    return;
  }

  dims0= dimsof(x);
  if (is_void(dims)) {dims= dims0; save, use, dims;}
  if (anyof(dims0!=dims)) error,"init with diff. dims";

  if (is_void(ljdir)) {ljdir= ljdir0; save, use, ljdir;}
  if (!is_void(ljdir0) && anyof(ljdir0!=ljdir)) 
    error,"init with diff. ljdir";

  if (is_void(rjdir)) {rjdir= rjdir0; save, use, rjdir;}
  if (!is_void(rjdir0) && anyof(rjdir0!=rjdir))
    error,"init with diff. rjdir";

  if (!is_void(inplace)) {
    if (inplace && !am_subroutine()) 
      error,"setup/plan was in-place, called as function.";
  } else {
    save, use, inplace=am_subroutine();
  }

  if (is_void(setup)) {
    if (usefftw) {
      if (inplace) {
        local p;
        fftw, dims, ljdir, rjdir, p;      
        setup= p;      
      } else {
        setup= fftw(dims, ljdir, rjdir);      
      }
    } else {
      setup= fft_setup(dims, ljdir, rjdir);
    }
  }

  if (inplace) {
    if (usefftw) {
      fftw, x, ljdir, rjdir, setup;
    } else {
      fft_inplace, x, ljdir, rjdir, setup=setup;
    }
  } else {
    if (usefftw) {
      return fftw(x, ljdir, rjdir, setup);
    } else {
      return fft(x, ljdir, rjdir, setup=setup);
    }
  }
}
offt= closure(offt, restore(tmp));
restore, scratch;

