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

scratch= save(scratch, tmp);
tmp= save(eval,add,reset,_cfg,_get);
func fftwrap (base, dims, ljdir, rjdir, usefftw=, inplace=, fcplx=)
/* DOCUMENT  oft= fftwrap();
             oft= fftwrap(dims,ljdir,rjdir,usefftw=,inplace=,fcplx=)
   create an fft operator with or without a first initialization.
   oft, add, dims,ljdir,rjdir,usefftw=,inplace=,fcplx=;
   oft, reset, dims,ljdir,rjdir,usefftw=,inplace=,fcplx=,index=;
   oft, eval,x,ljdir,rjdir,reset=,fcplx=;
   or
   y= oft(eval,x,ljdir,rjdir,reset=,fcplx=);


 * SEE ALSO: fftw, fft.
 */
{
  ob= base(:);

  cfgdb= save();

  if (!is_void(dims)) {
      cfg= ob(_cfg,dims,ljdir,rjdir,fcplx=fcplx,inplace=inplace,usefftw=usefftw);
      save, cfgdb, string(0), cfg;
  }

  save, ob, cfgdb;

  return ob;
}
func add (dims, ljdir, rjdir, usefftw=, inplace=, fcplx=)
{
    use_method,_get,dims,ljdir,rjdir,fcplx=fcplx,inplace=inplace,usefftw=usefftw;
}
func _cfg (dims,ljdir,rjdir,fcplx=,inplace=,usefftw=,verbose=)
{
  if (is_void(fcplx)) fcplx= 0;
  if (is_void(inplace)) inplace= 0;
  if (is_void(usefftw)) usefftw= 0;

  dirs= fft_dirs(dims(1),ljdir,rjdir);
  cfg= save(dims,dirs,fcplx,inplace,usefftw);

  if (fcplx==1 && usefftw!=1)
    error,"fcplx=1 only possible with fftw";

  if (usefftw==1) {
      if (fcplx==1 && dims(2)!=2)
          error,"fcplx=1 requires leading length-2 dimension";
      if (inplace==1) {
          local p;
          fftw, dims, ljdir, rjdir, p, fcplx=fcplx;
          save, cfg, setup=p;
          p= [];
      } else {
          save, cfg, setup= fftw(dims, ljdir, rjdir, fcplx=fcplx);
      }
      if (verbose)
          write,pr1(dims),pr1(dirs),fcplx,inplace, \
              format="FFTW plan: dims %s, dirs %s, fcplx %i, inplace %i\n";
  } else {
      save, cfg, setup= fft_setup(dims, ljdir, rjdir);
      if (verbose)
          write,pr1(dims),pr1(dirs),                    \
              format="FFT plan: dims %s, dirs %s\n";
  }

  return cfg;
}
func _get (dims,ljdir,rjdir,&i,fcplx=,inplace=,usefftw=,verbose=)
{
  use, cfgdb;

  if (is_void(fcplx)) fcplx= 0;
  if (is_void(inplace)) inplace= 0;
  if (is_void(usefftw)) usefftw= 0;

  dirs= fft_dirs(dims(1),ljdir,rjdir);
  wd= where(dirs);
  n= cfgdb(*);

  if (fcplx==1 && usefftw!=1)
    error,"fcplx=1 only possible with fftw";

  for (i=1;i<=n;i++) {
      o= cfgdb(noop(i));
      if ( usefftw==1 && o(usefftw)==1 && fcplx==o(fcplx) && \
           inplace==o(inplace) && dims(1)==o(dims,1) && \
           allof(dirs==o(dirs)) )
          return o;
      owd= where(o(dirs));
      if ( usefftw==0 && o(usefftw)==0 && numberof(owd)>=numberof(wd) && \
           allof((dims(wd+1)(-,)==o(dims,owd+1))(max,)) )
          return o;
  }
  cfg= use_method(_cfg,dims,ljdir,rjdir,fcplx=fcplx,inplace=inplace, \
                  usefftw=usefftw, verbose=verbose);
  save, cfgdb, string(0), cfg;
  i= cfgdb(*);

  return cfg;
}
func reset (dims,ljdir,rjdir,fcplx=,inplace=,usefftw=,index=)
{
    use, cfgdb;
  if (is_void(index))
      cfg= use_method(_get,dims,ljdir,rjdir,index,fcplx=fcplx, \
                      inplace=inplace,usefftw=usefftw);
  else
      cfg= cfgdb(noop(index));

  if (index>0) {
      if (cfg(usefftw)==1)
          fftw_clean, cfg(setup), fcplx=cfg(fcplx);
      cfgdb2= save();
      for (j=1;j<index;j++) save,cfgdb2,string(0),cfgdb(noop(j));
      for (j=index+1;j<=cfgdb(*);j++) save,cfgdb2,string(0),cfgdb(noop(j));
      cfgdb= cfgdb2;
  }
}
func eval (&x, ljdir, rjdir, &index, reset=, fcplx= ,usefftw=)
{
  inplace= am_subroutine();

  cfg= use_method(_get,dimsof(x),ljdir,rjdir,index,fcplx=fcplx, \
                  inplace=inplace,usefftw=usefftw,verbose=1);
  local setup;
  setup= cfg(setup);
  if (inplace==1) {
    if (usefftw==1) {
      fftw, x, ljdir, rjdir, setup, fcplx=fcplx, keep=1;
    } else {
      fft_inplace, x, ljdir, rjdir, setup=setup;
    }
  } else {
    if (usefftw==1) {
      return fftw(x, ljdir, rjdir, setup, fcplx=fcplx, keep=1);
    } else {
      return fft(x, ljdir, rjdir, setup=setup);
    }
  }

  if (reset==1)
      use_method, reset, index=index;
}
fftwrap= closure(fftwrap, restore(tmp));
restore, scratch;
