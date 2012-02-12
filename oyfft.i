plug_in, "yfftw";

/* WRAPPERS FOR FFT AND FFTW */
require, "yfftw.i";

/* Yeti is required for this package. */
if (noneof(is_func(h_new) == [1,2])) {
  include, "yeti.i", 1;
}

func op_fft (nil, ljdir=, rjdir=, dims=, usefftw=, inplace=, fcplx=)
/* DOCUMENT obj = op_fft(nil,ljdir=,rjdir=,dims=,usefftw=,fcplx=)
 *
 *   Return a new new linear operator object which can be used 
 *   to compute FFT by means of Swarztrauber's FFT or FFTW, 
 *   the "fastest FFT in the West".
 *
 *   Keyword DIMS can be used to pre-specify the dimension list of the
 *   arrays to be transformed.  If left unspecified, the actual
 *   dimension list will be initialized the first time the linear
 *   operator is applied.  In any cases, a given operator can only be
 *   used onto arrays with same dimension lists.
 *
 *   Examples:
 *
 *       OBJ = op_fft(ljdir, rjdir, usefftw=1)
 *
 *       OBJ.nevals = number of FFT computed so far by OBJ
 *
 *       OBJ(x)
 *   TEST
 *       See opcheck below
 * SEE ALSO: fftw, fft. 
 */
{
  if (! is_void(nil)) error, "no non-keyword argument allowed";

  this = h_new(ljdir=ljdir,rjdir=rjdir,dims=dims,usefftw=usefftw,\
               inplace=inplace,fcplx=fcplx);

  h_evaluator, this, "_opfft_eval";
  
  if ((!is_void(ljdir)||!is_void(rjdir))&&!is_void(dims)) {
    if (usefftw) {
      if (inplace) {
        local p;
        fftw, dims,ljdir,rjdir,p,fcplx=fcplx;
        h_set, this, setup=p;
        p= [];
      } else {
        h_set, this, setup= fftw(dims,ljdir,rjdir,fcplx=fcplx);
      }
    } else {
      if (fcplx==1)
        error,"use fftw for float complex.";
      h_set, this, setup= fft_setup(dims, ljdir, rjdir);
    }
  } else {
    h_set, this, setup=[];
  }

  return this;
}

func _opfft_eval (this, x, ljdir, rjdir)
{
  dx= dimsof(x);
  if (!is_void(ljdir)) {
    if (!is_void(this.ljdir)) {
      if (anyof(ljdir!=this.ljdir)) error,"init with diff. ljdir";
    } else {
      h_set,this,ljdir=ljdir;
    }
  }
  if (!is_void(rjdir)) {
    if (!is_void(this.rjdir)) {
      if (anyof(rjdir!=this.rjdir)) error,"init with diff. rjdir";
    } else {
      h_set,this,rjdir=rjdir;
    }
  }
  if (!is_void(this.dims)) {
    if (anyof(dx!=this.dims)) error,"init with diff. dims";
  } else {
    h_set,this,dims=dx;
  }
  if (!is_void(this.inplace)) {
    if (this.inplace && !am_subroutine()) error,"setup/plan was in-place";
  } else {
    h_set,this,inplace=am_subroutine();
  }
  if (is_void(this.setup)) {
    if (this.usefftw) {
      if (this.inplace) {
        local p;
        fftw,this.dims,this.ljdir,this.rjdir,p,fcplx=this.fcplx;      
        h_set,this,setup=p;      
      } else {
        h_set,this,setup=fftw(this.dims,this.ljdir,this.rjdir,\
              fcplx=this.fcplx);      
      }
    } else {
      h_set, this, setup=fft_setup(this.dims,this.ljdir,this.rjdir);
    }
  }
  if (this.inplace) {
    if (this.usefftw) {
      fftw,x,this.ljdir,this.rjdir,this.setup,fcplx=this.fcplx,keep=1;
    } else {
      fft_inplace,x,this.ljdir,this.rjdir,setup=this.setup;
    }
  } else {
    if (this.usefftw) {
      return fftw(x,this.ljdir,this.rjdir,this.setup,fcplx=this.fcplx,keep=1);
    } else {
      return fft(x,this.ljdir,this.rjdir,setup=this.setup);
    }
  }
}

func op_fft_reset (this)
{
  if (this.usefftw)
    fftw_clean, this.setup, fcplx=this.fcplx;
}

#if 0

func cplxgr (d)
{
  return random_n(d)+1i*random_n(d);
}

func cplx2f (z)
{
  d= dimsof(z);
  f= array(0.f,_(d(1)+1,2,d(2:)));
  f(1,..)= z.re;
  f(2,..)= z.im; 
  return f;
}

func opfftwcheck(void)
{
  d= [2,128,256];
  for (n=1,i=1;i<=d(1);i++) n*= d(i+1);

  oo= op_fft(dims=d,ljdir=1); 
  oow= op_fft(dims=d,ljdir=1,usefftw=1);
  ooi= op_fft(dims=d,ljdir=-1); 
  ooiw= op_fft(dims=d,ljdir=-1,usefftw=1);
  c= cplxgr(d); 
  cf= oo(c,1);
  cf= ooi(cf,-1)/n;
  cfw= oow(c,1);
  cfw= ooiw(cfw,-1)/n;
  op_fft_reset, oo;
  op_fft_reset, ooi;
  op_fft_reset, oow;
  op_fft_reset, ooiw;
  write,"\n double out-of-place";
  write,"fft  Mean, STD: ",(dd=abs(c-cf)(*))(avg),dd(rms);
  write,"fftw Mean, STD: ",(dd=abs(c-cfw)(*))(avg),dd(rms);


  oo= op_fft(dims=d,ljdir=1,inplace=1); 
  oow= op_fft(dims=d,ljdir=1,usefftw=1,inplace=1);
  ooi= op_fft(dims=d,ljdir=-1,inplace=1);
  ooiw= op_fft(dims=d,ljdir=-1,usefftw=1,inplace=1);
  c= cplxgr(d); 
  cf= c;
  oo, cf, 1;
  ooi, cf, -1;
  cf*= 1.0/n;
  cfw= c;
  oow, cfw, 1;
  ooiw, cfw, -1;
  cfw*= 1.0/n;
  op_fft_reset, oo;
  op_fft_reset, ooi;
  op_fft_reset, oow;
  op_fft_reset, ooiw;
  write,"\n double in-place";
  write,"fft  Mean, STD: ",(dd=abs(c-cf)(*))(avg),dd(rms);
  write,"fftw Mean, STD: ",(dd=abs(c-cfw)(*))(avg),dd(rms);

  d2= _(d(1)+1,2,d(2:));
  oow= op_fft(dims=d2,ljdir=1,usefftw=1,fcplx=1); 
  ooiw= op_fft(dims=d2,ljdir=-1,usefftw=1,fcplx=1);
  c= cplx2f(cplxgr(d)); 
  cfw= oow(c,1);
  cfw= ooiw(cfw,-1)/n;
  op_fft_reset, oow;
  op_fft_reset, ooiw;
  write,"\n float out-of-place";
  write,"fftw Mean, STD: ",(dd=abs(c-cfw)(*))(avg),dd(rms);


  oow= op_fft(dims=d2,ljdir=1,usefftw=1,inplace=1,fcplx=1); 
  ooiw= op_fft(dims=d2,ljdir=-1,usefftw=1,inplace=1,fcplx=1);
  c= cplx2f(cplxgr(d)); 
  cfw= c;
  oow, cfw, 1;
  ooiw, cfw, -1;
  cfw*= 1.0/n;
  op_fft_reset, oow;
  op_fft_reset, ooiw;
  write,"\n float in-place";
  write,"fftw Mean, STD: ",(dd=abs(c-cfw)(*))(avg),dd(rms);
}

func opfftwcheck2(void)
{
  d= [2,512,512];
  for(n=1,i=1;i<=d(1);i++)n*= d(i+1);
  d2= _(d(1)+1,2,d(2:))
  oow= op_fft(dims=d2,ljdir=1,usefftw=1,fcplx=1); 
  ooiw= op_fft(dims=d2,ljdir=-1,usefftw=1,fcplx=1);
  "init1"
  c= cplx2f(cplxgr(d)); 
  cfw= oow(c,1);
  cfw= ooiw(cfw,-1)/n;
  op_fft_reset, oow;
  op_fft_reset, ooiw;
  write,"\n float out-of-place";
  write,"fftw Mean, STD: ",(dd=abs(c-cfw)(*))(avg),dd(rms);


  oow= op_fft(dims=d2,ljdir=1,usefftw=1,inplace=1,fcplx=1); 
  ooiw= op_fft(dims=d2,ljdir=-1,usefftw=1,inplace=1,fcplx=1);
  c= cplx2f(cplxgr(d)); 
  cfw= c;
  oow, cfw, 1;
  ooiw, cfw, -1;
  cfw*= 1.0/n;
  op_fft_reset, oow;
  op_fft_reset, ooiw;
  write,"\n float in-place";
  write,"fftw Mean, STD: ",(dd=abs(c-cfw)(*))(avg),dd(rms);
}

#endif
