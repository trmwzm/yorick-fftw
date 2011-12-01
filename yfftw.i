plug_in, "yfftw";

require, "rollarr.i";

FFTW_NLIMIT_1D = 12;  /* how fast is your box */
FFTW_NLIMIT_2D = 12;
FFTW_TIMELIMIT = 3*60.0;   /* double seconds */

/* Yfftw needs its own copy of wisdom files found at locations below. */
/* System-wide files might be in /etc and could be copied from there. */
/* Fftw3 distribution comes with exec. fftw[f]-wisdom;  see man page. */

fftw_home= get_env("HOME");
//fftw_home= "/my/wisdom";

FFTW_WISDOM_FNM = fftw_home+"/etc/fftw/wisdom_yorick";
FFTWF_WISDOM_FNM = fftw_home+"/etc/fftw/wisdom_f_yorick";

func fftw_wisdom (void)
/* DOCUMENT func fftw_wisdom(void)
   this function will run each time this file is included.
   It reads wisdom files.  To create wisdom files in case they do not exist,
   call _init_fftwf_plans, this will create optimized FFTW plans and save 
   them in the wisdom files if wisdom files are not found.
   SEE ALSO: _init_fftwf_plans, fftw_init_wisdom, fftw,
 */
{
  write,"";
  wisdom_file = FFTWF_WISDOM_FNM;
  if (open(wisdom_file,"r",1) && _fftwfI(wisdom_file)) { 
    write,"fftwf wisdom file loaded: "+wisdom_file;
  } else { //file does not exist or couldn not import
    write,"FFTW3 [float] wisdom file not found: "+wisdom_file;
    write,"A write/readable copy of wisdom file required to store new plans.";
    write,"Run shell> fftwf-wisdom (in fftw distribution) to generate wisdom file.";
    write," see FFTWF_WISDOM_FNM file name definition in yfftw.i";
    write,"shell> fftwf-wisdom -v -c -n -t 1.5 -o $HOME/etc/fftw/wisdom_f";
    write,"File FFTWF_WISDOM_FNM can also be updated or created using:";
    write,"  yorick> fftwf_init_wisdom, [log2nmax_1D,log2nmax_2D,...];";
    write,""
  }
  wisdom_file = FFTW_WISDOM_FNM;
  if (open(wisdom_file,"r",1) && _fftwI(wisdom_file)) {
    write,"fftw wisdom file loaded: "+wisdom_file;
  } else { //file does not exist
    write,"FFTW3 [double] wisdom file not found: "+wisdom_file;
    write,"A write/readable copy of wisdom file required to store new plans.";
    write,"Run shell> fftw-wisdom (in fftw distribution) to generate wisdom file.";
    write," see FFTW_WISDOM_FNM file name in yfftw.i";
    write,"shell> fftw-wisdom -v -c -n -t 1.5 -o $HOME/etc/fftw/wisdom";
    write,"File FFTW_WISDOM_FNM can also be updated or created using:";
    write,"  yorick> fftw_init_wisdom, [log2nmax_1D,log2nmax_2D,...];";
    write,"";
  }
  _fftwT, double(FFTW_TIMELIMIT);
  write,"fftw maximum planning time [sec] (_fftwT) set to: ",double(FFTW_TIMELIMIT);
}

func fftw_init_wisdom (nlimit)
/* DOCUMENT func init_fftw_wisdom(nlimit)
   Optimization and caching of FFTW plans for n-equi-dimensional arrays. 
     nlimit==[log2_1d_nmax, log2_2d_nmax,....log2_md_nmax].
     nlimit defaults to [FFTW_NLIMIT_1D, FFTW_NLIMIT_2D]. 
   Plans are saved to the wisdom file. (See FFTW_WISDOM_FNM.)
   example: fftw_init_wisdom, [4,3];
     saves plans for fftw{array(complex,[1,4-or-8-or-16])}
                 and fftw{array(complex,[2,4-or-8,4-or-8])}
     and for in-place or out-of-place fft's.
   SEE ALSO: _init_fftw_plans,fftwf_init_wisdom
 */
{
  if (nlimit == []) nlimit=[FFTW_NLIMIT_1D, FFTW_NLIMIT_2D];
  
  for (m=1; m<=numberof(nlimit); m++) { /* need to plan dir= +/-1 ? */
    for (n=2;n<=nlimit(m);n++) {
      d= _(m,array(2^n,m));
      p= fftw(d,1);
      p= [0,0];fftw,d,1,,p;
    }
  }
  return _fftwO(FFTW_WISDOM_FNM);
}

func fftwf_init_wisdom (nlimit)
/* DOCUMENT fftwf_init_wisdom (nlimit)
   Optimization and caching of FFTW plans for n-equi-dimensional arrays. 
     nlimit==[log2_1d_nmax, log2_2d_nmax,....log2_md_nmax].
     nlimit defaults to [FFTW_NLIMIT_1D, FFTW_NLIMIT_2D]. 
   Plans are saved to the wisdom file. (See FFTW_WISDOM_FNM.)
   example: fftw_init_wisdom, [4,3];
       saves plans for fftw{array(complex,[1,4-or-8-or-16])}
               and fftw{array(complex,[2,4-or-8,4-or-8])}
   and for in-place or out-of-place fft's.
   SEE ALSO: 
 */
{
  if (nlimit == []) nlimit=[FFTW_NLIMIT_1D, FFTW_NLIMIT_2D];

  for (m=1; m<=numberof(nlimit); m++) {
    for (n=2;n<=nlimit(m);n++) {
      d= _(m+1,2,array(2^n,m));
      p= fftw(d,1);
      p= [0,0];fftw,d,1,,p;
    }
  }
  return _fftwfO(FFTWF_WISDOM_FNM);
}

func fftw (x, ljdir, rjdir, &plan, nosave=, fcplx=, keep=)
/* DOCUMENT y= fftw(x, ljdir, rjdir, &plan, nosave=, fcplx=, keep=)
   main yorick_fftw(3) interface, much of it similar to fft.i.
   returns Y the complex Fastest Fourier Transform of array X.
   If FCPLX==1 input arrays X==array(float,2,..) are treated 
   as float-complex.
   
   PLAN == array(long, 2) 
           where plan(1) is 0 if there is no dir=+1 transforms
                      or is a casted pointer to an fftw plan
           where plan(2) is 0 if there is no dir=-1 transforms
                      or is a casted pointer to an fftw plan
           *ATTN* plan differ for in- and out-of-place transforms!!

   USAGE TABLE:   OUT-of-place                      IN-place
     --planning and execution:.....................................
     y= fftw(x, dir)                     fftw, x, dir
     y= fftw(x, ljdir, rjdir)            fftw, x, ljdir, rjdir               
     --planning only:..............................................
     *attn* *plan is different out- or in-place*
     plan= fftw(dimsof(x), dir)          fftw, dimsof(x), dir, [], plan        
     plan= fftw(dimsof(x), ljdir, rjdir) fftw, dimsof(x), ljdir, rjdir, plan 
     --plan execution only:........................................
     y= fftw(x, dir,[], plan)            fftw, x, dir,[],plan
     y= fftw(x, ljdir, rjdir, plan)      fftw, x, ljdir, rjdir, plan         
     
   NOTE ABOUT USAGE TABLE ABOVE:
    + X is of type==complex OR array(float,[nd+1,2,d_1,...,d_nd])
    + right-column usage in-place
    + left-column usage out-of-place
    + setting plan=[0,0] on input forces planning
    + Since plans are archived as soon as they are optimized it is 
      NOT a good idea not to do repeated calls to planning. If it
      is required use NOSAVE==1 to turn off wisdom caching. 
                  
   DIRECTION determines which direction the transform is in --
     e.g.- from time to frequency or vice-versa -- see below:

   DIRECTION    meaning                              
   ---------    -------                              
       1    "forward" transform (coefficients of exp(+i * 2*pi*kl/N))    
        on every dimension of X                      
      -1    "backward" transform (coefficients of exp(-i * 2*pi*kl/N))   
        on every dimension of X                      
   [1,-1,1] forward transform on first and third dimensions of X,        
        backward transform on second dimension of X (any other       
        dimensions remain untransformed)                 
   [-1,0,0,1]   backward transform on first dimension of X, forward      
        transform on fourth dimension of X               
      etc.                                   

   The third positional argument, if present, allows the direction       
   of dimensions of X to be specified relative to the final dimension        
   of X, instead of relative to the first dimension of X.  In this       
   case, both LJDIR and RJDIR must be vectors of integers -- the         
   scalar form is illegal:                           

      LJDIR    RJDIR      meaning                        
      -----    -----      -------                        
      []       [1]   forward transform last dimension of X          
      [1]      []   forward transform first dimension of X         
      []       [-1,-1]   backward transform last two dimensions of X,       
                         leaving any other dimensions untransformed         
   [-1,0,0,1]  []   backward transform on first dimension of X,        
                    forward transform on fourth dimension of X         
      []       [-1,0,0,1]   backward transform 4th to last dimension of X,  
                            forward transform on last dimension of X       
      etc.                                   

   Note that the final element of RJDIR corresponds to the last dimension    
   of X, while the initial element of LJDIR corresponds to the first         
   dimension of X.                               

   The explicit meaning of "forward" transform -- the coefficients of        
   exp(+i * 2*pi*kl/N) -- is:                            

   for j=1,...,n ;   result(j)= sum from k=1,...,n of                                   
                                    x(k)*exp(-i*(j-1)*(k-1)*2*pi/n)              
         where i=sqrt(-1)                   

   Note that the result is unnormalized.  Applying the "backward"        
   transform to the result of a "forward" transform returns N times      
   the original vector of length N.  Equivalently, applying either       
   the "forward" or "backward" transform four times in succession        
   yields N^2 times the original vector of length N.                 

   Performing transforms requires PLANNING, which can be set up        
   beforehand by calling fftw on a dimsof(X) array.  This allow fftw
   to be called more than once with arrays X of the same shape with
   KEEP==1.

   SEE ALSO:  fft, fft_inplace, fft_setup, fftw[f]_init_wisdom, fftw_wisdom,
             _fftwfI, _fftwfO, _fftwfP, _fftwfE, _fftwfD, _fftwfS, _fftwfC,
             _fftwT,  _fftwI, _fftwO, _fftwP, _fftwE, _fftwD, _fftwS, _fftwC 
 */
{
  if (is_void(plan)) {
    plan= array(long,2);
    iplan= 1;
  } else if (numberof(plan)==2 && allof(plan==[0,0])) {
    iplan= 1;
  } else {
    iplan= 0;
  }

  if (structof(x)==long) {  /* planning X dimsof(x) if d1==2 assume fftwf */
    iplan= 1;
    plan= array(long,2);
    planex= 0;
    if (x(1)==0) return;
    if (x(1)>7) error,"numberof dimensions > 7 ?";
    dims= x;
    if (fcplx!=1) {
      x= array(complex, dims);
    } else {
      if (dims(2)!=2) 
        error,"fcplx==1: float arrays leading dim==2, but dimsof(x)(2)="+pr1(dims(2));
      x= array(float, dims);
    }
  } else {
    planex= 1;
  }

  /* form list of dimension lengths (dims) and corresponding list of
     transform directions (dirs) */
  dims= dimsof(x);
  if (structof(x)==complex) {
    dcplx= 1;
    if (fcplx==1)
      error,"fcplx==1, but double-complex array to transform";
  } else if (structof(x)==float) {
    if (fcplx!=1)
      error,"fcplx!=1, but float-complex array to transform, expecting fcplx==1";
    dcplx= 0;
    dims= _(dims(1)-1,dims(3:));
  } else {
    x+= 0i;
    dcplx= 1;
  } 

  /* get dims & strides associated with each dimension & fft dimension mask */
  ndims= dims(1);
  if (ndims<1) return;
  dims= dims(2:0);
  dirs= fft_dirs(ndims, ljdir, rjdir);
  stds= array(1, ndims);
  for (i=1 ; i<ndims ; i++) stds(i+1)= stds(i)*dims(i);
  //tops= numberof(x)/(stds*dims); /* tops == total number of vectors for each dim */

  if (am_subroutine()) {  /*inplace in==x out=xx */
    inplace= 1;
    if (planex==1) {
      eq_nocopy, xx, x;
    }
    if (iplan==1) {
      if (!planex) {
        eq_nocopy, xop, x;      
      } else {
        xop= array(structof(x), dimsof(x));
      }
      eq_nocopy, xip, xop;      
    }
  } else {
    inplace= 0;
    if (planex==1) {
      xx= array(structof(x), dimsof(x));
    }
    if (iplan==1) {
      if (!planex) {
        eq_nocopy, xop, x; 
        xip = array(structof(x), dimsof(x));
      } else {
        eq_nocopy, xop, xx; 
        xip = array(structof(x), dimsof(x));
      }
    }
  }
    
  for(idir=1;idir>-2;idir-=2) {  /* loop 1,-1 : plan(i) i=1: out= c_k[*eikl/N]  i=2:..(-i)...*/
    ip= 1+(idir<0);
    mksdir= dirs==idir;
    if (anyof(mksdir)) {
      if (iplan==1) {
        list= where(mksdir);    /* list of dimensions to be transformed */
        n= numberof(list);
        nlist= where(!mksdir);  /* list of dimensions to be looped over */
        dimlp= dims(nlist);
        stdlp= stds(nlist);
        if (is_void(dimlp)) dimlp= [0];
        if (is_void(stdlp)) stdlp= [0];
        if (dcplx==1) {
          plan(ip)= _fftwP(xip,xop,n,dims(list),stds(list),ndims-n,dimlp,stdlp,idir);
          if (nosave != 1) _fftwO,FFTW_WISDOM_FNM; 
        } else {
          plan(ip)= _fftwfP(xip,xop,n,dims(list),stds(list),ndims-n,dimlp,stdlp,idir);
          if (nosave != 1) _fftwfO,FFTWF_WISDOM_FNM;
        }
      }
      if (planex==1) {
        if (dcplx==1) {
          _fftwE, plan(ip), x, xx;
          if (iplan==1 && keep!=1) {_fftwD, plan(ip);plan(ip)= 0;};
        } else {
          _fftwfE, plan(ip), x, xx;
          if (iplan==1 && keep!=1) {_fftwfD, plan(ip);plan(ip)= 0;}
        }
        if (idir==1) x= xx;
      } 
    }
  }
  if (inplace!=1) {
    if (planex==0) {
      return plan;
    } else {
      return xx;
    }
  }
}

func fftw_clean (&p, fcplx=)
/* DOCUMENT fftw_clean, p, fcplx=;
   free-up memory space allocated by FFTW plan P.
   FCPLX=1 if the plan was for a float-complex transform.
   Once space is deallocated, returns a null plan P==[0,0];
 */
{
  if (fcplx!=1) 
    for (i=1;i<=numberof(p);i++) 
      if (p(i)) {_fftwD, p(i);p(i)= 0;}
  else
    for (i=1;i<=numberof(p);i++) 
      if (p(i)) {_fftwfD, p(i);p(i)= 0;}
}

/////////////////////////////////////////////////////////////////////////

extern  _fftwT
/* PROTOTYPE
   long _fftwT (double tlim)
*/
/* DOCUMENT long _fftwT (double tlim)
   set rough [T]time limit on planning
 */
extern  _fftwI
/* PROTOTYPE
   long _fftwI (string wisdom_file)
*/
/* DOCUMENT long _fftwI (string wisdom_file)
   [I]nput plans from wisdom file "wisdom_file"
 */
extern  _fftwO
/* PROTOTYPE
   long _fftwO (string wisdom_file)
*/
/* DOCUMENT long _fftwO (string wisdom_file)
   [0]utput plan(s) to wisdom file "wisdom_file"
 */
extern  _fftwP
/* PROTOTYPE
   long _fftwP (complex array in, complex array out, long ft_rank, long array ft_dims,
                long array ft_strides, long ft_loop, long array loop_dims,
            long array loop_strides, long dir)
*/
/* DOCUMENT long _fftwP (complex array in, complex array out, long ft_rank, long array ft_dims,
                         long array ft_strides, long ft_loop, long array loop_dims,
                     long array loop_strides, long dir)
   [P]lan multidimensional ffts for double-complex
   in-place in==out      ! fftw3 doc says in must be allocated and IS overwritten
   out-of-place in!=out  ! fftw3 doc says out must be allocated but is not overwritten
   
   ft-rank #dimensions of the fft (number of array dimensions to trasform)
   ft-dims array of the lengths of each of the array dimensions to transform
   ft-strides array of the strides of each of the dimensions to transform
   ft-loop #dimensions to loop over (number of array dimensions NOT to trasform)
   loop_dims array of the lengths of each of the array dimensions NOT to transform
   loop_strides array of the strides of each of the dimensions NOT to transform
   
   dir is in YORICK's convention.
 */
extern  _fftwE
/* PROTOTYPE
   void _fftwE (long plan, complex array in, complex array out)
*/
/* DOCUMENT _fftwE (long plan, complex array in, complex array out)
   [E]xecute plan ! data & fft geometry must be identical to plan
 */
extern  _fftwD
/* PROTOTYPE
   void _fftwD (long plan) 
*/
/* DOCUMENT  void _fftwD (long plan) 
   usage: _fftwD,p; p= [];
   [D]estroy plan
*/
extern  _fftwS
/* PROTOTYPE
   void _fftwS (long plan) 
*/
/* DOCUMENT  void _fftwS (long plan) 
   print plan to [S]creen (needs a newline) 
*/
extern  _fftwC
/* PROTOTYPE
   void _fftwC (void)
*/
/* DOCUMENT  void _fftwC (void)
   [C]learing plans 
*/
/*................................... id for float-complex ........................*/
extern _fftwfI
/* PROTOTYPE
   long _fftwfI (string wisdom_file)
*/
/* DOCUMENT long _fftwfI (string wisdom_file) 
   [I]nput plans from wisdom file "wisdom_file"
*/
extern _fftwfO
/* PROTOTYPE
   long _fftwfO (string wisdom_file)
*/
/* DOCUMENT  long _fftwfO (string wisdom_file)
   [0]utput plan(s) to wisdom file "wisdom_file"
*/
extern _fftwfP
/* PROTOTYPE
   long _fftwfP (float array in, float array out, long ft_rank, long array ft_dims,
                 long array ft_strides, long ft_loop, long array loop_dims,
             long array loop_strides, long dir)
*/
/* DOCUMENT long _fftwfP (real array in, real array out, long ft_rank, long array ft_dims,
                          long array ft_strides, long ft_loop, long array loop_dims,
                      long array loop_strides, long dir)
   [P]lan multidimensional ffts for float-complex
   in-place in==out      ! fftw3 doc says in must be allocated and IS overwritten
   out-of-place in!=out  ! fftw3 doc says out must be allocated but is not overwritten
   
   ft-rank #dimensions of the fft (number of array dimensions to trasform)
   ft-dims array of the lengths of each of the array dimensions to transform
   ft-strides array of the strides of each of the dimensions to transform
   ft-loop #dimensions to loop over (number of array dimensions NOT to trasform)
   loop_dims array of the lengths of each of the array dimensions NOT to transform
   loop_strides array of the strides of each of the dimensions NOT to transform
   
   dir is in YORICK's convention.
 */
extern _fftwfE
/* PROTOTYPE
   void _fftwfE (long plan, float array in, float array out)
*/
/* DOCUMENT  void _fftwfE (long plan, float array in, float array out)
   [E]xecute plan ! data & fft geometry must be identical to plan
*/
extern _fftwfD
/* PROTOTYPE
   void _fftwfD (long plan) 
*/
/* DOCUMENT  void _fftwfD (long plan) 
   usage: _fftwfD,p; p= [];
   [D]estroy plan
*/
extern _fftwfS
/* PROTOTYPE
   void _fftwfS (long plan) 
*/
/* DOCUMENT  void _fftwfS (long plan)
   print plan to [S]creen (needs a newline) 
*/
extern _fftwfC
/* PROTOTYPE
   void _fftwfC (void)  
*/
/* DOCUMENT  void _fftwfC (void) 
   [C]learing plans 
*/
//////////////////////////////////////////////

fftw_wisdom;
