plug_dir, ["."];

require, "yfftw.i";
require, "offt.i";

func _fftw_check (void)
{
  if (is_void(n)) n= 1024;
  
  write,""
  write,"Exec. this twice to see the effect of cached wisdom";
 
  write,""
  write,"Caching or retrieving optimization info for pow(2) 1&2-D fft's....patience";

  elapsed= elapsed0= split1= split2= split3= split4= split5= [0., 0., 0.];
  //check init function
  timer, elapsed;elapsed0= elapsed; 
    ck_init;
  timer, elapsed, split1;
  write,"Done caching wisdom";

  //various usages examples
  d= [1,128];
  fd= [2,2,128]; //watch it (no fun)
  ld= 1;
  //double /*---------------------------------------*/
  c= cplxgr(d);
  //out-of-place
  fw= fftw(c, ld);            /* plan and exec */
  plan= fftw(dimsof(c), ld);          /* plan */
  fw= fftw(c, ld, [], plan);  /* exec */
  fftw_clean, plan;
  //in-place
  fw= c;
  fftw, fw, ld;               /* plan and exec */
  fftw, dimsof(c), ld, [], plan;      /* plan */
  fw= c;
  fftw, fw, ld, [], plan;      /* exec */
  fftw_clean, plan;
  write,"Done double use cases";
  //float /*---------------------------------------*/
  c= cplx2f(c);
  //out-of-place
  fw= fftw(c, ld, fcplx=1);            /* plan and exec */
  plan= fftw(dimsof(c), ld, fcplx=1);          /* plan */
  fw= fftw(c, ld, [], plan, fcplx=1);  /* exec */
  fftw_clean, plan, fcplx=1;
  //in-place
  fw= c;
  fftw, fw, ld, fcplx=1;               /* plan and exec */
  fftw, dimsof(c), ld, [], plan, fcplx=1;      /* plan */
  fw= c;
  fftw, fw, ld, [], plan, fcplx=1;      /* exec */
  fftw_clean, plan, fcplx=1;
  write,"Done float use cases";

  timer, elapsed, split2;

  //many dimensions
  ck1,[3,32,17,13],[0,1,-1];
  ck1,[5,12,5,9,22,11],[1,1,0,0,-1];

  //Swarztrauber vs. fftw timing
  d= [2,512,512]
  c= cplxgr(d);
  ld= 1;
  write,"";
  write,"Timing comparison: fft(x,ld) with dimsof(x)="+pr1(d)+" and ld="+pr1(ld);

  //fftw
  plan= fftw(d, ld);          /* plan */
  for (i=0 ; i<32 ; i++) fw= fftw(c, ld, [], plan, keep=1);  /* exec */
  fftw_clean, plan;
  timer, elapsed, split3;
  write,"Fftw done:";

  //fft 
  ws= fft_setup(d);
  for (i=0 ; i<32 ; i++) fs= fft(c, ld, setup=ws);  /* exec */
  timer, elapsed, split4;
  write,"Swarztrauber done:";
 
  //testm.i fft tests
  write,""
  write, "testing fftw routines as in testm.i...";
  force2d = 0;
  fftw_test, 95;
  fftw_test, 32;
  fftw_test, 12;
  fftw_test, 1024;
  fftw_test, 4096;
  for (i=0 ; i<32 ; i++) fftw_test, 4096, nosave=1;
  force2d = 1;  fftw_test, 128;  force2d = 0;
  timer, elapsed, split5;

  write, "\n\n------------  yfftw tests complete ------------";
  write, "\n\n------------------ timing info ----------------";
  timer_print, "FFTW pow2 initialization   ", split1,\
               "1-d use tests              ", split2,\
               "FFTW plan and 32Xexec      ", split3,\
               "Swarztrauber setup+32Xexec ", split4,\
               "testm.i, but with FFTW     ", split5,\
               "-total                     ", elapsed-elapsed0;
  return;
}

func ck1(d,ld)
{
  write,""
  write,"fft(x,ld) with dimsof(x)="+pr1(d)+" and ld="+pr1(ld);
  
  c= cplxgr(d);
    
  fs= fft(c, ld);
  fw= fftw(c, ld);
  write,(fs.re-fw.re)(*)(ptp),(fs.im-fw.im)(*)(ptp),\
     format=" PTP-Diffs {Double_Swarztrauber - Double_fftw} [Re, IM]: %15.8lg %15.8lg \n";

  fc= cplx2f(c);
  fw= fftw(fc, ld, fcplx=1);
  fw= fw(1,..)+1i*fw(2,..);
  write,(fs.re-fw.re)(*)(ptp),(fs.im-fw.im)(*)(ptp),\
     format=" PTP-Diffs {Double_Swarztrauber - Float_fftw}  [Re, IM]: %15.8lg %15.8lg \n";
  return;
}

func ck_init (void)
{
  fftw_init_wisdom, [12,11];
  fftwf_init_wisdom, [12,11];
  return;
}

func fcabs(x) {return abs(x(1,..),x(2,..));}

func fftw_test(n, nosave=)
{
  index= (2*pi*indgen(0:n-1))/n;
  z= sin(index*3);
  fz= array(float,[2,2,n]);
  fz(1,..)= z;
  //
  zf= fftw(z, 1,nosave=nosave);
  fzf= fftw(fz, 1,nosave=nosave,fcplx=1);
  z3= z2= array(0i, n);
  z3(4)= -0.5i*n;
  z3(-2)= 0.5i*n;
  zb= fftw(z, -1,nosave=nosave);
  fzb= fftw(fz, -1,nosave=nosave,fcplx=1);
  if (max(abs(zf-z3))>1.e-12*n || max(abs(zb-conj(z3)))>1.e-12*n)
    write, "***WARNING*** failed double 1D fftw test";
  if (max(fcabs(fzf-cplx2f(z3)))>1.e-6*n || max(fcabs(fzb-cplx2f(conj(z3))))>1.e-6*n)
    write, "***WARNING*** failed float 1D fftw test";
  if (n<=96 || force2d) {
    z*= cos(index*2)(-,);
    fz*= float(cos(index*2))(-,-,);
    zf= fftw(z, [0, 1],nosave=nosave);
    fzf= fftw(fz, [0, 1],nosave=nosave,fcplx=1);
    z2(3)= z2(-1)= 0.5*n;
    zb= fftw(z, [], [-1, 0],nosave=nosave);
    fzb= fftw(fz, [], [-1, 0],nosave=nosave,fcplx=1);
    if (max(abs(zf-sin(index*3)*z2(-,)))>1.e-12*n ||
        max(abs(zb-conj(z3)*cos(index*2)(-,)))>1.e-12*n)
      error, "***WARNING*** failed first double 2D fftw test";
    if (max(fcabs(fzf-cplx2f(sin(index*3)*z2(-,))))>1.e-6*n ||
        max(fcabs(fzb-cplx2f(conj(z3)*cos(index*2)(-,))))>1.e-6*n)
      error, "***WARNING*** failed first float 2D fftw test";
    zf= fftw(z, 1,nosave=nosave);
    zb= fftw(z, -1,nosave=nosave);
    fzf= fftw(fz, 1,nosave=nosave,fcplx=1);
    fzb= fftw(fz, -1,nosave=nosave,fcplx=1);
    if (max(abs(zf-z3*z2(-,)))>1.e-12*n ||
        max(abs(zb-conj(z3)*z2(-,)))>1.e-12*n)
      error, "***WARNING*** failed second double 2D fftw test";
    if (max(fcabs(fzf-cplx2f(z3*z2(-,))))>1.e-5*n ||
        max(fcabs(fzb-cplx2f(conj(z3)*z2(-,))))>1.e-5*n)
      error, "***WARNING*** failed second float 2D fftw test";
  }
}

func cplxgr (d) {
  return random_n(d)+1i*random_n(d);
}

func cplx2f (z) {
  d= dimsof(z);
  f= array(0.f,_(d(1)+1,2,d(2:)));
  f(1,..)= z.re; f(2,..)= z.im; 
  return f;
}

func _offt_check (d,usefftw=)
{
  if (is_void(d)) d= [2,128,256];
  for (n=1,i=1;i<=d(1);i++)
    n*= d(i+1);

  oo= offt(usefftw=usefftw);   
  ooi= offt(usefftw=usefftw); 
  
  c= cplxgr(d); 
  cf= oo(c,1);
  cf= ooi(cf,-1)/n;
  oo, reset=1;
  ooi, reset=1;
 
  sf= (usefftw==1?" FFTW":string(0));

  write,"\n double out-of-place"+sf;
  write,"fft  Mean, STD: ",(dd=abs(c-cf)(*))(avg),dd(rms);

  oo= offt(usefftw=usefftw);
  ooi= offt(usefftw=usefftw);

  c= cplxgr(d); 
  cf= c;
  oo, cf, 1;
  ooi, cf, -1;
  cf*= 1.0/n;
  oo, reset=1;
  ooi, reset=1;

  write,"\n double in-place"+sf;
  write,"fft  Mean, STD: ",(dd=abs(c-cf)(*))(avg),dd(rms);

  if (usefftw==1) {
    d2= _(d(1)+1,2,d(2:));

    oow= offt(usefftw=1); 
    ooiw= offt(usefftw=1);

    c= cplx2f(cplxgr(d)); 
    cfw= oow(c,1,fcplx=1);
    cfw= ooiw(cfw,-1,fcplx=1)/n;
    oow, reset=1;
    ooiw, reset=1;
    write,"\n float out-of-place"+sf;
    write,"fftw Mean, STD: ",(dd=abs(c-cfw)(*))(avg),dd(rms);

    oow= offt(usefftw=1,fcplx=1); 
    ooiw= offt(usefftw=1,fcplx=1);
    c= cplx2f(cplxgr(d)); 
    cfw= c;
    oow, cfw, 1;
    ooiw, cfw, -1;
    cfw*= 1.0/n;
    oow, reset=1;
    ooiw, reset=1;
    write,"\n float in-place"+sf;
    write,"fftw Mean, STD: ",(dd=abs(c-cfw)(*))(avg),dd(rms);
  }
}

if (batch()) {
  _fftw_check;
  _offt_check;
  _offt_check,usefftw=1;
}

