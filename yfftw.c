#include <stdlib.h>
#include <fftw3.h>
#include <stdlib.h>

/* GENERAL ------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
long _fftwT (double tlim) /* set time limit on planning */
{
  fftw_set_timelimit(tlim);
  return (1);
}

/* DOUBLE ------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
long _fftwI (char *wisdom_file) /* import */
{
  FILE *fp;
  fp = fopen(wisdom_file, "r");
  if((fp = fopen(wisdom_file, "r"))==NULL) {
    printf("Error reading wisdom file\"%s\"\n",wisdom_file);
    fflush(stdout);
    exit(0);
  }
  fftw_import_wisdom_from_file(fp);
  fclose(fp);
  return (1);
}
/*--------------------------------------------------------------------------*/
long _fftwO (char *wisdom_file) /* export */
{
  FILE *fp;
  if((fp = fopen(wisdom_file, "w"))==NULL) {
    printf("Error creating wisdom file\"%s\"\n",wisdom_file);
    fflush(stdout);
    exit(0);
  }
  fftw_export_wisdom_to_file(fp);
  fflush(fp);
  fclose(fp);
  return (1);
}
/*--------------------------------------------------------------------------*/
fftw_plan _fftwP (fftw_complex *in,   /* plan */
	          fftw_complex *out,
	          long ft_rank, long *ft_dims, long *ft_strides, 
	          long ft_loop, long *loop_dims, long *loop_strides, 
	          long dir)                         /* forward (1) or reverse (-1) */
{
  /* Declarations */
  int i;
  int sign= -1*dir;
  fftw_plan p;
  long plan_mode;
  fftw_iodim *id_fft;
  fftw_iodim *id_loop;

  plan_mode = FFTW_MEASURE; /* FFTW_ESTIMATE FFTW_MEASURE FFTW_PATIENT FFTW_EXHAUSTIVE;*/

  id_fft = (fftw_iodim *)fftw_malloc(sizeof(fftw_iodim) * ft_rank);
  id_loop = (fftw_iodim *)fftw_malloc(sizeof(fftw_iodim) * ft_loop);

  for (i = 0 ; i < ft_rank ; i++) {
    id_fft[i].n = ft_dims[i];
    id_fft[i].is = ft_strides[i];
    id_fft[i].os = ft_strides[i];
  }
  for (i = 0 ; i < ft_loop ; i++) {
    id_loop[i].n = loop_dims[i];
    id_loop[i].is = loop_strides[i];
    id_loop[i].os = loop_strides[i];
  }
    
  if (in == out){ /*in place*/
    p = fftw_plan_guru_dft((int)ft_rank, id_fft, (int)ft_loop, id_loop, in, in, sign, plan_mode);
  }else{
    p = fftw_plan_guru_dft((int)ft_rank, id_fft, (int)ft_loop, id_loop, in, out, sign, plan_mode);
  }
  
  fftw_free(id_fft);
  fftw_free(id_loop);

  return p;
}
/*--------------------------------------------------------------------------*/
void _fftwE (fftw_plan p, fftw_complex *in, fftw_complex *out)  /*execute*/
{
  fftw_execute_dft(p, in, out);
  return;
}
/*--------------------------------------------------------------------------*/
void _fftwD (fftw_plan p)  /*destroy*/
{
  fftw_destroy_plan(p);
  return;
}

/*--------------------------------------------------------------------------*/
void _fftwS (fftw_plan p)  /*print S==screen ;)*/
{
  fftw_print_plan(p);
  return;
}

/*--------------------------------------------------------------------------*/
void _fftwC (void)  /*cleanup*/
{
  fftw_cleanup();
  return;
}

/* FLOAT ------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
long _fftwfI (char *wisdom_file) /* import */
{
  FILE *fp;
  fp = fopen(wisdom_file, "r");
  if((fp = fopen(wisdom_file, "r"))==NULL) {
    printf("Error reading wisdom file\"%s\"\n",wisdom_file);
    fflush(stdout);
    exit(0);
  }
  fftwf_import_wisdom_from_file(fp);
  return (1);
}
/*--------------------------------------------------------------------------*/
long _fftwfO (char *wisdom_file) /* export */
{
  FILE *fp;
  if((fp = fopen(wisdom_file, "w"))==NULL) {
    printf("Error creating wisdom file\"%s\"\n",wisdom_file);
    fflush(stdout);
    exit(0);
  }
  fftwf_export_wisdom_to_file(fp);
  fflush(fp);
  return (1);
}
/*--------------------------------------------------------------------------*/
fftwf_plan _fftwfP (fftwf_complex *in,   /* plan */
	          fftwf_complex *out,
	          long ft_rank, long *ft_dims, long *ft_strides, 
	          long ft_loop, long *loop_dims, long *loop_strides, 
	          long dir)                         /* forward (1) or reverse (-1) */
{
  /* Declarations */
  int i;
  int sign= -1*dir;
  fftwf_plan p;
  long plan_mode;
  fftwf_iodim *id_fft;
  fftwf_iodim *id_loop;

  plan_mode = FFTW_PATIENT; /* FFTW_ESTIMATE FFTW_MEASURE FFTW_PATIENT FFTW_EXHAUSTIVE;*/

  id_fft = (fftwf_iodim *)fftwf_malloc(sizeof(fftwf_iodim) * ft_rank);
  id_loop = (fftwf_iodim *)fftwf_malloc(sizeof(fftwf_iodim) * ft_loop);

  for (i = 0 ; i < ft_rank ; i++) {
    id_fft[i].n = ft_dims[i];
    id_fft[i].is = ft_strides[i];
    id_fft[i].os = ft_strides[i];
  }
  for (i = 0 ; i < ft_loop ; i++) {
    id_loop[i].n = loop_dims[i];
    id_loop[i].is = loop_strides[i];
    id_loop[i].os = loop_strides[i];
  }
    
  if (in == out){ /*in place*/
    p = fftwf_plan_guru_dft((int)ft_rank, id_fft, (int)ft_loop, id_loop, in, in, sign, plan_mode);
  }else{
    p = fftwf_plan_guru_dft((int)ft_rank, id_fft, (int)ft_loop, id_loop, in, out, sign, plan_mode);
  }
  
  fftwf_free(id_fft);
  fftwf_free(id_loop);

  return p;
}
/*--------------------------------------------------------------------------*/
void _fftwfE (fftwf_plan p, fftwf_complex *in, fftwf_complex *out)  /*execute*/
{
  fftwf_execute_dft(p, in, out);
  return;
}
/*--------------------------------------------------------------------------*/
void _fftwfD (fftwf_plan p)  /*destroy*/
{
  fftwf_destroy_plan(p);
  return;
}

/*--------------------------------------------------------------------------*/
void _fftwfS (fftwf_plan p)  /*print S==screen ;)*/
{
  fftwf_print_plan(p);
  return;
}

/*--------------------------------------------------------------------------*/
void _fftwfC (void)  /*cleanup*/
{
  fftwf_cleanup();
  return;
}
