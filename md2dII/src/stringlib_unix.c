readint_(a,len,ix)
 
char *a;
int *ix,*len;
{
  char b[20];
  int i;
 
  for (i=0; i < *len; i++)
    b[i] = a[i];
  b[*len] = 0;
  sscanf(b,"%d",ix);
}
 
readdouble_(a,len,x)
 
char *a;
int *len;
double *x;
{
  char b[20];
  int i;
 
  for (i=0; i < *len; i++)
    b[i] = a[i];
  b[*len] = 0;
  sscanf(b,"%lf",x);
}
 
writeint_(x,a,len)
 
int *x;
char *a;
int *len;
{
  char b[20];
  int i;
 
  sprintf(b,"%d",*x);
  *len = strlen(b);
  for (i=0; i < *len; i++)
    a[i] = b[i];
}                               
