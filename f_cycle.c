#include<stdio.h>
#include<stdlib.h>
#include<math.h>

void gauss_seidel(double v[],double f[], int l);
void w_jacobi(int M, double v[], double f[]);
void restriction_methods(int l, double f[], double r[]);
void prolongation_methods(int l, double v[] ,double c[]);
void v_cycle(int l,double v[][513], double f[][513]);
void FMG(int l,double v[][513], double f[][513]);
double compute_norm(int M,double f[]);

void main()
{
    int i,j,k,M,L,iteration,nu,method;
    FILE *f1 = fopen("input_parameters.txt","r");
    fscanf(f1,"grid_size = %d\n",&M);
    fscanf(f1,"nu = %d\n",&nu);
    fscanf(f1,"k = %d\n",&k);
    L = log(M)/log(2);
    //Declaring v anf f as 2-dimensional array where 1st row corresponds to level 1 and so on
    double h,v[L][M+1],f[L][M+1],residue[M+1],r_norm,PI,sigma;
    fscanf(f1,"sigma = %lf",&sigma);
    PI = 22.0/7;
    h = 1.0/M;
    iteration = 0;
    fclose(f1);
    //Initializing the v and f vectors
    for(j=0;j<L;j++)
    {
        for(i=0;i<=M;i++)
        {
            v[j][i] = 0.0;
            f[j][i] = (PI*PI*k*k + sigma)*sin(k*PI*h*i);
        }
    }
    printf("Choose the method to solve from follows (Enter 1,2 or 3):\n 1. V-cycle(using gauss-seidel relaxation) 2. Basic gauss-seidel 3. Weighted jacobi\n ");
    scanf("%d",&method);
    FILE *fw = fopen("output.txt","w");
    fprintf(fw,"Iteration\t%Residual_2norm\n");
    do
    {
        if(method == 1)
        {
            //Choose whether to carry out first iteration as F cycle or not
            if(iteration == 0)
            {
                FMG(9,v,f);
                //v_cycle(1,v,f);
            }
            else
            {
                v_cycle(1,v,f);
            }
            //Resetting the values of initial error at levels greater than 1(for coarser grids)
            for(j=1;j<L;j++)
            {
                for(i=0;i<=M;i++)
                {
                    v[j][i] = 0.0;
                }
            }
        }
        else if(method == 2)
        {
            gauss_seidel(v[0],f[0],1);
        }
        else if(method == 3)
        {
            w_jacobi(M, v[0],f[0]);
        }
        //Calculating Residue
        for(i=1;i<M;i++)
            {
                residue[i] = f[0][i] + ( (v[0][i+1] + v[0][i-1]) - (2+h*h)*v[0][i] )/(h*h);
            }
            residue[0] = 0.0;
            residue[M] = 0.0;
        //Calculating second norm of residual error
        r_norm = compute_norm(M,residue);
        iteration++;
        printf("iteration %d\t r_2norm %lf\n",iteration,r_norm);
        fprintf(fw,"%d\t%lf\n",iteration,r_norm);
    }while(r_norm>.000001); //termination condition
    for(i=0;i<=M;i++)
    {
        printf("v %lf \t f %lf \t res %lf \n",v[0][i],f[0][i],residue[i]);
    }
    fclose(fw);
}

double compute_norm(int M,double f[])
{
    //calculating e2 norm
    int i;
    double sum=0;
    for(i=0;i<=M;i++)
    {
        sum = sum + pow(f[i],2);
    }
    sum = sqrt(sum);
    return sum;
}

void prolongation_methods(int l, double v[] ,double c[])
{
    //using full weighting operator
    int m,n;
    m = 512/pow(2,l-1);
    //n = 2*m;

    for(int i=0;i<=m;i++)
    {
        if(i>0)
        {
            c[2*i-1] = (v[i-1] + v[i])/2.0;
        }
        c[2*i] = v[i];
    }
}

void w_jacobi(int M, double v[], double f[])
{
    int i;
    double v_old[M+1],h,beta;
    h = 1.0/M;
    beta = 1.0/(2 + h*h);
    //saving previous values
    for(i=1;i<M;i++)
    {
        v_old[i] = v[i];
    }
    for(i=1;i<M;i++)//computing error only at interior points
    {
        v[i] = (1.0/3)*v_old[i] + (2.0/3)*beta*(h*h*f[i] + (v_old[i-1] + v_old[i+1]) );
    }
}

void gauss_seidel(double v[],double f[], int l)
{
    int i,k,M,nu;
    M = 512/pow(2,l-1);
    double h,beta;
    h = 1.0/M;
    beta = 1.0/(2 + h*h);
    k = 0;
    nu = 2;
    do
    {
        for(i=1;i<M;i++)//computing error only at interior points
        {
            v[i] = h*h*beta*f[i] + beta*(v[i-1] + v[i+1]);
        }
     k=k+1;
    }while(k < nu);
}

void restriction_methods(int l, double f[], double r[])
{
    //using full weighting operator
    int m,n;
    m = 512/pow(2,l-1);
    n = m/2;
    for(int i=1;i<n;i++)
    {
        f[i] = ( r[2*i-1] + 2*r[2*i] + r[2*i+1] )/4;
    }
    f[0] = r[0];
    f[n] = r[m];
}

void v_cycle(int l,double v[][513], double f[][513])
{
    int i,M,stop;
    stop = l; //to control the stopping step in FMG cycle
    M = 512/pow(2,l-1);
    double residue[513],correction[513],h;
    do
    {
        /*
        printf("At level %d, f vector is :\n",l);
        for(i=0;i<=M;i++)
        {
            printf("%lf ",f[l-1][i]);
        }
        printf("\n");
        */
        h = 1.0/M;
        //Relaxation step
        gauss_seidel(v[l-1],f[l-1],l);
        //calculating residue
        for(i=1;i<M;i++)
        {
            residue[i] = f[l-1][i] + ( (v[l-1][i+1] + v[l-1][i-1]) - (2+h*h)*v[l-1][i] )/(h*h);
        }
        residue[0] = 0.0;
        residue[M] = 0.0;
        //restricting to coarser grid
        restriction_methods(l,f[l],residue);
        M = M/2;
        l++;
    }while(l<9);
    gauss_seidel(v[l-1],f[l-1],l);
    do
    {
        //initializing array to store interpolated error values
        for(i=0;i<=512;i++)
        {
            correction[i] = 0.0;
        }
        //carrying out interpolation
        prolongation_methods(l,v[l-1],correction);
        l--;
        for(i=0;i<=512;i++)
        {
            v[l-1][i] = v[l-1][i] + correction[i];
        }
        //relaxation step
        gauss_seidel(v[l-1],f[l-1],l);
    }while(l>stop);
}

void FMG(int l,double v[][513], double f[][513])
{
    int i,j;
    gauss_seidel(v[l-1],f[l-1],l);
    do
    {
        //carrying out interpolation
        prolongation_methods(l,v[l-1],v[l-2]);
        l--;
        //Resetting the values of initial error at levels greater than 1(for coarser grids)
        for(j=l;j<9;j++)
        {
            for(i=0;i<=512;i++)
            {
                v[j][i] = 0.0;
            }
        }
        //one iteration of v cycle from current level
        v_cycle(l,v,f);
    }while(l>1);
    /*
    printf("At level %d, v approximated is :\n",l);
    for(i=0;i<=512;i++)
    {
        printf("%lf ",v[l-1][i]);
    }
    printf("\n");
    */
}
