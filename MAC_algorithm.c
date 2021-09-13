#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#define tol pow(10,-5)
#define nx 129
#define ny 129
#define dt 0.001
#define xlen 1.0
#define ylen 1.0
#define Re 400.0

int main()
{
    
    double dx=xlen/(nx-1), dy=ylen/(ny-1), b=dx/dy;
    double u[nx][ny+1]={0.}, v[nx+1][ny]={0.}, p[nx+1][ny+1]={0.}, RHSU[nx][ny+1]={0.}, RHSV[nx+1][ny]={0.};
    
    double err_u;
    int it_p=0, it_u=0;

    //Initializing u
    for (int i=0;i<=(nx-1);i++)
    {
        u[i][ny-1]=u[i][ny]=1.;
    }
    FILE *fp;
    fp=fopen("u,v,p_grid.txt","w");
    FILE *fp1;
    fp1=fopen("u_central.txt","w");
    FILE *fp2;
    fp2=fopen("v_central.txt","w");
    
    //start of iteration for MAC
    do
    {
        
        //initializing RHSU
        it_u++;
        err_u=0.;
        for(int i=0;i<=(nx-1);i++){
            for(int j=0;j<=ny;j++){
                RHSU[i][j]=0.0;
            }
        }

        //initializing RHSV
        for (int i = 0; i <= nx; i++)
        {
            for (int j = 0; j <= (ny-1); j++)
            {
                RHSV[i][j]=0.0;
            }
        }
        
        double u2, u2E, uvne, uvse, uvnw, v2, v2N;
        //RHSU

        for(int i=1; i<=(nx-2); i++)
        {
            for(int j=1; j<=(ny-1); j++)
            {
                u2=pow((u[i][j]+u[i-1][j]),2)/4.0;
                u2E=pow((u[i+1][j]+u[i][j]),2)/4.0;
                uvne=(u[i][j]+u[i][j+1])*(v[i][j]+v[i+1][j])/4.;
                uvse=(u[i][j]+u[i][j-1])*(v[i][j-1]+v[i+1][j-1])/4.;
                RHSU[i][j]=u[i][j] - (dt/dx)*(u2E-u2) - (dt/dy)*(uvne-uvse) + (dt/(dx*dx*Re))*(u[i-1][j]-2.*u[i][j]+u[i+1][j])
                     + (dt/(Re*dy*dy))*(u[i][j-1]-2.0*u[i][j]+u[i][j+1]);
            }
        }

        //RHSV
        for(int i=1; i<=(nx-1); i++){
            for(int j=1; j<=(ny-2); j++){
                uvne=(u[i][j]+u[i][j+1])*(v[i][j]+v[i+1][j])/4.;
                uvnw=(u[i-1][j]+u[i-1][j+1])*(v[i][j]+v[i-1][j])/4.;
                v2=pow((v[i][j]+v[i][j-1]),2)/4.;
                v2N=pow((v[i][j+1]+v[i][j]),2)/4.;
                RHSV[i][j]=v[i][j] - (dt/dx)*(uvne-uvnw) - (dt/dy)*(v2N-v2) + (dt/(dx*dx*Re))*(v[i-1][j]-2.0*v[i][j]+v[i+1][j])
                     + (dt/(Re*dy*dy))*(v[i][j-1]-2.*v[i][j]+v[i][j+1]);
            }
        }

        double err_p;
        it_p=0;
        
        //pressure
        do
        {
            it_p++;
            err_p=0.0;
            for(int i=1; i<=(nx-1); i++)
            {
                for(int j=1; j<=(ny-1); j++)
                {                    
                    double temp=(0.5/(1.0+b*b))*((p[i+1][j]+p[i-1][j]+b*b*(p[i][j+1]+p[i][j-1])) - (dx*dx/dt)*((RHSU[i][j]-RHSU[i-1][j])/dx
                     + (RHSV[i][j]-RHSV[i][j-1])/dy));
                    err_p+=fabs(temp-p[i][j]);
                    p[i][j]=temp;
                }
            }

            //updating pressure fictitious points on top and bottom
            for(int i=1; i<=(nx-1); i++){
                p[i][ny]=p[i][ny-1];
                p[i][0]=p[i][1];
            }

            //updating pressure fictitious points on right and left
            for(int j=0; j<=ny; j++){
                p[0][j]=p[1][j];
                p[nx][j]=p[nx-1][j];
            }

        } while (err_p>tol);
       
        //u velocity
        for(int i=1; i<=(nx-2); i++)
        {
            for(int j=1; j<=(ny-1); j++)
            {
                double temp=-(dt/dx)*(p[i+1][j]-p[i][j]) + RHSU[i][j];
                err_u+=fabs(u[i][j]-temp);
                u[i][j]=temp;
            }
        }

        //update u fictitious points on top and bottom
        for(int i=0; i<=(nx-1); i++){
            u[i][0]=-u[i][1];   u[i][ny]=2.0-u[i][ny-1];
        }

        //v velocity
        for(int i=1; i<=(nx-1); i++){
            for(int j=1; j<=(ny-2); j++){
                v[i][j]=-(dt/dy)*(p[i][j+1]-p[i][j]) + RHSV[i][j];
            }
        }

        //update v fictitious points on right and left
        for(int j=0; j<=(ny-1); j++){
            v[0][j]=-v[1][j];   v[nx][j]=-v[nx-1][j];
        }

        //print u error at every 1000 iteration
        if(it_u % 1000==0 || it_u==1)
            printf("%d: %f\n",it_u,err_u);
    } while (err_u>tol);
    
    //interpolating u,v and p at grid points
    double xpos,ypos, ug[nx][ny], vg[nx][ny], pg[nx][ny];
    for(int i=0;i<nx;i++)
    {
        for(int j=0;j<ny;j++)
        {
            ug[i][j]=0.5*(u[i][j]+u[i][j+1]);
            vg[i][j]=0.5*(v[i][j]+v[i+1][j]);
            pg[i][j]=0.25*(p[i][j]+p[i+1][j]+p[i+1][j+1]+p[i][j+1]);
        }
    }

    //printing u, v and p in output file
    fprintf(fp,"x\ty\tu\tv\tp\n");
    for(int i=0;i<nx;i++)
    {
        for(int j=0;j<ny;j++)
        {
            xpos=i*dx; ypos=j*dy;
            fprintf(fp,"%f\t%f\t%5.8f\t%5.8f\t%5.8f\n",xpos,ypos,ug[i][j],vg[i][j],pg[i][j]);
        }
    }

    //printing u velocity along central vertical line
    fprintf(fp1,"u\ty\n");
    int x=(nx-1)/2;
    for(int j=0;j<ny;j++)
    {
        ypos=(double)j*dy;
        fprintf(fp1,"%5.8f\t%f\n",ug[x][j],ypos);
    }

    //printing v velocity along central horizontal line
    fprintf(fp2,"v\tx\n");
    int y=(ny-1)/2;
    for(int i=0;i<nx;i++)
    {
        xpos=(double)i*dx;
        fprintf(fp2,"%5.8f\t%f\n",vg[i][y],xpos);
    }
    return 0;
}