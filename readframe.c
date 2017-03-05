/*no top file readed in, so assign atomic information by input.inp and cg_type by cgtype.inp*/
#include<stdio.h>
#include<stdlib.h>
#include"xdrfile.h"
#include"xdrfile_trr.h"
//#define n_cg_sites 3456
//try to modify this define for general use
int main(int argc,char **argv)
{
    int i,j,k;
    //read number of atom
    int n_atom; //total number of atoms in the traj
    if(read_trr_natoms(argv[1],&n_atom)) {printf("cannot read trr file\n");return 1;}
    printf("%d atoms read\n",n_atom);
    //read snapshots
    XDRFILE* fp_trr=xdrfile_open(argv[1],"r");
    matrix box;
    rvec *x=(rvec *)malloc(n_atom*sizeof(rvec));
    rvec *v=(rvec *)malloc(n_atom*sizeof(rvec));
    rvec *f=(rvec *)malloc(n_atom*sizeof(rvec));
    int i_step,sstep;
    sscanf(argv[2],"%d",&sstep);
    float t,lambda;
    while(!read_trr(fp_trr,n_atom,&i_step,&t,&lambda,box,x,v,f)) //loop for snapshots
    {
    //    printf("%d step\t %f time\n",i_step,t); 
    if(i_step==sstep) 
for(i=0;i<n_atom;i++)
{
     for(j=0;j<3;j++)
{
    if(j==0) printf("%d\t",i+1);
    if(j==2)
        printf("%f\t",10*x[i][j]);
    else
        printf("%f\t",10*x[i][j]);
}
    for(j=0;j<3;j++)
{
    //if(j==0) printf("%d\t",i);
    if(j==2)
        printf("%f\n",10*v[i][j]);
    else
        printf("%f\t",10*v[i][j]);
}

    for(j=0;j<3;j++)
{
    //if(j==0) printf("%d\t",i);
    if(j==2)
        printf("%f\n",f[i][j]/41.84);
    else
        printf("%f\t",f[i][j]/41.84);
}
}
    }
    xdrfile_close(fp_trr);
	return 0;
}
