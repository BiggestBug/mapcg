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
    //read input.inp
    double *weights=malloc(n_atom*sizeof(double)); //allocate weight list
    int *cg_indice=malloc(n_atom*sizeof(int));//allocate cg_indice list
    //int *cg_types=malloc(n_atom*sizeof(int)); //allocate cg_types list
    /*weights[i] gives the weight in mapping of the ith atom, cg_indice[i] indicates ith atom is mapped into cg_indice[i]th cg site, cg_types[i] indicate the cg type of cg_indice[i]th cg site*/
    for(i=0;i<n_atom;i++) { weights[i]=0.0; cg_indice[i]=0; /*cg_types[i]=0;*/ } //initialization of the three list
    FILE *fp_input=fopen(argv[2],"r");
    if(fp_input==NULL) { printf("cannot read input file\n"); return 1;}
    char line[100];
    int n_cg_sites; fgets(line,100,fp_input); sscanf(line,"%d",&n_cg_sites);//read the number of cg_sites in total
//print for debug
printf("# of cgsites %d\n",n_cg_sites);
    double weight;
    while(!feof(fp_input))
    {
        fgets(line,100,fp_input);
        //printf("%s",line);
        //sscanf(line,"%d %lf %d %d\n",&i,&weight,&j,&k);
        sscanf(line,"%d %lf %d",&i,&weight,&j);
        weights[i-1]=weight;//assign weight of the ith atom in the system
        cg_indice[i-1]=j;
        /*in input.inp atom indice starts from 1 but weights and cg_indice starts from 0*/
        //cg_types[i]=k;
    }
    fclose(fp_input);
//print for debug 
//for(i=0;i<n_atom;i++) printf("%d th atom has weight %lf and belong to cg site %d\n ",i,weights[i],cg_indice[i]);//printf("%d %lf %d %d\n",i,weights[i],cg_indice[i],cg_types[i]);
    //read cgtypes.inp
    FILE *fp_cgtypes=fopen(argv[3],"r");
    fgets(line,100,fp_cgtypes);
    int n_cg_types; sscanf(line,"%d",&n_cg_types);//read first line in cgtypes.inp to get n_cg_types 
    int *cg_types=malloc(n_cg_sites*sizeof(int));
    while(!feof(fp_cgtypes))
    {
        fgets(line,100,fp_cgtypes);
        sscanf(line,"%d %d",&i,&j);
        cg_types[i-1]=j; //ith cg site has type cg_types[i]
        /*in cgtypes.inp indice of cg sites start from 1 but in cg_types start from 0*/
    }
//print for debug
printf("# of cg_types %d\n",n_cg_types);
//for(i=0;i<n_cg_sites;i++) printf("%d th cg site had type %d\n",i,cg_types[i]);
    fclose(fp_cgtypes);
    //read snapshots
    XDRFILE* fp_trr=xdrfile_open(argv[1],"r");
    FILE *fp_cgmap=fopen(argv[4],"w");//output lammps traj
    matrix box;
    rvec *x=(rvec *)malloc(n_atom*sizeof(rvec));
    rvec *v=(rvec *)malloc(n_atom*sizeof(rvec));
    int i_step;
    float t,lambda;
    rvec *xcg=(rvec *)malloc(n_cg_sites*sizeof(rvec));
    rvec *vcg=(rvec *)malloc(n_cg_sites*sizeof(rvec));
    while(!read_trr(fp_trr,n_atom,&i_step,&t,&lambda,box,x,v,NULL)) //loop for snapshots
    {
        for(i=0;i<n_cg_sites;i++)//initialization for xcg,fcg, important!
            for(j=0;j<3;j++)
            {
                xcg[i][j]=0.0; vcg[i][j]=0.0;
            }
        printf("%d step\t %f time\n",i_step,t); 
        for(i=0;i<n_atom;i++) //loop for atom, ith atom
            for(j=0;j<3;j++) //loop for xyz
            {

//print for debug
/*if(i_step==0&&i>106636&& i<=106742) 
{
    if(j==0) printf("%d\t",i);
    if(j==2)
        printf("%f\n",x[i][j]);
    else
        printf("%f\t",x[i][j]);
}*/


                    xcg[cg_indice[i]-1][j]+=10.0*weights[i]*x[i][j];//do mapping
                    vcg[cg_indice[i]-1][j]+=10.0*weights[i]*v[i][j];//do mapping
            }
        //print lammps head
        fprintf(fp_cgmap,"ITEM: TIMESTEP\n");
        fprintf(fp_cgmap,"%d\n",i_step);
        fprintf(fp_cgmap,"ITEM: NUMBER OF ATOMS\n");
        fprintf(fp_cgmap,"%d\n",n_cg_sites);
        fprintf(fp_cgmap,"ITEM: BOX BOUNDS pp pp pp\n");
        fprintf(fp_cgmap,"0.0 %f\n0.0 %f\n0.0 %f\n",10.0*box[0][0],10.0*box[1][1],10.0*box[2][2]);
        fprintf(fp_cgmap,"ITEM: ATOMS id type xu yu zu fx fy fz\n");
        for(i=0;i<n_cg_sites;i++) //loop for cg_sites
        {
            fprintf(fp_cgmap,"%d %d %f %f %f %f %f %f\n",i+1,cg_types[i],xcg[i][0],xcg[i][1],xcg[i][2],vcg[i][0],vcg[i][1],vcg[i][2]);
        }
    }
    xdrfile_close(fp_trr);
    printf("nojump is important!!\n");
	return 0;
}
