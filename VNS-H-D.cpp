
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <cstdlib>
#include <string.h>
#include <time.h>
#include <string>
#include <fstream>
#include <iostream>
#include <vector>
#include <algorithm>
#include <sstream>

#include <cstdlib>
#include "conio.h"
#include <time.h>
#include <fstream>
#include <iomanip>

using namespace std;

	int i,j,k,p;
	clock_t t1;

	int run;
	int Nv,Np;
	int *IDp,*IDv;

	double bestcost=0,cost,cost1;
	double bestt;
	int *usedPM,*PM;
	int * VM;
	
	
	int voisinage=0;
	int decomp;int recherche;
	int count_neig;
	int*disk;int *b;int *used_phys_disk,*y;
	int *tabp,*tabv;
	int * listeVm, *listePm;
	
	int *X,*X1,*Xbest;
	int **Y,**Y1,**Ybest;
	int *n_vm_par_pm,*n_vm_par_pm1,*n_vm_par_pm_best;
	double *res_CPU_phy,*res_CPU_phy1,*res_CPU_phy_best;
	double *res_MEM_phy,*res_MEM_phy1,*res_MEM_phy_best;
	double **res_DISK_phy,**res_DISK_phy1,**res_DISK_phy_best;
	
	int* R;
	int c[2];int currentshaking;
	int v1,v2,p1,p2;
	string instance;
	int index;
	int cpt=0;
	int*yy1,*yy2;

	float initcost;

	double CPU_phy[15]={8,8,8,8,16,16,16,16,16,32,48,64,80,120,120};
	double MEM_phy[15]={16,32,64,64,32,64,128,256,256,256,512,1024,2048,4096,4096};
	int nDISK_phy[15]={1,1,2,4,2,4,4,8,16,4,8,4,16,4,24};
	double COST_phy[15]={100,120,200,300,600,700,900,1500,1800,2500,3500,5000,7000,9000,1200};

	double DISK_phy[15][24]={{256,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
							{512,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
							{512,512,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
							{512,512,512,512,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
							{512,512,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
							{512,512,512,512,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
							{1000,1000,1000,1000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
							{1000,1000,1000,1000,1000,1000,1000,1000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
							{512,512,512,512,512,512,512,512,512,512,512,512,512,512,512,512,0,0,0,0,0,0,0,0},
							{1000,1000,1000,1000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
							{1000,1000,1000,1000,1000,1000,1000,1000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
							{1000,1000,1000,1000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
							{1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,0,0,0,0,0,0,0,0},
							{1000,1000,1000,1000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
							{1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600},
	};
							

	double CPU_vir[18]={1,2,4,8,2,4,8,16,32,2,4,8,16,32,4,8,16,32};
	double MEM_vir[18]={3.75,7.5,15,30,3.75,7.5,15,30,60,15.25,30.5,61,122,244,30.5,61,122,244};
	double REVENE_vir[18]={0.067,0.133,0.266,0.532,0.105,0.21,0.42,0.84,1.68,0.166,0.333,0.665,1.33,2.66,0.853,1.705,3.41,6.82};
	int nDISK_vir[18]={1,1,2,2,2,2,2,2,2,1,1,1,1,2,1,2,4,8};
	double DISK_vir[18][8]={
		{4,0,0,0,0,0,0,0},{32,0,0,0,0,0,0,0},{40,40,0,0,0,0,0,0},
		{80,80,0,0,0,0,0,0},{16,16,0,0,0,0,0,0},{40,40,0,0,0,0,0,0},
		{80,80,0,0,0,0,0,0},{160,160,0,0,0,0,0,0},{320,320,0,0,0,0,0,0},
		{32,0,0,0,0,0,0,0},{80,0,0,0,0,0,0,0},{160,0,0,0,0,0,0,0},
		{320,0,0,0,0,0,0,0},{320,320,0,0,0,0,0,0},{800,0,0,0,0,0,0,0},
		{800,800,0,0,0,0,0,0},{800,800,800,800,0,0,0,0},{800,800,800,800,800,800,800,800}
	};
	
	



double CPUp(int pm){	return CPU_phy[IDp[pm]];}
double COSTp(int pm){return COST_phy[IDp[pm]];}
double COSTUNITp(int pm){return (COST_phy[IDp[pm]]/CPUp(pm));}
double MEMp(int pm){	return MEM_phy[IDp[pm]];}
int nDISKp(int pm){	return nDISK_phy[IDp[pm]];}
double DISKp(int pm, int nd){ return DISK_phy[IDp[pm]][nd];}

double CPUv(int vm){	return CPU_vir[IDv[vm]];}
double MEMv(int vm){	return MEM_vir[IDv[vm]];}
double REVENUEv(int vm){return REVENE_vir[IDv[vm]];}
int nDISKv(int vm){	return nDISK_vir[IDv[vm]];}
double DISKv(int vm,int nd){ return DISK_vir[IDv[vm]][nd];}


double utilization(int pm){

	double x1,x2,x3=1;
	x1=(res_CPU_phy[pm])/CPUp(pm);
	x2=(res_MEM_phy[pm])/MEMp(pm);

	return x1*x2*x3*1000;
}
double utilization_add_vm(int vm,int pm){
	
	double x1,x2,x3=1;
	x1=(res_CPU_phy[pm]-CPUv(vm))/CPUp(pm);
	x2=(res_MEM_phy[pm]-MEMv(vm))/MEMp(pm);
	return x1*x2*x3*1000;
}
double utilization_permute(int pm,int vm1,int vm2){
	double x1,x2,x3=1;
	x1=(res_CPU_phy[pm]+CPUv(vm1)-CPUv(vm2))/CPUp(pm);
	x2=(res_MEM_phy[pm]+MEMv(vm1)-MEMv(vm2))/MEMp(pm);
	return x1*x2*x3*1000;
}

void merge_VM(int,int,int);
void merge_sort_VM(int low,int high)
{
 int mid;
 if(low<high)
 {
  mid = low + (high-low)/2; //This avoids overflow when low, high are too large
  merge_sort_VM(low,mid);
  merge_sort_VM(mid+1,high);
  merge_VM(low,mid,high);
 }
}
void merge_VM(int low,int mid,int high)
{
 int h,i,j,k;

 h=low;
 i=low;
 j=mid+1;

 while((h<=mid)&&(j<=high))
 {
	 if(REVENUEv(VM[h])>=REVENUEv(VM[j]))
  {
   b[i]=VM[h];
   h++;
  }
  else
  {
   b[i]=VM[j];
   j++;
 }
  i++;
 }
 if(h>mid)
 {
  for(k=j;k<=high;k++)
  {
   b[i]=VM[k];
   i++;
  }
 }
 else
 {
  for(k=h;k<=mid;k++)
  {
   b[i]=VM[k];
   i++;
  }
 }
 for(k=low;k<=high;k++) VM[k]=b[k];

}

void merge_cost_PM(int low,int mid,int high)
{
 int h,i,j,k;
 
 h=low;
 i=low;
 j=mid+1;

 while((h<=mid)&&(j<=high)) {
	 if(COSTp(PM[h])<=COSTp(PM[j])){
	   b[i]=PM[h];
	   h++;
	  }
	  else{
	   b[i]=PM[j];
	   j++;
	 }
  i++;
 }
 if(h>mid)
  for(k=j;k<=high;k++){
   b[i]=PM[k];
   i++;
  }
 else
  for(k=h;k<=mid;k++){
   b[i]=PM[k];
   i++;
  }
 
 for(k=low;k<=high;k++) 
	 PM[k]=b[k];
 
}
void merge_sort_cost_PM(int low,int high)
{
 int mid;
 if(low<high){
  mid = low + (high-low)/2; //This avoids overflow when low, high are too large
  merge_sort_cost_PM(low,mid);
  merge_sort_cost_PM(mid+1,high);
  merge_cost_PM(low,mid,high);
 }
}

void merge_disk(int,int,int,int);
void merge_sort_disk(int low,int high,int pm)
{
 int mid;
 if(low<high)
 {
  mid = low + (high-low)/2; //This avoids overflow when low, high are too large
  merge_sort_disk(low,mid,pm);
  merge_sort_disk(mid+1,high,pm);
  merge_disk(low,mid,high,pm);
 }
}
void merge_disk(int low,int mid,int high,int pm)
{
 int h,i,j,k;
  

 h=low;
 i=low;
 j=mid+1;
 if(b==NULL) 
    {
            printf("b valeur inexistante\n");
    }
 while((h<=mid)&&(j<=high))
 {
	 if(res_DISK_phy[pm][disk[h]]>=res_DISK_phy[pm][disk[j]])
		  {
		   b[i]=disk[h];
		   h++;
		  }
	  else
		  {
		   b[i]=disk[j];
		   j++;
		 }
      i++;
 }
 if(h>mid)
 {
  for(k=j;k<=high;k++)
  {
   b[i]=disk[k];
   i++;
  }
 }
 else
 {
  for(k=h;k<=mid;k++)
  {
   b[i]=disk[k];
   i++;
  }
 }
 for(k=low;k<=high;k++) disk[k]=b[k];
 
}

void merge_PM(int low,int mid,int high)
{
 int h,i,j,k;
 h=low;
 i=low;
 j=mid+1;

 while((h<=mid)&&(j<=high)) {
if(COSTp(listePm[h])>=COSTp(listePm[j])){
	   b[i]=listePm[h];
	   h++;
	  }
	  else{
	   b[i]=listePm[j];
	   j++;
	 }
  i++;
 }
 if(h>mid)
  for(k=j;k<=high;k++){
   b[i]=listePm[k];
   i++;
  }
 else
  for(k=h;k<=mid;k++){
   b[i]=listePm[k];
   i++;
  }
 
 for(k=low;k<=high;k++) 
	 listePm[k]=b[k];
}

void merge_sort_PM(int low,int high)
{
 int mid;
 if(low<high){
  mid = low + (high-low)/2; //This avoids overflow when low, high are too large
  merge_sort_PM(low,mid);
  merge_sort_PM(mid+1,high);
  merge_PM(low,mid,high);
 }
}

void merge_utilization_PM(int low,int mid,int high)
{
 int h,i,j,k;

 h=low;
 i=low;
 j=mid+1;

 while((h<=mid)&&(j<=high)) {
	 if(utilization(listePm[h])>=utilization(listePm[j])){
	   b[i]=listePm[h];
	   h++;
	  }
	  else{
	   b[i]=listePm[j];
	   j++;
	 }
  i++;
 }
 if(h>mid)
  for(k=j;k<=high;k++){
   b[i]=listePm[k];
   i++;
  }
 else
  for(k=h;k<=mid;k++){
   b[i]=listePm[k];
   i++;
  }
 
 for(k=low;k<=high;k++) 
	 listePm[k]=b[k];

}
void merge_sort_utilization_PM(int low,int high)
{
 int mid;
 if(low<high){
  mid = low + (high-low)/2; //This avoids overflow when low, high are too large
  merge_sort_utilization_PM(low,mid);
  merge_sort_utilization_PM(mid+1,high);
  merge_utilization_PM(low,mid,high);
 }
}

double calcul_cos_angle(int pm){
double a=0,b=0,s=0;
if(n_vm_par_pm[pm]==0) return 1;
	s=(CPUp(pm)-res_CPU_phy[pm])/CPUp(pm);
	a+=s;
	b+=s*s;
	s=(MEMp(pm)-res_MEM_phy[pm])/MEMp(pm);
	a+=s;
	b+=s*s;
	return((a)/(sqrt(b)*1,414));
}

double calcul_cos_angle_permute(int pm,int vm1,int vm2){
	//vm1 is embedded on pm
	if(n_vm_par_pm[pm]==0) return 1;
	double a=0,b=0,s=0;
	s=(CPUp(pm)-res_CPU_phy[pm]-CPUv(vm1)+CPUv(vm2))/CPUp(pm);
	a+=s;
	b+=s*s;
	s=(MEMp(pm)-res_MEM_phy[pm]-MEMv(vm1)+MEMv(vm2))/MEMp(pm);
	a+=s;
	b+=s*s;

	return(a/(sqrt(b)*1,414));

}


void lecture_data_center(string fichier){
		
	fstream fdc(fichier, ios::in); 
	if(!fdc) printf("Error opening Data_center.txt file...\n");
		
	fdc >> Np;
	
	IDp=(int*) malloc (sizeof(int)*Np);

	for(i=0; i<Np;i++){
		fdc>>IDp[i];
		
	}

	fdc.close();
}
void lecture_requette(string fichier){
	
	fstream fr(fichier, ios::in);  
	if(!fr) printf("Error opening %s.txt file...\n",fichier);
	
	fr>>Nv;

	IDv=(int*) malloc (sizeof(int)*Nv);

	for(i=0; i<Nv;i++){
		fr>>IDv[i];
	}

	fr.close();	
}
int* embedd_PM_VM(int pm, int vm){

	if((res_CPU_phy[pm]<CPUv(vm))||(res_MEM_phy[pm]<MEMv(vm))||(nDISKp(pm)<nDISKv(vm))) return nullptr;
	else{
		
		for(k=0;k<nDISKp(pm);k++) {disk[k]=k;used_phys_disk[k]=0;}

		int n=0;
		merge_sort_disk(0,nDISKp(pm)-1,pm);

		for(p=0;p<nDISKv(vm);p++){
			
			for(k=0;k<nDISKp(pm);k++){
				if((used_phys_disk[disk[k]]==0)&&(res_DISK_phy[pm][disk[k]]>=DISKv(vm,p))){
					y[p]=disk[k];
					used_phys_disk[disk[k]]=1;
					n++;
					
					break;
				}
			}
			
		}
		
		if(n<nDISKv(vm)) { return nullptr;}
		else{
			res_CPU_phy[pm]-=CPUv(vm);
			if(res_CPU_phy[pm]<0) printf("..cpu..\n");
			res_MEM_phy[pm]-=MEMv(vm);
			X[vm]=pm;
			
			for(int p=0;p<nDISKv(vm);p++){
				Y[vm][p]=y[p];
				res_DISK_phy[pm][y[p]]-=DISKv(vm,p);
			}		
			
			return y;
		}
	
	}
	
	return nullptr;
}
bool prim_test_embed_PM_VM(int i , int j){
	if((res_CPU_phy[i]<CPUv(j))||(res_MEM_phy[i]<MEMv(j))||(nDISKp(i)<nDISKv(j))){ return false;}
	else return true;
}
int* test_embedd_PM_VM(int pm, int vm){

	if((res_CPU_phy[pm]<CPUv(vm))||(res_MEM_phy[pm]<MEMv(vm))||(nDISKp(pm)<nDISKv(vm))){ return nullptr;}
	else{
		
		for(k=0;k<nDISKp(pm);k++) {disk[k]=k;used_phys_disk[k]=0;}
		
		int n=0;
		merge_sort_disk(0,nDISKp(pm)-1,pm);
		for(p=0;p<nDISKv(vm);p++){
			for(k=0;k<nDISKp(pm);k++){
				
				if((used_phys_disk[disk[k]]==0)&&(res_DISK_phy[pm][disk[k]]>=DISKv(vm,p))){
					
					used_phys_disk[disk[k]]=1;
					n++;
					y[p]=disk[k];
					break;
				}
			}
			
		}

		if(n<nDISKv(vm)) return nullptr;
		else return y;
		
	
	}
	
	return nullptr;
}
void remove_PM_VM(int pm, int vm, int*y){
	if(X[vm]==pm){
	res_CPU_phy[pm]+=CPUv(vm);
	res_MEM_phy[pm]+=MEMv(vm);
	for(int p=0;p<nDISKv(vm);p++)
		res_DISK_phy[pm][y[p]]+=DISKv(vm,p);	
	}else printf("\n------->erreur remove\n");
		
}
void add_PM_VM(int pm, int vm, int*y){
		res_CPU_phy[pm]-=CPUv(vm);
		res_MEM_phy[pm]-=MEMv(vm);
		X[vm]=pm;
			
		for(int p=0;p<nDISKv(vm);p++){
				Y[vm][p]=y[p];
				res_DISK_phy[pm][y[p]]-=DISKv(vm,p);
		}	
		
	
}


double calcul_cost(double c,int pm1,int pm2){
	//remove vm from pm1
	//add vm to pm2
	if(n_vm_par_pm[pm2]==0) c=c+COSTp(pm2);
	if(n_vm_par_pm[pm1]-1==0) c=c-COSTp(pm1);
	
	return c;
}
double calcul_cost_remove(double c,int pm){
	//remove vm from pm
	if(n_vm_par_pm[pm]-1==0) c=c-COSTp(pm);
	return c;
}
double calcul_cost_add(double c,int pm){
	//add vm to pm
	if(n_vm_par_pm[pm]==0) c=c+COSTp(pm);
	return c;
}

bool test_move_PM(int pm1,int pm2){

	bool pos=true;
	for(int k=0;k<Nv;k++){
		if(X[k]==pm1){
			int *y21=test_embedd_PM_VM(pm2,k);
			
			if(y21==nullptr){ pos=false;break;}
			else{
				res_CPU_phy[pm2]-=CPUv(k);
				res_MEM_phy[pm2]-=MEMv(k);
			
				for(int p=0;p<nDISKv(k);p++){
						res_DISK_phy[pm2][y21[p]]-=DISKv(k,p);
				}
			}

		}
	}
	res_CPU_phy[pm2]=CPUp(pm2);
	res_MEM_phy[pm2]=MEMp(pm2);
	for(int p=0;p<nDISKp(pm2);p++){
		res_DISK_phy[pm2][p]=DISKp(pm2,p);
	}
	return pos;
}

//----------------------------------------------------------------
int count_neighboor_structure(){
	int cn=decomp;
	int x=1;
	while(x<decomp){
		for(int cpt=0;cpt<decomp-x;cpt++)
			cn++;
		x++;
	}

	
	return cn;

}
void generate_table(){
	tabv=(int*) calloc ((decomp+1),sizeof(int));
	tabp=(int*) calloc ((decomp+1),sizeof(int));
	
	for(int cpt=0;cpt<decomp+1;cpt++){
		tabv[cpt]=(cpt*Nv)/decomp;
		tabp[cpt]=(cpt*Np)/decomp;

	}tabp[decomp]=Np-1;
	tabv[decomp]=Nv;

}
void init_recherche(){

	yy1=(int*) calloc(8,sizeof(int));
	yy2=(int*) calloc(8,sizeof(int));
	R= (int*) calloc (Nv,sizeof(int));
	c[0]=1;c[1]=1;
	count_neig=count_neighboor_structure();
	generate_table();

	
	listeVm=(int*) calloc (Nv,sizeof(int));
	listePm=(int*) calloc (Np,sizeof(int));
	for(i=0;i<Np;i++) 
		listePm[i]=i;

	n_vm_par_pm=(int*) calloc (Np,sizeof(int));
	for(k=0;k<Nv;k++){
		n_vm_par_pm[X[k]]++;
	}

	Xbest=(int*) calloc (Nv,sizeof(int));

	Ybest=(int**) calloc (Nv,sizeof(int*));
	for(k=0;k<Nv;k++){
		Ybest[k]=(int*) calloc (nDISKv(k),sizeof(int));
	}
	

	for (int p=0;p<Nv;p++)
		Xbest[p]=X[p];

	for (int k=0;k<Nv;k++){
		for (int p=0;p<nDISKv(k);p++)
		Ybest[k][p]=Y[k][p];
	}

	n_vm_par_pm_best=(int*) calloc (Np,sizeof(int));
	for(int k=0;k<Np;k++){
		n_vm_par_pm_best[k]=n_vm_par_pm[k];
	}

	
	
	cost=bestcost;	
	res_CPU_phy_best=(double*) malloc (sizeof(double)*Np);
	res_MEM_phy_best=(double*) malloc (sizeof(double)*Np);
	res_DISK_phy_best=(double**) malloc (sizeof(double*)*Np);
	for(i=0;i<Np;i++)
		res_DISK_phy_best[i]=(double*) malloc (sizeof(double)*nDISKp(i));
	for(i=0;i<Np;i++){
		res_CPU_phy_best[i]=res_CPU_phy[i];
		res_MEM_phy_best[i]=res_MEM_phy[i];
		for(int k=0;k<nDISKp(i);k++)
			res_DISK_phy_best[i][k]=res_DISK_phy[i][k];
	}
	
}
void init_solution(){
	
	
	b=(int*) calloc(Np+Nv,sizeof(int));
	used_phys_disk=(int*) calloc (24,sizeof(int));
	disk=(int*) calloc (24,sizeof(int));
	y=(int*) calloc (8,sizeof(int));

	X=(int*) calloc (Nv,sizeof(int));

	Y=(int**) calloc (Nv,sizeof(int*));
	for(k=0;k<Nv;k++){
		Y[k]=(int*) calloc (8,sizeof(int));
	}
		
	res_CPU_phy=(double*) malloc (sizeof(double)*Np);
	res_MEM_phy=(double*) malloc (sizeof(double)*Np);
	res_DISK_phy=(double**) malloc (sizeof(double*)*Np);
	for(i=0;i<Np;i++)
		res_DISK_phy[i]=(double*) malloc (sizeof(double)*nDISKp(i));
	
	usedPM=(int*) calloc (Np,sizeof(int));
	PM=(int*) malloc (sizeof(int)*Np);
	VM=(int*) calloc (Nv,sizeof(int));
	for(i=0;i<Np;i++) 
		PM[i]=i;
	for(k=0;k<Nv;k++)
		VM[k]=k;
	
	for(i=0;i<Np;i++){
		res_CPU_phy[i]=CPUp(i);
		res_MEM_phy[i]=MEMp(i);
		for(k=0;k<nDISKp(i);k++)
			res_DISK_phy[i][k]=DISKp(i,k);
	}
	

}
	
void calcul_bornes(int val,int r){
	int x=0;
		int va=val;
		while(va-decomp>=1){
			va=va-decomp;
			x++;
		}
	if(r==1){
		
		p1=tabp[decomp-x-1];p2=tabp[decomp-x];
		v1=tabv[va-1];v2=tabv[va];
	}else if(r==2){
		v1=tabv[decomp-va];v2=tabv[decomp-va+1];
		p1=tabp[decomp-va];p2=tabp[decomp-va+1];
	}else if(r==3){
		p1=tabp[decomp-va];p2=tabp[decomp-va+1];
	}
}

bool test_permute_PM_VM(int vm1,int vm2){
	
	int pm1=X[vm1];
	int pm2=X[vm2];
	
	for(k=0;k<nDISKv(vm1);k++) yy1[k]=Y[vm1][k];
	for(k=0;k<nDISKv(vm2);k++) yy2[k]=Y[vm2][k];
	

	remove_PM_VM(pm1,vm1,yy1);
	remove_PM_VM(pm2,vm2,yy2);

	int *z1= test_embedd_PM_VM(pm1,vm2);
	int *z2= test_embedd_PM_VM(pm2,vm1);

	add_PM_VM(pm1,vm1,yy1);
	add_PM_VM(pm2,vm2,yy2);

	if((z1==nullptr)||(z2==nullptr)) return false;
	else{
		return true;
	}
	

}


int nbR;

void recherche_local_N1(){
cpt=0;
		
	for(i=0;i<Np;i++) 
	listePm[i]=i;
	
	merge_sort_PM(0,Np-1);
	for(i=0;i<Np;i++)
		for(k=0;k<Nv;k++)
			if(X[k]==listePm[i]){listeVm[cpt]=k;cpt++;}
		
	int index=0;
	double c1;
	int pm1;
	int vm1;
	
	for(int i=v1;i<v2;i++){
		index=0;
		
		double a1=utilization(X[listeVm[i]]);
		for(int j=p2;j>=p1;j--){
			
			if((test_embedd_PM_VM(listePm[j],listeVm[i])!=nullptr)&&(listePm[j]!=X[listeVm[i]])){
				double a11=utilization_add_vm(listeVm[i],listePm[j]);
				
				//if(a11<a1)				
					if((index==0)||(a11<c1)){
						c1=a11;
						pm1=j;
						vm1=i;
						index++;
					}
			}	
		}
			
	
		if(index!=0){
			cost=calcul_cost_add(cost,listePm[pm1]);
			cost=calcul_cost_remove(cost,X[listeVm[vm1]]);
			
			n_vm_par_pm[X[listeVm[vm1]]]--;
			n_vm_par_pm[listePm[pm1]]++;
		
			remove_PM_VM(X[listeVm[vm1]],listeVm[vm1],Y[listeVm[vm1]]);
			embedd_PM_VM(listePm[pm1],listeVm[vm1]);
			
		}
	}	

	
}
void recherche_local_N2(){

	for(i=0;i<Np;i++) 
		listePm[i]=i;
	
	cpt=0;

	for(i=p2;i>=p1;i--)
		for(k=0;k<Nv;k++)
			if(X[k]==listePm[i]){listeVm[cpt]=k;cpt++;}

	
	double c1;
	int vm1,vm2;

	for(int i=0;i<cpt;i++){
		index=0;
		double a1=calcul_cos_angle(X[listeVm[i]]);
		for(int j=0;j<cpt;j++){
			double a2=calcul_cos_angle(X[listeVm[j]]);
			if((X[listeVm[j]]!=X[listeVm[i]])){
				
			if(test_permute_PM_VM(listeVm[j],listeVm[i])){
				
				double a11=calcul_cos_angle_permute(X[listeVm[i]],listeVm[i],listeVm[j]);
				
				if((a11>a1)&&(a11>a2)){
					
					if((index==0)||(a11-a1>c1)){
						
						c1=a11-a1;
						vm1=listeVm[j];
						vm2=listeVm[i];
						index++;
					}
				}
			}
			}
			
		}
	
		if(index!=0){
			remove_PM_VM(X[vm2],vm2,Y[vm2]);
			remove_PM_VM(X[vm1],vm1,Y[vm1]);
			
			int pm1=X[vm1];
			int pm2=X[vm2];
			
			embedd_PM_VM(pm1,vm2);
			embedd_PM_VM(pm2,vm1);
			
		}
	}
	
	
	
}
void recherche_local_N3(){
	//remove PM far from being densly packed
	for(i=0;i<Np;i++) 
		listePm[i]=i;
	merge_sort_utilization_PM(0,Np-1);

	cpt=0;
	for(i=0;i<Np;i++)
		for(k=0;k<Nv;k++)
			if(X[k]==listePm[i]){listeVm[cpt]=k;cpt++;}
	merge_sort_PM(0,Np-1);

	double c1;
	int pm1,vm1;
	for(int i=0;i<0.2*Nv;i++){
		
		index=0;
		for(int j=p2;j>=p1;j--){//listePm[j] destination (target)
			
				if((test_embedd_PM_VM(listePm[j],listeVm[i])!=nullptr)&&(listePm[j]!=X[listeVm[i]])){					
					if(n_vm_par_pm[listePm[j]]>0){
							vm1=listeVm[i];
							pm1=listePm[j];
							index++;
							break;
					}
				}
			
			
		}
		if(index!=0){
			
			cost=calcul_cost_add(cost,pm1);
			cost=calcul_cost_remove(cost,X[vm1]);
			
			n_vm_par_pm[X[vm1]]--;
			n_vm_par_pm[pm1]++;
		
			remove_PM_VM(X[vm1],vm1,Y[vm1]);
			embedd_PM_VM(pm1,vm1);
			
		}



	}
	
}
void recherche_local_N4(){
	
	//tries to change the corresponding pm types in order to use cheaper pm if possible
	int cpt=0;
	for(int k=0;k<Np;k++){
		if(n_vm_par_pm[k]>0){
			listePm[cpt]=k;
			cpt++;
		}
	}
	merge_sort_PM(0,cpt-1);
	
	for(int i=cpt-1;i>=0;i--){
		
		int c1=0;
		for(int k=0;k<Np;k++){
			if(n_vm_par_pm[k]==0){
				PM[c1]=k;
				c1++;
			}
		}
		merge_sort_cost_PM(0,c1-1);

		for(int j=0;j<c1;j++){
			if((COSTp(PM[j])<COSTp(listePm[i]))&&(CPUp(PM[j])>=CPUp(listePm[i])-res_CPU_phy[listePm[i]])&&(MEMp(PM[j])>=MEMp(listePm[i])-res_MEM_phy[listePm[i]])){
				
				bool pos= test_move_PM(listePm[i],PM[j]);
				if(pos==true){ 
					for(int k=0;k<Nv;k++){
						if(X[k]==listePm[i]){
							remove_PM_VM(listePm[i],k,Y[k]);
							embedd_PM_VM(PM[j],k);
						}
					}
					n_vm_par_pm[PM[j]]=n_vm_par_pm[listePm[i]];
					n_vm_par_pm[listePm[i]]=0;
					
					cost-=COSTp(listePm[i]);
					cost+=COSTp(PM[j]);
					break;
				}
			}
		}
	}


}

void no_update_best(){
	cost=bestcost;
	for (int p=0;p<Nv;p++)
		X[p]=Xbest[p];

	for (int k=0;k<Nv;k++){
		for (int p=0;p<nDISKv(k);p++)
			Y[k][p]=Ybest[k][p];
	}
	for(int k=0;k<Np;k++)
		n_vm_par_pm[k]=n_vm_par_pm_best[k];
			
	for(int p=0;p<Np;p++){
		res_CPU_phy[p]=res_CPU_phy_best[p];
		res_MEM_phy[p]=res_MEM_phy_best[p];
		for(int k=0;k<nDISKp(p);k++)
			res_DISK_phy[p][k]=res_DISK_phy_best[p][k];
	}
			
	
}
void update_best(){
	
			
	printf("\n--------------UPDATE------------\n");	
			
	bestcost=cost;
	for (int p=0;p<Nv;p++)
		Xbest[p]=X[p];

	for (int k=0;k<Nv;k++)
		for (int p=0;p<nDISKv(k);p++)
			Ybest[k][p]=Y[k][p];
	for(int k=0;k<Np;k++)
		n_vm_par_pm_best[k]=n_vm_par_pm[k];
			
	for(int p=0;p<Np;p++){
		res_CPU_phy_best[p]=res_CPU_phy[p];
		res_MEM_phy_best[p]=res_MEM_phy[p];
		for(int k=0;k<nDISKp(p);k++)
			res_DISK_phy_best[p][k]=res_DISK_phy[p][k];
	}
			
}

void update1(){
	cost1=cost;
	for (int p=0;p<Nv;p++)
		X1[p]=X[p];

	for (int k=0;k<Nv;k++){
		for (int p=0;p<nDISKv(k);p++)
			Y1[k][p]=Y[k][p];
	}
	for(int k=0;k<Np;k++)
		n_vm_par_pm1[k]=n_vm_par_pm[k];
			
	for(int p=0;p<Np;p++){
		res_CPU_phy1[p]=res_CPU_phy[p];
		res_MEM_phy1[p]=res_MEM_phy[p];
		for(int k=0;k<nDISKp(p);k++)
			res_DISK_phy1[p][k]=res_DISK_phy[p][k];
	}
		
}
void no_update1(){
	cost=cost1;
	for (int p=0;p<Nv;p++)
		X[p]=X1[p];

	for (int k=0;k<Nv;k++){
		for (int p=0;p<nDISKv(k);p++)
			Y[k][p]=Y1[k][p];
	}
	for(int k=0;k<Np;k++)
		n_vm_par_pm[k]=n_vm_par_pm1[k];
			
	for(int p=0;p<Np;p++){
		res_CPU_phy[p]=res_CPU_phy1[p];
		res_MEM_phy[p]=res_MEM_phy1[p];
		for(int k=0;k<nDISKp(p);k++)
			res_DISK_phy[p][k]=res_DISK_phy1[p][k];
	}
	
}

int* removal1(){
	bool removed;
	int r;
	for(int i=0;i<nbR;i++){
		do{
			r=rand()%(Nv);
			removed=false;
			for(int j=0; j<i;j++){
				if(R[j]==r) removed=true;
			}
		}while(removed==true);
		
		remove_PM_VM(X[r],r,Y[r]);
		cost=calcul_cost_remove(cost,X[r]);
		n_vm_par_pm[X[r]]--;	
		R[i]=r;		
	}
	return R;
}
bool insertion(int * R){
	bool embedd;int cpt=0;
	for(int i=0;i<Np;i++){
		if(n_vm_par_pm[i]==0) {PM[cpt]=i; cpt++;}
	}

	for(int i=0;i<nbR;i++){
		embedd=false;
		for(int j=0; j<cpt;j++){
			if(test_embedd_PM_VM(PM[j],R[i])!=nullptr){
				embedd=true;
				
				cost=calcul_cost_add(cost,PM[j]);
				n_vm_par_pm[PM[j]]++;
				embedd_PM_VM(PM[j],R[i]);
				break;
			}
		
		}
		if(embedd==false) {			
			int t;int ccc=0;
			do{
				t= rand()%Np;
				ccc++;
			}while((test_embedd_PM_VM(t,R[i])==nullptr)&&(ccc<Np*10));
			
			if(test_embedd_PM_VM(t,R[i])!=nullptr){			
				cost=calcul_cost_add(cost,t);
				n_vm_par_pm[t]++;
				embedd_PM_VM(t,R[i]);
			}else{
				no_update_best();
				return false;
			}
		}
	}
	return true;

}
void shaking(int ck)
{
	int N=0;
	for(int k=0;k<Np;k++)
		if(n_vm_par_pm[k]!=0) N++;
	
	if(ck==1)
	{
		int nb=N*0.15;
		int tab[200];
		for(k=0;k<Np;k++)
			PM[k]=k;
		for(int k=1;k<=nb;k++)
		{
			int nb_vm=0,nb_vm2=0;
			//select random used PM
			int r,r2;
			do
			{
				r= rand() % (Np);
			}while(n_vm_par_pm[r]==0) ;
			//select target unused PM
			int t;
			do
			{
				t= rand() % (Np);
			}while(n_vm_par_pm[t]!=0) ;
			//remove VMs from PM r
			for(int v=0;v<Nv;v++)
			{
				if(X[v]==r)
				{ 
					remove_PM_VM(X[v],v,Y[v]);
					cost=calcul_cost_remove(cost,r);
					n_vm_par_pm[r]--;
					tab[nb_vm]=v;
					nb_vm++;
				}
			}			
			//try to insert VMs on target PM
			for(int v=0;v<nb_vm;v++)
			{
				if(embedd_PM_VM(t,tab[v])!=nullptr)
				{
					//printf("migrate vm %d from %d to pm %d\n",tab[v],r,t );
					cost=calcul_cost_add(cost,t);
					n_vm_par_pm[t]++;
				}
				else
				{
					VM[nb_vm2]=tab[v];
					nb_vm2++;
				}
			}
			//insert lefted VMs on used PMs
			bool embed=false;
			if(nb_vm2>0)
			{
				merge_sort_VM(0,nb_vm2-1);
				for(int v=0;v<nb_vm2;v++)
				{
					embed=false;
					for(int j=0;j<Np;j++)
					{
						if(n_vm_par_pm[PM[j]]!=0)
							if(embedd_PM_VM(PM[j],VM[v])!=nullptr)
							{
								//printf("migrate vm %d from %d to pm %d\n",VM[v],r,PM[j]);
								cost=calcul_cost_add(cost,PM[j]);
								n_vm_par_pm[PM[j]]++;
								embed=true;
								break;
							}
					}
					while(embed==false)
					{
						for(int j=0;j<Np;j++)
						{
							if(n_vm_par_pm[PM[j]]==0)
							if(embedd_PM_VM(PM[j],VM[v])!=nullptr)
							{
								//printf("migrate vm %d from %d to pm %d\n",VM[v],r,PM[j]);
								cost=calcul_cost_add(cost,PM[j]);
								n_vm_par_pm[PM[j]]++;
								embed=true;
								break;
							}
						}
					}					
				}
			}
		}
	}
	else if(ck==2)
	{
		int nb=N*0.05;
		int tab[1000];
		int nb_vm=0;
		//select random PM to deleat
		for(int k=1;k<=nb;k++)
		{
			int r;
			do
			{
				r= rand() % (Np);
			}while(n_vm_par_pm[r]==0) ;
			for(int v=0;v<Nv;v++){
				if(X[v]==r) 
				{
					//printf("remove vm %d from %d\n",v,X[v],t);
					remove_PM_VM(X[v],v,Y[v]);
					cost=calcul_cost_remove(cost,X[v]);
					n_vm_par_pm[X[v]]--;
					tab[nb_vm]=v;
					nb_vm++;
				}
			}
		}
		//assign the VMs to new PMs randomly
		for(int v=0;v<nb_vm;v++)
		{
			//select target unused PM
			int in=0;int t;
			do
			{
				t= rand() % (Np);in++;
				if (n_vm_par_pm[t]==0)
				{
					if(embedd_PM_VM(t,tab[v])!=nullptr)
					{
						//printf("migrate vm %d from %d to pm %d\n",tab[v],X[tab[v]],t);
						cost=calcul_cost_add(cost,t);
						n_vm_par_pm[t]++;
						break;
					}
				}
			}while(in<Np*50);
			if(in>=Np*50)
			{ 
				bool embed=false;
				for(int t=0;t<Np;t++)
				{
					if (n_vm_par_pm[t]!=0)
					if(embedd_PM_VM(t,tab[v])!=nullptr)
					{
						//printf("migrate vm %d from %d to pm %d\n",tab[v],X[tab[v]],t);
						cost=calcul_cost_add(cost,t);
						n_vm_par_pm[t]++;
						embed=true;
						break;
					}
				}
				if(embed==false) 
				{
					no_update_best();
					shaking(4);
					currentshaking=4;
					break;
				}				
			}			
		}
	}
	else
	{
		//select random vm to remove
		nbR=0.3*Nv;
		int* rem=removal1();
		insertion(rem);
	}
}
void choose_shaking()
{
	int x= rand()%(c[0]+c[1]+c[2]);
	if(x<=c[0])
    {
        shaking(1);
        currentshaking=0;
    }
	else if(x<=c[0]+c[1])
    {
        shaking(2);
        currentshaking=1;
    }
	else
    {
        shaking(3);
        currentshaking=2;
    }
}

void initial(){
	
	cost1=bestcost;
	n_vm_par_pm1=(int*) calloc (Np,sizeof(int));
	
	X1=(int*) calloc (Nv,sizeof(int));

	Y1=(int**) calloc (Nv,sizeof(int*));
	for(k=0;k<Nv;k++){
		Y1[k]=(int*) calloc (nDISKv(k),sizeof(int));
	}
	
	res_CPU_phy1=(double*) malloc (sizeof(double)*Np);
	res_MEM_phy1=(double*) malloc (sizeof(double)*Np);
	res_DISK_phy1=(double**) malloc (sizeof(double*)*Np);
	for(i=0;i<Np;i++)
		res_DISK_phy1[i]=(double*) malloc (sizeof(double)*nDISKp(i));
	
	
	
	for (int p=0;p<Nv;p++)
		X1[p]=Xbest[p];

	for (int k=0;k<Nv;k++){
		for (int p=0;p<nDISKv(k);p++)
			Y1[k][p]=Ybest[k][p];
	}
	for(int k=0;k<Np;k++)
		n_vm_par_pm1[k]=n_vm_par_pm_best[k];
			
	for(int p=0;p<Np;p++){
		res_CPU_phy1[p]=res_CPU_phy_best[p];
		res_MEM_phy1[p]=res_MEM_phy_best[p];
		for(int k=0;k<nDISKp(p);k++)
			res_DISK_phy1[p][k]=res_DISK_phy_best[p][k];
	}
		
}

bool FFD(){
	//sort pm list
	merge_sort_cost_PM(0,Np-1);
	//sort vm list
	merge_sort_VM(0,Nv-1);

	bestcost=0;
	bool embed;

	for(j=0;j<Nv;j++){
		embed=false;
			
		//try to embedd vm on a used pm
		for(i=0;i<Np;i++){
			if (usedPM[i]==1)
				if(embedd_PM_VM(PM[i],VM[j])!=nullptr){
					embed=true;
					break;
				}
		}
		//try to embed vm on an unsed pm
		if(embed==false){
				
			for(i=0;i<Np;i++){
				if (usedPM[i]==0)
					if(embedd_PM_VM(PM[i],VM[j])!=nullptr){
						embed=true;
						bestcost=bestcost+COSTp(PM[i]);
						usedPM[i]=1;
						break;
					}
			}
		}
		if(embed==false){ printf("\n%d)NO embedd VM(%d) diskv:%d!",j, VM[j],nDISKv(VM[j]));break;}
		
	}
	return embed;
}


bool verification_vesabilite(){
	for(int i=0;i<Np;i++){
		double c=0,m=0;
		for(int j=0;j<Nv;j++){
			if(X[j]==i){
				c=c+CPUv(j);
				m=m+MEMv(j);
			}
		}
		if((c>CPUp(i))||(m>MEMp(i))) return false;
	}
	return true;


}

void main()
{

	fstream fl("instances_liste_large.txt", ios::in);
	if(!fl)
        printf("Error opening file instances_liste.txt...\n");
	
	ofstream ffres("res.xls", ios::out | ios::app );
	
	
	while(!fl.eof()){
	fl>>instance;
	lecture_data_center("all_instances_large/data_center_"+instance);
	lecture_requette("all_instances_large/"+instance);
	for(run=1;run<=5;run++){
	decomp=20;
	
	
	t1 = clock();
	init_solution();
	printf("begin resolution...\n");
	bool initialSol=FFD();
	
		
	double t = (double)(clock() - t1)/(double)CLK_TCK ;
	double t_best=t;

	ffres<<instance<<"\t"<<decomp<<"\t"<<bestcost<<"\t";
	
	
	if(initialSol) {	


	init_recherche();initial();		
	int iter=0,iter_best=0,iter_tot=0;
	int val=1;
	srand(1000*run);
	
	update_best();
	
	do{
		iter++;
		iter_tot++;
	
		if(iter_tot!=1)		choose_shaking();
		
		printf("%d) best=%.0f  shaked=%.0f   ",iter,bestcost,cost);
		update1();
		
		recherche=1; val=1;
		do{
			calcul_bornes(val,recherche);
		
			recherche_local_N1();
			if(cost<cost1) update1();
			else no_update1();
			val++;
		}while(val<=decomp*decomp);
		
	
	
		recherche=2; val=1;
		do{
			calcul_bornes(val,recherche);
			recherche_local_N2();
			val++;
		}while(val<=decomp);
		update1();
		
		recherche=3; val=1;
		do{
			calcul_bornes(val,recherche);
			recherche_local_N3();
			if(cost<cost1) update1();
			else no_update1();
			val++;
		}while(val<=decomp);
	
		
		recherche_local_N4();
		update1();
			

		if(cost<bestcost) {
			update_best();
			bestt = (double)(clock() - t1)/(double)CLK_TCK ;
			if(iter_tot>0) c[currentshaking]++;
			iter_best=iter_tot;
			iter=0;
		}
		else no_update_best();
	
	
	
	}while(iter<5);
	

	
	t = (double)(clock() - t1)/(double)CLK_TCK ;		
	
	ffres<<"\t"<<bestcost<<"\t"<<bestt<<"\t"<<iter_best<<"\t"<<t<<"\t"<<iter_tot<<endl;


	for(i=0;i<Np;i++){
		free(res_DISK_phy[i]);free(res_DISK_phy1[i]);free(res_DISK_phy_best[i]);
	}
	free(res_DISK_phy);free(res_DISK_phy1);free(res_DISK_phy_best);
	for(i=0;i<Nv;i++){
		free(Y[i]);free(Y1[i]);free(Ybest[i]);
	}
	
	free(res_MEM_phy);free(res_MEM_phy1);free(res_MEM_phy_best);
	free(yy1);free(yy2);

	free(n_vm_par_pm);free(n_vm_par_pm1);free(n_vm_par_pm_best);
	free(R);
	free(tabp);free(tabv);
	free(X);free(Xbest);free(X1);
	free(listeVm);free(listePm);free(PM);free(VM);free(usedPM);


	free(b);free(used_phys_disk);free(disk);free(y);


	}
	}
	free(IDv); free(IDp);

	}

	ffres.close();

	system("pause");
	
	
}
