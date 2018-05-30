///////////////Librerias//////////////////     
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "genetic.h"
/////////////////////////////////////////


///////////////////////////Variables Globales//////////////////////////////////////
//Los rangos para cruza y mutacion solo para x2 y x3 entonces seran (10, 10*Num_var)
#define puntocruza_mutacion 10
/*Numero de variables*/
#define Num_var	3
/*Number of bits used to log2(lsup-linf) represent an individual */	
#define codesize 10*Num_var
int popsize;
int Globalpopsize;
char nombreArchivo[100];
char comando[300];
double varibles_interval[21]={-0.5, -0.45, -0.4, -0.35, -0.3, -0.25, -0.2, -0.15, -0.1, -0.05, 0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5};
///////////////////////////////////////////////////////////////////////////////////

////////////////////////Ficheros Globales///////////////////////////////////////////
FILE *ficheroGrafica;
FILE *ficheroReporte;
FILE *fp_aux;
////////////////////////////////////////////////////////////////////////////////////


void show_pop();
int main(int argc , char *argv[])
{
	struct individual *temp;
	if(argc != 7)
	{
		printf("Faltan argumentos\n");
		printf("Sintaxis: ./prog , popsize , mutation_rate , crossover_rate , num_generation , seed(0,1) , output_file\n");
	}
	else
	{
		//Captura de variables pasadas por parámetros
		//cls();
		Globalpopsize=0;				
		popsize=atoi(argv[1]);	//Size of Population		
		Pm = atof(argv[2]);		//Mutation rate		
		Pc=atof(argv[3]);		//Crossover rate		
		Gmax = atoi(argv[4]);	//Number of generation		
		Rseed = atof(argv[5]);	//Seed
		strcpy(nombreArchivo, argv[6]);		
		strcat( nombreArchivo, ".txt" );
		/////////////////////////////////////////////
		ficheroGrafica = fopen("grafica.txt", "w" );
		ficheroReporte = fopen(nombreArchivo, "w" );
		/////////////Memory Allocation///////////////
		memory_allocation();	
		poblacion_evaluada();
		/////////////////////////////////////////////
		/* Initialize random number generator */	
		randomize();
		//////////////////////////////////////////////////////
		/* Set up variables to 0 */		
		nmutation = 0;
		ncross = 0;
		bestfit.fitness = 0;
		bestfit.valorFuncion = 0;
		bestfit.frecuencia = 0;
 		bestfit.amplitud = 0;
		bestfit.generation = 0;
		gen =0;		
		//////////////////////////////////////////////////////
		/* Initialize the populations and report statistics */			
		initialize_pop();
		statistics(oldpop);
		initial_report();
		//show_pop();		
		//////////////////////////////////////////////////////
		/* Iterate the given number of generations */
		for (gen=0; gen<Gmax; gen++)
		{
			show_pop();
			/* Create a new generation */
			generation();
			/* Compute fitness statistics on new populations */
			statistics(newpop);
			/* Report results for new generation */
			report();
			/* Advance the generation */
			temp = oldpop;
			oldpop = newpop;
			newpop = temp;
			//free(temp);			
		}
		////////////////////////////////////////////////////////
		show_pop();
		fclose(ficheroGrafica);
		fclose(ficheroReporte);
		printf("\n\n");
		free_all();		
	}	
}

//////////////////////////Allocate Memory for structures////////////////////////////////////
/* Allocate memory for populations */
void memory_allocation()
{
	unsigned numbytes,numbytes2;	
	int i;
	/* Allocate memory for old and new populations of individuals */
	numbytes = popsize*sizeof(struct individual); //Neccesary bytes
	numbytes2 = 500*sizeof(struct poblacion_); //Neccesary bytes 2	
	if ((oldpop = (struct individual *) malloc(numbytes)) == NULL)
	{
		nomemory("old population");
	}
	if ((newpop = (struct individual *) malloc(numbytes)) == NULL)
	{
		nomemory("new population");
	}
	if ((globalpop = (struct poblacion_ *) malloc(numbytes2)) == NULL)
	{
		nomemory("new population");
	}

	/* Allocate memory for chromosome strings in populations */
	/* Inicializamos el tam. del cromosoma para las 4 variables*/
	for (i=0; i< popsize; i++)
	{
		if ((oldpop[i].chrom = (unsigned *) malloc(codesize*sizeof(unsigned))) == NULL)
		{
			nomemory("old population chromosomes");
		}
		if ((newpop[i].chrom = (unsigned *) malloc(codesize*sizeof(unsigned))) == NULL)
		{
			nomemory("new population chromosomes");
		}
		if ((bestfit.chrom = (unsigned *) malloc(codesize*sizeof(unsigned))) == NULL)
		{
			nomemory("bestfit chromosomes");
		}
	}
	memory_for_selection();
}
/* Allocate memory for performing tournament selection */
void memory_for_selection()
{
	unsigned numbytes;	
	numbytes = popsize*sizeof(int);
	if ((tournlist = (int *) malloc(numbytes)) == NULL)
		nomemory("tournament list");
	/* Use binary tournament selection */
	tournsize = 2;  
}

/* Free memory for population structures*/
void free_all()
{
	int i;
	printf("Liberando globalpop\n");
	free(globalpop);
	printf("Liberando cromosomas\n");
	for (i=0; i < popsize; i++)
	{
		printf("Liberando oldpop %d\n",i );
		free(oldpop[i].chrom);
		printf("Liberando newpop %d\n",i );
		free(newpop[i].chrom);
	}
	printf("Liberando oldpop\n");
	free(oldpop);
	printf("Liberando newpop\n");
	free(newpop);	
	printf("Liberando Torneo\n");
	free_selection();
}

/* Free memory for tournment */
void free_selection()
{
	free(tournlist);
}

void nomemory(char *string)
{
	printf("ERROR!! --> malloc: out of memory making %s\n",string);
	exit(1);
}
///////////////////////////////////////////////////////////////////////////////////////////////

//////////////////////////Randomizers//////////////////////////////////////////////////////////
/* Initialize random numbers batch */
void randomize()
{
    int j1;
    for(j1=0; j1<=54; j1++)
      oldrand[j1] = 0.0;
    jrand=0;
    warmup_random(Rseed);
}

void warmup_random(float random_seed)
{
    int j1, ii;
    double new_random, prev_random;
    oldrand[54] = random_seed;
    new_random = 0.000000001;
    prev_random = random_seed;
    for(j1 = 1 ; j1 <= 54; j1++)
    {
		ii = (21*j1)%54;
		oldrand[ii] = new_random;
		new_random = prev_random-new_random;
		if(new_random<0.0) new_random = new_random + 1.0;
		prev_random = oldrand[ii];
    }
    advance_random();
    advance_random();
    advance_random();
    jrand = 0;
}

void advance_random()
{
    int j1;
    double new_random;
    for(j1 = 0; j1 < 24; j1++)
    {
		new_random = oldrand[j1] - oldrand[j1+31];
		if(new_random < 0.0) new_random = new_random + 1.0;
			oldrand[j1] = new_random;		
    }
    for(j1 = 24; j1 < 55; j1++)
    {
		new_random = oldrand [j1] - oldrand [j1-24];
		if(new_random < 0.0) new_random = new_random + 1.0;
			oldrand[j1] = new_random;
    }
}

//////////////////////////////////////////////////////////////////////////

/*Valores de la funcion ya evaluados*/
void poblacion_evaluada()
{
	FILE *file = fopen("poblacion_.dat", "r");
	struct poblacion_ auxiliar;
	printf("Leyendo archivo\n");   
	while(fread(&auxiliar, sizeof(auxiliar), 1, file))
	{		
    	globalpop[Globalpopsize].x1=auxiliar.x1;
		globalpop[Globalpopsize].x2=auxiliar.x2;
		globalpop[Globalpopsize].x3=auxiliar.x3;
		globalpop[Globalpopsize].valorFuncion = auxiliar.valorFuncion;
		globalpop[Globalpopsize].frecuencia = auxiliar.frecuencia;
		globalpop[Globalpopsize].amplitud = auxiliar.amplitud;
		Globalpopsize++;
		printf("Global Pop size is now:%d\n",Globalpopsize );		
	}	
    fclose(file);      
} 

/**Inicializacion de la Poblacion */
void initialize_pop()
{
	int j,k,i=0;
	double mayor=0;
	for (j=0; j < popsize; j++)//Itera sobre tam de Poblacion
	{
		for (k=0; k < codesize; k++)
		{
			//printf("Valor j:%d \n",j);
			/*if(j==0) // metemos la solucion 0.2 0.0 0.0
			{				
				//printf("Valor k:%d \n",k);
				if(k==0)
				{
					printf("Primer individuo\n");
					//0.2 ---111110100---	
					oldpop[j].chrom[0]=0;
					oldpop[j].chrom[1]=0;
					oldpop[j].chrom[2]=1;
					oldpop[j].chrom[3]=0;
					oldpop[j].chrom[4]=1;
					oldpop[j].chrom[5]=1;
					oldpop[j].chrom[6]=1;
					oldpop[j].chrom[7]=1;
					oldpop[j].chrom[8]=1;
					oldpop[j].chrom[9]=0;
					k=9;				
				}
				else
				{	
					if( k==19 || k==29)
					{
						oldpop[j].chrom[k] = 0;
					}
					else
					{
						oldpop[j].chrom[k] = 1;
					}
				}
			}
			else
			{
				oldpop[j].chrom[k]=flip(0.5);
			}*/	
			oldpop[j].chrom[k]=flip(0.5);
		}
		//Decode binary chain to real number		
		oldpop[j].x1=decode(oldpop[j].chrom,0,codesize/Num_var);
		oldpop[j].x2=decode(oldpop[j].chrom,1,codesize/Num_var);
		oldpop[j].x3=decode(oldpop[j].chrom,2,codesize/Num_var);
		printf(">>>>> Evaluando Individuo:%d <<<<\n",j);		
		objfunc(&(oldpop[j]));
		oldpop[j].parent[0] = 0;        /* Initialize parent info */
		oldpop[j].parent[1] = 0;
		oldpop[j].xsite1 = 0;           /* Initialize crossover sites */
		oldpop[j].xsite2 = 0;
		oldpop[j].placemut = 0;         /* Initialize mutation place */
		oldpop[j].mutation = 0;         /* Initialize mutation indicator */
	}
}

/* Decode the chromosome string into its corresponding decimal value */
/*Decodificacion le mando los extremos de genotipo*/
double decode(unsigned *chrom,int i,int length)
{
	int j,m;
	double accum=0.0;
	for (j=(i*codesize/ Num_var),m=0; j<length*(i+1); j++)
	{
		if (chrom[j]==1)
		{
			accum+=pow(2.0,(double)m);
		}
		m++;
	}
	return accum;
}


void show_pop()
{
	FILE *generaciones;
	generaciones = fopen("R_Generaciones.txt", "a" );
	int i,n=0;
	for (i=0;i<popsize;i++)
	{
		printf("[x1:%f ,",oldpop[i].x1);
		printf("x2:%f ,",oldpop[i].x2);
		printf("x3:%f ]\n",oldpop[i].x3);
		printf("Fitness:%f \n",oldpop[i].fitness);
		printf("Amplitud:%f \n",oldpop[i].amplitud);
		printf("Frecuencia:%f \n",oldpop[i].frecuencia);
		fprintf( generaciones, "I[%i]-----> %f\n",i,oldpop[i].fitness);
		for(n=0;n<codesize;n++)
		{
			printf("%d ",oldpop[i].chrom[n]);
			fprintf(generaciones,"%d ",oldpop[i].chrom[n]);
		}
		printf("\n");
		fprintf(generaciones,"\n");
	}
	fclose(generaciones);
}

int flip(float prob)
{    
    if(randomperc() <= prob)
		return(1);
    else
		return(0);
}


/* Fetch a single random number between 0.0 and 1.0 - Subtractive Method */
/* See Knuth, D. (1969), v. 2 for details */
/* name changed from random() to avoid library conflicts on some machines*/
float randomperc()
{
    jrand++;
    if(jrand >= 55)
    {
		jrand = 1;
		advance_random();
    }
    return((float) oldrand[jrand]);
}

/////////////////////////////Evaluacion Poblacion/////////////////////////////////
void objfunc(struct individual *ind)
{
	double x1,x2,x3;
	float linf =-0.5,linf_x1=-0.1;
	float lsup=0.5,lsup_x1=0.5;
	int l=10;
	ind->fitness=0.0;
	char c[10];
	int indice =-1;
	/*Codifico mis variables*/	
	x1 = round_number(linf_x1+((lsup_x1-linf_x1)/(pow(2,l)-1))*ind->x1);
	ind->x1=x1;
	x2 = round_number(linf+((lsup-linf)/(pow(2,l)-1))*ind->x2);
	ind->x2=x2;
	x3 = round_number(linf+((lsup-linf)/(pow(2,l)-1))*ind->x3);
	ind->x3=x3;
	//Buscamos que ya no esté evaluada la solición
	if(Globalpopsize>0)
	{
		indice = buscar(ind);
	}
	/*Script que usa FDTD para el decaimiento de la energia*/
	if(indice==-1)
	{
		//meep shift1=0.1 shift2=0.4 shift3=0.5 funcion/l3defect.ctl | tee a.out
		//mpirun.openmpi -np 2 meep-mpi-default shift1=0.1 shift2=0.4 shift3=0.5 funcion/l3defect.ctl | tee a.out
		//strcpy(comando,"mpirun.openmpi -np 3 meep-mpi-default shift1=");
		strcpy(comando,"meep shift1=");
		sprintf(c, "%g", x1);strcat(comando,c);
		strcat(comando," shift2=");
		sprintf(c, "%g", x2);strcat(comando,c);
		strcat(comando," shift3=");
		sprintf(c, "%g", x3);strcat(comando,c);
		strcat(comando," funcion/l3defect.ctl > funcion/a.out");
		printf("\n\n%s\n",comando);
		system(comando);	

		//PARA EL FACTOR FRECUENCIA AL FINAL ES -f 1
		printf("Obteniendo Frecuencia\n");	
		strcpy(comando,"grep harminv funcion/a.out | cut -d "); 
		strcat(comando,",");
		strcat(comando," -f 2,4,5 | sort -r -g -k 3 -t , | sed -n 1,1p | cut -d ");
		strcat(comando,",");
		strcat(comando," -f 1");
		fp_aux = popen(comando, "r");
			fscanf(fp_aux, "%lf", &ind->frecuencia);
			printf("Frecuencia:%lf \n",ind->frecuencia );
		pclose(fp_aux);
	
		//PARA EL FACTOR Q AL FINAL ES -f 2  lo guardaremos en el ind->valorfunción 
		printf("Obteniendo Q\n");	
		strcpy(comando,"grep harminv funcion/a.out | cut -d ");
		strcat(comando,",");
		strcat(comando," -f 2,4,5 | sort -r -g -k 3 -t , | sed -n 1,1p | cut -d ");
		strcat(comando,",");
		strcat(comando," -f 2");
		fp_aux = popen(comando, "r");
			fscanf(fp_aux, "%lf", &ind->fitness);
			printf("Fitness:%lf \n",ind->fitness );
		pclose(fp_aux);
     
		//PARA LA AMPLITUD AL FINAL ES -f 3
		//grep harminv a.out | cut -d "," -f 2,4,5 | sort -r -g -k 3 -t , | sed -n 1,1p | cut -d "," -f 2
		printf("Obteniendo Amplitud\n");	
		strcpy(comando,"grep harminv funcion/a.out | cut -d "); 
		strcat(comando,",");
		strcat(comando," -f 2,4,5 | sort -r -g -k 3 -t , | sed -n 1,1p | cut -d ");
		strcat(comando,",");
		strcat(comando," -f 3");
		fp_aux = popen(comando, "r"); 
			fscanf(fp_aux, "%lf", &ind->amplitud);
			printf("Amplitud:%lf \n\n",ind->amplitud );
		pclose(fp_aux);
		//En este caso el valor de la función el fitness es igual y es dado por el factor Q
		ind->valorFuncion = ind->fitness;
		add_individual(ind);
	}
	else
	{
		printf("Individuo ya evaluado\n");
		ind->valorFuncion = globalpop[indice].valorFuncion;
		ind->fitness = globalpop[indice].valorFuncion;
		ind->frecuencia = globalpop[indice].frecuencia;
		ind->amplitud = globalpop[indice].amplitud;
		printf("Frecuencia:%lf \n",ind->frecuencia );
		printf("Fitness:%lf \n",ind->fitness );
		printf("Amplitud:%lf \n\n",ind->amplitud );
		
	}		
}

double round_number(double number)
{
	int i=0;
	for (i=0;i<21;i++)
		if (number<=varibles_interval[i])
			break;
	return varibles_interval[i];
}

void add_individual(struct individual * ind)
{
	if(Globalpopsize<(500))
	{
		//Agregar Individuo a Arreglo en Memoria
		globalpop[Globalpopsize].x1=ind->x1;
		globalpop[Globalpopsize].x2=ind->x2;
		globalpop[Globalpopsize].x3=ind->x3;
		globalpop[Globalpopsize].valorFuncion = ind->valorFuncion;
		globalpop[Globalpopsize].valorFuncion = ind->fitness;
		globalpop[Globalpopsize].frecuencia = ind->frecuencia;
		globalpop[Globalpopsize].amplitud = ind->amplitud;		
		//Agregar Individuo a Archivo
		printf("Agregare:\n");
		printf("x1:%f,x2:%f,x3:%f\n",ind->x1,ind->x2,ind->x3 );
		FILE *file = fopen("poblacion_.dat", "a");
		struct poblacion_ individual;
		individual.x1=ind->x1;
		individual.x2=ind->x2;
		individual.x3=ind->x3;
		individual.frecuencia=ind->frecuencia;
		individual.valorFuncion=ind->fitness;
		individual.amplitud=ind->amplitud;    
			fwrite(&individual, sizeof(individual), 1, file);			
    	fclose(file);  
    	Globalpopsize++;
	}
}

int buscar (struct individual * ind)
{
	for (int i=0; i < Globalpopsize; i++) 
	{
		if ((globalpop[i].x1==ind->x1)&&(globalpop[i].x2==ind->x2)&&(globalpop[i].x3==ind->x3))
			{
				return i;
			}
	}
	return -1;
}

void statistics(struct individual *pop)
{
	int i,j;
	int mejor=0;
	int peor=0; 
	sumfitness = 0.0;	
	min = pop[0].fitness;
	max = pop[0].fitness;
	localmax = pop[0].chrom;
	/* Loop for max, min, sumfitness */
	for (j=0; j < popsize; j++)
	{
		sumfitness = sumfitness + pop[j].fitness;            /* Accumulate */
		if (pop[j].fitness > max)
		{ 
			/* New maximum  fitness*/
			max = pop[j].fitness; 
			/* Store string */
			localmax = pop[j].chrom; 
			/* New maximum  fitness index*/
		    mejor=j;
		}
		if (pop[j].fitness < min)
		{
			/* New minimim fitness */
			min = pop[j].fitness;
			/* New minimum index*/
			peor=j;
		}
		/* Define new global best-fit individual */
		if (pop[j].fitness  > bestfit.fitness)
		{
			for (i=0; i < codesize; i++)
			{
		  		bestfit.chrom[i] = pop[j].chrom[i];
		  	}			
			bestfit.fitness = pop[j].fitness;
			bestfit.valorFuncion = pop[j].valorFuncion;
			bestfit.frecuencia = pop[j].frecuencia;
	 		bestfit.amplitud = pop[j].amplitud;
			bestfit.x1=pop[j].x1;
			bestfit.x2=pop[j].x2;  
			bestfit.x3=pop[j].x3;   
	        bestfit.xsite1=pop[j].xsite1;             
	        bestfit.xsite2=pop[j].xsite2;             
	        bestfit.placemut=pop[j].placemut;		 
	        bestfit.mutation=pop[j].mutation;
		  	for (i=0; i < 2; i++)
		  	{
		  		bestfit.parent[i] =pop[j].parent[i];
		  	}		  	        
            bestfit.ve=pop[j].ve;
		  	bestfit.generation = gen;
		}
	}
	/*Elitismo saco el peor e ingreso al mejor*/
	if( pop[mejor].fitness!=bestfit.fitness)
	{		
		for (i=0; i < codesize; i++)
		{
			 pop[peor].chrom[i]=bestfit.chrom[i];	
		}		
	 	pop[peor].fitness=bestfit.fitness ; 
	 	pop[peor].valorFuncion=bestfit.valorFuncion ;
		pop[peor].frecuencia=bestfit.frecuencia ;
        pop[peor].amplitud=bestfit.amplitud;
	 	pop[peor].x1=bestfit.x1;
	 	pop[peor].x2=bestfit.x2;
		pop[peor].x3=bestfit.x3;    
        pop[peor].xsite1=bestfit.xsite1;           
        pop[peor].xsite2=bestfit.xsite2;            
        pop[peor].placemut=bestfit.placemut;		 
        pop[peor].mutation=bestfit.mutation;
		for (i=0; i < 2; i++)
		{
	  		pop[peor].parent[i]=bestfit.parent[i]; 
	  	}       
        pop[peor].ve=bestfit.ve;
	}
	/* Calculate average */
	avg = sumfitness/popsize;
	fprintf( ficheroGrafica, "%d\t%.10f\n",gen,bestfit.valorFuncion);
}

void initial_report()
{
	printf("\n Parameters used with the Genetic Algorithm \n");
 	fprintf( ficheroReporte,"\n Parameters used with the Genetic Algorithm \n");
	printf(" Total population size          =       %d\n",popsize);
	fprintf( ficheroReporte," Total population size          =       %d\n",popsize);
	printf(" Chromosome length              =       %d\n",codesize);
	fprintf( ficheroReporte," Chromosome length              =       %d\n",codesize);
	printf(" Maximum number of generations  =       %d\n",Gmax);
	fprintf( ficheroReporte," Maximum number of generations  =       %d\n",Gmax);
	printf(" Crossover probability          =       %f\n",Pc);
	fprintf( ficheroReporte," Crossover probability          =       %f\n",Pc);
	printf(" Mutation probability           =       %f\n",Pm);
	fprintf( ficheroReporte," Mutation probability           =       %f\n",Pm);
	fprintf( ficheroReporte," Seed           =       %f\n",Rseed);
	printf("\n\n");
	fprintf( ficheroReporte,"\n\n");
}


void generation()
{
  int mate1, mate2, jcross1, jcross2, j = 0;
  int ban=0;
  /* perform any preselection actions necessary before generation */
  int i=0;
  preselect(); 
  /* select, crossover, and mutation */
  do
    {
    	printf("Selection:\n");
      	/* pick a pair of mates */
     	 mate1 = selection();
     	 mate2 = selection();
     	 printf("->Mate1:%d \n",mate1);
     	 printf("->Mate2:%d \n",mate2);
     	 /* En esta parte aplique el elitismo lo cruzo con el mismo asi genero el mismo indiviuo*/
      	if(oldpop[mate1].fitness==bestfit.fitness && ban==0)
	  	{
			mate2=mate1;
			ban++;
		}
		else
		{
			if(oldpop[mate2].fitness==bestfit.fitness && ban==0)
			{
				mate1=mate2;
				ban++;
			}
		}
		//printf("->Mate1_new:%d \n",mate1);
     	//printf("->Mate2_new:%d \n",mate2);
	    /* Crossover and mutation */
	    crossover(oldpop[mate1].chrom, oldpop[mate2].chrom,newpop[j].chrom,newpop[j+1].chrom,&jcross1,&jcross2);
	    /* Dejo una copia intacta del mejor indiviuo y la otro la muto*/
	    if(ban!=1)
		{		 
			mutation(&(newpop[j]));
		   	mutation(&(newpop[j+1]));
	    }
		else
		{
			ban++;
		}
	    /* Decode string, evaluate fitness, & record */
	    /* parentage date on both children */
	    newpop[j].x1=decode(newpop[j].chrom, 0,codesize/Num_var);
	    newpop[j+1].x1=decode(newpop[j+1].chrom, 0,codesize/Num_var);
	    newpop[j].x2=decode(newpop[j].chrom, 1,codesize/Num_var);
	    newpop[j+1].x2=decode(newpop[j+1].chrom, 1,codesize/Num_var);
	    newpop[j].x3=decode(newpop[j].chrom, 2,codesize/Num_var);
	    newpop[j+1].x3=decode(newpop[j+1].chrom, 2,codesize/Num_var);
	    objfunc(&(newpop[j]));
	    newpop[j].parent[0] = mate1+1;
	    newpop[j].xsite1 = jcross1;
	    newpop[j].xsite2 = jcross2;
	    newpop[j].parent[1] = mate2+1;
	    objfunc(&(newpop[j+1]));
	    newpop[j+1].parent[0] = mate1+1;
	    newpop[j+1].xsite1 = jcross1;
	    newpop[j+1].xsite2 = jcross2;
	    newpop[j+1].parent[1] = mate2+1;
	    /* Increment population index */
	    printf("Added childs %d and %d\n",j,j+1);
	    j = j + 2;	    
    }while(j < (popsize-1));
}

void preselect()
{
	reset();
	tournpos = 0;
}

/* Proceso de seleccion mediante Ruleta*/
int selection()
{
	int pick, winner, i;
	double acumulado,ra;
	/* If remaining members not enough for a tournament, then reset list */
	if ((popsize - tournpos) < tournsize)
	{
		reset();
		tournpos = 0;
	}
	// Selecciona por ruleta 
	ra=rndreal(0,popsize-1);
	printf("Random selection number:%f\n",ra);
	acumulado=0;
	for (i=0; i < popsize; i++)
	{
		acumulado+=oldpop[tournlist[i]].ve;
		if(acumulado>=ra)
		{
			winner= tournlist[i];
			break;
		}
	}
	/* Update tournament position */
	tournpos += tournsize;
	return(winner);
}

/* Shuffle the tournament list at random */
/* aqui se barejan los individuos*/
void reset()
{
	int i, rand1, rand2, temp;
	double suma=0;
	double media=0;
	/* Aqui le pones el valor esperado */
	for (i=0; i < popsize; i++) 
	{
		suma+=oldpop[i].fitness;
	}
	media= suma/popsize;
	printf("Media:%f\n",media);	
	for (i=0; i < popsize; i++) 
	{
		oldpop[i].ve=oldpop[i].fitness/media;
		printf("Ve[%i]:%f\n",i,oldpop[i].ve);	
	}
	for (i=0; i < popsize; i++) 
	{
		tournlist[i]=i;
	}
	for (i=0; i < popsize; i++)
	{
		rand1=rnd(i,popsize-1);
		rand2=rnd(i,popsize-1);
		temp=tournlist[rand1];
		tournlist[rand1]=tournlist[rand2];
		tournlist[rand2]=temp;
	}
	printf("Jugadores del Torneo Barajeados\n");
	for (i=0; i < popsize; i++)
	{
		printf("T[%d]=%d\n",i,tournlist[i]);		
	}
}

int rnd(int low, int high)
{
    int i;
    float randomperc();
    if(low >= high)
    {
		i = low;
	}
    else
    {
		i = (randomperc() * (high - low + 1)) + low;
		if(i > high)
		{
		 	i = high;
		}
    }
    return(i);
}

/* real random number between specified limits */
float rndreal(float lo ,float hi)
{
    return((randomperc() * (hi - lo)) + lo);
}



/* La cruza realizada es la cruza simple de un solo punto */

void crossover (unsigned *parent1, unsigned *parent2,unsigned *child1, unsigned *child2,int *jcross1, int *jcross2)
{
    int j;
    unsigned temp;
    /* Do crossover with probability Pc */
    if(flip(Pc))
    {
    	printf("Cruza Exitosa\n");
    	 /* escogo un punto aleatorio para la cruza */
		*jcross1 = rnd(puntocruza_mutacion,(codesize - 1));         
		printf("Cruza en punto:%d\n",jcross1[0]);
		ncross++;
		/*La cruza de un punto*/
		for (j=(codesize-1); j>(codesize - *jcross1); j--)
		{
			child2[j]=parent1[j];
			child1[j]=parent2[j];			
		}
		for (j=j;j>=0; j--)
		{
			child1[j]=parent1[j];
			child2[j]=parent2[j];
		}	
    }
	else 
	{
		printf("Sin Cruza\n");
		for (j=0; j<=(codesize-1); j++)
		{
		  child1[j]=parent1[j];
		  child2[j]=parent2[j];
		}
		*jcross1=0; *jcross2=0;
	}
	*jcross2=*jcross1;
	int n;
	printf("Child1\n");
	for (n=0;n<codesize; n++)
		{
			printf("%u ",child1[n]);
		}
	printf("\nChild2\n");
	for (n=0;n<codesize; n++)
		{
			printf("%u ",child2[n]);
		}
}

void mutation(struct individual *ind)
{
	int placemut,k;
	int result;
	/* mutacion uniforme con el flip a cada alelo y emepzando de las variables x2 en adelante*/
	for(k=puntocruza_mutacion;k<codesize;k++)
	{
		if(flip(Pm))
		{
			if ((ind->chrom[k])==0)
				ind->chrom[k]=1;
			else 
				ind->chrom[k]=0;	
			nmutation++;
		}
	}
}

void report()
{	
	printf("\nGeneracion # %d Estadisticas Acumuladas: \n",gen);
	fprintf( ficheroReporte,"\nGeneracion # %d Estadisticas Acumuladas: \n",gen);
	printf("Total Cruzas = %d, Total Mutacion = %d\n", ncross,nmutation);
	fprintf( ficheroReporte,"Total Cruzas = %d, Total Mutacion = %d\n", ncross,nmutation);
	printf("Aptitud minima = %.10f   Maxima = %.10f   Media = %.10f   Suma de Aptitudes = %.10f\n",min,max,avg,sumfitness);
	fprintf( ficheroReporte,"Aptitud minima = %.10f   Maxima = %.10f   Media = %.10f   Suma de Aptitudes = %.10f\n",min,max,avg,sumfitness);
	printf("Mejor individuo de esta generacion= ");
	fprintf( ficheroReporte,"Mejor individuo de esta generacion= ");
	writechrom(localmax); printf("\n");fprintf( ficheroReporte,"\n");
	printf("El mejor individuo Global hasta hora, Generacion # %d:\n",bestfit.generation);
	fprintf( ficheroReporte,"El mejor individuo Global hasta hora, Generacion # %d:\n",bestfit.generation);
	printf("Cadena Binaria = ");	
	fprintf( ficheroReporte,"Cadena Binaria = ");
	writechrom((&bestfit)->chrom);
	printf("    Aptitud = %.10f\n", bestfit.fitness);
	fprintf( ficheroReporte,"    Aptitud = %.10f\n", bestfit.fitness);
	printf("     x1: %.2f \t x2: %.2f \t x3: %.2f\n", bestfit.x1, bestfit.x2, bestfit.x3);
	fprintf( ficheroReporte,"     x1: %.2f \t x2: %.2f \t x3: %.2f\n",  bestfit.x1, bestfit.x2, bestfit.x3);
	printf("    Valor de la funcion = %.10f\n", bestfit.valorFuncion);
	fprintf( ficheroReporte,"    Valor de la funcion = %.10f\n", bestfit.valorFuncion);
	printf("    Valor de la frecuencia = %.10f\n", bestfit.frecuencia);
	fprintf( ficheroReporte,"    Valor de la frecuencia = %.10f\n", bestfit.frecuencia);
	printf("    Valor de la amplitud = %.10f\n", bestfit.amplitud);
	fprintf( ficheroReporte,"    Valor de la amplitud = %.10f\n", bestfit.amplitud);
 }

void writechrom(unsigned *chrom)
{
	int j;
	for (j=(codesize-1); j>=0; j--)
	{
		if (chrom[j]==0)
		{
			 printf("0");
			 fprintf( ficheroReporte,"0");
		}
		else
	    {
	    	printf("1");
	    	fprintf( ficheroReporte,"1");
	    }
	}
}
