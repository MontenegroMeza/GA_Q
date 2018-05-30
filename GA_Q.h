/* Global structures and variables */

struct individual 
{
    unsigned *chrom;             /*      Chromosome string for the individual */
    double x1,x2,x3;		     /*      Value of the decoded string */
    double   fitness;            /*      Fitness of the individual */
    int      xsite1;             /*      Crossover site 1 at mating */
    int      xsite2;             /*      Crossover site 2 at mating */
    int      placemut;		     /*      Mutation place */
    int      mutation;		     /*      Mutation indicator */
    int      parent[2];          /*      Who the parents of offspring were */
    double   ve;        
    double valorFuncion;            
    double frecuencia;	         /*      Para le problema de cavidades de nanoestrcuturas */
    double amplitud;		     /*      Para le problema de cavidades de nanoestrcuturas */
   };

struct bestever
{
    unsigned *chrom;             /*      chromosome string for the individual */
    double x1,x2,x3;		 /*	          value of the decoded string */
    double   fitness;            /* 	            fitness of the individual */
    int      xsite1;             /* 	           crossover site 1 at mating */
    int      xsite2;             /* 	           crossover site 2 at mating */
    int      placemut;		 /*			       mutation place */
    int      mutation;		 /*		           mutation indicator */
    int      parent[2];          /*         who the parents of offspring were */
    double   ve;
    double valorFuncion;
    double frecuencia;		/*        Para le problema de cavidades de nanoestrcuturas */
    double amplitud;		/*        Para le problema de cavidades de nanoestrcuturas */
    int      generation;    /*                   generation which produced it */
};

struct poblacion_
{
   
    double x1,x2,x3;		 /*	          value of the decoded string */
    double frecuencia;
    double valorFuncion;		/*        Para le problema de cavidades de nanoestrcuturas */
    double amplitud;		/*        Para le problema de cavidades de nanoestrcuturas */
};

/* Functions prototypes */
void poblacion_evaluada();
double round_number(double number);
int buscar (struct individual *);
void memory_allocation();
void memory_for_selection();
void free_all();
void free_selection();
void nomemory(char *);
void initialize_pop();
void initial_report();
void statistics(struct individual *);
void warmup_random(float);
float rndreal(float,float);
int rnd(int,int);
float randomperc();
double randomnormaldeviate();
void randomize();
double noise(double,double);
void initrandomnormaldeviate();
int flip(float);
void advance_random();
void generation();
int selection();
void reset();
void preselect();
void mutation(struct individual *);
void crossover(unsigned *,unsigned *,unsigned *,unsigned *,int *, int *);
void report();
void writepop();
void writechrom(unsigned *);
void cls();
double decode(unsigned *,int,int);
void objfunc(struct individual *);
void add_individual(struct individual *);

static int *tournlist, tournpos, tournsize; /* Tournment list, position in list */
struct individual *oldpop;                    /* last generation of individuals */
struct individual *newpop;                    /* next generation of individuals */
struct poblacion_ *globalpop;                    /* next generation of individuals */
struct bestever bestfit;                           /* fittest individual so far */
double sumfitness;                      /* summed fitness for entire population */
double max;                                    /* maximum fitness of population */
unsigned *localmax;		   /* String corresponding to the local maximum */
double avg;                                    /* average fitness of population */
double min;                                    /* minimum fitness of population */
float  Pc;                                          /* probability of crossover */
float  Pm;                                           /* probability of mutation */
int    gen;                                        /* current generation number */
int    Gmax;                                   /* maximum number of generations */
int    nmutation;                                        /* number of mutations */
int    ncross;                                          /* number of crossovers */
float  Rseed;			       		         /* Random numbers seed */
double oldrand[55];                               /* Array of 55 random numbers */
int jrand;                                             /* current random number */
double rndx2;                                /* used with random normal deviate */
int rndcalcflag;                             /* used with random normal deviate */

