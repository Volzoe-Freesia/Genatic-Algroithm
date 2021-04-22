# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <time.h>

/* Change any of these parameters to match your needs */
#define POPSIZE 900 /* ÖÖÈº´óÐ¡ */
#define MAXGENS 40
/* ×î¶àÉú³ÉµÄ´úÊý */
#define NVARS 2 /* ±äÁ¿µÄ¸öÊý */
#define PXOVER 0.8 /* ½»²æµÄ¸ÅÂÊ */
#define PMUTATION 0.16 /* Í»±äµÄ¸ÅÂÊ */
#define TRUE 1
#define FALSE 0
#define branch_no 8



int generation; /* µ±Ç°µÄ´úÊý */
int cur_best; /* ×î¼ÑµÄ¸öÌå */
int all_branchfit[branch_no];   /*Õû¸öÖÖÈºµÄ·ÖÖ§¸²¸ÇÇé¿ö*/
int randint = 0;

FILE* galog; /* an output file */

/* genotype (GT), a member of the population */
struct genotype
{
	double gene[NVARS]; /* a string of variables */
	int branchfit[branch_no]; /* ·ÖÖ§¸²¸ÇÇé¿ö */
	double fitness; /* GT's fitness */

	double rfitness; /* relative fitness */
	double cfitness; /* cumulative fitness */
};
struct genotype population[POPSIZE + 1]; /* population */
struct genotype newpopulation[POPSIZE + 1]; /* new population; */

/************±»²â³ÌÐò***********/
void testedPro(double a, double b, int* fi)
{
    if(a<100)
    {
        fi[0]=1;
        goto label;
    }
    label:
     if(a<90)
    {
        fi[1]=1;
        goto label1;
    }
    label1:
     if(a<80)
    {
        fi[2]=1;
        goto label3;
    }
    label3:
     if(a<70)
    {
        fi[3]=1;
        goto label4;
    }
    label4:
     if(b<60)
    {
        fi[4]=1;
        goto label5;
    }
    label5:
     if(b<50)
    {
        fi[5]=1;
        goto label6;
    }
    label6:
     if(b<40)
    {
        fi[6]=1;
        goto label7;
    }
    label7:
    fi[7]=1;
}


/***********************************************************/
/* Random value generator: Generates a value within bounds */
/***********************************************************/
double randval(void)
{
	double val;
	randint++;

	val = (double)(rand() % 100 - 10);
	return (val);
}


//ÖÖÈº³õÊ¼»¯
void initialize(void)
{

	int i, j, z;
	/* initialize variables within the bounds */
	srand((unsigned)time(NULL));
	for (i = 0; i < NVARS; i++)
	{
		for (j = 0; j < POPSIZE; j++)
		{
			population[j].fitness = 0;
			population[j].rfitness = 0;
			population[j].cfitness = 0;
			population[j].gene[i] = randval();
			for (z = 0; z < branch_no; z++)
				population[j].branchfit[z] = 0;
		}
	}
}

/*************************************************************/
/* Evaluation function: This takes a user defined function. */
/* Each time this is changed, the code has to be recompiled. */
/* The current function is:  */
/*************************************************************/
void evaluate(void)
{
	int mem, bno;
	for (mem = 0; mem < POPSIZE; mem++)
	{
		testedPro(population[mem].gene[0], population[mem].gene[1], population[mem].branchfit);
	}

	for (bno = 0; bno < branch_no; bno++)
	{
		for (mem = 0; mem < POPSIZE; mem++)
		{
			all_branchfit[bno] += population[mem].branchfit[bno];
		}
	}
	for (mem = 0; mem < POPSIZE; mem++)
	{
		for (bno = 0; bno < branch_no; bno++)
		{
			if (all_branchfit[bno] != 0)
				population[mem].fitness += population[mem].branchfit[bno] / (all_branchfit[bno]);
		}
	}
}


/**************************************************************/
/* Selection function: Standard proportional selection for */
/* maximization problems incorporating elitist model - makes */
/* sure that the best member survives */
/**************************************************************/
void select(void)
{
	int mem, i, j;
	double sum = 0;
	double p;
	/* find total fitness of the population */
	for (mem = 0; mem < POPSIZE; mem++)
	{
		sum += population[mem].fitness;
	}
	/* calculate relative fitness */
	for (mem = 0; mem < POPSIZE; mem++)
	{
		population[mem].rfitness = population[mem].fitness / sum;
	}
	population[0].cfitness = population[0].rfitness;
	/* calculate cumulative fitness */
	for (mem = 1; mem < POPSIZE; mem++)
	{
		population[mem].cfitness = population[mem - 1].cfitness + population[mem].rfitness;
	}
	/* finally select survivors using cumulative fitness. */
	for (i = 0; i < POPSIZE; i++)
	{
		p = rand() % 1000 / 1000.0;
		if (p < population[0].cfitness)
			newpopulation[i] = population[0];
		else
		{
			for (j = 0; j < POPSIZE; j++)
				if (p >= population[j].cfitness && p < population[j + 1].cfitness)
					newpopulation[i] = population[j + 1];
		}
	}
	/* once a new population is created, copy it back */
	for (i = 0; i < POPSIZE; i++)
		population[i] = newpopulation[i];
}

/*************************************************************/
/* Swap: A swap procedure that helps in swapping 2 variables */
/*************************************************************/
void swap(double* x, double* y)
{
	double temp;
	temp = *x;
	*x = *y;
	*y = temp;
}

/**************************************************************/
/* Crossover: performs crossover of the two selected parents. */
/**************************************************************/
void Xover(int one, int two)
{
	int i;
	int point; /* crossover point */
	/* select crossover point */
	if (NVARS > 1)
	{
		if (NVARS == 2)
			point = 1;
		else
			point = (rand() % (NVARS - 1)) + 1;
		for (i = 0; i < point; i++)
		{
			swap(&population[one].gene[i], &population[two].gene[i]);
		}

	}
}


/***************************************************************/
/* Crossover selection: selects two parents that take part in */
/* the crossover. Implements a single point crossover */
/***************************************************************/
void crossover(void)
{
	int mem, one;
	int first = 0; /* count of the number of members chosen */
	double x;
	for (mem = 0; mem < POPSIZE; ++mem)
	{
		x = rand() % 1000 / 1000.0;
		if (x < PXOVER)
		{
			++first;
			if (first % 2 == 0)
				Xover(one, mem);
			else
				one = mem;
		}
	}
}

/**************************************************************/
/* Mutation: Random uniform mutation. A variable selected for */
/* mutation is replaced by a random value between lower and */
/* upper bounds of this variable */
/**************************************************************/
void mutate(void)
{
	int i, j;
	double x;
	srand((unsigned)time(NULL));
	for (i = 0; i < POPSIZE; i++)
		for (j = 0; j < NVARS; j++)
		{
			x = rand() % 1000 / 1000.0;
			if (x < PMUTATION)
			{
				/* find the bounds on the variable to be mutated */
				population[i].gene[j] = randval();
			}
		}
}


void setallbf()
{
	int bno;
	for (bno = 0; bno < branch_no; bno++)
		all_branchfit[bno] = 0;
}

void setbf()
{
	int bno, mem;
	for (mem = 0; mem < POPSIZE; mem++)
	{
		for (bno = 0; bno < branch_no; bno++)
			population[mem].branchfit[bno] = 0;
	}
}

/**************************************************************/
/* Main function: Each generation involves selecting the best */
/* members, performing crossover & mutation and then */
/* evaluating the resulting population, until the terminating */
/* condition is satisfied */
/**************************************************************/
int main(void)
{
    int begin=clock();
	int end;
	int i, mem, bno;
	if ((galog = fopen("galog.txt", "w")) == NULL)
	{
		exit(1);
	}
	generation = 0;
	initialize();
	evaluate();

	while (generation < MAXGENS)
	{
		generation++;
		fprintf(galog, "%d:\n", generation);
		for (bno = 0; bno < branch_no; bno++)
			fprintf(galog, " all_branchfit(%d) = %d\n", bno, all_branchfit[bno]);

		setallbf();
		select();
		crossover();
		mutate();
		setbf();
		evaluate();
	}

	fprintf(galog, "Simulation completed\n");
	fprintf(galog, "Result: ");
	for (mem = 0; mem < POPSIZE; mem++)
	{
		fprintf(galog, "%d:\n", mem);
		for (i = 0; i < NVARS; i++)
		{
			fprintf(galog, "%var(%d) = %3.3f  ", i, population[mem].gene[i]);

		}
		/*for(bno = 0;bno < branch_no;bno++)
		{
			fprintf (galog,"\n %d : branchfit(%d) = %d",mem,bno,population[mem].branchfit[branch_no]);
		}*/
		fprintf(galog, "\n");
	}

	for (bno = 0; bno < branch_no; bno++)
		fprintf(galog, "\n all_branchfit(%d) = %d", bno, all_branchfit[bno]);
	fclose(galog);
	printf("Success!\n");
	end=clock();
	for (mem = 0; mem < POPSIZE; mem++)
	{
		fprintf(galog, "%d:\n", mem);
		for (i = 0; i < NVARS; i++)
		{
			fprintf(galog, "%var(%d) = %3.3f  ", i, population[mem].gene[i]);
            printf( "%var(%d) = %3.3f \n ", i, population[mem].gene[i]);

		}
		/*for(bno = 0;bno < branch_no;bno++)
		{
			fprintf (galog,"\n %d : branchfit(%d) = %d",mem,bno,population[mem].branchfit[branch_no]);
		}*/
		fprintf(galog, "\n");
	}

	int time;
	time=end-begin;
	printf("running time"":""%d millisecond \n",time);

	return 0;
}
