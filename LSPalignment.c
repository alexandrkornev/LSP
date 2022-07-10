#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#define MAXLEN 2048 
#define MAXNEI 1024
double THETA_CUTOFF = 30.0;
double DELTA_ALPHA = 0.4;
double DELTA_BETA = 0.75;
struct VERTEX
{
  char res;
  char chain;
  char pdbnum[10];
  int  neinum;
  int  vertex[MAXNEI];
  char neires[MAXNEI];
  char neichain[MAXNEI];
  char neipdbnum[MAXNEI][10];
  double AA[MAXNEI];
  double AB[MAXNEI];
  double BA[MAXNEI];
  double alpha[MAXNEI];
};
struct EDGE
{
  int graph;
  int v1;
  int v2;
  double length;
};
struct DOUBLE_EDGE
{
  int v11;
  int v12;
  int v21;
  int v22;
  double weight;
};
struct PROTEIN
{
	char res;
	char chain;
	char pdbnum[10];
	int degree;
};
int readgraph(FILE*, struct VERTEX*);
int read_edges(struct EDGE*, struct VERTEX*, int);
int compare_edge(struct EDGE*,struct EDGE*);
int StripLeadingSpaces(char*,char*);
void Copy(struct EDGE*, struct EDGE *);
int find_neighbor_num(struct VERTEX*,int,int);
int the_same(struct VERTEX,struct VERTEX);
int part_of_a_triangle(struct DOUBLE_EDGE*,int,int,int);
int main(int argc, char *argv[]){
	char *file1;
	char *file2;
	char output_file[256];
	int i,j,i11,i12,i21,i22;
	int g1_len,g2_len;
	int e1_len,e2_len;
	int de_index;
	int graph1, graph2;
	int across;
	int neighbor_num1,neighbor_num2;
	int seq_length;
	double **adjacency_matrix;
	double delta_CA;
	double delta_CB1;
	double delta_CB2;
	double delta_alpha;
	FILE *FILE1;
	FILE *FILE2;
	FILE *OUTPUT;
	struct VERTEX *g1;
	struct VERTEX *g2;
	struct EDGE *e1;
	struct EDGE *e2;
	struct EDGE *pool;
	struct DOUBLE_EDGE *de;
	struct PROTEIN *protein;
	printf("Graph similarity search by Edge Comparison and Combination. MD version. v.0.96 (November 21 2019)\n"
		"Usage: MDLSP_W.exe -1graph1 -2graph2\n"
		"Optional parameters:\n"
		"\t-a<RealNumber> tolerance for Ca-Ca distance comparison (Angstroms). Default = 0.4A\n"
		"\t-b<RealNumber> tolerance for Ca-Cb distance comparison (Angstroms). Default = 0.75A\n"
		"\t-g<RealNumber> tolerance for dihedral angle comparison (degrees). Default = 30 deg\n"
		"\nWritten by Alexandr Kornev (2003-2019)\n");
	if (argc < 3){
		printf("\n\tWrong parameters.\n");
		return 0;
	}
	while ((argc > 1) && (argv[1][0] == '-')){
		switch (argv[1][1]){
			case '1':
			  file1=&argv[1][2];
			  break;
			case '2':
			  file2=&argv[1][2];
			  break;
			case 'a':
			  DELTA_ALPHA  = atof(&argv[1][2]);
			  break;
			case 'b':
			  DELTA_BETA  = atof(&argv[1][2]);
			  break;
			case 'g':
			  THETA_CUTOFF = atof(&argv[1][2]);
			  break;
			default:
			  printf("Bad option %s\n", argv[1]);
			  return 0;
		}
		++argv;
		--argc;
	}
	if ((FILE1=fopen(file1,"r")) == NULL ){
	  printf("\n Can't open %s file.\n", file1);
	  return 0;
	}
	if ((FILE2=fopen(file2,"r")) == NULL ){

	  printf("\n Can't open %s file.\n", file2);
	  return 0;
	}

	if((g1 = (struct VERTEX *) calloc (MAXLEN, sizeof(struct VERTEX))) == NULL){
		printf("Can't allocate memory for g1\n");
	return 0;
	}

	if((g2 = (struct VERTEX *) calloc (MAXLEN, sizeof(struct VERTEX))) == NULL)
	{
		printf("Can't allocate memory for g2\n");
		return 0;
	}
	g1_len = readgraph(FILE1,g1);
	g2_len = readgraph(FILE2,g2);
	if(g1_len!=g2_len) {
		printf("Warning: Graphs have different sizes: %i vs. %i\n", g1_len, g2_len);
		return 0;
	} else {
		seq_length=g1_len;
	}
	if((protein = (struct PROTEIN *) calloc (seq_length, sizeof(struct PROTEIN))) == NULL)
		{
			printf("Can't allocate memory for protein\n");
			return 0;
		}
	if((adjacency_matrix = (double **) calloc(seq_length,sizeof(double *))) == NULL){
	    printf("Can't allocate memory for adjacency matrix\n");
	    return 0;
	}
	for (i = 0; i < seq_length; i++) {
		if((adjacency_matrix[i] = (double *) calloc(seq_length,sizeof(double))) == NULL){
			printf("Can't allocate memory for substitution matrix\n");
			return 0;
	    }
	}
	if((e1 = (struct EDGE *) calloc (g1_len*MAXNEI, sizeof(struct EDGE))) == NULL)
	{
		printf("Can't allocate memory for e1\n");
		return 0;
	}
	if((e2 = (struct EDGE *) calloc (g2_len*MAXNEI, sizeof(struct EDGE))) == NULL)
	{
		printf("Can't allocate memory for e1\n");
		return 0;
	}
	e1_len = read_edges(e1,g1,1);
	e2_len = read_edges(e2,g2,2);

	if((pool = (struct EDGE *) calloc ((e1_len+e2_len), sizeof(struct EDGE))) == NULL)
	{
		printf("Can't allocate memory for pool\n");
		return 0;
	}
	if((de = (struct DOUBLE_EDGE *) calloc (50*(e1_len+e2_len), sizeof(struct DOUBLE_EDGE))) == NULL)
	{
		printf("Can't allocate memory for de\n");
		return 0;
	}
	for(i=0;i<e1_len;i++)
		Copy(&pool[i],&e1[i]);
	for(i=e1_len,j=0;i<(e1_len+e2_len);i++,j++)
		Copy(&pool[i],&e2[j]);
	qsort(pool,
		  e1_len+e2_len,
		  sizeof(struct EDGE),
		  (void *) compare_edge);
	de_index=0;
	for(i=0;i<(e1_len+e2_len);i++)
	{ 
	  j=i+1;
	  while(fabs(pool[i].length-pool[j].length)<DELTA_ALPHA)
	  {
		  delta_CA = fabs(pool[i].length-pool[j].length);
		  graph1 = pool[i].graph;
		  graph2 = pool[j].graph;
		  across = 0;
		  if(graph1 != graph2)
		  {	
			i11 = pool[i].v1;
			i12 = pool[i].v2;
			i21 = pool[j].v1;
			i22 = pool[j].v2;
			if (graph1 == 1)
			{
			 if
			  (  
			   (the_same(g1[i11],g2[i21]) && the_same(g1[i12],g2[i22]))
			   || 
			   (across = (the_same(g1[i11],g2[i22]) && the_same(g1[i12],g2[i21])))
			  )
			  {
				neighbor_num1 = find_neighbor_num(g1,i11,i12);
				neighbor_num2 = find_neighbor_num(g2,i21,i22);
				if(
					(fabs(g1[i11].AB[neighbor_num1] - g2[i21].AB[neighbor_num2]) < DELTA_BETA) &&
					(fabs(g1[i11].BA[neighbor_num1] - g2[i21].BA[neighbor_num2]) < DELTA_BETA) &&
					(fabs(g1[i11].alpha[neighbor_num1] - g2[i21].alpha[neighbor_num2]) < THETA_CUTOFF)
				  )
				{
				  delta_CB1 = fabs(g1[i11].AB[neighbor_num1] - g2[i21].AB[neighbor_num2]);
				  delta_CB2 = fabs(g1[i11].BA[neighbor_num1] - g2[i21].BA[neighbor_num2]);
				  delta_alpha = fabs(g1[i11].alpha[neighbor_num1] - g2[i21].alpha[neighbor_num2]);
				  de[de_index].v11 = i11;
				  de[de_index].v12 = i12;
				  de[de_index].v21 = i21;
				  de[de_index].v22 = i22;
				  de[de_index].weight = 0.25 *
						  	  	 (
						  	  	  (1-delta_CA/DELTA_ALPHA)+
								  (1-delta_CB1/DELTA_BETA)+
								  (1-delta_CB2/DELTA_BETA)+
								  (1-delta_alpha/THETA_CUTOFF)
								 );

				  if(across)
				  {
					de[de_index].v21 = i22;
					de[de_index].v22 = i21;
				  }
				  de_index++;
				}
			  }
			}
			else
			{
			  if
				(
				(the_same(g2[i11],g1[i21]) && the_same(g2[i12],g1[i22]))
				||
				(across = (the_same(g2[i11],g1[i22]) && the_same(g2[i12],g1[i21])))
				)
			  {
				neighbor_num1 = find_neighbor_num(g2,i11,i12);
				neighbor_num2 = find_neighbor_num(g1,i21,i22);
				if(
					(fabs(g2[i11].AB[neighbor_num1] - g1[i21].AB[neighbor_num2]) < DELTA_BETA) &&
					(fabs(g2[i11].BA[neighbor_num1] - g1[i21].BA[neighbor_num2]) < DELTA_BETA) &&
					(fabs(g2[i11].alpha[neighbor_num1] - g1[i21].alpha[neighbor_num2]) < THETA_CUTOFF)
				  )
				{
				  delta_CB1 = fabs(g2[i11].AB[neighbor_num1] - g1[i21].AB[neighbor_num2]);
				  delta_CB2 = fabs(g2[i11].BA[neighbor_num1] - g1[i21].BA[neighbor_num2]);
				  delta_alpha = fabs(g2[i11].alpha[neighbor_num1] - g1[i21].alpha[neighbor_num2]);
				  de[de_index].v11 = i21;
				  de[de_index].v12 = i22;
				  de[de_index].v21 = i11;
				  de[de_index].v22 = i12;
				  de[de_index].weight = 0.25 *
								 (
								  (1-delta_CA/DELTA_ALPHA)+
								  (1-delta_CB1/DELTA_BETA)+
								  (1-delta_CB2/DELTA_BETA)+
								  (1-delta_alpha/THETA_CUTOFF)
								 );
				  if(across)
				  {
					de[de_index].v21 = i12;
					de[de_index].v22 = i11;
				  }
				  de_index++;
				}
			  }
			}
		  }
		  j++;
	  }
	}
	free(e1);
	free(e2);
	for(i=0;i<seq_length;i++){
		protein[i].res = g1[i].res;
		protein[i].chain = g1[i].chain;
		strcpy(protein[i].pdbnum,g1[i].pdbnum);
		protein[i].degree =0;
	}
	for(i=0;i<seq_length;i++){
		for(j=0;j<seq_length;j++){
			adjacency_matrix[i][j]=0;
		}
	}
	for(i=0;i<de_index;i++){
		for(j=0;j<de_index;j++){
			if(part_of_a_triangle(de,de_index,i,j))
			{
				adjacency_matrix[de[i].v11][de[i].v12]=de[i].weight;
				adjacency_matrix[de[i].v12][de[i].v11]=de[i].weight;
				adjacency_matrix[de[j].v11][de[j].v12]=de[j].weight;
				adjacency_matrix[de[j].v12][de[j].v11]=de[j].weight;
			}
		}
	}
	strcpy(output_file,file1);
	strcat(output_file,file2);
	strcat(output_file,".csv");
	if ((OUTPUT=fopen(output_file,"w")) == NULL ){
		printf("\nCan't open %s file.\n", output_file);
		return 0;
	}
	for(i=0;i<seq_length-1;i++){
		for(j=0;j<seq_length-1;j++){
			fprintf(OUTPUT,"%4.3f,",adjacency_matrix[i][j]);
		}
		fprintf(OUTPUT,"%4.3f\n",adjacency_matrix[i][seq_length-1]);
	}
	for(i=0;i<seq_length-1;i++){
		fprintf(OUTPUT,"%4.3f,",adjacency_matrix[seq_length-1][i]);
	}
	fprintf(OUTPUT,"%4.3f",adjacency_matrix[seq_length-1][seq_length-1]);

	fclose(OUTPUT);
	return 1;
}
int part_of_a_triangle(struct DOUBLE_EDGE * array, int length, int i, int j)
{
  int k;
  if
  (   
    (array[i].v11 == array[j].v11) &&
    (array[i].v21 == array[j].v21)
  ){
  	for (k = 0; k < length; k++) {
  		if
  		(
  			(array[k].v11 == array[i].v12) && //right from "i" is on the left "k"
  			(array[k].v21 == array[i].v22) &&
  			(array[k].v12 == array[j].v12) && //right from "j" is on the right "k"
  			(array[k].v22 == array[j].v22)
  		)
	    return 1;
  		if
  		(
  			(array[k].v12 == array[i].v12) && //right from "i" is on the right "k"
  			(array[k].v22 == array[i].v22) &&
  			(array[k].v11 == array[j].v12) && //right from "j" is on the left "k"
  			(array[k].v21 == array[j].v22)
  		)
	    return 1;
		}
  }

  if
  ( 
    (array[i].v12 == array[j].v12) &&
    (array[i].v22 == array[j].v22)
  ){
  	for (k = 0; k < length; k++) {
  		if
  		(
  			(array[k].v11 == array[i].v11) &&
  			(array[k].v21 == array[i].v21) &&
  			(array[k].v12 == array[j].v11) &&
  			(array[k].v22 == array[j].v21)
  		)
  		return 1;
  		if
  		(
  			(array[k].v12 == array[i].v11) && //right from "i" is on the right "k"
  			(array[k].v22 == array[i].v21) &&
  			(array[k].v11 == array[j].v11) && //right from "j" is on the left "k"
  			(array[k].v21 == array[j].v21)
  		)
  		return 1;
		}
  }
  return 0;
}
int find_neighbor_num(struct VERTEX *g, int n1, int n2)
{ 
  int i;

  for(i=0;i<g[n1].neinum;i++)
  { 
    if(g[n1].vertex[i] == n2)
    {
      return i;
    }
  }
  printf("Can't find number for %c-%s as a neighbor of %c-%s\n", g[n2].res, g[n2].pdbnum, g[n1].res, g[n1].pdbnum);
  return 0;
}


int the_same(struct VERTEX gr1, struct VERTEX gr2){

	if( strcmp(gr1.pdbnum,gr2.pdbnum) == 0
			&&
			gr1.chain == gr2.chain
			&&
			gr1.res == gr2.res){
		return 1;
	}
	return 0;
}


void Copy(struct EDGE *dest, struct EDGE *src)
{
  dest->graph=src->graph;
  dest->v1=src->v1;
  dest->v2=src->v2;
  dest->length=src->length;
  if (src->graph > 2)
      	printf("kuku\n");

}
int read_edges(struct EDGE *e, struct VERTEX *g, int graph_num)
{
  int i=0;
  int j;
  int edge_num=0;

  while (g[i].res != '\0')
  {
    j=0;
    for(j=0;j<g[i].neinum;j++)
    {
      e[edge_num].graph = graph_num;
      e[edge_num].v1 = i;
      e[edge_num].v2 = g[i].vertex[j];
      e[edge_num].length = g[i].AA[j];
      edge_num++;
    }
    i++;
  }
  return edge_num;
}
int compare_edge(struct EDGE * eA,struct EDGE * eB)
{
  if( (double) eA->length > (double) eB->length) return 1;
  else if ((double) eA->length < (double) eB->length) return -1;
  else return 0;
}
int StripLeadingSpaces (char * source, char * target){
	int i,s;
	s = strlen(source);
	if(s<1) return 0;
	for (i=0;i<s;i++){
		if(source[i]!=' '){
			strcpy(target, (char *)(source+i));
			return 1;
		}
	}
	return 1;
}
int readgraph(FILE *fpt, struct VERTEX *g)
{
  char *buffer;
  char *ptr;
  char tempstr[256];
  char tempstr2[256];
  int i,j;
  int number_of_vertices=0;
  buffer = (char *)calloc(MAXNEI*36,sizeof(char));
  i=0;
  while(1)
  {
    fgets (buffer, MAXNEI*36, fpt);
    j=0;
    if (feof(fpt))
      break;
    if
    (
	  (strlen(buffer)<3) ||
	  (buffer[0] == '#')
	  )
	    continue;
    else
	  {
    	number_of_vertices++;
  	  strtok(buffer," ");
  	  strcpy(tempstr,(char *)strtok(NULL,"-"));
  	  StripLeadingSpaces(tempstr,tempstr2);
  	  g[i].res = tempstr2[0];
  	  g[i].chain = tempstr2[2];
  	  strcpy(g[i].pdbnum,(char *) strtok(NULL," :"));
  	  while ((ptr = (char *)strtok(NULL," :()-\n")) != NULL)
 	    {	
	      g[i].vertex[j] = atoi(ptr) - 1;
	      strcpy(tempstr,(char *) strtok(NULL,"-"));
	      g[i].neires[j] = tempstr[0];
	      g[i].neichain[j] = tempstr[2];
	      strcpy(g[j].neipdbnum,(char *) strtok(NULL," "));
 	      strcpy(tempstr,(char *)strtok(NULL," ()\n"));
 	      g[i].AA[j] = atof(tempstr);
 	      strcpy(tempstr,(char *)strtok(NULL," ()\n"));
 	      g[i].AB[j] = atof(tempstr);
 	      strcpy(tempstr,(char *)strtok(NULL," ()\n"));
 	      g[i].BA[j] = atof(tempstr);
 	      strcpy(tempstr,(char *)strtok(NULL," ()\n"));
 	      g[i].alpha[j] = atof(tempstr);
 	      g[i].neinum++;
 	      j++;
 	    }
  	  i++;
  	}
  }
  free(buffer);
  return number_of_vertices;
}
