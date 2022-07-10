#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <sys/types.h>
#include <dirent.h>
#include <math.h>
#define MAX_NUM_OF_ATOMS 65536
#define MAX_NUM_OF_RESIDUES 2048
#define MAX_TYPES_OF_RESIDUES 256
#define PI 3.14159265
double DISTANCE;
struct VECTOR_DESCRIPTOR
	{
		char record_name[8];
		char three_letter_name [8];
		char one_letter_name;
		char atom_name_S[8];
		char atom_name_E[8];
	};
struct PDB_RECORD
{
	char record_name[8];
	int atom_number;
	char atom_name[8];
	char altloc;
	char residue_name[8];
	char chain_ID;
	int residue_seq_number;
	double x;
	double y;
	double z;
};
struct VERTEX
{
	char res;
	char chain;
	int sequence_number;
	double xS;
	double yS;
	double zS;
	double xE;
	double yE;
	double zE;
	int * neighbor;
	double * AA;
	double * AB;
	double * BA;
	double * alpha;
};
struct VECTOR
{
	double x;
	double y;
	double z;
};
struct PDB_RECORD * pdb_records;
struct VERTEX * g;
struct VECTOR_DESCRIPTOR * vector_descriptions;
int number_of_atoms;
int number_of_vectors;
int number_of_vertices;
int main(int argc, char *argv[] ) {
	int i,j,vn;
	int current_resnum;
	int temp_int;
	int nei_num;
	double d;
	double norm1,norm2;
	char input_filename[256];
	char vectors_filename[256];
	char chain[256];
	char diststr[256];
	char buffer[256];
	char tag[256];
	char temp_str[256];
	char temp_char;
	struct VECTOR Ca1Cb1;
	struct VECTOR Ca1Ca2;
	struct VECTOR Ca2Cb2;
	struct VECTOR CProduct1;
	struct VECTOR CProduct2;
	FILE * residue_table;
	FILE * pdbfile;
	FILE * graph;
	void read_from(char *,char *,int,int);
	char * strtoupper (char *);
	void trim_right(char * );
	int find_vector_start(int);
	int find_vector_end(int);
	char long_to_short(char *);
	double norm(struct VECTOR *);
	void CrossProduct(struct VECTOR *, struct VECTOR *, struct VECTOR *);
	double DotProduct(struct VECTOR *, struct VECTOR *);
	double Distance(double,double,double,double,double,double);
	printf("Graph of residues represented by vectors. v.3.1 (August 15 2019)\n"
			 "\tUsage: GraphGenerator.exe -fPdbfile -cChains -dDistance -rResidues\n"
			 "\tWhere Distance is a limit for the neighborhood (Angstroms)\n"
			 "\tPdbfile - filename for the protein coordinates\n"
			 "\tResidues - file containing residues representation\n"
			 "Written by Alexandr Kornev (2003-2019)\n");
	fflush( stdout );
	if (argc < 3)
	{
		printf("\n\tWrong parameters.\n");
	    return 0;
	}
	strcpy(chain," ");
	while ((argc > 1) && (argv[1][0] == '-'))
	{
		switch (argv[1][1])
	  	{
	    case 'f':
	      strcpy(input_filename,&argv[1][2]);
	      if ((pdbfile=fopen(input_filename,"r")) == NULL )
	      {
	    	printf("\n Can't open %s file.\n", &argv[1][2]);
	    	return 0;
	      }
	      break;
	    case 'r':
	      strcpy(vectors_filename,&argv[1][2]);
		  if ((residue_table=fopen(vectors_filename,"r")) == NULL )
		  {
			printf("\n Can't open %s file.\n", &argv[1][2]);
			return 0;
		  }
		  break;
	    case 'd':
	  	  DISTANCE = atof(&argv[1][2]);
	  	  strcpy(diststr,&argv[1][2]);
	  	  break;
	    default:
	  	  fprintf(stderr,"Bad option %s\n", argv[1]);
	  	  return 0;
	  	}
	    ++argv;
	    --argc;
	}
	if ((vector_descriptions = (struct VECTOR_DESCRIPTOR *)calloc(MAX_TYPES_OF_RESIDUES,sizeof(struct VECTOR_DESCRIPTOR )))==NULL)
	{
		printf("Can't allocate memory for vector descriptions\n");
	    return 0;
	}
	i=0;
	while (!feof(residue_table))
	{
		fgets (buffer,256,residue_table);
	    if
		(
	    		(buffer[0] != '#') &&
				(!feof(residue_table)) &&
				strlen(buffer)>1
		)
	    {
	        strcpy(vector_descriptions[i].record_name, (char *) strtok(buffer," ;\n\r"));
	        strcpy(vector_descriptions[i].three_letter_name, (char *) strtok(NULL," ;\n\r"));
	        strcpy(temp_str, (char *) strtok(NULL," ;"));
	        vector_descriptions[i].one_letter_name = temp_str[0];
	        strcpy(vector_descriptions[i].atom_name_S,(char *) strtok(NULL," ;\n\r"));
	        strcpy(vector_descriptions[i].atom_name_E,(char *) strtok(NULL," ;\n\r"));
	        i++; /* next residue */
	    }
	}
	number_of_vectors=i;
	if ((pdb_records = (struct PDB_RECORD *)calloc(MAX_NUM_OF_ATOMS,sizeof(struct PDB_RECORD )))==NULL)
	{
		printf("Can't allocate memory for pdb records\n");
		return 0;
	}

	i=0;
	while (!feof(pdbfile))
	{
	    fgets (buffer,256,pdbfile);
	    if(feof(pdbfile)) break;
	    read_from(buffer,tag,0,6);
	    if
		(
				strcmp(tag,"ATOM  ") == 0 ||
				strcmp(tag,"HETATM") == 0
		)
	    {
	    	strcpy(pdb_records[i].record_name, tag);
	    	read_from(buffer, temp_str, 6, 5);
	    	pdb_records[i].atom_number = atoi(temp_str);
	    	read_from(buffer, temp_str, 13, 4);
	    	trim_right(temp_str);
	    	strcpy(pdb_records[i].atom_name,temp_str);
	    	read_from(buffer, temp_str, 17, 3);
	    	pdb_records[i].altloc=buffer[16];
	    	strcpy(pdb_records[i].residue_name, temp_str);
	    	pdb_records[i].chain_ID = buffer[21];
	    	read_from(buffer, temp_str, 22, 4);
	    	pdb_records[i].residue_seq_number=atoi(temp_str);
	    	read_from(buffer, temp_str, 30, 8);
	    	pdb_records[i].x = atof(temp_str);
	    	read_from(buffer, temp_str, 38, 8);
	    	pdb_records[i].y = atof(temp_str);
	    	read_from(buffer, temp_str, 46, 8);
	    	pdb_records[i].z = atof(temp_str);
	    	i++;
	    }
	}
	number_of_atoms=i;
	if ((g = (struct VERTEX *)calloc(MAX_NUM_OF_RESIDUES,sizeof(struct VERTEX )))==NULL)
	{
		printf("Can't allocate memory for graph\n");
		return 0;
	}
	current_resnum = pdb_records[0].residue_seq_number;
	number_of_vertices =0;
	for(i=0;i<number_of_atoms;i++)
	{
		temp_char = long_to_short(pdb_records[i].residue_name);
		if (
			(i==0
			||
			(current_resnum != pdb_records[i].residue_seq_number))
			&&									
			temp_char != 'X')							
		{
			current_resnum=pdb_records[i].residue_seq_number;
			if(temp_char == 'X') break;
			g[number_of_vertices].sequence_number = pdb_records[i].residue_seq_number;
			g[number_of_vertices].chain = pdb_records[i].chain_ID;
			g[number_of_vertices].res = long_to_short(pdb_records[i].residue_name);
			temp_int = find_vector_start(i);
			g[number_of_vertices].xS=pdb_records[temp_int].x;
			g[number_of_vertices].yS=pdb_records[temp_int].y;
			g[number_of_vertices].zS=pdb_records[temp_int].z;
			temp_int = find_vector_end(i);
			g[number_of_vertices].xE=pdb_records[temp_int].x;
			g[number_of_vertices].yE=pdb_records[temp_int].y;
			g[number_of_vertices].zE=pdb_records[temp_int].z;
			number_of_vertices++;
		}
	}
	for (i=0;i<number_of_vertices;i++)
	{
		if ((g[i].neighbor = (int *)calloc(number_of_vertices,sizeof(int)))==NULL)
		{
			printf("Can't allocate memory for neighbors\n");
			return 0;
		}
		if ((g[i].AA = (double *)calloc(number_of_vertices,sizeof(double)))==NULL)
		{
			printf("Can't allocate memory for AA\n");
			return 0;
		}
		if ((g[i].AB = (double *)calloc(number_of_vertices,sizeof(double)))==NULL)
		{
			printf("Can't allocate memory for AB\n");
			return 0;
		}
		if ((g[i].BA = (double *)calloc(number_of_vertices,sizeof(double)))==NULL)
		{
			printf("Can't allocate memory for BA\n");
			return 0;
		}
		if ((g[i].alpha = (double *)calloc(number_of_vertices,sizeof(double)))==NULL)
		{
			printf("Can't allocate memory for alpha\n");
			return 0;
		}
	}
	for(i=0;i<number_of_vertices;i++)
	{
		nei_num=0;
	    for(j=i;j<number_of_vertices;j++)
		{
	    	if
		    (
		     (i!=j) &&
		     ((d = Distance(g[i].xS,g[i].yS,g[i].zS,g[j].xS,g[j].yS,g[j].zS))<DISTANCE)
		    )
		    {
	    		g[i].neighbor[nei_num] = j+1;
	    		g[i].AA[nei_num] = d;
	    		g[i].AB[nei_num]=Distance(g[i].xS,g[i].yS,g[i].zS,g[j].xE,g[j].yE,g[j].zE);
	    		g[i].BA[nei_num]=Distance(g[i].xE,g[i].yE,g[i].zE,g[j].xS,g[j].yS,g[j].zS);
	    		Ca1Cb1.x = g[i].xE - g[i].xS;
	    		Ca1Cb1.y = g[i].yE - g[i].yS;
	    		Ca1Cb1.z = g[i].zE - g[i].zS;
	    		Ca1Ca2.x = g[j].xS - g[i].xS;
	    		Ca1Ca2.y = g[j].yS - g[i].yS;
	    		Ca1Ca2.z = g[j].zS - g[i].zS;
	    		CrossProduct(&Ca1Ca2,&Ca1Cb1,&CProduct1);
	    		Ca2Cb2.x = g[j].xE - Ca1Ca2.x;
	    		Ca2Cb2.y = g[j].yE - Ca1Ca2.y;
	    		Ca2Cb2.z = g[j].zE - Ca1Ca2.z;
	    		Ca2Cb2.x = Ca2Cb2.x - g[i].xS;
	    		Ca2Cb2.y = Ca2Cb2.y - g[i].yS;
	    		Ca2Cb2.z = Ca2Cb2.z - g[i].zS;
	    		CrossProduct(&Ca1Ca2,&Ca2Cb2,&CProduct2);
	    		norm1 = norm(&CProduct1);
	    		norm2 = norm(&CProduct2);
	    		g[i].alpha[nei_num] = (acos(DotProduct(&CProduct1,&CProduct2)/(norm1*norm2))/PI)*180;
	    		if (DotProduct(&CProduct1,&Ca2Cb2) < 0)
	    			g[i].alpha[nei_num] *= -1;
	    		nei_num++;
		    }
		  }
	  }
	  strcat(input_filename,".gr");
	  graph = fopen(input_filename,"w");
	  fprintf(graph,"#Adjacency list for %s. Residue vectors defined in %s."
		  "Cut-off for CA-CA distance %3.1fA.\n"
		  "#Format: Vertex_number Residue(Chain)-PDB_number "
		  ": Neighbor vertex_number (AA AB BA theta)\n",
		  input_filename, vectors_filename, DISTANCE);
	  for (i=0;i<number_of_vertices;i++)
	  {
	    fprintf(graph,"%-4i %c(%c)-%-5i:",i+1,g[i].res,g[i].chain,g[i].sequence_number);
	    /* its neighbors */
	    j=0;
	    while((vn = g[i].neighbor[j]) != 0)
		  {
		    fprintf(graph,"%4i %c(%c)-%-5i (%5.2f %5.2f %5.2f %6.1f)",
			  vn,g[vn-1].res,g[vn-1].chain,g[vn-1].sequence_number,g[i].AA[j],g[i].AB[j],g[i].BA[j],g[i].alpha[j]);
		    j++;
		   }
	     fprintf (graph,"\n");
	   }
	   fclose (graph);
	return 1;
}
char * strtoupper (char * string)
{
  int i,length;

  length=strlen(string);
  for(i=0;i<length;i++)
    string[i] = toupper(string[i]);
  return string;
}
void read_from (char* source, char* target, int offset, int length)
{
  int i;
  for (i=0;i<length;i++)
	target[i] = source[i+offset];
  target[length]=0;
}
void trim_right(char * string){
	int i;

	i = strlen (string) - 1;
	while (string[i]==' '){
		string[i] = '\0';
		i--;
	}
}
char long_to_short(char * long_name )
{
	int i;

	for (i=0;i<number_of_vectors;i++)
	{
		if (strcmp(vector_descriptions[i].three_letter_name,long_name)==0)
			return vector_descriptions[i].one_letter_name;
	}
	printf("Can't find residue %s in the table, skipping...\n", long_name);
	return 'X';
}
int find_vector_start(int atom_num)
{
	int i;
	char atom_type[8];

	strcpy(atom_type,"NF");
	i=0;
	while(i<number_of_vectors)
	{
		if(strcmp(vector_descriptions[i].three_letter_name,pdb_records[atom_num].residue_name)==0)
		{
			strcpy(atom_type,vector_descriptions[i].atom_name_S);
			break;
		}
		i++;
	}
	if(strcmp(atom_type, "NF")==0) printf("Residue for atom %i was not found\n",atom_num);

	i=atom_num;
	while(i<number_of_atoms)
	{
		if( strcmp (pdb_records[i].atom_name, atom_type) == 0)
			return i;
		i++;
	}
	printf("Starting atom for %i was not found\n", atom_num);
	return 0;
}
int find_vector_end(int atom_num)
{
	int i;
	char atom_type[8];

	strcpy(atom_type,"NF");
	i=0;
	while(i<number_of_vectors)
	{
		if(strcmp(vector_descriptions[i].three_letter_name,pdb_records[atom_num].residue_name)==0)
		{
			strcpy(atom_type,vector_descriptions[i].atom_name_E);
			break;
		}
		i++;
	}
	if(strcmp(atom_type, "NF")==0) printf("Residue for atom %i was not found\n",atom_num);
	i=atom_num;
	while(i<number_of_atoms)
	{
		if( strcmp (pdb_records[i].atom_name, atom_type) == 0)
			return i;
		i++;
	}
	printf("Ending atom for %i was not found\n", atom_num);
	return 0;
}
double norm(struct VECTOR * v)
{
  return (sqrt(v->x*v->x + v->y*v->y + v->z*v->z));
}
void CrossProduct(struct VECTOR *v1, struct VECTOR *v2, struct VECTOR *result)
{
  result->x = v1->y*v2->z - v1->z*v2->y;
  result->y = v1->z*v2->x - v1->x*v2->z;
  result->z = v1->x*v2->y - v1->y*v2->x;
}
double DotProduct(struct VECTOR *v1, struct VECTOR *v2)
{
  return (v1->x*v2->x + v1->y*v2->y + v1->z*v2->z);
}

double Distance(double x1,double y1,double z1,double x2,double y2,double z2)
{
  double x,y,z,d;

  x = x1 - x2;
  y = y1 - y2;
  z = z1 - z2;
  d = sqrt(x*x + y*y + z*z);
  return d;
}
