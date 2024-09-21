#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>



void setup(char *filename);
void setup2(char *filename);
int rank1(int i, char c,int choice);
int rank(char c, int i, int choice);
int count(char *query);
int find_mem(int x, int min_len, int len, char* query);
void forward_backward(char *query);
void set_interval(int* s, int*e, char c, int choice);
void reverse(char* filename);
void generate_bwt(char* filename, char* output);

char *bitvector_a;
char *bitvector_c;
char *bitvector_g;
char *bitvector_T;

uint16_t *miniheaders_a;
uint16_t *miniheaders_c;
uint16_t *miniheaders_g;
uint16_t *miniheaders_t;

uint64_t *macroheaders_a;
uint64_t *macroheaders_c;
uint64_t *macroheaders_g;
uint64_t *macroheaders_t;

int pref[1000000][4];
int dollarPos;
int C1[4], C2[4];
uint64_t masks[64];
int size = 1000000; //********* smake sure to change this to 1e6  after testing *****************
int MIN_LEN = 40;


char T[1000000];

int suffixCmp(const void *a, const void *b) {
	return(strcmp(&T[*(int *) a], &T[*(int *) b]));
}

void print_bitvectors(){
  printf("\n");

 for(int v = 0; v<4;v++){
   char* bitvector;
   
   if(v == 0)
	 bitvector = bitvector_a;
   else if(v == 1)
	 bitvector = bitvector_c;
   else if(v == 2)
	 bitvector = bitvector_g;
   else
	 bitvector = bitvector_T;

   printf("bit vector %d: ",v);
   
  for(int i = 0;i<10;i++)
	for(int j = 0;j<7;j++){
	  if((bitvector[i] & (1<<j)) != 0 ){
		printf("1");
	  }else
		printf("0");
	}
	  printf("\n");
  }
}

int main(int argc, char *argv[]) {
    bitvector_a = (char *) malloc(size / 8);
    bitvector_c = (char *) malloc(size / 8);
    bitvector_g = (char *) malloc(size / 8);
    bitvector_T = (char *) malloc(size / 8);
	
	miniheaders_a = (uint16_t *) malloc(1000000 / 32);
	miniheaders_c = (uint16_t *) malloc(1000000 / 32);
	miniheaders_g = (uint16_t *) malloc(1000000 / 32);
	miniheaders_t = (uint16_t *) malloc(1000000 / 32);

	macroheaders_a = (uint64_t *) malloc((1065535 / 65536) * 8);
	macroheaders_c = (uint64_t *) malloc((1065535 / 65536) * 8);
	macroheaders_g = (uint64_t *) malloc((1065535 / 65536) * 8);
	macroheaders_t = (uint64_t *) malloc((1065535 / 65536) * 8);

	generate_bwt(argv[1], "DNA_BWT.txt"); // use the DNA text file to build the BWT. This is utilized in the index used for backward extension

	printf("done with generating the BWT for DNA.txt file.....\n");

	reverse(argv[1]); // reverse the DNA text file

	printf("done with reversing the DNA.txt file.....\n");

	generate_bwt("DNAR.txt", "DNAR_BWT.txt"); // generate the BWT for the reverse DNA text file. This is utilized in the index used for forward extension

	printf("done with generating the BWT for reversed DNA text.....\n");

	
	setup("DNA_BWT.txt");
	setup2("DNAR_BWT.txt");

	printf("set up done...\n\n");


	//print_bitvectors(); ******************************comment this while testing for the original text*********************************************
	
	C1[0] = rank('A',size-1,0);
	C1[1] = rank('C',size-1,0);
	C1[2] = rank('G',size-1,0);

	C2[0] = rank('A',size-1,1);
	C2[1] = rank('C',size-1,1);
	C2[2] = rank('G',size-1,1);
	
	char query[101];

	while (scanf("%100s", &query) == 1) {
	  //printf("%i occurrences of %s.\n", count(query), query);
	  printf("The mems are as follows:\n");
	  forward_backward(query);
	  printf("\n");
	  
	}

	return 0;
}

void generate_bwt(char* filename, char* output){
    int SA[size];

	
    FILE *inputFile = fopen(filename, "r");
    FILE *outputFile = fopen(output, "w");
	
	char ch;
	int i = 0,length = 0;
	
	while((ch = fgetc(inputFile)) != EOF){
	  T[i] = ch;
	  SA[i] = i;
	  i++; length++;
	}
  		
	qsort(SA, i, sizeof(int), suffixCmp);

	for ( i = 0; i < length; i++) {
	  fputc(T[(SA[i] + (length-1)) % length], outputFile);
	}

	fflush(outputFile);
	fclose(inputFile);
	fclose(outputFile);
}

void reverse(char* filename){
  char c[size];

  FILE *inputFile = fopen(filename,"r");
  FILE *outputFile = fopen("DNAR.txt","w");

  char ch;
  int cnt = 0;
  while((ch = fgetc(inputFile)) != EOF){
	c[cnt] = ch;
	cnt++;
  }
  cnt--;
  for(int i = 0;i<cnt/2;i++){
	char temp = c[cnt-i-1];
	c[cnt-i-1] = c[i];
	c[i] = temp;
  }
  
  c[cnt+1] = '$';
  
  for(int i = 0;i<cnt+1;i++){
	fputc(c[i],outputFile);
  }

  
	fflush(outputFile);
	fclose(inputFile);
	fclose(outputFile);
}

void setup(char *filename) {
			   
	FILE *file = fopen(filename, "r");
	int currentRank[4] = {0};	
	int miniRank[4] = {0};

	memset(bitvector_a, 0, 1000000 / 8);
	memset(bitvector_c, 0, 1000000 / 8);
    memset(bitvector_g, 0, 1000000 / 8);
	memset(bitvector_T, 0, 1000000 / 8);


	for (int i = 0; i < size; i++) {

		if (i % 64 == 0) {
		  miniheaders_a[i/64] = (uint16_t) (miniRank[0]);
		  miniheaders_c[i/64] = (uint16_t) (miniRank[1]);
		  miniheaders_g[i/64] = (uint16_t) (miniRank[2]);
		  miniheaders_t[i/64] = (uint16_t) (miniRank[3]);
		}

		if (i % 65536 == 0) {
		    miniRank[0] = 0; miniRank[1] = 0; miniRank[2] = 0; miniRank[3] = 0;
			macroheaders_a[i / 65536] = (uint64_t) currentRank[0];
			macroheaders_c[i / 65536] = (uint64_t) currentRank[1];
			macroheaders_g[i / 65536] = (uint64_t) currentRank[2];
			macroheaders_t[i / 65536] = (uint64_t) currentRank[3];
		}

		char c = (char) fgetc(file);
		
		switch (c) {
            case 'A':
                currentRank[0]++;
                bitvector_a[i / 8] |= (1 << (i % 8));
				//pref[i][0]++;
                miniRank[0]++;
                break;
            case 'C':
                currentRank[1]++;
                bitvector_c[i / 8] |= (1 << (i % 8));
				//pref[i][1]++;
                miniRank[1]++;
                break;
            case 'G':
                currentRank[2]++;
                bitvector_g[i / 8] |= (1 << (i % 8));
				//pref[i][2]++;
                miniRank[2]++;
                break;
            case 'T':
                currentRank[3]++;
                bitvector_T[i / 8] |= (1 << (i % 8));
				//pref[i][3]++;
                miniRank[3]++;
                break;
            case '$':
                dollarPos = i;
                break;
        }


	}
	

	fclose(file);

	memset(&masks[63], 255, 8);

	for (int i = 62; i >= 0; i--) {
		memcpy(&masks[i], &masks[i + 1], 8);
		char * charMask = (char *) &masks[i];
		charMask[(i + 1) / 8] = charMask[(i + 1) / 8] & ~(1 << (i + 1) % 8);
	}

	return;
}

void setup2(char* filename){
  FILE *file = fopen(filename,"r");
  memset(pref, 0, sizeof(pref));

  //using the prefix array instead of bitvectors. 
  for(int i = 0;i<size;i++){
	char c = fgetc(file);
		switch (c) {
            case 'A':
				pref[i][0]++;
                break;
            case 'C':
				pref[i][1]++;
                break;
            case 'G':
				pref[i][2]++;
                break;
            case 'T':
				pref[i][3]++;
                break;
            case '$':
                dollarPos = i;
                break;
        }
	for(int j = 0;j<4;j++){
		  if(i == 0)
			continue;
		  pref[i][j] += pref[i-1][j];
	 }
	
  }
  
}

int rank1(int i, char c, int choice) {

  int sum = 0;
  if(choice == 0){
	switch (c) {
            case 'A':
                sum = (int) miniheaders_a[i / 64] + (int) macroheaders_a[i / 65536];
			    sum += __builtin_popcountll(*(uint64_t *) &bitvector_a[(i / 64) * 8] & masks[i % 64]);
                break;
            case 'C':
                sum = (int) miniheaders_c[i / 64] + (int) macroheaders_c[i / 65536];
			    sum += __builtin_popcountll(*(uint64_t *) &bitvector_c[(i / 64) * 8] & masks[i % 64]);

                break;
            case 'G':
                sum = (int) miniheaders_g[i / 64] + (int) macroheaders_g[i / 65536];
			    sum += __builtin_popcountll(*(uint64_t *) &bitvector_g[(i / 64) * 8] & masks[i % 64]);
                break;
            case 'T':
                sum = (int) miniheaders_t[i / 64] + (int) macroheaders_t[i / 65536];
				sum += __builtin_popcountll(*(uint64_t *) &bitvector_T[(i / 64) * 8] & masks[i % 64]);
                break;
        }
  }
  else{
	switch (c) {
        case 'A':
		  sum = pref[i][0];
            break;
        case 'C':
		  sum = pref[i][1];
            break;
        case 'G':
		  sum = pref[i][2];
            break;
        case 'T':
          sum = pref[i][3];
            break;
	 }
  }
  
	return(sum);
}

int rank(char c, int i, int choice) {
	if (i < 0) {
		return(0);
	}

	int rank1i = rank1(i, c, choice);

	return rank1i;	
}

void forward_backward(char *query){
  int x = 0;
  int len = strlen(query);
  
  do {
	x = find_mem(x, MIN_LEN ,len,query);
  } while (x < len);
}


int find_mem(int x, int min_len, int len, char* query)
{
  int i,j;
  int  s = 0;
  int e = size-1;

  
  if (x + min_len - 1 >= len) return len; //making sure we have enough room to extend

  for (i = x; i < len && s<=e ; i++) { // forward extension: finds the end of the mem starting at x
	  set_interval(&s,&e, query[i], 1);
  }
  
  if(s>e) // move i back to the position where mismatch occurred if any
	i--;
  
  if(i-x >= min_len){ // prints the mem if it satisfies the min length
   for(int z = x;z<i;z++){
	  printf("%c",query[z]);
	}
	printf("\n");
  }
  
  s = 0,e = size-1;
  for(j = i; j>=x && i<len && s<=e; j--) { // backward extension: finds the beginning of the next mem
	set_interval(&s,&e, query[j],0);
  }

  return j+2; //starting of the mem because j is decreemented even after e>s
}

void set_interval(int* s, int* e, char c, int choice){
  int *C = choice == 0?C1:C2;

  	  switch (c) {
            case 'A':
			  *s = rank('A', *s - 1, choice) + 1;
              *e = rank('A', *e, choice);
              break;
            case 'C':
			  *s = C[0]+rank('C', *s - 1, choice) + 1;
			  *e = C[0]+rank('C', *e, choice);
              break;
            case 'G':
			  *s = C[0]+C[1]+rank('G', *s - 1, choice) + 1;
              *e = C[0]+C[1]+rank('G', *e, choice);
              break;
            case 'T':
			  *s = C[0]+C[1]+C[2]+rank('T', *s - 1, choice) + 1;
			  *e = C[0]+C[1]+C[2]+rank('T', *e, choice);
              break;
        }
}
/*
int count(char *query) {
	int m = strlen(query);

	int s = 0;
	int e = size-1;

	for (int i = 0; i <= m-1 && s <= e; i++) {
	  switch (query[i]) {
            case 'A':
			  s = rank('A', s - 1,0) + 1;
                e = rank('A', e,0);
                break;
            case 'C':
			  s = C[0]+rank('C', s - 1,0) + 1;
			  e = C[0]+rank('C', e,0);
                break;
            case 'G':
			  s = C[0]+C[1]+rank('G', s - 1,0) + 1;
                e = C[0]+C[1]+rank('G', e,0);
                break;
            case 'T':
			  s = C[0]+C[1]+C[2]+rank('T', s - 1,0) + 1;
			  e = C[0]+C[1]+C[2]+rank('T', e,0);
                break;
        }
	}

	return(e - s + 1);
}
*/
