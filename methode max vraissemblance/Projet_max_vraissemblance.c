#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <float.h>

//#define N 100
#define R 500
#define sigma_m 5 //entre 0.2 et 5
#define sigma_y 5 //entre 0.2 et 5


double gauss(void)
{
  static double x0, x1;
  static unsigned int nb_ready;

  double u, v, w, z;

  if (nb_ready == 0)
  {
    u = ((double) rand()) / RAND_MAX;
    v = ((double) rand()) / RAND_MAX;
    w = sqrt(-2.0 * log(u));
    z = 2 * M_PI * v;
    x0 = w * cos(z);
    x1 = w * sin(z);

    nb_ready = 1;
    return x0;
  }
  else
  {
    nb_ready = 0;
    return x1;
  }
};



int main()
{
	FILE *fichier = NULL;
	fichier = fopen("échantillon1.txt","w");
	unsigned int N;
	srand(time(NULL));
	double reussite;
	long double pr_hm_hy_sachant_k;
	unsigned int hw[256] =
	{
	  0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4, 1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
	  1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
	  1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
	  2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
	  1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
	  2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
	  2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
	  3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7, 4, 5, 5, 6, 5, 6, 6, 7, 5, 6, 6, 7, 6, 7, 7, 8
	};
	unsigned char SBOX[256] =
	{
	   0x63, 0x7C, 0x77, 0x7B, 0xF2, 0x6B, 0x6F, 0xC5, 0x30, 0x01, 0x67, 0x2B, 0xFE, 0xD7, 0xAB, 0x76,
	   0xCA, 0x82, 0xC9, 0x7D, 0xFA, 0x59, 0x47, 0xF0, 0xAD, 0xD4, 0xA2, 0xAF, 0x9C, 0xA4, 0x72, 0xC0,
	   0xB7, 0xFD, 0x93, 0x26, 0x36, 0x3F, 0xF7, 0xCC, 0x34, 0xA5, 0xE5, 0xF1, 0x71, 0xD8, 0x31, 0x15,
	   0x04, 0xC7, 0x23, 0xC3, 0x18, 0x96, 0x05, 0x9A, 0x07, 0x12, 0x80, 0xE2, 0xEB, 0x27, 0xB2, 0x75,
	   0x09, 0x83, 0x2C, 0x1A, 0x1B, 0x6E, 0x5A, 0xA0, 0x52, 0x3B, 0xD6, 0xB3, 0x29, 0xE3, 0x2F, 0x84,
	   0x53, 0xD1, 0x00, 0xED, 0x20, 0xFC, 0xB1, 0x5B, 0x6A, 0xCB, 0xBE, 0x39, 0x4A, 0x4C, 0x58, 0xCF,
	   0xD0, 0xEF, 0xAA, 0xFB, 0x43, 0x4D, 0x33, 0x85, 0x45, 0xF9, 0x02, 0x7F, 0x50, 0x3C, 0x9F, 0xA8,
	   0x51, 0xA3, 0x40, 0x8F, 0x92, 0x9D, 0x38, 0xF5, 0xBC, 0xB6, 0xDA, 0x21, 0x10, 0xFF, 0xF3, 0xD2,
	   0xCD, 0x0C, 0x13, 0xEC, 0x5F, 0x97, 0x44, 0x17, 0xC4, 0xA7, 0x7E, 0x3D, 0x64, 0x5D, 0x19, 0x73,
	   0x60, 0x81, 0x4F, 0xDC, 0x22, 0x2A, 0x90, 0x88, 0x46, 0xEE, 0xB8, 0x14, 0xDE, 0x5E, 0x0B, 0xDB,
	   0xE0, 0x32, 0x3A, 0x0A, 0x49, 0x06, 0x24, 0x5C, 0xC2, 0xD3, 0xAC, 0x62, 0x91, 0x95, 0xE4, 0x79,
	   0xE7, 0xC8, 0x37, 0x6D, 0x8D, 0xD5, 0x4E, 0xA9, 0x6C, 0x56, 0xF4, 0xEA, 0x65, 0x7A, 0xAE, 0x08,
	   0xBA, 0x78, 0x25, 0x2E, 0x1C, 0xA6, 0xB4, 0xC6, 0xE8, 0xDD, 0x74, 0x1F, 0x4B, 0xBD, 0x8B, 0x8A,
	   0x70, 0x3E, 0xB5, 0x66, 0x48, 0x03, 0xF6, 0x0E, 0x61, 0x35, 0x57, 0xB9, 0x86, 0xC1, 0x1D, 0x9E,
	   0xE1, 0xF8, 0x98, 0x11, 0x69, 0xD9, 0x8E, 0x94, 0x9B, 0x1E, 0x87, 0xE9, 0xCE, 0x55, 0x28, 0xDF,
	   0x8C, 0xA1, 0x89, 0x0D, 0xBF, 0xE6, 0x42, 0x68, 0x41, 0x99, 0x2D, 0x0F, 0xB0, 0x54, 0xBB, 0x16
	};
	

	long double M[256][9][9];
	//initialisation des modèles théoriques
	for(unsigned int i = 0; i<9 ; i++)
	{
		for(unsigned int j = 0; j<9; j++)
		{
			for(unsigned int m = 0; m<256; m++)
			{
				M[m][i][j] = 0;
			}
		}
	}


	//calcul des modeles pour une clé fixée
	for(unsigned int k=0; k<256; k++)
	{
		for(unsigned int m = 0; m<256; m++)
		{
			M[k][hw[m]][hw[SBOX[m^k]]]++;
		}
		//normalisation
		for(unsigned int i = 0; i<9; i++)
		{
			for(unsigned int j = 0; j<9; j++)
			{
				M[k][i][j] = M[k][i][j]/256;
			}
		}
	}

	double lm;
	double ly;
	unsigned int nombre_succes[91];
	//double constante = (1/(sigma_m*sqrt(2*M_PI))) * (1/(sigma_y*sqrt(2*M_PI)));
	for(unsigned int i=0; i<91; i++)
	{
		nombre_succes[i] = 0;
	}

	for(unsigned int tour =0; tour<R; tour++)
	{
		//initialisation du tableau de probabilité pour chaque clé
		long double Pr[256];
		for(unsigned int i =0; i<256; i++)
		{
			Pr[i] = 1000000;
		}

		unsigned int key = rand() % 256;
		for(N=1; N<=1000; N++)
		{
			//calcul des courants bruités + de la probabilité de clé au fur et à mesure
			unsigned int m = rand()%256;
			lm = hw[m] + sigma_m*gauss();
			ly = hw[SBOX[m^key]] + sigma_y*gauss();
			for(unsigned int k=0; k<256; k++)
			{
				pr_hm_hy_sachant_k = (long double)0;
				for(unsigned int i=0; i<9; i++)
				{
					for(unsigned int j=0; j<9; j++)
					{
						// ici on devrait multiplier par constante mais celle-ci ne dépend pas de la clé et on perd juste un coefficient de proportionnalité qui est le même pour toutes les probabilités
						pr_hm_hy_sachant_k += M[k][i][j] * exp(-0.5 * ((lm - i)/sigma_m) * ((lm - i)/sigma_m)) * exp(-0.5 * ((ly - j)/sigma_y) * ((ly - j)/sigma_y));				
					}
				}
				Pr[k] = Pr[k] * pr_hm_hy_sachant_k * 5;
			}
			

			if((N>=100) & ((N%10) == 0))
			{
				//calcul de la clé avec le maximum de vraissemblance
				unsigned int cle_suppose = 0;
				long double proba = Pr[0];
				for(unsigned int k=1; k<256; k++)
				{
					if(Pr[k]-proba > 0)
					{
						proba = Pr[k];
						cle_suppose = k;
					}
				}
				if(cle_suppose == key)
				{
					nombre_succes[(N-100)/10]++;
				}
			}
		}
		printf("%d\n",tour);
	}
	for(unsigned int i=0; i<91; i++)
	{
		printf("nombre de succés : %d sur %d attaques\nPour %d messages par attaques\n",nombre_succes[i],R,(100+10*i));
		reussite = (double)nombre_succes[i]/(double)R;
		printf("Le taux de succés est donc de %f\n\n",reussite);
		if(fichier != NULL)
		{
			fprintf(fichier,"%d %f\n",(100+10*i),reussite);
		}
	}
	fclose(fichier);
	return(0);
};