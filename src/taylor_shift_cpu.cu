/* This file is part of the CUMODP library

    CUMODP is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    CUMODP is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with CUMODP.  If not, see <http://www.gnu.org/licenses/>.

    Copyright: Sardar Anisul Haque <shaque4@uwo.ca>
               Xin Li <xli.software@gmail.com>
               Farnam Mansouri <mansouri.farnam@gmail.com>
               Davood Mohajerani <dmohajer@uwo.ca>
               Marc Moreno Maza  <moreno@csd.uwo.ca>
               Wei Pan <wei.pan@intel.com>
               Ning Xie <nxie6@csd.uwo.ca>
*/



#include "taylor_shift_conf.h"
#include "taylor_shift_cpu.h"

/** error message if there is a lack of arguments to make the program
* @m : number of arguments
**/
void error_message(sfixn m)
{
  if (m < 3)
    {
      printf("********** ERROR, not enough arguments ! **********\nThe program works wi\
th the following parameters:\n\n");
      printf("1st parameter : file containing coefficients of the polynomial you want t\
o consider.\n");
      printf("2nd parameter : prime number p.\n");

      exit(1);
    }
}

/** error message if there is a lack of arguments to make the program
* @m : number of arguments
**/
void error_message2(sfixn m)
{
  if (m < 2)
    {
      printf("********** ERROR, not enough arguments ! **********\nThe program works wi\
th the following parameters:\n\n");
      printf("1st parameter : file containing coefficients of the polynomial you want t\
o consider.\n");
      exit(1);
    }
}


/** computes the nomber of blocks needed in computation in parallel
* @n : number of instructions to do in parallel
**/
sfixn number_of_blocks(sfixn n)
{
  sfixn res;
  res = n/NB_THREADS;
  if ( n % NB_THREADS != 0)
    res++;
  return res;
}


/** stocks a file in an array
* @filename : name of the file
* @n : number of instructions to do in parallel
* @a : number of instructions to do in parallel
**/
void stock_file_in_array(char* filename, sfixn n, sfixn* & a)
{
  ifstream data_file;
  sfixn i;
  data_file.open(filename);

  if (! data_file.is_open())
    {
      printf("\n Error while reading the file %s. Please check if it exists !\n", filen\
ame);
      exit(1);
    }

  a = (sfixn*) malloc (n*sizeof(sfixn));

  for (i=0; i<n; i++)
    data_file >> a[i];

  data_file.close();
}


// stockes the array of Newton's coefficients in a file
void stock_array_in_file(const char *name_file, sfixn *T, sfixn size)
{
  sfixn i;
  FILE* file = NULL;

  file = fopen(name_file, "w+");
  if (file == NULL)
  {
    printf("error when opening the file !\n");
    exit(1);
  }

  // writting the file
  fprintf(file, "%d", T[0]);
  for (i=1; i<size; i++)
    fprintf(file, "\n%d", T[i]);
  fclose(file);
}


// computes the number of lines of a file
sfixn size_file(char* filename)
{
  sfixn size = 0;
  ifstream in(filename);
  std::string line;

  while(std::getline(in, line))
    size++;
  in.close();

  return size;
}


// display of an array
void display_array(sfixn *T, sfixn size)
{
  sfixn k;
  printf("[ ");
  for (k=0; k<size; k++)
    printf("%d ", T[k]);
  printf("] \n");
}


// addition of two arrays
void add_arrays(sfixn *res, sfixn *T1, sfixn *T2, sfixn size, sfixn p)
{
  sfixn i;
  for (i=0; i<size; i++)
    res[i] = (T1[i] + T2[i]) % p;
}


// Horner's method to compute g(x) = f(x+1) (equivalent to Shaw & Traub's method for a=1)
void horner_shift_CPU(sfixn *Polynomial, sfixn *Polynomial_shift, sfixn n, sfixn p)
{
  sfixn i;
  sfixn *temp;
  temp = (sfixn*) calloc (n, sizeof(sfixn));

  Polynomial_shift[0] = Polynomial[n-1];

  for (i=1; i<n; i++)
  {
    memcpy(temp+1, Polynomial_shift, i*sizeof(sfixn));
    add_arrays(Polynomial_shift, Polynomial_shift, temp, n, p);
    Polynomial_shift[0] = (Polynomial_shift[0] + Polynomial[n-1-i]) % p;
  }

  free(temp);
}
