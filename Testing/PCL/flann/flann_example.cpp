// Copyright (c) 2012-2019 University of Lyon and CNRS (France).
// All rights reserved.
//
// This file is part of MEPP2; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of
// the License, or (at your option) any later version.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
#include <stdio.h>
#include <stdlib.h>

#if defined _MSC_VER
#pragma warning(disable : 4267 4244) // MT
#endif

#include <flann/flann.h>

float *
read_points(const char *filename, int rows, int cols)
{
  float *data;
  float *p;
  FILE *fin;
  int i, j;

  fin = fopen(filename, "r");
  if(!fin)
  {
    printf("Cannot open input file.\n");
    exit(1);
  }

  data = (float *)malloc(rows * cols * sizeof(float));
  if(!data)
  {
    printf("Cannot allocate memory.\n");
    exit(1);
  }
  p = data;

  int r = 0;
  for(i = 0; i < rows; ++i)
  {
    for(j = 0; j < cols; ++j)
    {
      r = fscanf(fin, "%g ", p);
      p++;
    }
  }

  fclose(fin);

  return data;
}

void
write_results(const char *filename, int *data, int rows, int cols)
{
  FILE *fout;
  int *p;
  int i, j;

  fout = fopen(filename, "w");
  if(!fout)
  {
    printf("Cannot open output file.\n");
    exit(1);
  }

  p = data;
  for(i = 0; i < rows; ++i)
  {
    for(j = 0; j < cols; ++j)
    {
      fprintf(fout, "%d ", *p);
      p++;
    }
    fprintf(fout, "\n");
  }
  fclose(fout);
}


int
main(int argc, char **argv)
{
  if(argc < 3)
  {
    printf("Usage: ./flann_example dataset.dat testset.dat\n");
    exit(EXIT_FAILURE);
  }

  float *dataset;
  float *testset;
  int nn;
  int *result;
  float *dists;
  struct FLANNParameters p;
  float speedup;
  flann_index_t index_id;

  int rows = 9000;
  int cols = 128;
  int tcount = 1000;

  /*
   * The files dataset.dat and testset.dat can be downloaded from:
   * http://people.cs.ubc.ca/~mariusm/uploads/FLANN/datasets/dataset.dat
   * http://people.cs.ubc.ca/~mariusm/uploads/FLANN/datasets/testset.dat
   */
  printf("Reading input data file.\n");
  dataset = read_points(argv[1], rows, cols);
  printf("Reading test data file.\n");
  testset = read_points(argv[2], tcount, cols);

  nn = 3;
  result = (int *)malloc(tcount * nn * sizeof(int));
  dists = (float *)malloc(tcount * nn * sizeof(float));

  p = DEFAULT_FLANN_PARAMETERS;
  p.algorithm = FLANN_INDEX_KDTREE;
  p.trees = 8;
  p.log_level = FLANN_LOG_INFO;
  p.checks = 64;

  printf("Computing index.\n");
  index_id = flann_build_index(dataset, rows, cols, &speedup, &p);
  flann_find_nearest_neighbors_index(
      index_id, testset, tcount, result, dists, nn, &p);

  write_results("results.dat", result, tcount, nn);

  flann_free_index(index_id, &p);
  free(dataset);
  free(testset);
  free(result);
  free(dists);

  return 0;
}
