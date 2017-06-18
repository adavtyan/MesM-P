/* This script converts nomral or in-plane vector file to a trajectory, 
 * where a particle pair indicates the direction of each vector.
 * It requires an vector file as an input obtained using one of the folling 
 * two custum dump commands, for normal and in-plane vectors, respectively.
 *
 * dump  nd mem custom 1000 dump.normal_vecs id type c_vecs[1] c_vecs[2] c_vecs[3] c_vecs[4] c_vecs[5] c_vecs[6]
 * dump  pd mem custom 1000 dump.inplane_vecs id type c_vecs[7] c_vecs[8] c_vecs[9] c_vecs[10] c_vecs[11] c_vecs[12]
 *
 * Additionally, sorting the enteries using the following dump_modify command is recommended.
 * dump_modify  [ID] sort id
 *
 * The output is a lammpstrj file specified as the second argumnet.
 *
 * Author: Aram Davtyan
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

typedef struct Data {
        int id;
        int ty;
	double x[3];
} Data;

typedef struct Frame {
	Frame() : data(NULL) {}
	int step;
	double time;
	double box[3][2];
	Data *data;
} Frame;

int main(int argc, char* argv[])
{
	if (argc<3 || argc>6) {
		printf("\n%s lammps_dump_traj output.dat [offset] [freq] [tstep]\n\n", argv[0]);
		return -1;
	}

	char *inp_name, *out_name;
	int freq = 1;
	int offset  = 0;
	double tstep = 0.001;
	inp_name = new char[strlen(argv[1])+1];
	out_name = new char[strlen(argv[2])+1];
	if (argc>=5) offset = atoi(argv[3]);
	if (argc>=6) freq = atoi(argv[4]);
	if (argc>=7) tstep = atof(argv[5]);

	strcpy(inp_name, argv[1]);
	strcpy(out_name, argv[2]);

	FILE * input = fopen(inp_name, "r");
	if (!input) {
		printf("\nERROR: Input file not found!\n\n");
		return -1;
	}

        const int block_size = 1000;
	int cur_size = block_size;
	Frame *frames = (Frame*)calloc(cur_size, sizeof(Frame));

	int i, j, k, atom_id, atom_type, nstr, j_st, aid;
	double x[6];
	double *x1, *x2;
        int i_ind_col=-1, i_tycol=-1, i_xcol=-1;
        int i_vecs_col[6] = {-1, -1, -1, -1, -1, -1};
	int n_vecs_col=0;
	bool bcols = false;
	int n_frame = 0;
        int i_frame = 0;
        int step = 0;
        int n_data = 0;
	int i_st = 0;
	double A[3][2], prd[3]; // Box boundaries and sizes
	char line[1024], buff[1024], item[128], *str[30], *str2;
	while ( fgets(line, sizeof(line), input) != NULL ) {
		if (strncmp(line, "ITEM:", 5)==0) {
			strcpy(item, line+6);
			
			if (strncmp(item, "TIMESTEP", 8)==0) {
				fgets(buff, sizeof(buff), input);
				step = atoi(buff);
			} else if (strncmp(item, "BOX BOUNDS", 10)==0) {
				for (i=0; i<3; i++) {
					fgets(buff, sizeof(buff), input);

					str[0] = strtok(buff," \t\n");
					str[1] = strtok(NULL," \t\n");

					A[i][0] = atof(str[0]);
					A[i][1] = atof(str[1]);
				}
			} else if (strncmp(item, "NUMBER OF ATOMS", 15)==0) {
				fgets(buff, sizeof(buff), input);
				n_data = atoi(buff);
			} else if (strncmp(item, "ATOMS", 5)==0) {
				if (step>=offset && i_frame%freq==0) {
					if (n_frame%1000==0) printf("Processing frame # %d\n", step);
					if (!bcols) {
						nstr = 0;
						n_vecs_col = 0;
						str2 = strtok(line+12," \t\n");
						while ( str2!=NULL ) {
							if (strcmp(str2, "id")==0) {
								i_ind_col = nstr;
							} else if (strcmp(str2, "type")==0) {
								i_tycol = nstr;
							} else if (strcmp(str2, "xcm")==0) {
								i_xcol = nstr;
							} else if (strncmp(str2, "c_vecs", 6)==0) {
								if (n_vecs_col>=6) {
									printf("Error: wrong input data format!\n");
									exit(0);
								}
								i_vecs_col[n_vecs_col] = nstr;
								n_vecs_col++;
							}

							nstr++;
							str2 = strtok(NULL," \t\n");
						}

						if (i_ind_col<0 || i_tycol<0 || n_vecs_col!=6 ) {
							printf("ID, type, or vecs colume not found!\n");
							exit(0);
						}

						bcols = true;
					}

					if (n_frame>=cur_size) {
						cur_size += block_size;
						frames = (Frame*)realloc(frames, cur_size*sizeof(Frame));
					}
					frames[n_frame].data = (Data*)malloc(2*n_data*sizeof(Data));

					frames[n_frame].step = step;
					frames[n_frame].time = tstep*step;

					for (i=0;i<3;i++) {
						frames[n_frame].box[i][0] = A[i][0];
						frames[n_frame].box[i][1] = A[i][1];
					}

					// Read coordinates
					for (i=0;i<n_data;++i) {
						fgets(buff, sizeof(buff), input);
						
						// Split the line string and read the first 5 values 
						nstr = 0;
						str[nstr] = strtok(buff," \t\n");
						while ( str[nstr]!=NULL ) {
							nstr++;
							if (nstr>=30) break;
							str[nstr] = strtok(NULL," \t\n");
						}

						atom_id = atoi(str[i_ind_col]);
						atom_type = atoi(str[i_tycol]);
						x[0] =  atof(str[i_vecs_col[0]]);
						x[1] =  atof(str[i_vecs_col[1]]);
						x[2] =  atof(str[i_vecs_col[2]]);
						x[3] =  atof(str[i_vecs_col[3]]);
						x[4] =  atof(str[i_vecs_col[4]]);
						x[5] =  atof(str[i_vecs_col[5]]);

						j = 2*i;
						aid = 2*(atom_id-1)+1;

						frames[n_frame].data[j].id = aid;
						frames[n_frame].data[j].ty = atom_type;
						frames[n_frame].data[j].x[0] = x[0];
						frames[n_frame].data[j].x[1] = x[1];
						frames[n_frame].data[j].x[2] = x[2];

						j++;
						aid++;

						frames[n_frame].data[j].id = aid;
						frames[n_frame].data[j].ty = atom_type+1;
						frames[n_frame].data[j].x[0] = x[3];
						frames[n_frame].data[j].x[1] = x[4];
						frames[n_frame].data[j].x[2] = x[5];
					}

					n_frame++;
				}
				if (step>=offset) i_frame++;
			}
		}
	}
	fclose(input);

	printf("Number of frames %d\n", n_frame);
	printf("Number of atoms %d\n\n", n_data);

	if (n_frame<=0) {
		printf("No frame was processed!\n");
		exit(0);
	}

	FILE *out = fopen(out_name, "w");
	for (i=0;i<n_frame;i++) {
		fprintf(out, "ITEM: TIMESTEP\n");
		fprintf(out, "%d\n", frames[i].step);
		fprintf(out, "ITEM: NUMBER OF ATOMS\n");
		fprintf(out, "%d\n", 2*n_data);
		fprintf(out, "ITEM: BOX BOUNDS pp pp pp\n");
		fprintf(out, "%f %f\n", frames[i].box[0][0], frames[i].box[0][1]);
		fprintf(out, "%f %f\n", frames[i].box[1][0], frames[i].box[1][1]);
		fprintf(out, "%f %f\n", frames[i].box[2][0], frames[i].box[2][1]);
		fprintf(out, "ITEM: ATOMS id type x y z\n");
		for (j=0;j<2*n_data;j++) {
			Data &a_d = frames[i].data[j];
			fprintf(out, "%d %d %f %f %f\n", a_d.id, a_d.ty, a_d.x[0], a_d.x[1], a_d.x[2]);
		}
	}
	fclose(out);

	for (i=0;i<n_frame;++i) {
		if (frames[i].data!=NULL) free(frames[i].data);
	}
	free(frames);

	return 0;
}
