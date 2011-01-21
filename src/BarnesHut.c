/*
 * burnesHut.c
 *
 *  Created on: 13/gen/2011
 *      Author: Claudio e Pino
 */
#include <stdio.h>
#include <stdlib.h>
#include "BarnesHut.h"
#include <string.h>
#include <math.h>
#include "genlib.h"
#include <papi.h>

//node_t* null_childs[8] = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL };

node_t** bodies;
double diameter;
double seed = 1234567890;
double center[3];
int curr = 0;
node_t *root = NULL;

//double random_generator(double min, double max);
void create_bodies();
void compute_center_and_diameter();
void insert(node_t* sub_root, node_t* node, double r);
node_t* new_node(double mass, double pos[3], double acc[3], double vel[3],
		int type);
void compute_center_of_mass(node_t* node);
void compute_force(node_t* root, node_t* body, double diameter, int where);
void recourse_force(node_t* root, node_t* body, double dsq);
void advance(node_t* body);
void deallocate_tree(node_t* node);

int main(int argc, char* argv[]) {

	FILE *inputf;
	FILE *outputf;
	long_long clockStart, clockEnd;
	long_long values[2];
	int event_set = PAPI_NULL;
	int events[2] = { PAPI_L2_TCM, PAPI_TOT_CYC };

	char* inputfile = argv[1];
	inputf = fopen(inputfile, "r");

	if (inputf == NULL) {
		printf("impossibile leggere da file");
		exit(1);
	}

	fscanf(inputf, "%d", &nbodies);
	fscanf(inputf, "%d", &steps);
	fscanf(inputf, "%lf", &dt);
	fscanf(inputf, "%lf", &eps);
	fscanf(inputf, "%lf", &tol);

	fclose(inputf);

	bodies = malloc(nbodies * sizeof(node_t*));
	create_bodies();

	dthf = 0.5 * dt;
	epssq = eps * eps;
	itolsq = 1.0 / (tol * tol);

	if (PAPI_library_init(PAPI_VER_CURRENT) != PAPI_VER_CURRENT) {
		printf("Inizializzazione della libreria di papi fallita \n");
		exit(1);
	}

	if (PAPI_create_eventset(&event_set) != PAPI_OK) {
		printf("E' andata a male la creazione dell'eventSet \n");
		exit(1);
	}

	if (PAPI_add_events(event_set, events, 2) != PAPI_OK) {
		printf("E' andata a male l'aggiunta degli eventi\n");
		exit(1);
	}

	int step = 0;
	root = NULL;
	PAPI_start(event_set);
	clockStart = PAPI_get_real_usec();
	for (step = 0; step < steps; step++) {
		compute_center_and_diameter();

		root = malloc(sizeof(struct node_t)); // "new" is like "malloc"

		root->type = 1;
		root->mass = 0.0;
		root->pos[0] = center[0];
		root->pos[1] = center[1];
		root->pos[2] = center[2];
		root->cell.childs[0] = NULL;
		root->cell.childs[1] = NULL;
		root->cell.childs[2] = NULL;
		root->cell.childs[3] = NULL;
		root->cell.childs[4] = NULL;
		root->cell.childs[5] = NULL;
		root->cell.childs[6] = NULL;
		root->cell.childs[7] = NULL;

		double radius = diameter * 0.5;

		int i = 0;
		for (i = 0; i < nbodies; i++) {
			insert(&(*root), &(*bodies[i]), radius); // questo è il modo per passare i dati per riferimento... cioè mandare l'indirizzo della struttura puntata dal puntatore
		}
		curr = 0;
		compute_center_of_mass(&(*root));

		for (i = 0; i < nbodies; i++) {
			compute_force(&(*root), &(*bodies[i]), diameter, step);
			advance(&(*bodies[i]));
		}
		//		for (i = 0; i < nbodies; i++) {
		//		}

		deallocate_tree(root);

		//		int p = 0;
		//		for (p = 0; p < nbodies; p++)
		//			printf("%lf, %lf, %lf \n", bodies[p]->pos[0], bodies[p]->pos[1],
		//					bodies[p]->pos[2]);
		//		printf("*************************************** \n");
	}
	clockEnd = PAPI_get_real_usec();
	PAPI_stop(event, values);
	int i = 0;
	outputf = fopen("output", "w");
	fprintf(outputf, "Tempo di esecuzione: %lld \n", clockEnd - clockStart);
	fprintf(outputf, "Cache miss L2: %lld \n", values[0]);
	fprintf(outputf, "# cicli: %lld \n", values[1]);
	for (i = 0; i < nbodies; i++) {
		fprintf(outputf, "%lf, %lf, %lf \n", bodies[i]->pos[0],
				bodies[i]->pos[1], bodies[i]->pos[2]);
	}

	fflush(outputf);
	fclose(outputf);
	printf("Esecuzione completata\n");

	return 0;
}

/**
 * Basic uniform random generator: Minimal Standard in Park and
 Miller (1988): "Random Number Generators: Good Ones Are Hard to
 Find", Comm. of the ACM, 31, 1192-1201.
 Parameters: m = 2^31-1, a=48271.
 */

//double random_generator(double min, double max) {
//	int m = INT_MAX;
//	int a = 48271;
//	double q = m / a;
//	double r = m % a;
//
//	double k = seed / q;
//	seed = a * (seed - k * q) - r * k;
//	if (seed < 1)
//		seed += m;
//	r = seed / m;
//	return r * (max - min) + min;
//
//
//
//}

void create_bodies() {
	FILE* inputdataf = fopen("inputData", "r");
	int i = 0;
	char line[100];
	char* token;
	char* end;
	double d = 0.0;
	//	const gsl_rng_type *T;
	//	gsl_rng *r;
	//	gsl_rng_env_setup();
	//
	//	T = gsl_rng_default;
	//	r = gsl_rng_alloc(T);

	for (i = 0; i < nbodies; i++) {

		//		printf("i = %d \n", i);
		bodies[i] = malloc(sizeof(node_t));

		//		double mass = gsl_rng_uniform(r);
		//		pos[0] = gsl_rng_uniform(r);
		//		pos[1] = gsl_rng_uniform(r);
		//		pos[2] = gsl_rng_uniform(r);
		//		vel[0] = gsl_rng_uniform(r);
		//		vel[1] = gsl_rng_uniform(r);
		//		vel[2] = gsl_rng_uniform(r);


		fgets(line, 100, inputdataf);

		//		printf("%s\n", line);

		token = strtok(line, " ");

		d = strtold(token, &end);
		bodies[i]->mass = d;
		//		printf("%lf\n",d);
		token = strtok(NULL, " ");
		d = strtold(token, &end);
		bodies[i]->pos[0] = d;
		//		printf("%lf\n",d);
		token = strtok(NULL, " ");
		d = strtold(token, &end);
		bodies[i]->pos[1] = d;
		//		printf("%lf\n",d);
		token = strtok(NULL, " ");
		d = strtold(token, &end);
		bodies[i]->pos[2] = d;
		//		printf("%lf\n",d);
		token = strtok(NULL, " ");
		d = strtold(token, &end);
		bodies[i]->cell.leaf.vel[0] = d;
		//		printf("%lf\n",d);
		token = strtok(NULL, " ");
		d = strtold(token, &end);
		bodies[i]->cell.leaf.vel[1] = d;
		//		printf("%lf\n",d);
		token = strtok(NULL, " ");
		d = strtold(token, &end);
		bodies[i]->cell.leaf.vel[2] = d;
		//		printf("%lf\n",d);


		bodies[i]->cell.leaf.acc[0] = 0.0;
		bodies[i]->cell.leaf.acc[1] = 0.0;
		bodies[i]->cell.leaf.acc[2] = 0.0;

	}

	fclose(inputdataf);
}
void compute_center_and_diameter() {
	double min[3] = { 1.0e90, 1.0e90, 1.0e90 }, max[3] = { -1.0e90, -1.0e90,
			-1.0e90 }, pos[3];

	int i = 0;
	for (i = 0; i < nbodies; i++) {
		pos[0] = bodies[i]->pos[0];
		pos[1] = bodies[i]->pos[1];
		pos[2] = bodies[i]->pos[2];

		if (min[0] > pos[0])
			min[0] = pos[0];

		if (min[1] > pos[1])
			min[1] = pos[1];

		if (min[2] > pos[2])
			min[2] = pos[2];

		if (max[0] < pos[0])
			max[0] = pos[0];

		if (max[1] < pos[1])
			max[1] = pos[1];

		if (max[2] < pos[2])
			max[2] = pos[2];

	}
	diameter = max[0] - min[0];
	if (diameter < (max[1] - min[1]))
		diameter = (max[1] - min[1]);

	if (diameter < (max[2] - min[2]))
		diameter = (max[2] - min[2]);

	center[0] = (max[0] + min[0]) * 0.5;
	center[1] = (max[1] + min[1]) * 0.5;
	center[2] = (max[2] + min[2]) * 0.5;
}

void insert(node_t* sub_root, node_t* node, double r) {

	bool finished = FALSE;

	do {
		double x = 0.0, y = 0.0, z = 0.0;
		int i = 0;
		if (sub_root->pos[0] < node->pos[0]) {
			i = 1;
			x = r;
		}

		if (sub_root->pos[1] < node->pos[1]) {
			i += 2;
			y = r;
		}

		if (sub_root->pos[2] < node->pos[2]) {
			i += 4;
			z = r;
		}

		if (sub_root->cell.childs[i] == NULL) {
			sub_root->cell.childs[i] = node;
			finished = TRUE;
		} else if (sub_root->cell.childs[i]->type == 1) {
			//				insert(&(*sub_root->cell.childs[i]), &(*node), 0.5 * r);
			r *= 0.5;
			sub_root = sub_root->cell.childs[i];
		} else {
			double rh = 0.5 * r;
			double position[3] = { sub_root->pos[0] - rh + x, sub_root->pos[1]
					- rh + y, sub_root->pos[2] - rh + z };
			node_t *cell = new_node(0.0, position, NULL, NULL, 1);
			//				insert(&(*cell), &(*node), rh);
			int k = 0;
			if (cell->pos[0] < node->pos[0]) {
				k = 1;
			}

			if (cell->pos[1] < node->pos[1]) {
				k += 2;
			}

			if (cell->pos[2] < node->pos[2]) {
				k += 4;
			}

			cell->cell.childs[k] = node;
			node_t* tmp = &(*sub_root->cell.childs[i]);
			sub_root->cell.childs[i] = cell;
			sub_root = cell;
			node = tmp;
			r = rh;
			//			insert(&(*cell), &(*sub_root->cell.childs[i]), rh);
		}
	} while (!finished);
}

node_t* new_node(double mass, double pos[3], double acc[3], double vel[3],
		int type) {
	node_t* node = malloc(sizeof(node_t)); // "new" is like "malloc"
	node->type = type;
	node->mass = mass;
	node->pos[0] = pos[0];
	node->pos[1] = pos[1];
	node->pos[2] = pos[2];
	if (type == 0) {
		//		node->cell.leaf->acc = acc;
		node->cell.leaf.acc[0] = acc[0];
		node->cell.leaf.acc[1] = acc[1];
		node->cell.leaf.acc[2] = acc[2];
		node->cell.leaf.vel[0] = vel[0];
		node->cell.leaf.vel[1] = vel[1];
		node->cell.leaf.vel[2] = vel[2];
	} else if (type == 1) {
		node->cell.childs[0] = NULL;
		node->cell.childs[1] = NULL;
		node->cell.childs[2] = NULL;
		node->cell.childs[3] = NULL;
		node->cell.childs[4] = NULL;
		node->cell.childs[5] = NULL;
		node->cell.childs[6] = NULL;
		node->cell.childs[7] = NULL;
	}

	return (node);
}
void compute_center_of_mass(node_t* node) {
	double m = 0, p[3] = { 0.0, 0.0, 0.0 };
	//		node_t* firstChild = node->cell.internal_node.child0;
	//		node_t** childs = &firstChild;

	node_t* ch;

	int j = 0;
	node->mass = 0.0;
	int i;
	for (i = 0; i < 8; i++) {
		ch = &(*node->cell.childs[i]);
		if (ch != NULL) {
			node->cell.childs[i] = NULL;
			node->cell.childs[j++] = ch;

			if (ch->type == 0) {
				bodies[curr++] = ch;
			} else {
				compute_center_of_mass(&(*ch));
			}

			m = ch->mass;
			node->mass += m;
			p[0] += ch->pos[0] * m;
			p[1] += ch->pos[1] * m;
			p[2] += ch->pos[2] * m;
		}
	}

	m = 1.0 / node->mass;
	node->pos[0] = p[0] * m;
	node->pos[1] = p[1] * m;
	node->pos[2] = p[2] * m;
}

void compute_force(node_t* root, node_t* body, double diameter, int where) {
	double a[3] = { body->cell.leaf.acc[0], body->cell.leaf.acc[1],
			body->cell.leaf.acc[2] };

	body->cell.leaf.acc[0] = 0.0;
	body->cell.leaf.acc[1] = 0.0;
	body->cell.leaf.acc[2] = 0.0;

	recourse_force(&(*root), &(*body), diameter * diameter * itolsq);

	if (where > 0) {

		body->cell.leaf.vel[0] += (body->cell.leaf.acc[0] - a[0]) * dthf;
		body->cell.leaf.vel[1] += (body->cell.leaf.acc[1] - a[1]) * dthf;
		body->cell.leaf.vel[2] += (body->cell.leaf.acc[2] - a[2]) * dthf;
	}

}

void recourse_force(node_t* root, node_t* body, double dsq) {

	double dr[3], drsq, nphi, scale, idr;

	dr[0] = root->pos[0] - body->pos[0];
	dr[1] = root->pos[1] - body->pos[1];
	dr[2] = root->pos[2] - body->pos[2];

	drsq = dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2];

	if (drsq < dsq) {
		if (root->type == 1) {
			dsq *= 0.25;
			if (root->cell.childs[0] != NULL) {
				recourse_force(&(*root->cell.childs[0]), &(*body), dsq);
				if (root->cell.childs[1] != NULL) {
					recourse_force(&(*root->cell.childs[1]), &(*body), dsq);
					if (root->cell.childs[2] != NULL) {
						recourse_force(&(*root->cell.childs[2]), &(*body), dsq);
						if (root->cell.childs[3] != NULL) {
							recourse_force(&(*root->cell.childs[3]), &(*body),
									dsq);
							if (root->cell.childs[4] != NULL) {
								recourse_force(&(*root->cell.childs[4]),
										&(*body), dsq);
								if (root->cell.childs[5] != NULL) {
									recourse_force(&(*root->cell.childs[5]),
											&(*body), dsq);
									if (root->cell.childs[6] != NULL) {
										recourse_force(
												&(*root->cell.childs[6]),
												&(*body), dsq);
										if (root->cell.childs[7] != NULL) {
											recourse_force(
													&(*root->cell.childs[7]),
													&(*body), dsq);
										}
									}
								}
							}
						}
					}
				}
			}
		} else {
			if (root != body) {
				drsq += epssq;
				idr = 1 / sqrt(drsq);
				nphi = root->mass * idr;
				scale = nphi * idr * idr;
				body->cell.leaf.acc[0] += (dr[0] * scale);
				body->cell.leaf.acc[1] += (dr[1] * scale);
				body->cell.leaf.acc[2] += (dr[2] * scale);
			}
		}
	} else { // il nodo è abbastanza distante, non andiamo più in profondità... non ne vale la pena!
		drsq += epssq;
		idr = 1 / sqrt(drsq);
		nphi = root->mass * idr;
		scale = nphi * idr * idr;
		body->cell.leaf.acc[0] += (dr[0] * scale);
		body->cell.leaf.acc[1] += (dr[1] * scale);
		body->cell.leaf.acc[2] += (dr[2] * scale);
	}
}

void advance(node_t* body) {

	double dvel[3], velh[3];

	dvel[0] = body->cell.leaf.acc[0] * dthf;
	dvel[1] = body->cell.leaf.acc[1] * dthf;
	dvel[2] = body->cell.leaf.acc[2] * dthf;

	velh[0] = body->cell.leaf.vel[0] + dvel[0];
	velh[1] = body->cell.leaf.vel[1] + dvel[1];
	velh[2] = body->cell.leaf.vel[2] + dvel[2];

	body->pos[0] += velh[0] * dt;
	body->pos[1] += velh[1] * dt;
	body->pos[2] += velh[2] * dt;

	body->cell.leaf.vel[0] = velh[0] + dvel[0];
	body->cell.leaf.vel[1] = velh[1] + dvel[1];
	body->cell.leaf.vel[2] = velh[2] + dvel[2];

}

void deallocate_tree(node_t* node) {

	if ((node == NULL) || (node->type == 0))
		return;
	else {
		deallocate_tree(node->cell.childs[0]);
		deallocate_tree(node->cell.childs[1]);
		deallocate_tree(node->cell.childs[2]);
		deallocate_tree(node->cell.childs[3]);
		deallocate_tree(node->cell.childs[4]);
		deallocate_tree(node->cell.childs[5]);
		deallocate_tree(node->cell.childs[6]);
		deallocate_tree(node->cell.childs[7]);
		free(node);
	}
}
