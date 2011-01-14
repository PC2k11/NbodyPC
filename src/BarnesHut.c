/*
 * burnesHut.c
 *
 *  Created on: 13/gen/2011
 *      Author: Claudio e Pino
 */
#include <stdio.h>
#include <stdlib.h>
#include "BarnesHut.h"
#include <limits.h>
#include <string.h>
#include <math.h>

node_t* null_childs[8] = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL };

node_t* bodies[];
double diameter;
double center[3];
int curr = 0;
node_t *root;

void create_bodies();
void compute_center_and_diameter();
void insert(node_t* sub_root, node_t* node, double r);
node_t* new_node(double mass, double pos[3], double acc[3], double vel[3],
		node_t* childs[8], int type);
void compute_center_of_mass(node_t* node);
void compute_force(node_t* root, node_t* body, double diameter, int where);
void recourse_force(node_t* root, node_t* body, double dsq);
void advance(node_t* body);

int main(int argc, char* argv[]) {

	FILE *inputf;
	FILE *outputf;

	char* inputfile = argv[1];
	inputf = fopen(inputfile, "r");

	if (inputf == NULL) {
		printf("impossibile leggere ad file");
		exit(1);
	}

	fscanf(inputf, "%d", &nbodies);
	fscanf(inputf, "%d", &steps);
	fscanf(inputf, "%lf", &dt);
	fscanf(inputf, "%lf", &eps);
	fscanf(inputf, "%lf", &tol);

	fclose(inputf);
	dthf = 0.5 * dt;
	epssq = eps * eps;
	itolsq = 1.0 / (tol * tol);

	*bodies = (node_t*) malloc(nbodies * sizeof(node_t));

	create_bodies();

	int step = 0;
	//TODO: far partire papi
	for (step = 0; step < steps; step++) {
		compute_center_and_diameter();

		//		root = new_node(0.0, center, NULL, NULL, null_childs, 1);
		root = (node_t*) malloc(sizeof(struct node_t*)); // "new" is like "malloc"

		root->type = 1;
		root->mass = 0.0;
		root->pos[0] = center[0];
		root->pos[1] = center[1];
		root->pos[2] = center[2];
		root->cell.internal_node.child0 = NULL;
		root->cell.internal_node.child1 = NULL;
		root->cell.internal_node.child2 = NULL;
		root->cell.internal_node.child3 = NULL;
		root->cell.internal_node.child4 = NULL;
		root->cell.internal_node.child5 = NULL;
		root->cell.internal_node.child6 = NULL;
		root->cell.internal_node.child7 = NULL;

		double radius = diameter * 0.5;
		int i = 0;
		for (i = 0; i < nbodies; i++) {
			insert(root, bodies[i], radius);
		}

		curr = 0;
		compute_center_of_mass(root);

		for (i = 0; i < nbodies; i++) {
			compute_force(root, bodies[i], diameter, step);
			advance(bodies[i]);
		}

		free(root);

	}
	int i = 0;
	outputf = fopen("output", "w");
	for (i = 0; i < nbodies; i++) {
		fprintf(outputf, "%lf, %lf, %lf \n", bodies[i]->pos[0],
				bodies[i]->pos[1], bodies[i]->pos[2]);
		fflush(outputf);
	}

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

double random_generator(int seed, double min, double max) {
	int m = INT_MAX;
	int a = 48271;
	double q = m / a;
	double r = m % a;

	double k = seed / q;
	seed = a * (seed - k * q) - r * k;
	if (seed < 1)
		seed += m;
	r = seed / m;
	return r * (max - min) + min;

}

void create_bodies() {
	int i = 0;
	for (i = 0; i < nbodies; i++) {
		double pos[3] = { 0.0, 0.0, 0.0 };
		double vel[3] = { 0.0, 0.0, 0.0 };
		double acc[3] = { 0.0, 0.0, 0.0 };
		double mass = random_generator(rand(), 0.0, drand48() + 0.1);
		pos[0] = random_generator(rand(), 0.0, drand48() + 0.1);
		pos[1] = random_generator(rand(), 0.0, drand48() + 0.1);
		pos[2] = random_generator(rand(), 0.0, drand48() + 0.1);
		vel[0] = random_generator(rand(), 0.0, drand48() + 0.1);
		vel[1] = random_generator(rand(), 0.0, drand48() + 0.1);
		vel[2] = random_generator(rand(), 0.0, drand48() + 0.1);

		node_t *body = new_node(mass, pos, acc, vel, null_childs, 0);
		bodies[i] = body;
	}
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

		diameter = max[0] - min[0];
		if (diameter < (max[1] - min[1]))
			diameter = (max[1] - min[1]);

		if (diameter < (max[2] - min[2]))
			diameter = (max[2] - min[2]);

		center[0] = (max[0] + min[0]) * 0.5;
		center[1] = (max[1] + min[2]) * 0.5;
		center[1] = (max[1] + min[2]) * 0.5;
	}
}

void insert(node_t* sub_root, node_t* node, double r) {

	int i = 0;
	double x = 0.0, y = 0.0, z = 0.0;

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

	//creiamo un array di nodi che puntano ai childs del nodo passato come param
	node_t *childs[8] = { sub_root->cell.internal_node.child0,
			sub_root->cell.internal_node.child1,
			sub_root->cell.internal_node.child2,
			sub_root->cell.internal_node.child3,
			sub_root->cell.internal_node.child4,
			sub_root->cell.internal_node.child5,
			sub_root->cell.internal_node.child6,
			sub_root->cell.internal_node.child7 };
	if (childs[i] == NULL) {
		childs[i] = node;
	} else if (node->type == 1) {
		insert(childs[i], node, 0.5 * r);
	} else {
		double rh = 0.5 * r;
		double position[3] = { sub_root->pos[0] - rh + x, sub_root->pos[1] - rh
				+ y, sub_root->pos[2] - rh + z };
		node_t *cell = new_node(0.0, position, NULL, NULL, null_childs, 1);
		insert(cell, node, rh);
		insert(cell, childs[i], rh);
		childs[i] = cell;
	}
}

node_t* new_node(double mass, double pos[3], double acc[3], double vel[3],
		node_t* childs[8], int type) {
	struct node_t* node = (node_t*) malloc(sizeof(struct node_t*)); // "new" is like "malloc"

	node->type = type;
	if (type == 0) {
		//		node->cell.leaf->acc = acc;
		node->cell.leaf.acc[0] = acc[0];
		node->cell.leaf.acc[1] = acc[1];
		node->cell.leaf.acc[2] = acc[2];
		node->cell.leaf.vel[0] = vel[0];
		node->cell.leaf.vel[1] = vel[1];
		node->cell.leaf.vel[2] = vel[2];
	} else if (type == 2) {
		node->mass = mass;
		node->pos[0] = pos[0];
		node->pos[1] = pos[1];
		node->pos[2] = pos[2];
		node->cell.internal_node.child0 = childs[0];
		node->cell.internal_node.child1 = childs[1];
		node->cell.internal_node.child2 = childs[2];
		node->cell.internal_node.child3 = childs[3];
		node->cell.internal_node.child4 = childs[4];
		node->cell.internal_node.child5 = childs[5];
		node->cell.internal_node.child6 = childs[6];
		node->cell.internal_node.child7 = childs[7];
	}

	return (node);
}
void compute_center_of_mass(node_t* node) {
	double m, p[3] = { 0.0, 0.0, 0.0 };
	node_t *childs[8] = { node->cell.internal_node.child0,
			node->cell.internal_node.child1, node->cell.internal_node.child2,
			node->cell.internal_node.child3, node->cell.internal_node.child4,
			node->cell.internal_node.child5, node->cell.internal_node.child6,
			node->cell.internal_node.child7 };
	node_t* ch;

	int j = 0;
	node->mass = 0.0;
	int i;
	for (i = 0; i < 8; i++) {
		ch = childs[i];
		if (ch != NULL) {
			childs[i] = NULL;
			childs[j++] = ch;

			if (ch->type == 0) {
				bodies[curr++] = ch;
			} else {
				compute_center_of_mass(ch);
			}

			m = ch->mass;
			node->mass += m;
			p[0] = ch->pos[0] * m;
			p[1] = ch->pos[1] * m;
			p[2] = ch->pos[2] * m;
		}
	}

	m = 1.0 / node->mass;
	node->pos[0] = p[0] * m;
	node->pos[1] = p[1] * m;
	node->pos[2] = p[2] * m;
}

void compute_force(node_t* root, node_t* body, double diameter, int where) {
	double a[3] = { body->cell.leaf.acc[0], body->cell.leaf.acc[0],
			body->cell.leaf.acc[0] };

	body->cell.leaf.acc[0] = 0.0;
	body->cell.leaf.acc[1] = 0.0;
	body->cell.leaf.acc[2] = 0.0;

	recourse_force(root, body, diameter * diameter * itolsq);

	if (where > 0) {

		body->cell.leaf.vel[0] += (body->cell.leaf.acc[0] - a[0]) * dthf;
		body->cell.leaf.vel[1] += (body->cell.leaf.acc[1] - a[1]) * dthf;
		body->cell.leaf.vel[2] += (body->cell.leaf.acc[2] - a[2]) * dthf;
	}

}

void recourse_force(node_t* root, node_t* body, double dsq) {

	double dr[3], drsq, nphi, scale, idr;

	dr[0] = root->pos[0] - body->pos[0];
	dr[0] = root->pos[1] - body->pos[1];
	dr[0] = root->pos[2] - body->pos[2];

	drsq = dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2];

	if (drsq < dsq) {
		if (root->type == 1) {
			dsq *= 0.25;
			if (root->cell.internal_node.child0 != NULL) {
				recourse_force(root->cell.internal_node.child0, body, dsq);
				if (root->cell.internal_node.child1 != NULL) {
					recourse_force(root->cell.internal_node.child1, body, dsq);
					if (root->cell.internal_node.child2 != NULL) {
						recourse_force(root->cell.internal_node.child2, body,
								dsq);
						if (root->cell.internal_node.child3 != NULL) {
							recourse_force(root->cell.internal_node.child3,
									body, dsq);
							if (root->cell.internal_node.child4 != NULL) {
								recourse_force(root->cell.internal_node.child4,
										body, dsq);
								if (root->cell.internal_node.child5 != NULL) {
									recourse_force(
											root->cell.internal_node.child5,
											body, dsq);
									if (root->cell.internal_node.child6 != NULL) {
										recourse_force(
												root->cell.internal_node.child6,
												body, dsq);
										if (root->cell.internal_node.child7
												!= NULL) {
											recourse_force(
													root->cell.internal_node.child7,
													body, dsq);
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

	body->pos[0] += velh[0] + dt;
	body->pos[1] += velh[1] + dt;
	body->pos[2] += velh[2] + dt;

	body->cell.leaf.vel[0] = velh[0] + dvel[0];
	body->cell.leaf.vel[1] = velh[1] + dvel[1];
	body->cell.leaf.vel[2] = velh[2] + dvel[2];

}
