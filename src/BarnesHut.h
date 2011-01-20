/*
 * burnesHut.h
 *
 *  Created on: 13/gen/2011
 *      Author: Claudio e Pino
 */

#ifndef BarnesHut_h_
#define BarnesHut_h_

int nbodies = 1000;
long double dt = 0.01;
int steps = 25;
long double eps = 0.0;
long double tol = 0.5;

long double dthf = 0.0;
long double epssq = 0.0;
long double itolsq = 0.0;

typedef struct leaf_t{
	long double vel[3];
	long double acc[3];
} leaf_t;

//typedef struct internal_node_t{
//	struct node_t *child0;
//	struct node_t *child1;
//	struct node_t *child2;
//	struct node_t *child3;
//	struct node_t *child4;
//	struct node_t *child5;
//	struct node_t *child6;
//	struct node_t *child7;
//} internal_node_t;

typedef union cell{
	struct leaf_t leaf;
	struct node_t* childs[8];
} cell;

// type definisce il tipo espresso dalla union: 0 = body_t, 1 = internal_node_t
typedef struct node_t{
	int type;
	long double mass;
	long double pos[3];
	cell cell;
} node_t;

#endif /* BURNESHUT_H_ */

