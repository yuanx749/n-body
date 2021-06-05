#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <omp.h>
#include "graphics/graphics.h"

const double e0 = 1e-3;

typedef struct particle {
    double x;
    double y;
    double ux;
    double uy;
    double ax;
    double ay;
    double m;
} particle_t;

typedef struct node {
    double m; // mass of box
    double x; // center of mass
    double y; // center of mass
    double w; // width of box
    double cx; // center of box
    double cy; // center of box
    particle_t *p;
    struct node *q1;
    struct node *q2;
    struct node *q3;
    struct node *q4;
} node_t;

static double get_wall_seconds() {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    double seconds = tv.tv_sec + (double)tv.tv_usec / 1000000;
    return seconds;
}

node_t * create_node(double w, double cx, double cy) {
    node_t *node = (node_t *)malloc(sizeof(node_t));
    if (node == NULL) {
        printf("Out of memory.\n");  
        exit(1);
    }
    node->m = 0;
    node->x = 0;
    node->y = 0;
    node->w = w;
    node->cx = cx;
    node->cy = cy;
    node->p = NULL;
    node->q1 = NULL;
    node->q2 = NULL;
    node->q3 = NULL;
    node->q4 = NULL;
    return node;
}

void insert_quadtree(node_t *node, particle_t *p) {
    if (p->x < 0 || p->x > 1 || p->y < 0 || p->y > 1) {
        printf("Particle moves outside\n");
        exit(1);
    }
    if (node->m > 0) {
        // box already has particle
        if (node->p != NULL && p->x == node->p->x && p->y == node->p->y) {
            printf("Two particles are at the same position\n");
            exit(1);
        }
        node->m += p->m;
        double w2 = node->w / 4;
        double w = w2 * 2; // partition width
        // find quadrant for new particle
        if (p->x >= node->cx && p->y >= node->cy) {
            if (node->q1 == NULL) {
                node->q1 = create_node(w, node->cx+w2, node->cy+w2);
            }
            insert_quadtree(node->q1, p);
        } else if (p->x < node->cx && p->y >= node->cy) {
            if (node->q2 == NULL) {
                node->q2 = create_node(w, node->cx-w2, node->cy+w2);
            }
            insert_quadtree(node->q2, p);
        } else if (p->x < node->cx && p->y < node->cy) {
            if (node->q3 == NULL) {
                node->q3 = create_node(w, node->cx-w2, node->cy-w2);
            }
            insert_quadtree(node->q3, p);
        } else if (p->x >= node->cx && p->y < node->cy) {
            if (node->q4 == NULL) {
                node->q4 = create_node(w, node->cx+w2, node->cy-w2);
            }
            insert_quadtree(node->q4, p);
        }
        // move existing particle to new leaf
        if (node->p != NULL) {
            if (node->p->x >= node->cx && node->p->y >= node->cy) {
                if (node->q1 == NULL) {
                    node->q1 = create_node(w, node->cx+w2, node->cy+w2);
                }
                insert_quadtree(node->q1, node->p);
            } else if (node->p->x < node->cx && node->p->y >= node->cy) {
                if (node->q2 == NULL) {
                    node->q2 = create_node(w, node->cx-w2, node->cy+w2);
                }
                insert_quadtree(node->q2, node->p);
            } else if (node->p->x < node->cx && node->p->y < node->cy) {
                if (node->q3 == NULL) {
                    node->q3 = create_node(w, node->cx-w2, node->cy-w2);
                }
                insert_quadtree(node->q3, node->p);
            } else if (node->p->x >= node->cx && node->p->y < node->cy) {
                if (node->q4 == NULL) {
                    node->q4 = create_node(w, node->cx+w2, node->cy-w2);
                }
                insert_quadtree(node->q4, node->p);
            }
            node->p = NULL; // only store particle in leaf
        }
    } else {
        // new leaf
        node->m = p->m;
        node->p = p;
    }
}

void reset_quadtree(node_t *node) {
    if (node->q1) {
        reset_quadtree(node->q1);
    }
    if (node->q2) {
        reset_quadtree(node->q2);
    }
    if (node->q3) {
        reset_quadtree(node->q3);
    }
    if (node->q4) {
        reset_quadtree(node->q4);
    }
    node->m = 0;
    node->p = NULL;
}

int prune_quadtree(node_t *node) {
    if (node == NULL) {
        return 0;
    }
    if (prune_quadtree(node->q1)) {
        node->q1 = NULL;
    }
    if (prune_quadtree(node->q2)) {
        node->q2 = NULL;
    }
    if (prune_quadtree(node->q3)) {
        node->q3 = NULL;
    }
    if (prune_quadtree(node->q4)) {
        node->q4 = NULL;
    }
    if (node->m == 0) {
        free(node);
        return 1;
    } else {
        return 0;
    }
}

void compute_center(node_t *node) {
    if (node == NULL) {
        return;
    }
    if (node->p != NULL) {
        // leaf
        node->x = node->p->x;
        node->y = node->p->y;
    } else {
        double sum_mx = 0, sum_my = 0;
        compute_center(node->q1);
        compute_center(node->q2);
        compute_center(node->q3);
        compute_center(node->q4);
        if (node->q1) {
            sum_mx += node->q1->m * node->q1->x;
            sum_my += node->q1->m * node->q1->y;
        }
        if (node->q2) {
            sum_mx += node->q2->m * node->q2->x;
            sum_my += node->q2->m * node->q2->y;
        }
        if (node->q3) {
            sum_mx += node->q3->m * node->q3->x;
            sum_my += node->q3->m * node->q3->y;
        }
        if (node->q4) {
            sum_mx += node->q4->m * node->q4->x;
            sum_my += node->q4->m * node->q4->y;
        }
        node->x = sum_mx / node->m;
        node->y = sum_my / node->m;
    }
}

void compute_acceleration(node_t *node, particle_t *p, double theta_max, double G) {
    if (node == NULL) {
        return;
    }
    if (node->p == NULL) {
        // internal node
        double rcx = p->x - node->cx;
        double rcy = p->y - node->cy;
        double distance = sqrt(rcx*rcx+rcy*rcy);
        if (node->w <= theta_max * distance) {
            // Newton's law of gravitation
            double rx = p->x - node->x;
            double ry = p->y - node->y;
            double r = 1.0 / (sqrt(rx*rx+ry*ry) + e0);
            double r_cubic = r*r*r;
            // calculate acceleration
            p->ax += -G * node->m * rx * r_cubic;  
            p->ay += -G * node->m * ry * r_cubic;
        } else {
            compute_acceleration(node->q1, p, theta_max, G);
            compute_acceleration(node->q2, p, theta_max, G);
            compute_acceleration(node->q3, p, theta_max, G);
            compute_acceleration(node->q4, p, theta_max, G);
        }
    } else {
        // leaf
        // Newton's law of gravitation
        double rx = p->x - node->x;
        double ry = p->y - node->y;
        double r = 1.0 / (sqrt(rx*rx+ry*ry) + e0);
        double r_cubic = r*r*r;
        // calculate acceleration
        p->ax += -G * node->m * rx * r_cubic;  
        p->ay += -G * node->m * ry * r_cubic;
    }
}

void move_particle(double delta_t, particle_t *p) {
    // calculate velocity
    p->ux += delta_t*p->ax;
    p->uy += delta_t*p->ay;
    // calculate position
    p->x += delta_t*p->ux;
    p->y += delta_t*p->uy;
    // reset acceleration
    p->ax = 0;
    p->ay = 0;
}

int main(int argc, char* argv[]) {
    // read input
    if (argc != 8) {
        printf("Give 7 input args: ./nbody N filename nsteps delta_t theta_max graphics n_threads\n");
        return 1;
    }
    const int N = atoi(argv[1]);
    const int nsteps = atoi(argv[3]);
    const double delta_t = atof(argv[4]);
    const double theta_max = atof(argv[5]);
    const int graphics = atoi(argv[6]);
    const int n_threads = atoi(argv[7]);
    const int file_size = 6*N*sizeof(double);
    const double G = 100.0/N;

    // read file
    double *buffer;
    buffer = (double *)malloc(file_size);
    FILE *inFile;
    long lSize;
    size_t result;
    inFile = fopen(argv[2], "rb");
    if (inFile == NULL) {fputs("File error\n", stderr); return 1;}
    fseek(inFile, 0, SEEK_END);
    lSize = ftell(inFile);
    if (lSize != file_size) {fputs("File size does not match N\n", stderr); return 1;}
    rewind(inFile);
    if (buffer == NULL) {fputs("Memory error\n", stderr); return 1;}
    result = fread(buffer, 1, lSize, inFile);
    if (result != lSize) {fputs("Reading error\n", stderr); return 1;}
    fclose(inFile);

    // initialize particles
    particle_t *particles = (particle_t *)malloc(N*sizeof(particle_t));
    for (int i = 0; i < N; i++) {
        particles[i].x  = buffer[6*i];
        particles[i].y  = buffer[6*i+1];
        particles[i].ux = buffer[6*i+3];
        particles[i].uy = buffer[6*i+4];
        particles[i].ax = 0;
        particles[i].ay = 0;
        particles[i].m = buffer[6*i+2];
    }
    
    // iterate
    double start_time = 0;
    int i;
    // initialize root box
    if (graphics) {
        InitializeGraphics(argv[0], 800, 800);
        SetCAxes(0, 1);
    }

    node_t *root = create_node(1.0, 0.5, 0.5);

    #pragma omp parallel private(i) num_threads(n_threads)
    {
    for (int step = 0; step < nsteps; step++) {
        #pragma omp single
        {
        if (graphics) {
            ClearScreen();
            for (i = 0; i < N; i++) {
                DrawCircle(particles[i].x, particles[i].y, 1, 1, 0.0025, 0);
            }
            Refresh();
            usleep(1000);
        }

        if (step == 0) {
            start_time = get_wall_seconds();
        }
        for (i = 0; i < N; i++) {
            insert_quadtree(root, &particles[i]);
        }
        // free memory of empty nodes
        prune_quadtree(root);
        compute_center(root);
        }
        #pragma omp for schedule(static)
        for (i = 0; i < N; i++) {
            compute_acceleration(root, &particles[i], theta_max, G);
        }
        #pragma omp for schedule(static)
        for (i = 0; i < N; i++) {
            move_particle(delta_t, &particles[i]);
        }
        // reset all the nodes to empty
        #pragma omp single
        reset_quadtree(root);
    }
    }

    // free quadtree
    prune_quadtree(root);
    printf("%d %.3f\n", N, get_wall_seconds() - start_time);

    if (graphics) {
        FlushDisplay();
        CloseDisplay();
    }

    // write buffer
    for (int i = 0; i < N; i++) {
        buffer[6*i]   = particles[i].x;
        buffer[6*i+1] = particles[i].y;
        buffer[6*i+3] = particles[i].ux;
        buffer[6*i+4] = particles[i].uy;
    }
    free(particles);

    // write file
    FILE *outFile;
    outFile = fopen("result.gal", "wb");
    fwrite(buffer, sizeof(double), 6*N, outFile);
    fclose(outFile);
    free(buffer);

    return 0;
}
