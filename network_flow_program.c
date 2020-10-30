/**
 * NAME         :   Debajyoti Dasgupta
 * ROLL NO      :   18CS30051
 * COURSE       :   ALGORITHMS - II
 * TOPIC        :   NETWORK FLOW ALGORITHM
 * LANGUAGE     :   C
 * DATE         :   15 . 10 . 2020
 * ASSIGNMENT   :   1
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// constant used as infinity in the program
const int INF = 1e9+7;

//======================    DATA TYPE DEFINATIONS   =========================

/**
 * struct for storing the information of the edges
 * this data type will form the individual element
 * of the linked list of the edges
 */ 
struct EDGE{
    int y;              // endpoint of the edge(x,y), i.e. edge from x to y, x->y
    int c;              // stores the capacity of this edge
    int f;              // stores the flow value to be assigned to this edge
    struct EDGE *next;  // pointer to the next edge of the linked list
};
typedef struct EDGE EDGE;

/**
 * struct for storing the information about the vertices
 * of the graph. Each vertex wil contain a linked list of 
 * the outgoing edges for that vertex. 
 */
typedef struct {
    int x;      // represents the ID of the current vertex
    int n;      // stores the need value of this vertex
    EDGE *p;    // points to an edge node which is first node of the adjacency list
} VERTEX ;

/**
 * the following data structure will store the directed graph.
 * this is the core structure which will have a list of vertex
 * inside it and each vertex will further have a linked list 
 * for each outgoing edge. 
 */
typedef struct {
    int V;      //  stores the number of vertices in the graph
    int E;      //  stores the number of edges in the graph
    VERTEX *H;  //  pointer to an array of vertex nodes
} GRAPH ;

/**
 * The following data structure will hold the the pair of 
 * the current node and the pointer to the next node. This 
 * data structure will help in the implementation of the 
 * queue in the BFS for max flow
 */
struct pair{                                // structure required for implementing the queue
    int node;                               // current node to be processed
    struct pair *next;                      // pointer for making a linked list
};
typedef struct pair pair;

//======================    HELPER FUNCTIONS DEFINATION   ========================

// Find minimum of two numbers
int MIN(int a, int b){
    return a<=b ? a : b;
}

// Find maximum of two numbers
int MAX(int a, int b){
    return a>=b ? a : b;
}

/**
 * Helper function that reads the data from the file name 
 * mentioned in the parameter fname and createss a graph 
 * according to the input
 * 
 * Return 
 * ------
 * GRAPH * : a pointewr of the type GRAPH that will contain 
 *          the entire Directed Graph read in the input
 */ 
GRAPH *ReadGraph(char *fname){
    FILE *f = fopen(fname, "r");                                // open file with the standard input stream

    GRAPH *head = (GRAPH*)malloc(sizeof(GRAPH));
    
    // read the number of nodes and edges
    int n, m;
    fscanf(f, "%d%d", &n, &m);
    head->V = n;                                                // initialize number of vertices with the input n
    head->E = m;                                                // initialize number of edges with the input m
    head->H = (VERTEX*)malloc((head->V+1) * sizeof(VERTEX));    // create an array of (number of) verices

    // read in the need valuesfor each vertex
    for (int i = 1; i <= head->V; i++){
        head->H[i].x = i;                                       // initialize the id of the vertex
        head->H[i].p = NULL;                                    // initialize the linked list as NULL List
        fscanf(f, "%d", &head->H[i].n);                         // read the need value from the input stream 
    }

    int x, y, c;
    // read and insert a new edge for every edge read
    for (int i = 0; i < head->E; i++){
        fscanf(f, "%d%d%d", &x, &y, &c);
        EDGE *cur = (EDGE*)malloc(sizeof(EDGE));                // make a new EDGE object
        cur->y = y;                                             // set the end point of the edge
        cur->c = c;                                             // initialize the capacity
        cur->f = 0;                                             // initialize the flow
        if(head->H[x].p) cur->next = head->H[x].p;
        else cur->next = NULL;
        head->H[x].p = cur;                                     // Insert the new element in the front of the linked list
    }

    fclose(f);
    return head;
}

/**
 * Helper function to print the Directed Graph ( that will be 
 * passed as an input parameter ) according to the format mentioned
 * N -> (x, c, f) -> (x, c, f) -> ....
 */

void PrintGraph(GRAPH G){
    printf("\n");
    for (int i = 1; i <= G.V; i++) {
        printf("%d", G.H[i].x);    // print the id of the current node that is being traversed

        // iteratively taverse through each node of the
        // linked list and print their contents

        for (EDGE *cur = G.H[i].p; cur; cur=cur->next){
            printf(" -> ");
            printf("(%d, %d, %d)", cur->y, cur->c, cur->f);
        }   
        printf("\n");
    }
}

/**
 * Following function is a search function that checks whether the 
 * current edge from _x to _y is added or not. The function returns
 * a pointer to the edge if the edge _x -> _y is already exists or 
 * else NULL pointer.
 * 
 * Return :
 * ------
 * EDGE* : pointer to the found edge else NULL
 */ 

EDGE* does_edge_exist(GRAPH *g, int _x, int _y){
    for (EDGE *i = g->H[_x].p; i ; i = i->next){
        if(i->y == _y)return i;                         // if found an edge _x -> _y return true
    }
    return NULL;                                           // return false if the edge has not been added yet
}

/**
 * The following function will construct and return a residual graph of the 
 * current network flow graph (that will be passed as a parameter)  
 * 
 *  #########################################
 *  #                                       #
 *  #       For Residual Graph              #
 *  #       ------------------              #
 *  #   1. Ignore the flow in the structure #
 *  #       (Not Required)                  #
 *  #   2. Capacity parameter denotes the   #
 *  #       residual Capacity               #
 *  #                                       #
 *  #########################################
 * 
 * Return 
 * ------
 * GRAPH* : a pointer to a GRAPH type object which contains 
 *          the residual graph 
 */

GRAPH *ConstructResidualGraph(GRAPH *G){
    GRAPH *R = (GRAPH*)malloc(sizeof(GRAPH));                   // Initialize the residual Graph
    R->V = G->V;
    R->H = (VERTEX*)malloc((R->V+1)*sizeof(VERTEX));
    int edge = 0;
    for(int i=1;i<=R->V;i++)R->H[i].p=NULL;

    for (int i = 1; i <= G->V; i++){
        R->H[i].x = i;
        R->H[i].n = 0;
        for (EDGE* j = G->H[i].p; j; j = j->next) {
            if(does_edge_exist(R, i, j->y))continue;            // if edge is already added no need to process it again     
            
            int c1, c2, f1, f2;
            c1 = j->c;
            f1 = j->f;
            EDGE *rev = does_edge_exist(G, j->y, i);            // reverse edge of the current edge
            if(rev){
                c2 = rev->c;
                f2 = rev->f;
            } else c2=0, f2=0;
            
            // f1 - f2 gives the net flow
            if(c1-(f1-f2) > 0){                                 // Add the edge only if residual capacity is greater than 0
                edge++;
                EDGE *new = (EDGE*)malloc(sizeof(EDGE));
                new->y = j->y;                                  // end point remains the same
                new->c = c1-(f1-f2);                            // set the residual capacity
                new->f = 0;                                     // no need of the flow
                if(R->H[i].p)new->next = R->H[i].p;             // set the next pointer
                else new->next = NULL;
                R->H[i].p = new;
            }

            if(c2-(f2-f1) > 0){                                 // Add the edge only if residual capacity is greater than 0
                edge++;
                EDGE *new = (EDGE*)malloc(sizeof(EDGE));
                new->y = i;                                     // end point remains the same
                new->c = c2-(f2-f1);                            // set the residual capacity
                new->f = 0;                                     // no need of the flow
                if(R->H[j->y].p)new->next = R->H[j->y].p;       // set the next pointer
                else new->next = NULL;
                R->H[j->y].p = new;
            }            
        }
    }
    R->E = edge;
    return R;
}

/**
 * The following function will do a breadth first search for 
 * the smallest augmenting path with the maximumal residual
 * capacity and return the flow so obtained.
 * 
 * Return 
 * ------
 * int : the flow for the required augmenting path
 */

int BFS(GRAPH *G, int s, int t, int *p){
    int *flow = (int*)malloc((G->V+1) * sizeof(int));   // array to store the maximal residual capacity to an element along the path discovered
    int *dist = (int*)malloc((G->V+1) * sizeof(int));   // array to store the distance of the node frome the root node
    int *vis = (int*)malloc((G->V+1) * sizeof(int));    // array to store whether the current node is visited or not
    for (int i = 1; i <= G->V; i++) {
        flow[i] = 0;
        dist[i] = -1;
        vis[i] = 0;
        p[i] = -1;
    }
    pair *front, *back;                                 // queue implementation sing a linked list
    
    //-----------INITIALISATION------------
    dist[s] = 0;
    flow[s] = INF;
    front = (pair*)malloc(sizeof(pair));
    front->node = s;
    front->next = NULL;
    back = front;

    int found = 0;                                      // flag that checks whether augmnting path is found or not
    int aug_path_flow = 0;                              // stores the flow of the maximaum residual capacity path
    while(front){                             // loop for implementing the breadth first search
        pair *cur = front;
        vis[cur->node]=1;

        for (EDGE *i = G->H[cur->node].p; i; i=i->next){

            if(!vis[i->y] && (dist[i->y]==-1 || dist[i->y]==1+dist[cur->node]) && MIN(flow[cur->node],i->c)>=flow[i->y] && i->c>0){
                if(p[i->y]==-1){                                    // make a new node in queue since node has not been yet discovered
                    pair *tmp = (pair*)malloc(sizeof(pair));

                    tmp->next = NULL;
                    tmp->node = i->y;

                    back->next = tmp;
                    back = back->next;
                } 

                flow[i->y] = MIN(flow[cur->node],i->c);             // Update the flow
                p[i->y] = cur->node;                                // update the parent
                dist[i->y] = dist[cur->node]+1;                     // update the distance from the root

                if(i->y==t){                                        // if the destination is found then store the flow
                    found = 1;
                    aug_path_flow = MAX(aug_path_flow, flow[i->y]);
                }
            }
        }   
        front = front->next;
    }
    return found ? aug_path_flow : 0;                               // if an augmented path is discovered then return the flow else return 0
}

/**
 * The following function is used to compute the 
 * maximum flow between the source and the sink
 * in the given Graph which will be given as input 
 */

void ComputeMaxFlow(GRAPH *G, int s, int t){
    int flow = 0;                                       // total flow of all augmenting path till now 
    int new_flow;                                       // flow of the current augmenting path
    int ctr = 1;                                        // counts the number of the augmenting path
    int *par = (int*)malloc((G->V+1)*sizeof(int));               // storing the parent of each node in the current augmenting path

    struct node{                                        // local stucture for making a stack for the augmenting path
        int v;                                          // node id
        struct node *nxt;                               // pointer to next node
    } *head;
    head = NULL;

    GRAPH *R = ConstructResidualGraph(G);
    while(new_flow = BFS(R, s, t, par)){
        flow += new_flow;                               // add the new flow value to total flow
        printf("\n\
                 ================     AUGMENTING PATH %d :     ================ \n\n\
                 FLOW (of the augmented path): %d\n\
                 TOTAL FLOW:                   %d\n\
                 FLOW PATH:  ", ctr++, new_flow, flow);      

        int cur = t;
        while(cur != s){                                    // loop to augment the flow in the flow graph and the residual graph
            // augment the flow returned from the BFS
            EDGE *a = does_edge_exist(G, par[cur], cur);    // get the forward edge
            EDGE *b = does_edge_exist(G, cur, par[cur]);    // get the backward edge
            if(b==NULL || b->f==0){
                a->f += new_flow;
            } else {
                if(b->f >= new_flow){
                    b->f -= new_flow;
                } else {
                    a->f = new_flow - b->f;
                    b->f = 0;
                }
            }
            EDGE *ra = does_edge_exist(R, par[cur], cur);                       // Forward  edge in the residual graph
            EDGE *rb = does_edge_exist(R, cur, par[cur]);                       // backward edge in the residual graph 

            ra->c -= new_flow;
            if(rb==NULL){                                                       // if backward edge didnt exist , then create one
                rb = (EDGE*)malloc(sizeof(EDGE));
                rb->next = R->H[cur].p;
                rb->f = 0;
                rb->c = 0;
                rb->y = par[cur];
                R->H[cur].p = rb;
            }
            rb->c += new_flow;                                                  // append the forward flow

            if(ra->c==0){                                                       // if the edge has 0 residual capacity the remove it from the graph
                if(R->H[par[cur]].p->y == ra->y) R->H[par[cur]].p = NULL;
                else {
                    for (EDGE *i = R->H[par[cur]].p; i->next; i=i->next){
                        if(i->next->y==ra->y){
                            i->next = i->next->next;
                            break;
                        }
                    }
                }
            }

            cur = par[cur];
            struct node *new = (struct node*)malloc(sizeof(struct node));           // make a new node to push in the stack containing the path
            new->nxt = head;                                                        // initialize the head
            new->v = cur;                                                           // initialize the node id
            head = new;
        }  
        for (struct node *i = head; i; i=i->nxt){                                   // printing the current augmented poth
            printf("%d -> ", i->v);
        }
        printf("%d\n",t);
        head = NULL;
    }
    printf("\n\
            -------------------------------------\n\
            - NO MORE AUGMENTING PATH REMAINING -\n\
            - ALGORITHM COMPLETED SUCCESSFULLY  -\n\
            -------------------------------------\n\n\
            MAXFLOW OF THE NETWORK GRAPH : %d UNITS\n\n", flow);
}

/**
 * The following function is used to verify whether the 
 * flow assignment in the graph is a valid assignment by 
 * comparing it to the need of the vertex and returns
 * true or false based on that
 * 
 * Return
 * ------
 * 0 -> Invalid assignment
 * 1 -> Valid assignment
 */
int valid_assignment_of_flow(GRAPH *G){
    int *net_flow = (int*)malloc((G->V+1)*sizeof(int));
    for (int i = 1; i <= G->V; i++) net_flow[i] = 0;
    for (int i = 1; i <= G->V; i++) {
        for (EDGE *j = G->H[i].p; j; j=j->next){
            net_flow[i] -= j->f;
            net_flow[j->y] += j->f;
        }
    }
    for (int i = 1; i <= G->V; i++){
        if(net_flow[i]!=G->H[i].n)return 0;
    }
    return 1;
}

/**
 * the following function returns a flow assignment if 
 * a valid flow assignment exsts according to the needs 
 * of the vertex.
 * 
 * For a achieving this we first covert out problem 
 * to a max flow problem by appending 2 more vertices 
 * 
 * - > SUPER SOURCE -> V+1 th vertex
 * - > SUPER SINK   -> v+2 th vertex 
 * 
 * then add an edge from the super source to the node
 * if the need of the node is negative (i.e. it is a 
 * producer) and add an edge from th node to super sink 
 * if the need is positive (i.e. the node is a consumer)
 * 
 * the capacity of the above vertices will be equal to
 * the absolute value of the need of the vertex
 */

void NeedBasedFlow(GRAPH *G){
    int v = G->V;
    int e = G->E;

    G->H = (VERTEX*)realloc(G->H,((v+2)+1)*sizeof(VERTEX));     // make space for two new vertices+
    G->V+=2;

    //----------INITIALIZE SUPER SOURCE-------------
    {
        G->H[v+1].n=-INF;
        G->H[v+1].x=v+1;
        G->H[v+1].p=NULL;
    }
    //----------INITIALIZE SUPER SINK-------------
    {
        G->H[v+2].n=-INF;
        G->H[v+2].x=v+2;
        G->H[v+2].p=NULL;
    }

    // v+1->source , v+2->sink
    for (int i = 1; i <= v; i++){
        if(G->H[i].n > 0){
            
            // Add Edges for the consumer nodes
            EDGE *tmp = (EDGE*)malloc(sizeof(EDGE));
            tmp->y = v+2;
            tmp->next = G->H[i].p;
            tmp->f = 0;
            tmp->c = G->H[i].n;
            G->H[i].p = tmp;
            G->E++;
        
        } else if(G->H[i].n < 0){
            
            // add edges for the producer nodes
            EDGE *tmp = (EDGE*)malloc(sizeof(EDGE));
            tmp->y = i;
            tmp->next = G->H[v+1].p;
            tmp->f = 0;
            tmp->c = -1*G->H[i].n;
            G->H[v+1].p = tmp;
            G->E++;
        
        }
    }
    
    // compute the flow assignment for the newly created graph
    ComputeMaxFlow(G, v+1, v+2);

    // restore the graph components
    G->V = v;                                                           // restore the number of vertex
    G->E = e;                                                           // restore the number of edges
    G->H = (VERTEX*)realloc(G->H,(v+1)*sizeof(VERTEX));                 // reallocate the initial space of the number vertex
    for (int i = 1; i <= v; i++){
        if(G->H[i].n > 0)
            G->H[i].p = G->H[i].p->next;
    }

    // check if the flow assignment is a valid assignment or not
    if(valid_assignment_of_flow(G)){
        printf("\n\n\
                '########:::'#######:::'######:::'######::'####:'########::'##:::::::'########:\n\
                 ##.... ##:'##.... ##:'##... ##:'##... ##:. ##:: ##.... ##: ##::::::: ##.....::\n\
                 ##:::: ##: ##:::: ##: ##:::..:: ##:::..::: ##:: ##:::: ##: ##::::::: ##:::::::\n\
                 ########:: ##:::: ##:. ######::. ######::: ##:: ########:: ##::::::: ######:::\n\
                 ##.....::: ##:::: ##::..... ##::..... ##:: ##:: ##.... ##: ##::::::: ##...::::\n\
                 ##:::::::: ##:::: ##:'##::: ##:'##::: ##:: ##:: ##:::: ##: ##::::::: ##:::::::\n\
                 ##::::::::. #######::. ######::. ######::'####: ########:: ########: ########:\n\
                ..::::::::::.......::::......::::......:::....::........:::........::........::\n\
                \n\
                #############################################\n\
                #                                           #\n\
                #       FLOW ASSIGNMENT IS POSSIBLE         #\n\
                #         ACCORDING TO THE NEEDS            #\n\
                #                                           #\n\
                #############################################\n");
    } else {
        printf("\n\n\
                '####:'##::: ##:'##::::'##::::'###::::'##:::::::'####:'########::\n\
                . ##:: ###:: ##: ##:::: ##:::'## ##::: ##:::::::. ##:: ##.... ##:\n\
                : ##:: ####: ##: ##:::: ##::'##:. ##:: ##:::::::: ##:: ##:::: ##:\n\
                : ##:: ## ## ##: ##:::: ##:'##:::. ##: ##:::::::: ##:: ##:::: ##:\n\
                : ##:: ##. ####:. ##:: ##:: #########: ##:::::::: ##:: ##:::: ##:\n\
                : ##:: ##:. ###::. ## ##::: ##.... ##: ##:::::::: ##:: ##:::: ##:\n\
                '####: ##::. ##:::. ###:::: ##:::: ##: ########:'####: ########::\n\
                ....::..::::..:::::...:::::..:::::..::........::....::........:::\n\
                \n\
                #############################################\n\
                #                                           #\n\
                #     FLOW ASSIGNMENT IS NOT POSSIBLE       #\n\
                #         ACCORDING TO THE NEEDS            #\n\
                #                                           #\n\
                #         NEEDS ARE NOT SATISFIED           #\n\
                #    IGNORE THE FLOW ASSIGNMENT BELO        #\n\
                #                                           #\n\
                #############################################\n");
    }
}

int main(int argc, char **argv){
    char file[100];
    int s, t;
    printf("ENTER THE FILE NAME: ");
    scanf("%s",file);

    GRAPH *G = ReadGraph(file);
    PrintGraph(*G);

    printf("Enter Two number corresponding to SOURCE and SINK\n");
    printf("SOURCE:  ");
    scanf("%d",&s);
    printf("SINK:    ");
    scanf("%d",&t);
    ComputeMaxFlow(G, s, t);
    
    PrintGraph(*G);
    
    G = ReadGraph(file);
    NeedBasedFlow(G);
    PrintGraph(*G);

    return 0;
}