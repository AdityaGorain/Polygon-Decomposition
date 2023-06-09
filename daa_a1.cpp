#include <bits/stdc++.h> // includes a large collection of standard C++ libraries
#include <chrono>
#include <algorithm>
#include <climits>
#include <new>
using namespace std;

/// These are the classes for the data structure
class Vertex;
class Edge;
class Face;
vector<pair<Vertex, Vertex>> diagonals;
vector<int> diagonals1;
vector<vector<pair<Vertex,Vertex>>> di;
vector<vector<int>> di1;

/// It is the class for the vertex, which stores its x and y coordinates and edge index
class Vertex
{
public:
    /// x-coordinate of the vertex
    double x;
    /// y-coordinate of the vertex
    double y;
    /// Index of any of the edges that connect to the vertex
    int edge; 

    /// It initializes the vertex's x and y coordinates and sets the edge index to -1
    /// @param x x-coordinate
    /// @param y y-coordinate
    Vertex(double x, double y) : x(x), y(y), edge(-1) {} 
    bool operator==(Vertex const& obj){
        if(x == obj.x && y == obj.y)
            return true;
        else
            return false;
    }
};

/// It is the class for the edge, which stores information on its origin, twin, previous and next edges, and right face
class Edge
{
public:
    /// Index of the vertex at which the edge starts
    int origin; 
    int destination;
     /// Index of the opposite edge
    int twin; 
    /// Index of the previous edge
    int prev;  
    /// Index of the next edge
    int next;  
    /// Index of the face on the right side of the edge
    int right; 
    /// It initializes the edge's origin vertex and sets its twin, previous, next, and right indices to -1
    Edge(int origin, int destination) : origin(origin), destination(destination), prev(-1), next(-1), right(-1) {} 
};

/// It is the class for the face, which stores the index of one of its edges
class Face
{
public:
    /// Index of one of the edges that make up the face
    int edge; 
    /// It initializes the face's edge index
    Face(int edge) : edge(edge) {} 
};

/// It is the class for the DCEL data structure, which stores the vertices, edges, and faces
class DCEL
{
public:
    /// vector that stores the vertices
    vector<Vertex> vertices; 
    /// vector that stores the edges
    vector<Edge> edges;    
    /// vector that stores the faces  
    vector<Face> faces;      
     map<pair<int, int>, int> vertex_index;

    /// Function to add a vertex and return its index
    int addVertex(double x, double y) 
    {
        vertices.push_back(Vertex(x, y)); // create a new vertex and add it to the vertices vector
        vertex_index[{x, y}] = vertices.size() - 1;
        return vertices.size() - 1;       // return the index of the vertex that was just added
    }

    /// Function to add an edge
    void addEdge(int origin, int destination) 
    {
        edges.push_back(Edge(origin, destination)); // create a new edge with the given origin vertex and add it to the edges vector
        int e1 = edges.size() - 1;                  // store the index of the new edge
        edges.push_back(Edge(destination, origin)); // create a new edge with the given destination vertex and add it to the edges vector
        int e2 = edges.size() - 1;                  // store the index of the new edge
        edges[e1].twin = e2;                        // set the twin index of the first edge to the index of the second edge
        edges[e2].twin = e1;                        // set the twin index of the second edge to the index of the first edge
        vertices[origin].edge = e1;                 // set the edge index of the origin vertex to the index of the first edge
        vertices[destination].edge = e2;            // set the edge index of the destination vertex to the index of the second edge
    }

    /// To connect all the edges that are in DCEL
    void connectAllEdges()
    {
        if (edges.size() < 2)
        {
            cout << "Error: not enough edges to connect" << endl; // A minimum of 2 edges is needed to be able to create a DCEL
            return;
        }

        /// Connect the even-numbered edges to their next and previous edges.
        /// Following even-numbered edges you follow a clockwise traversal.
        for (int i = 0; i < edges.size(); i += 2)
        {
            if (i + 2 < edges.size())
                edges[i].next = i + 2;
            if (i - 2 >= 0)
                edges[i].prev = i - 2;
        }

        /// Connect the odd-numbered edges to their next and previous edges.
        /// Following odd-numbered edges you follow an anti-clockwise traversal.
        for (int i = 1; i < edges.size(); i += 2)
        {
            if (i + 2 < edges.size())
                edges[i].prev = i + 2;
            if (i - 2 >= 0)
                edges[i].next = i - 2;
        }

        // Connect the first and last edges to their neighbors
        edges[edges.size() - 2].next = 0;
        edges[0].prev = edges.size() - 2;
        edges[edges.size() - 1].prev = 1;
        edges[1].next = edges.size() - 1;
    }

    /// Add a new face given an edge index
    void addFace(int edge_index)
    {
        // Check if the edge already has a face
        if (edges[edge_index].right != -1)
        {
            cout << "Face already exists for " << printVertex(vertices[edges[edge_index].origin].x, vertices[edges[edge_index].origin].y) << endl;
            return;
        }

        // Create a new face and add it to the list of faces
        Face f1(edge_index);
        faces.push_back(f1);

        // Set the right face of all edges in the new face
        int i = edge_index;
        do
        {
            edges[i].right = faces.size() - 1;
            i = edges[i].next;
        } while (i != edge_index);
    }

    /// @brief To find next vertex in DCEL
    /// @param index 
    /// @return next vertex
    int nextVertex(int index)
    {
        return edges[2 * index].destination;
    }

    /// @brief To find next vertex in DCEL
    /// @param index 
    /// @return previous vertex
    int prevVertex(int index)
    {
        int i = index;
        int prev = -1;
        do
        {
            prev = i;
            i = nextVertex(i);
        } while (i != index);

        return prev;
    }

    /// Insertion of edge as the final edge inside the DCEL
    void insertEdge(int from, int to) 
    {
        edges[2 * from].destination = to;
        edges[2 * from].next = 2 * to;
        edges[2 * from].twin = 2 * from + 1;

        edges[2 * from + 1].origin = to;
        edges[2 * from].prev = 2 * to;
        edges[2 * from].twin = 2 * from;

        edges[2 * to].prev = 2 * from;
    }

    /// Print all vertices in the DCEL
    void printVertices()
    {
        if (vertices.size() == 0)
            return;
        cout << "Vertices: " << endl;
        cout << "Count: " << vertices.size() << endl;
        for (auto v : vertices)
        {
            cout << "(" << v.x << ", " << v.y << ")" << endl;
        }
        cout << endl;
    }

    /// Print all edges in the DCEL
    void printEdges()
    {
        if (edges.size() == 0)
            return;

        cout << "Edges: " << endl;
        cout << "Count: " << edges.size() << endl;

        // Print each edge and its twin
        for (int i = 0; i < edges.size(); i += 2)
        {
            Edge e = edges[i];
            cout << printVertex(vertices[e.origin].x, vertices[e.origin].y);
            cout << " <=> ";
            cout << printVertex(vertices[edges[e.twin].origin].x, vertices[edges[e.twin].origin].y) << endl;
        }
        cout << endl;
    }

    /// Traverse all edges in the DCEL starting from a given vertex
    /// @param startVertex It is the start vertex from where the polygon starts
    void traverseDCEL(int startVertex)
    {
        cout << "Traversing edges starting from"
             << " (" << vertices[startVertex].x << ", " << vertices[startVertex].y << ")" << endl;

        int i = vertices[startVertex].edge;
        do
        {
            cout << printVertex(vertices[edges[i].origin].x, vertices[edges[i].origin].y) << endl;
            i = edges[i].next;
        } while (i != vertices[startVertex].edge);
        cout << endl;
    }

    /// Print all faces in the DCEL
    void printFaces()
    {
        if (faces.size() == 0)
            return;
        cout << "Faces: " << endl;
        cout << "Count: " << faces.size() << endl;

        // Print each face as a sequence of vertices
        for (auto f : faces)
        {
            cout << printVertex(vertices[edges[f.edge].origin].x, vertices[edges[f.edge].origin].y) << endl;
        }
        cout << endl;
    }

    /// Print the entire DCEL
    void print()
    {
        cout << endl;
        cout << printVertex(vertices[edges[0].origin].x, vertices[edges[0].origin].y);
        for (int i = 0; i < edges.size(); i += 2)
        {
            Edge e = edges[i];
            cout << " <=> ";
            cout << printVertex(vertices[edges[e.next].origin].x, vertices[edges[e.next].origin].y);
        }
        cout << endl;
    }

    void read()
    {
        int n;
        cout << "Enter Number of Vertices" << endl;
        cin >> n;

        for (int i = 0; i < n; i++)
        {
            int x, y;
            cout << "Enter Vertices in the order in which they are connected (Clockwise or Anticlockwise)" << endl;
            cin >> x;
            cin >> y;
            addVertex(x, y);
        }

        for (int i = 0; i < n; i++)
        {
            if (i + 1 < n)
                addEdge(i, i + 1);
            else
                addEdge(i, 0);
        }

        connectAllEdges();
        addFace(0);
        addFace(1);
    }

    /// This function generates a DCEL from a list of vertices
    void generateDCEL(vector<Vertex> L)
    {
        vertices.clear();
        edges.clear();
        faces.clear();
        for (Vertex v : L)
        {
            addVertex(v.x, v.y);
        }

        for (int i = 0; i < vertices.size(); i++)
        {
            if (i + 1 < vertices.size())
                addEdge(i, i + 1);
            else
                addEdge(i, 0);
        }

        connectAllEdges();

        addFace(0);
        addFace(1);
    }

private:
    /// Helper function to print a vertex as a string
    string printVertex(float a, float b)
    {
        return "(" + to_string(a) + ", " + to_string(b) + ")";
    }
};

/// Finding general dot product
float dot_product(Vertex a, Vertex b, Vertex c)
{
    float ba_x = a.x - b.x;
    float ba_y = a.y - b.y;
    float bc_x = c.x - b.x;
    float bc_y = c.y - b.y;
    return ba_x * bc_x + ba_y * bc_y;
}

/// Finding general cross product
float cross_product(Vertex a, Vertex b, Vertex c)
{
    float ba_x = a.x - b.x;
    float ba_y = a.y - b.y;
    float bc_x = c.x - b.x;
    float bc_y = c.y - b.y;
    return ba_x * bc_y - ba_y * bc_x;
}

/// This function checks whether if points are collinear
bool areCollinear(Vertex a, Vertex b, Vertex c)
{
    return (c.y - b.y) * (b.x - a.x) == (b.y - a.y) * (c.x - b.x);
}

/// This function checks if the angle between three vertices is reflex
bool isReflex(Vertex a, Vertex b, Vertex c)
{
    if (areCollinear(a, b, c))
        return false;
    float dot = dot_product(a, b, c);
    float cross = cross_product(a, b, c);
    if (cross < 0)
    {
        return true;
    }
    else if (cross == 0 && dot < 0)
    {
        return true;
    }
    else
    {
        return false;
    }
}

/// It checks whether if a point lies inside a polygon
/// @return true if intersection is odd, else false
bool isVertexInPolygon(Vertex point, vector<Vertex> polygon) {
    int intersections = 0;
    int n = polygon.size();

    for (int i = 0; i < n; i++) {
        Vertex p1 = polygon[i];
        Vertex p2 = polygon[(i+1)%n];

        // Check if the point lies on the edge
        if ((point.y == p1.y) && (point.x >= min(p1.x, p2.x)) && (point.x <= max(p1.x, p2.x))) {
            return true;
        }

        // Check if the ray intersects with the edge
        if ((point.y > min(p1.y, p2.y)) && (point.y <= max(p1.y, p2.y)) && (point.x <= max(p1.x, p2.x)) && (p1.y != p2.y)) {
            double x_intersection = (point.y - p1.y) * (p2.x - p1.x) / (p2.y - p1.y) + p1.x;
            if (p1.x == p2.x || point.x <= x_intersection) {
                intersections++;
            }
        }
    }

    // Return true if the number of intersections is odd, false otherwise
    return (intersections % 2 == 1);
}

/// This function checks whether if point is inside a rectangle
bool insideRectangle(int min_x, int max_x, int min_y, int max_y, int x, int y)
{
    if (x > min_x && x < max_x && y > min_y && y < max_y)
        return true;
    else
        return false;
}

/// This function backtracks when a vertex is making the polygon convcave
vector<Vertex> removeInternalVertices(DCEL P, vector<Vertex> L, int &start, int end)
{

    /// We generate a rectangle by finding out min max of vertices in L
    double min_x = INT_MAX;
    double max_x = INT_MIN;
    double min_y = INT_MAX;
    double max_y = INT_MIN;

    int i = start;
    do
    {
        Vertex t = P.vertices[i];
        for (int j = 0; j < L.size(); j++)
        {
            min_x = min(L[j].x, min_x);
            max_x = max(L[j].x, max_x);
            min_y = min(L[j].y, min_y);
            max_y = max(L[j].y, max_y);
        }
        Vertex previous = P.vertices[P.prevVertex(i)];
        Vertex next = P.vertices[P.nextVertex(i)];
        if (isReflex(previous, t, next) && insideRectangle(min_x, max_x, min_y, max_y, t.x, t.y) 
        && isVertexInPolygon(t, L)
        )
        {
            L.pop_back();
            start = P.prevVertex(start);
            // start = P.prevVertex(start);
            // i = start;
            continue;
        }
        i = P.nextVertex(i);
    } while (i != end);

    return L;
}


/// @brief It breaks the polygon into different parts
/// @param P It is the polygon stored in DCEL
/// @param v It is the starting vertex
/// @return Different partitions of the polygon
vector<DCEL> decompose(DCEL P, int v)
{

    vector<DCEL> OS;
    vector<Vertex> L;
    DCEL D;

    int index = v;                // to traverse the DCEL P
    int size = P.vertices.size(); // storing the size as it will keep changing
    int nextEdgeAt = v;
    while (size > 2)
    {                  // while there are more than 2 vertices in the polygon
        int count = 0; // to keep count of vertices traversed
        int start = index;
        while (count < size)
        { // while all vertices in P are not traversed
            if (L.size() > 1)
            {
                Vertex v_curr = L[L.size() - 1];   // Latest vertex that we pushed
                Vertex v_next = P.vertices[index]; // The one which we want to add
                Vertex v_prev = L[L.size() - 2];   // Added previously
                Vertex v1 = L[0];                  // First Vertex
                Vertex v2 = L[1];                  // Second Vertex

                if (isReflex(v_prev, v_curr, v_next) || isReflex(v_curr, v_next, v1) || isReflex(v_next, v1, v2)) // Main condition to check if the polygon we are trying to form using L has any notches
                    break;
                else
                    L.push_back(v_next);
            }
            else
            {
                L.push_back(P.vertices[index]);
            }

            count++;
            index = P.nextVertex(index);
        }

        if (L.size() < 3)
            return {};

        int end = P.prevVertex(start);
        L = removeInternalVertices(P, L, index, end); // We pass index because we need to go back everytime we pop an element in L
        if (L.size() > 2)
        {
            D.generateDCEL(L); // Generate a DCEL out of the vertex list L
        }
        else
            return {};

        index = P.prevVertex(index);

        P.insertEdge(nextEdgeAt, index);
        nextEdgeAt = index;
        size -= L.size() - 2; // Reduce size of P after decomposition
        L.clear();
        OS.push_back(D);
    }

    return OS;
}

/// This function extracts digonals from a partitioned DCEL pool
vector<pair<Vertex, Vertex>> extractDiagonals(vector<DCEL> D, int j){
    vector<pair<Vertex, Vertex>> ans;
    for(int i = 0 ; i < D.size()- 1; i++){
        ans.push_back({D[i].vertices[D[i].vertices.size()-1], D[i].vertices[0]});
        di1[j].push_back(i);
    }

    return ans;
}


vector<vector<pair<Vertex,Vertex>>> diagonalSet(vector<vector<DCEL>> OS){
    vector<vector<pair<Vertex,Vertex>>> ans;
    for(int i = 0; i < OS.size(); i++){
        di1.push_back({});
        ans.push_back(extractDiagonals(OS[i], i));

    }
        return ans;
}


/// This function takes a 2D vector of DCEL (Doubly Connected Edge List) objects as input
/// and fills a diagonal set (stored in the `di` and `di1` vectors) for each DCEL object
void fillDiagonalSet(vector<vector<DCEL>> OS){

    for(int i=0; i<OS.size(); i++){ // iterate over each DCEL object
        
        for(int j=0; j < OS[i].size(); j++){ // iterate over each face of the current DCEL object
            
            for(int k = 0; k< OS[i][j].vertices.size(); k++){ // iterate over each vertex of the current face

                Vertex x = OS[i][j].vertices[k]; // current vertex
                Vertex y(0,0); // next vertex
                if(k == OS[i][j].vertices.size()-1)  y = OS[i][j].vertices[0]; // if last vertex, wrap around to the first vertex
                else { y = OS[i][j].vertices[k+1];}
                
                for(int l=0; l<di[i].size(); l++){ // iterate over each existing diagonal of the DCEL object
                    Vertex xa = di[i][l].second; // start vertex of diagonal
                    Vertex ya = di[i][l].first; // end vertex of diagonal

                    // if the current and next vertices of the current face match the start and end vertices of an existing diagonal,
                    // then add the diagonal to the diagonal set for this face
                    if(xa.x == x.x && xa.y == x.y && ya.x == y.x && ya.y == y.y)
                    {
                        di[i].push_back({xa,ya});
                        di1[i].push_back(j);
                    }
                }


            }
            
        }
    }
}



vector<vector<DCEL>> minimizePartitions(DCEL P)
{
    int cardinality = INT_MAX; // initialize the minimum cardinality to a very large number
    vector<vector<DCEL>> OS; // initialize the output 2D vector of DCEL objects
    for (int i = 0; i < P.vertices.size(); i++) // iterate over each vertex of the input DCEL object
    {   
        vector<DCEL> D = decompose(P, i); // decompose the input DCEL object with the current vertex as the starting point
        if (D.size() != 0) // if the decomposition is not empty
        {
            if (D.size() == cardinality) // if the size of the decomposition is the same as the current minimum cardinality
            {
                OS.push_back(D); // add the decomposition to the output vector
            }
            else if (D.size() < cardinality) // if the size of the decomposition is smaller than the current minimum cardinality
            {
                cardinality = D.size(); // update the minimum cardinality
                OS.clear(); // clear the output vector
                OS.push_back(D); // add the decomposition to the output vector
            }
        }
    }
    return OS; // return the output vector of DCEL objects
}

/// checks if point is concave
bool isConvex(DCEL D, Vertex v)
{
    Vertex a = D.vertices[D.prevVertex(D.vertex_index[{v.x, v.y}])];
    Vertex c = D.vertices[D.nextVertex(D.vertex_index[{v.x, v.y}])];

    return !isReflex(c, v, a);
}

/// returns angle
double polar_angle(Vertex p, Vertex centroid) {
    double dx = p.x - centroid.x;
    double dy = p.y - centroid.y;
    return atan2(dy, dx);
}

/// generation of centroid
Vertex centroid(vector<Vertex>& points) {
    int n = points.size();
    double sum_x = 0, sum_y = 0;
    for (int i = 0; i < n; i++) {
        sum_x += points[i].x;
        sum_y += points[i].y;
    }
    return {sum_x / n, sum_y / n};
}

/// This function generates more polygons by merging two existing polygons
void generateMorePolygons(vector<DCEL> &D, int j, int u, vector<DCEL> &morePolygons, Vertex vs, Vertex vt){
    
    // create a new DCEL object and a vector of vertices to hold the new polygon
    DCEL d;
    vector<Vertex> V;

    // get the number of vertices in the two polygons to merge
    int n1 = D[j].vertices.size();
    int n2 = D[u].vertices.size();

    // loop through the vertices of the first polygon and add them to the new polygon
    for(int i=0; i<n1; i++)
    {
        V.push_back(D[j].vertices[i]);
        
        // optional: print the coordinates of each vertex in the first polygon
        // cout << D[j].vertices[i].x << "," << D[j].vertices[i].y << endl;
    }

    // loop through the vertices of the second polygon and add them to the new polygon,
    // except for the start and end points (vs and vt) of the common edge
    for(int i=0; i<n2; i++)
    {   
        if((D[u].vertices[i].x == vs.x && D[u].vertices[i].y == vs.y) || (D[u].vertices[i].x == vt.x && D[u].vertices[i].y == vt.y)){}
        else {
            V.push_back(D[u].vertices[i]);
        }
    }

    // find the centroid of the new polygon
    Vertex centroid_point = centroid(V);
    
    // sort the vertices of the new polygon in counterclockwise order around the centroid
    sort(V.begin(), V.end(), [&](Vertex p1, Vertex p2) {
        return polar_angle(p1, centroid_point) < polar_angle(p2, centroid_point);
    });

    // generate a new DCEL object for the new polygon and add it to the vector of DCEL objects
    d.generateDCEL(V);
    D.push_back(d);
}


void merge(vector<DCEL>& D, vector<pair<Vertex, Vertex>> LLE, vector<int> LLE1, DCEL d)
{   
    vector<DCEL> D1; // Create a new vector to hold merged polygons

    int m = LLE.size()/2;
    int NP = m + 1;
    vector<DCEL> morePolygons; // Create a new vector to hold polygons generated during merging process

    // Initialize arrays and maps
    vector<int> LDP;
    for (int i = 0; i < NP; i++)
    {
        LDP.push_back(true);
    }

    vector<int> LUP;
    for (int i = 0; i < NP; i++)
    {
        LUP.push_back(i);
    }

    unordered_map<int, int> marker;
    for (int i = 0; i < LLE.size(); i++)
    {
        marker[i] = 0;
    }

    std::map<pair<int, int>, vector<pair<int, Vertex>>> LPv;

    // Fill in LPv with edge endpoints for each polygon
    for (int i = 0; i <=NP-1 ; i++)
    {
        for (int j = 0; j < D[i].edges.size(); j=j+2) // iterate over edges
        {
            Vertex start = D[i].vertices[D[i].edges[j].origin];
            Vertex end = D[i].vertices[D[i].edges[(D[i].edges[j].twin)].origin];

            int flag = 0;

            // Check if edge is already marked as used
            for (int k = 0; k < LLE.size(); k++)
            {
                if (marker[k] == 0 &&((LLE[k].first.x == start.x && LLE[k].first.y == start.y && LLE[k].second.x == end.x && LLE[k].second.y == end.y) ))
                {
                    marker[k] = 1;
                    LPv[{start.x, start.y}].push_back({i, end});
                    flag = 1;
                }
                if ( marker[k] == 1 && ((LLE[k].first.x == start.x && LLE[k].first.y == start.y && LLE[k].second.x == end.x && LLE[k].second.y == end.y)) || (LLE[k].second.x == start.x && LLE[k].second.y == start.y && LLE[k].first.x == end.x && LLE[k].first.y == end.y) )
                {
                    flag = 1;
                }
            }
            if (flag == 0)
                LPv[{start.x, start.y}].push_back({i, end});

        }
    }

    // Merge polygons

    for (int j = 0; j < m; j++)
    {
        
        Vertex vs = LLE[j].first;
        Vertex vt = LLE[j].second;


        
        if ((LPv[{vs.x, vs.y}].size() > 2 && LPv[{vt.x, vt.y}].size() > 2) || (LPv[{vs.x, vs.y}].size() > 2 && isConvex(d, vt)) || (LPv[{vt.x, vt.y}].size() > 2 && isConvex(d, vs)) || (isConvex(d, vs) && isConvex(d, vt)))
        {
            int j11 = LLE1[j];
            Vertex j2 = vt;
            Vertex i2 = vs;
            Vertex j3 = D[j11].vertices[D[j11].nextVertex(D[j11].vertex_index[{vt.x, vt.y}])];
            Vertex i1 = D[j11].vertices[D[j11].prevVertex(D[j11].vertex_index[{vs.x, vs.y}])];

            int u;
            for (int i = 0; i < LPv[{vt.x, vt.y}].size(); i++)
            {
                if (LPv[{vt.x, vt.y}][i].second.x == vs.x && LPv[{vt.x, vt.y}][i].second.y == vs.y)
                {
                    u = LPv[{vt.x, vt.y}][i].first;
                    break;
                }
            }

            Vertex j1 = D[u].vertices[D[u].prevVertex(D[u].vertex_index[{vt.x, vt.y}])];
            Vertex i3 = D[u].vertices[D[u].nextVertex(D[u].vertex_index[{vs.x, vs.y}])];

            // if condition to recognise and remove a diagonal
            if (!isReflex(i1, i2, i3) && !isReflex(j1, j2, j3))
            {
                NP++;
                generateMorePolygons(D, LUP[u], LUP[j11], morePolygons, vs, vt);

                LDP[j11] = false;
                LDP[u] = false;
                
                int j111  = LUP[j11];
                int u111 = LUP[u];

                LUP[j11] = NP-1;
                LUP[u] = NP-1;

                for (int h = 0; h < NP -1; h++)
                {
                    if (LUP[h] == j111 || LUP[h] == u111)
                        LUP[h] = NP-1;
                }
                LUP.push_back(NP-1);
            }
        }
    }


    for(int i=0; i<LUP.size(); i++)
    {
        if(i == LUP[i]) D1.push_back(D[i]);
    }
    D = D1;

}

/// @brief It is the main function
/// @return zero
int main()
{
    
    auto start = chrono::high_resolution_clock::now();
    DCEL dcel;

    // dcel.generateDCEL({Vertex(0,0), Vertex(-0.5,0.5), Vertex(-1,1), Vertex(-0.5,1.5), Vertex(-1,3), Vertex(-0.5,2.5), Vertex(0,3), Vertex(0.5,2.5), Vertex(1,3), Vertex(1.5,2.5), Vertex(2,3), Vertex(2.5,2.5), Vertex(2,2), Vertex(2.5,1.5), Vertex(2,1), Vertex(2.5,0.5), Vertex(2,0), Vertex(1.5,0.5), Vertex(1,0), Vertex(0.5,0.5)});
    dcel.generateDCEL({Vertex(0,0), Vertex(0,3), Vertex(1,3), Vertex(2,2), Vertex(3,3), Vertex(4,3), Vertex(5,2), Vertex(6,3), Vertex(7,3), Vertex(8,2), Vertex(9,3), Vertex(10,3), Vertex(10,0), Vertex(9,0), Vertex(8,1), Vertex(7,0), Vertex(6,0), Vertex(5,1), Vertex(4,0), Vertex(3,0), Vertex(2,1), Vertex(1,0)});
    vector<vector<DCEL>> OS = minimizePartitions(dcel); // A vector of the partitions created with the minimum possible cardinality
    di = diagonalSet(OS);   
    fillDiagonalSet(OS);    // getting diagonals

    vector<int> sizeOfOS;
    for (int i=0; i< OS.size(); i++)
    {   
        vector<DCEL> D = OS[i];
        merge(D, di[i], di1[i], dcel); // merge function
        OS[i] = D;
        sizeOfOS.push_back(D.size());

            
    }


    int mini = INT_MAX;
    for (int i=0; i< OS.size(); i++)
    {   
        if(mini > OS[i].size()) {mini = OS[i].size();} 
    }

    for(int i=0; i< OS.size(); i++)
    {
        vector<DCEL> D = OS[i];
        if(OS[i].size() == mini){
        cout << "================================================";
        for (DCEL d : D)
        { 
            d.print();
        }
        }
    }
    

    if(OS.size() == 0) cout << "The given polygon cannot be divided using MP1 algorithm." << endl;
    cout << "================================================" << endl;

    auto stop = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<std::chrono::microseconds>(stop-start);
    cout << duration.count() << endl;
    return 0;
    
}
