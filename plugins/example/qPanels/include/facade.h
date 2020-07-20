#ifndef FACADE_H
#define FACADE_H

#include <vector>
#include <cmath>
#include <fstream>

#define FACADE_H

using namespace std;

// orientation for supporting areas
// ---------- horizontal
// | | | vertical
// | | |
// | | |
enum class orientation {
    vertical,
    horizontal
};

struct Point {
    double x;
    double y;

    bool operator == (const Point& p) const {
        return x == p.x && y == p.y;
    }
    
    bool operator < (const Point& p) const {
        return y < p.y || (y == p.y && x < p.x);
    }
};

struct Rectangle {
    double x[2];
    double y[2];

    double area() {
        return (y[1] - y[0]) * (x[1] - x[0]);
    }

    bool operator < (const Rectangle& p) const {
        return x[0] < p.x[0] || (x[0] == p.x[0] && y[0] < p.y[0]);
    }
    bool operator == (const Rectangle& p) const {
        return x[0] == p.x[0] && y[0] == p.y[0] && x[1] == p.x[1] && y[1] == p.y[1];
    }
    
    bool overlap(const Rectangle& r) {
        // if either rectangle is to the left of the other, or if either rectangle is on top of the other
        if (x[0] >= r.x[1] || r.x[0] >= x[1] || y[0] >= r.y[1] || r.y[0] >= y[1]) {
            return false;
        }
        return true;
    } 
    
    bool overlapStrict(const Rectangle& r) {
        if (x[0] > r.x[1] || r.x[0] > x[1] || y[0] > r.y[1] || r.y[0] > y[1]) {
            return false;
        }
        return true;
    }

    bool overlap(const Rectangle& r, double delta) {
        if (x[0] >= r.x[1] + delta || r.x[0] - delta >= x[1] || y[0] >= r.y[1] + delta || r.y[0] - delta >= y[1]) {
            return false;
        }
        return true;
    }

    bool overlap(Point p, double delta) {
        if (y[0] - delta <= p.y && p.y < y[1] + delta && x[0] - delta <= p.x && p.x < x[1]) {
            return true;
        }
        if (x[0] - delta <= p.x && p.x < x[1] + delta && y[0] - delta <= p.y && p.y < y[1]) {
            return true;
        }
        return false;
    }

    void print() {
        printf("(%f, %f) width: %f height: %f\n", x[0], y[0], x[1] - x[0], y[1] - y[0]);
    }

};

struct Mesh {
    double unit;
    double *mesh = NULL;
    int m; // rows = y
    int n; // columns = x

    // ~Mesh() {
    //     if (mesh) free(mesh);
    //     mesh = NULL;
    // }

    Mesh() {
        unit = 0;
        mesh = NULL;
        m = 0;
        n = 0;
    }

    Mesh(const Mesh& mesh) {
        this->unit = mesh.unit;
        this->m = mesh.m;
        this->n = mesh.n;
        this->mesh = (double *) calloc(this->m*this->n, sizeof(double));
        memcpy(this->mesh, mesh.mesh, m*n*sizeof(double));
    }

    double getDist(double x, double y) {
        int d = (floor(x / unit) + floor(y / unit) * n);
        return mesh[d];
    }
    double getDist(int i , int j) {
        return mesh[i*n + j];
    }

    void setVal(double x, double y, double val) {
        int d = (int) (floor(x / unit) + floor(y / unit) * n);
        if (d < 0 || d >= m*n) return;
        mesh[d] = val;
        return;
    }

    bool loadFromFile(string path) {
        ifstream inFile;
        inFile.open(path);

        if (!inFile) {
            printf("Error opening file.\n");
            return false;
        } 

        double* temp = NULL;
        try {
            inFile >> unit;
            inFile >> m >> n;
            temp = (double *) calloc(m*n, sizeof(double));
            if (!temp) {
                throw;
            }
            for (int i=0; i < m*n; i++) {
                inFile >> temp[i];
            }
        } catch (exception& e) {
            if (temp != NULL) free(temp);
            printf("An exception occurred\n");
            return false;
        }

        mesh = temp;

        inFile.close();
        return true;
    }
};

struct SupportingArea : public Rectangle {

    orientation dir;
    double load;

    bool operator < (const SupportingArea& sa) const {
        if (sa.dir == orientation::vertical && dir == orientation::horizontal) {
            return true;
        } else if (dir == sa.dir) {
            if (dir == orientation::horizontal) {
                return (y[0] < sa.y[0]) || (y[0] == sa.y[0] && x[0] < sa.x[0]);
            } else {
                return (x[0] < sa.x[0]) || (x[0] == sa.x[0] && y[0] < sa.y[0]);
            }
        } else {
            return false;
        } 
    }
};



struct Frame : public Rectangle {
    int marks; 
};

struct Panel : public Rectangle{

    vector<Frame> covered;
    vector<SupportingArea*> sa_covered_h;
    vector<SupportingArea*> sa_covered_v;
    SupportingArea* current_sa;
};

class FacadeSolver {
  double h_fac;
  double w_fac;

  double max_dist = 0.8;

  double min_height;
  double min_width;
  double max_height;
  double max_width;
  
  double delta; 

  double epsilon = 10e-7;

  Mesh mesh; 

  vector<Frame> frames;
  vector<Panel> solution;
  vector<Point> origins;

  vector<SupportingArea> horizontal;
  vector<SupportingArea> vertical;

  bool show = true;

  public:
    bool solve(void);
    void printSolution();
    FacadeSolver();
    void toggleGraphics();
    vector<Panel> getSolution();
    void init(double, double, double, double, double, double, double, vector<Frame>, vector<SupportingArea>, vector<SupportingArea>, Mesh);
  private:
    void adjustForFrames(Panel&);
    void reduceDimensions(Panel&, double, double);
    bool checkDimensions(Panel);
    void getEndPoint(Panel&);
    void adjustForDepthScan(Panel&);
    void populateCoveredSA(Panel&);
    void unpopulateCoveredSA(Panel &p);
    bool checkMesh(Panel&);
    void shrinkToMesh(Panel&);
    void drawState(int);
    
}; // class FacadeSolver

#endif // FACADE_H