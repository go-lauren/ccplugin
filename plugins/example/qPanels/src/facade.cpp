#include "facade.h"
#include "tools.h"
#include <iostream>
#include <algorithm>
#include <cfloat>

FacadeSolver::FacadeSolver() {
    h_fac = 10;
    w_fac = 10;
    Point op = {0, 0};
    origins.push_back(op);
}

void FacadeSolver::init(double w_fac, double h_fac, double min_height, double min_width, double max_height, double max_width, double delta, vector<Frame> frames,vector<SupportingArea> vertical, vector<SupportingArea> horizontal, Mesh mesh) {
    this->w_fac = w_fac; 
    this->h_fac = h_fac;
    this->min_height = min_height;
    this->min_width = min_width;
    this->max_height = max_height;
    this->max_width = max_width;
    this->delta = delta;
    this->frames.insert(this->frames.end(), frames.begin(), frames.end());
    this->vertical.insert(this->vertical.end(), vertical.begin(), vertical.end());
    sort(this->vertical.begin(), this->vertical.end());
    this->horizontal.insert(this->horizontal.end(), horizontal.begin(), horizontal.end());
    sort(this->horizontal.begin(), this->horizontal.end());
    this->mesh = Mesh(mesh);
    // remove values inside frames
    for (vector<Frame>::iterator it = frames.begin(); it != frames.end(); it++) {
        Frame f = *it;
        for (double x = f.x[0] - delta; x <= f.x[1]- delta; x += mesh.unit) {
            for (double y = f.y[0] - delta; y <= f.y[1] - delta; y += mesh.unit) {
                mesh.setVal(x, y, -1);
            }
        }
    }
}

void FacadeSolver::toggleGraphics() {
    show = !show;
}

void FacadeSolver::reduceDimensions(Panel &p, double upper_bound_y, double upper_bound_x) {
    SupportingArea * sa_h = NULL;
    SupportingArea * sa_v = NULL;

    // look for horizontal supporting areas to reduce to
    for (int i = p.sa_covered_h.size() -1; i >= 0; i--) {
        SupportingArea* sa = p.sa_covered_h[i];
        if (upper_bound_y >= sa->y[0]) {
            sa_h = sa;
            break;
        } 
    }
    // for vertical supporting areas, height is limited by the
    double vert_upper_bound = h_fac;
    if (p.sa_covered_v.size() > 1) {
        SupportingArea * sa_left_vertical = p.sa_covered_v.front();  
        // guarantee that a supporting area covers the left side of the panel    
        if (sa_left_vertical->x[0] <= p.x[0] && p.x[0] <= sa_left_vertical->x[1]) {
            vert_upper_bound = sa_left_vertical->y[1];
            // look for supporting areas for the right side of the panel
            for (int i = p.sa_covered_v.size() -1; i >= 0; i--) {
                SupportingArea* sa = p.sa_covered_v[i];
                if (upper_bound_x > sa->x[0]) {
                    sa_v = sa;
                    break;
                } 
            }
        }
    }
    
    if (sa_h == NULL && sa_v == NULL) {
        // could not find a supporting area to reduce to
        p.x[1] = p.x[0];
        p.y[1] = p.y[0];
        return;
    } 

    SupportingArea * new_sa;
    if (sa_h == NULL || sa_v == NULL) {
        new_sa = (sa_h != NULL) ? sa_h : sa_v;
    } else {
        // calculate cost of reducing to either a vertical or horizontal supporting area
        // TODO: maybe midpoint isn't the best matrix
        double new_p_y_1 = min(upper_bound_y, 0.5*(sa_h->y[0] + sa_h->y[1]));
        double new_p_x_1 = min(upper_bound_x, sa_h->x[1]);
        double cost_of_sa_h = p.area() - (new_p_y_1 - p.y[0]) * (new_p_x_1 - p.x[0]);
        new_p_x_1 = min(upper_bound_x, 0.5*(sa_v->x[0] + sa_v->x[1]));
        new_p_y_1 = min(min(vert_upper_bound, upper_bound_y), sa_v->y[1]);
        double cost_of_sa_v = p.area() - (new_p_y_1 - p.y[0]) * (new_p_x_1 - p.x[0]);

        new_sa = (cost_of_sa_h > cost_of_sa_v) ? sa_v : sa_h;
    }

    if (new_sa->dir == orientation::horizontal) {
        p.y[1] = min(upper_bound_y, 0.5*(new_sa->y[0] + new_sa->y[1]));
        p.x[1] = min(upper_bound_x, new_sa->x[1]);
    } else {
        p.y[1] = min(min(vert_upper_bound, upper_bound_y), sa_v->y[1]);
        p.x[1] = min(upper_bound_x, 0.5*(sa_v->x[0] + sa_v->x[1]));
    }
    p.current_sa = new_sa;
    unpopulateCoveredSA(p);
}

/* Remove supporting areas no longer covered by the panel after reducing dimensions */
void FacadeSolver::unpopulateCoveredSA(Panel &p) {
    for (int i = p.sa_covered_h.size() -1; i >= 0; i--) {
        SupportingArea * sa = p.sa_covered_h[i];
        if (sa->y[0] > p.y[1]) {
            p.sa_covered_h.erase(p.sa_covered_h.begin() + i);
        }
    }

    for (int i = p.sa_covered_v.size() - 1; i >= 0; i--) {
        SupportingArea *sa = p.sa_covered_v[i];
        if (sa->x[0] > p.x[1]) {
            p.sa_covered_v.erase(p.sa_covered_v.begin() + i);
        }
    }
}

/* Add all the supporting areas that are overlapped by the panel */
void FacadeSolver::populateCoveredSA(Panel& p) {
    for (int i = 0; i < horizontal.size(); i++) {
        SupportingArea sa = horizontal[i];
        if (p.y[0] <= sa.y[0] && sa.y[0] <= p.y[1] && sa.x[0] <= p.x[0] && p.x[0] <= sa.x[1] ) {
            p.sa_covered_h.push_back(&horizontal[i]);
        } else if (sa.y[0] <= p.y[0] && p.y[0] <= sa.y[1] && sa.x[0] <= p.x[0] && p.x[0] <= sa.x[1] ) {
            p.sa_covered_h.push_back(&horizontal[i]);
        }
    }

    for (int i = 0; i < vertical.size(); i++) {
        SupportingArea sa = vertical[i];
        if (p.x[0] <= sa.x[0] && sa.x[0] <= p.x[1] && sa.y[0] <= p.y[0] && p.y[0] <= sa.y[1]) {
            p.sa_covered_v.push_back(&vertical[i]);
        } else if (sa.x[0] <= p.x[0] && p.x[0] <= sa.x[1] && sa.y[0] <= p.y[0] && p.y[0] <= sa.y[1]) {
            p.sa_covered_v.push_back(&vertical[i]);
        }
    }
}

/* Adjust panel dimensions to avoid frames */
void FacadeSolver::adjustForFrames(Panel &p) {
 
    // leave room for the next panel
    if (p.y[1] < h_fac && h_fac - p.y[1] - mesh.unit < min_height) {
        p.y[1] = h_fac - min_height - mesh.unit;
    }
    if (p.x[1] < w_fac && w_fac - p.x[1] - mesh.unit < min_width) {
        p.x[1] = w_fac - min_width - mesh.unit;
    }

    // find all frames currently overlapping with the panel
    vector<Frame*> stack;
    for (vector<Frame>::iterator it = frames.begin(); it != frames.end(); it++) {
        Frame* f = &*it;
        if (p.overlap(*f, delta)) {
            stack.push_back(f);
        }
    }

    // for each frame in the stack, check if it either overlaps or if it is covered by the panel
    for(vector<Frame*>::iterator it = stack.begin(); it != stack.end(); it++) {
        Frame * f = *it;
        if (p.y[1] - f->y[1] < delta) {
            reduceDimensions(p, f->y[0] - delta - mesh.unit, p.x[1]);
        } else {
            f->marks++;
        }
    }
    for(vector<Frame*>::iterator it = stack.begin(); it != stack.end(); it++) {
        Frame* f = *it;
        if (p.x[1] - f->x[1] < delta) {
            reduceDimensions(p, p.y[1], f->x[0] - delta - mesh.unit);
        } else {
            f->marks++;
        }
    }

    // remove all frames covered by a panel so we don't have to check again
    vector<Frame>::iterator it = frames.begin();
    while (it != frames.end()) {
        Frame *f = &*it;
        if (f->marks == 2) {
            f->marks = 0;
            p.covered.push_back(*f);
            it = frames.erase(it);
        } else {
            f->marks = 0;
            it++;
        }
    }
};

/* Checks if a panel has valid dimensions */
bool FacadeSolver::checkDimensions(Panel p) {
    double width = p.x[1] - p.x[0];
    double height = p.y[1] - p.y[0];
    return width >= min_width && width <= max_width && height >= min_height && height <= max_height;
};

/* Returns an endpoint for a panel, guaranteeing that it does not overlap with any other panels and that it does not exceed width or height of facade*/
void FacadeSolver::getEndPoint(Panel &p) {
    p.x[1] = (p.x[0] + max_width >= w_fac) ? w_fac : p.x[0] + max_width;
    p.y[1] = (p.y[0] + max_height >= h_fac) ? h_fac : p.y[0] + max_height;
    for (vector<Panel>::iterator it = solution.begin(); it != solution.end(); it++) {
        Panel other_panel = *it;
        // check if overlaps other_panel
        if (other_panel.x[0] < p.x[1] && p.x[1] < other_panel.x[1] &&
        other_panel.y[0] < p.y[1] && p.y[1] < other_panel.y[1])
        {
            // TODO: change depending if u want skinner or longer panels...?
            p.x[1] = other_panel.x[0];
        } 
    }
    populateCoveredSA(p);
};

/* Prints the solution */
void FacadeSolver::printSolution() {
    sort(solution.begin(), solution.end());
    for (int i = 0; i < solution.size(); i++) {
        Panel p = solution[i];
        printf("Panel %d: (%f, %f) width: %f height: %f\n", i, p.x[0], p.y[0], p.x[1] - p.x[0], p.y[1] - p.y[0]);
    }
}

/* Checks to see if a panel is within MAX_DIST of the mesh */
bool FacadeSolver::checkMesh(Panel& p) {
    double max_val = -1;
    double min_val = DBL_MAX;
    for (double j = p.y[0]; j < p.y[1]; j += mesh.unit) {
        for (double i = p.x[0]; i < p.x[1]; i += mesh.unit) {
            double d = mesh.getDist(i, j);
            max_val = max(max_val, d);
            min_val = (d < 0 || min_val < 0) ? max(d, min_val) : min(d, min_val);
        }
    }
    return max_val - min_val <= max_dist;
}

/* Reduces panel dimensions till it fits the mesh */
void FacadeSolver::shrinkToMesh(Panel& p) {
    int m = (int) ceil((p.y[1] - p.y[0])/mesh.unit) +1;
    int n = (int) ceil((p.x[1] - p.x[0])/mesh.unit) +1;
    // double cost[m][n];
    // double max_val[m][n];
    // double min_val[m][n];
    
    vector<vector<double>> cost;
    cost.resize(m, vector<double>(n, 0));
    vector<vector<double>> max_val;
    max_val.resize(m, vector<double>(n, 0));
    vector<vector<double>> min_val;
    min_val.resize(m, vector<double>(n, 0));
    
    int offset_i = (int) ceil(p.y[0] / mesh.unit);
    int offset_j = (int) ceil(p.x[0] / mesh.unit); 
    cout << offset_i << " " << offset_j << endl;
    for (int i = 0; i < m; i++) 
    {
        double y = p.y[0] + i * mesh.unit;
        for (int j = 0; j < n; j++)
        {
            double x = p.x[0] + j * mesh.unit;
            double d = mesh.getDist(i + offset_i, j + offset_j);
            double current_max = d;
            double current_min = d;
            if (j != 0)
            {
                current_max = max(current_max, max_val[i][j-1]);
                current_min = (current_min < 0 || min_val[i][j-1] < 0) ? max(current_min, min_val[i][j-1]) : min(current_min, min_val[i][j-1]);
            }
            if (i != 0)
            {
                current_max = max(current_max, max_val[i-1][j]);
                current_min = (current_min < 0 || min_val[i][j-1]) ? max(min_val[i-1][j], current_min) : min(current_min, min_val[i-1][j]);
            }
           
            max_val[i][j] = current_max;
            min_val[i][j] = current_min;
            cost[i][j] = (current_max - current_min > max_dist || x - p.x[0] < min_width || y - p.y[0] < min_height) ? -1 : (x - p.x[0]) * (y - p.y[0]);
            // printf("(i,j): (%d %d), (x, y): (%f %f), max: %f min: %f val: %f %f\n", i, j, x, y, current_max, current_min, d, cost[i][j]);
        }
    }
    
    writeArrayToFile(cost, m, n, "cost.txt");

    while (cost[max((int) floor((p.y[1] - p.y[0])/mesh.unit)-1,0) ][max((int) floor((p.x[1] - p.x[0])/mesh.unit)-1,0)] < 0)
    {
        SupportingArea* current_sa = p.current_sa;
        if (!p.sa_covered_h.empty()) {
            SupportingArea* sa = p.sa_covered_h.back();
            if (sa->x[0] <= p.x[1] && p.x[1] <= sa->x[1] && sa->y[0] <= p.y[1] && p.y[1] <= sa->y[1]) {
                current_sa = sa;
            }
        } 
        if (!p.sa_covered_v.empty()){
            SupportingArea* sa = p.sa_covered_v.back();
                if (sa->x[0] <= p.x[1] && p.x[1] <= sa->x[1] && sa->y[0] <= p.y[1] && p.y[1] <= sa->y[1]) {
                    current_sa = sa;
                }
        }
    
        int best_coords[2] = {-1, -1};
        double best = -1;
        if (current_sa->dir == orientation::horizontal) {
            for (int i = 0; i < n; i++) {
                if (best < cost[m-1][i]) {
                    best = cost[m-1][i];
                    best_coords[0] = m-1;
                    best_coords[1] = i;
                }
            }
        }  else {
            for (int i = m / 4 * 3; i < m; i++) {
                if (best < cost[i][n-1]) {
                    best = cost[i][n-1];
                    best_coords[0] = i;
                    best_coords[1] = n-1;
                }
            }
        }
        if (best > 0 + epsilon) {
            reduceDimensions(p, best_coords[0]*mesh.unit + p.y[0], best_coords[1]* mesh.unit + p.x[0]);
        } else {
            double new_x = (current_sa->dir == orientation::horizontal) ? p.x[1] : current_sa->x[0] - epsilon;
            double new_y = (current_sa->dir == orientation::horizontal) ? current_sa->y[0] - epsilon : p.y[1];
            reduceDimensions(p, new_y, new_x);
        }
        m = (int) ceil((p.y[1] - p.y[0])/mesh.unit) +1;
        n = (int) ceil((p.x[1] - p.x[0])/mesh.unit) +1;
    }
}

/* Main solving algorithm */
bool FacadeSolver::solve(void) {
    if (origins.empty()) {
        return true;
    } 
    Panel p;
    p.x[0] = origins.front().x;
    p.y[0]= origins.front().y;
    origins.erase(origins.begin());
    getEndPoint(p);
    reduceDimensions(p, p.y[1], p.x[1]);
    shrinkToMesh(p);
    while (1) 
    {
        adjustForFrames(p);
        if (!checkDimensions(p)) {
            origins.push_back({p.x[0], p.y[0]});
            return false;
        }
        // TODO: weight restrictions

        /* Manage origin points */
        vector<Point> removed;
        // remove origin points that may be covered by the addition of this panel
        vector<Point>::iterator it = origins.begin();
        while (it != origins.end()) {
            Point op = *it;
            if (p.overlap(op, mesh.unit)) {
                removed.push_back(op);
                it = origins.erase(it);
            } else {
                it++;
            }
        }

        // add new origin points (upper left corner of panel, bottom right corner)
        Point op1 = {p.x[1] + mesh.unit, p.y[0]} ;
        bool covered1 = false;
        Point op2 = {p.x[0] , p.y[1] + mesh.unit};
        bool covered2 = false;

        // make sure origin points are not at th edge
        if (op1.x >= w_fac || op1.y >= h_fac) covered1 = true;
        if (op2.x >= w_fac || op2.y >= h_fac) covered2 = true;

        // check if these new points are covered by existing solutions
        for(vector<Panel>::iterator it = solution.begin(); it!= solution.end() && (!covered1 || !covered2); it++) {
            Panel p = *it;
            if (!covered1 && p.overlap(op1, mesh.unit)) {
                covered1 = true;
            }
            if (!covered2 && p.overlap(op2, mesh.unit)) {
                covered2 = true;
            }
        }

        // don't add if they are dupes
        for(vector<Point>::iterator it = origins.begin(); it != origins.end() && (!covered1 || !covered2); it++) {
            Point op = *it;
            if (!covered1 && op == op1) covered1 = true;
            if (!covered2 && op == op2) covered2 = true;
        }

        if (!covered1) origins.push_back(op1);
        if (!covered2) origins.push_back(op2);
        sort(origins.begin(), origins.end());

        solution.push_back(p);
        if (show) drawState(500);
        p.print();
        bool next = solve();
        
        if (next) {
            return true;
        } else {
            vector<Panel>::iterator rm = remove(solution.begin(), solution.end(), p);
            solution.erase(rm);

            origins.insert(origins.end(), removed.begin(), removed.end());
            printf("ORIGIN: %f %f ", p.x[0], p.y[0]);
            // removed added origin points
            if (!covered1) {
                vector<Point>::iterator rm = remove(origins.begin(), origins.end(), op1);
                origins.erase(rm);
            } 
            if (!covered2) {
                vector<Point>::iterator rm = remove(origins.begin(), origins.end(), op2);
                origins.erase(rm);
            }
            // add back frames
            frames.insert(frames.end(), p.covered.begin(), p.covered.end());
            
            // TODO: shrink by some factor?
            if (p.current_sa->dir == orientation::horizontal && p.x[1] - p.x[0] > min_width + epsilon) {
                printf("CONDITION1\n");
                reduceDimensions(p, p.y[1], max(0.9*p.x[1] + 0.1* p.x[0], p.x[0] + min_width));
            } else if (p.current_sa->dir == orientation::vertical && p.y[1] - p.y[0] > min_height) {
                printf("CONDITION2\n");
                reduceDimensions(p, max(p.y[0] + max_height, 0.1*p.y[0] + 0.9*p.y[1]), p.x[1]);
            } else {
                printf("CONDITION3\n");
                reduceDimensions(p, p.current_sa->y[0] - epsilon, p.x[1]);
            }
        }
    }
};

void FacadeSolver::drawState(int ms) {
    // r.drawState(solution, frames, vertical, horizontal, w_fac, h_fac, ms);
}
