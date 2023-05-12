#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <random>

#define X1 (double)(-1-sqrt(5))/2
#define X2 (double)(-1+sqrt(5))/2
#define X0 X1 + 0.1
#define Z 0.0

double V(double X) {
    return (3*pow(X, 4) + 4*pow(X, 3) - 6*X*X)/12;
}

double Vprime(double X) {
    return X*(X*X - 1 + X);
}

double Vseconde(double X) {
    return 3*X*X + 2*X - 1;
}

double MC_naif(double x0, double h, double epsilon, int N) {
    // variables
    int E = 0;
    double sigma = sqrt(2*epsilon*h);
    double Xn = x0;


    std::random_device rd{};
    std::mt19937 gen{rd()};
    std::normal_distribution<> d{0, 1};

    for (int i = 0; i < N; i++) {
        Xn = x0;
        while (Xn >= X1 && Xn <= X2) {
            Xn = Xn - Vprime(Xn)*h + sigma*d(gen);
        }
        if (Xn > X2) {
            E++;
        }
    }

    return (double)E/N;
}

std::vector<double> simulation_trajectoire(double x0, double h, double sigma) {

    std::random_device rd{};
    std::mt19937 gen{rd()};
    std::normal_distribution<> d{0, 1};

    std::vector<double> X;

    X.push_back(x0);
    double Xn = x0;

    while (Xn > X1 && Xn < X2) {
        Xn = Xn - Vprime(Xn)*h + sigma*d(gen);
        X.push_back(Xn);
    }

    return X;
}

double K(std::vector<std::vector<double>> trajectoires) {
    // renvoie le plus petit maximum des trajectoires
    std::vector<double> maxs;
    for (auto trajectoire : trajectoires) {
        maxs.push_back(*std::max_element(trajectoire.begin(), trajectoire.end()));
    }
    return *std::min_element(maxs.begin(), maxs.end());
}

std::vector<unsigned long> I(std::vector<std::vector<double>> trajectoires, double K_value) {
    // renvoie les indices des trajectoires dont le maximum est K_value
    std::vector<unsigned long> indices;
    for (unsigned long i = 0; i < trajectoires.size(); i++) {
        if (*max_element(trajectoires[i].begin(), trajectoires[i].end()) == K_value) {
            indices.push_back(i);
        }
    }
    return indices;
}

double estimateur_AMS(double x0, double h, double epsilon, int M) {
    // variables
    const double sigma = sqrt(2*epsilon*h);
    std::vector<std::vector<double>> trajectoires;
    double K_q = 0.0;
    std::vector<unsigned long> I_q;
    double p_estim = 1.0;

    std::vector<double> first_half;
    std::vector<double> second_half;

    // simulation initiale
    for (int i = 0; i < M; i++) {
        trajectoires.push_back(simulation_trajectoire(x0, h, sigma));
    }
    K_q = K(trajectoires);
    I_q = I(trajectoires, K_q);
    p_estim *= 1 - static_cast<double>(I_q.size()) / M;

    while (K_q < X2) {
        // on selectionne les indices entre 0 et M-1 non presents dans I_q
        std::vector<unsigned long> indices;
        for (int i = 0; i < M; i++) {
            if (std::find(I_q.begin(), I_q.end(), i) == I_q.end()) {
                indices.push_back(i);
            }
        }
        for (auto m : I_q) {
            // on prend un indice au hasard dans indices
            unsigned long p = indices[std::rand() % indices.size()];
            // on trouve le premier instant ou la trajectoire p depasse K_q
            int t = 0;
            while (trajectoires[p][t] <= K_q) t++;
            // on selectionne la trajectoire p jusqu'a l'instant t
            first_half.clear();
            first_half.insert(first_half.end(), trajectoires[p].begin(), trajectoires[p].begin() + t);
            // on simule la trajectoire m a partir de l'instant t
            second_half.clear();
            second_half = simulation_trajectoire(trajectoires[p][t], h, sigma);
            // on remplace la trajectoire m par la concatenation des deux moitiees
            trajectoires[m].clear();    
            trajectoires[m].insert(trajectoires[m].end(), first_half.begin(), first_half.end());
            trajectoires[m].insert(trajectoires[m].end(), second_half.begin(), second_half.end());
        }

        K_q = K(trajectoires);
        I_q = I(trajectoires, K_q);
        p_estim *= 1 - static_cast<double>(I_q.size()) / M;
    }

    return p_estim;
}

double T12_naif(double h, double epsilon, int N) {
    // variables
    double time = 0.0;
    double sigma = sqrt(2*epsilon*h);
    double Xn = X1;

    std::random_device rd{};
    std::mt19937 gen{rd()};
    std::normal_distribution<> d{0, 1};

    for (int i = 0; i < N; i++) {
        Xn = X1;
        while (Xn < X2) {
            Xn = Xn - Vprime(Xn)*h + sigma*d(gen);
            time += h;
        }
    }

    return time/N;
}

double approx_T12(double epsilon) {
    return 2*M_PI*exp((V(Z) - V(X1))/epsilon)/sqrt(-Vseconde(Z)*Vseconde(X1));
}

struct T02_result {
    double T02;
    double p_estim;
};

double T02_naif(double h, double epsilon, int N) {
    // variables
    int global_time = 0;
    int local_time = 0;
    double sigma = sqrt(2*epsilon*h);
    double Xn = X0;

    std::random_device rd{};
    std::mt19937 gen{rd()};
    std::normal_distribution<> d{0, 1};

    for (int i = 0; i < N; i++) {
        Xn = X0;
        while (Xn < X2) {
            Xn = Xn - Vprime(Xn)*h + sigma*d(gen);
            local_time++;
            if (Xn < X1) {
                local_time = 0;
                Xn = X0;
            }
        }
        global_time += local_time;
        local_time = 0;
    }

    return (double)global_time/N;
}

T02_result T02_AMS(double h, double epsilon, int M) {
    // variables
    const double sigma = sqrt(2*epsilon*h);
    std::vector<std::vector<double>> trajectoires;
    double K_q = 0.0;
    std::vector<unsigned long> I_q;
    double p_estim = 1.0;
    std::vector<double> longueurs;

    std::vector<double> first_half;
    std::vector<double> second_half;

    // simulation initiale
    for (int i = 0; i < M; i++) {
        trajectoires.push_back(simulation_trajectoire(X0, h, sigma));
    }
    K_q = K(trajectoires);
    I_q = I(trajectoires, K_q);
    p_estim *= 1 - static_cast<double>(I_q.size()) / M;

    while (K_q < X2) {
        // on selectionne les indices entre 0 et M-1 non presents dans I_q
        std::vector<unsigned long> indices;
        for (int i = 0; i < M; i++) {
            if (std::find(I_q.begin(), I_q.end(), i) == I_q.end()) {
                indices.push_back(i);
            }
        }
        for (auto m : I_q) {
            // on prend un indice au hasard dans indices
            unsigned long p = indices[std::rand() % indices.size()];
            // on trouve le premier instant ou la trajectoire p depasse K_q
            int t = 0;
            while (trajectoires[p][t] < K_q) t++;
            // on selectionne la trajectoire p jusqu'a l'instant t
            first_half.clear();
            first_half.insert(first_half.end(), trajectoires[p].begin(), trajectoires[p].begin() + t);
            // on simule la trajectoire m a partir de l'instant t
            second_half.clear();
            second_half = simulation_trajectoire(trajectoires[p][t], h, sigma);
            // on remplace la trajectoire m par la concatenation des deux moitiees
            trajectoires[m].clear();    
            trajectoires[m].insert(trajectoires[m].end(), first_half.begin(), first_half.end());
            trajectoires[m].insert(trajectoires[m].end(), second_half.begin(), second_half.end());
        }

        K_q = K(trajectoires);
        I_q = I(trajectoires, K_q);
        p_estim *= 1 - static_cast<double>(I_q.size()) / M;
    }

    for (auto t : trajectoires) {
        longueurs.push_back(t.size()-1);
    }
    // calcul de la moyenne des longueurs
    double moyenne = std::accumulate(longueurs.begin(), longueurs.end(), 0.0) / longueurs.size();

    T02_result result;
    result.T02 = moyenne;
    result.p_estim = p_estim;

    return result;
}

int main() {
    std::cout << approx_T12(0.1) << std::endl;
    std::cout << T12_naif(0.001, 0.1, 1000000) << std::endl;
    return 0;
}
