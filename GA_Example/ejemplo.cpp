#include <bits/stdc++.h>

#define CRUCE_1PUNTO 1
#define CRUCE_2PUNTOS 2
#define CRUCE_UNIFORME 3

using namespace std;

class Cromosoma{
public:
    string str;
    double fitness;
};
bool cmp(Cromosoma c1, Cromosoma c2){
    return c1.fitness < c2.fitness;
}

class Variable{
public:
    double upperBound;
    double lowerBound;
    int digit;
};

double randomDouble(double a, double b) {
    double random = ((double) rand()) / (double) RAND_MAX;
    double diff = b - a;
    double r = random * diff;
    return a + r;
}

double f(vector<double> x);

void cruce1punto(string &str1, string &str2);
void cruce2puntos(string &str1, string &str2);
void cruceUniforme(string &str1, string &str2);
void mutacion1bit(string &str);

void inicializa(vector<Cromosoma> &Gt, vector<Variable> v, int populationSize);
void evaluar(vector<Cromosoma> &Gt, vector<Variable> v);
void seleccionar(vector<Cromosoma> &Gt);
void operar(vector<Cromosoma> &Gt, double pCruce, double pMutacion, int tipoCruce);

string codificar(double x, Variable v);
string codificar(vector<double> x, vector<Variable> v);
double decodificar(string str, Variable v);
vector<double> decodificar(string str, vector<Variable> v);

int main(int argc, char* argv[]){
    if(argc != 6){
        printf("Parametros:\n");
        printf("  numero_poblacion_inicial: entero positivo\n");
        printf("  cantidad_generaciones: entero positivo\n");
        printf("  tipo_cruce: \n");
        printf("              1: 1 punto\n");
        printf("              2: 2 puntos\n");
        printf("              3: uniforme \n");
        printf("  probabilidad_cruce: real [0, 1]\n");
        printf("  probabilidad_mutacion: real [0, 1]\n");
        return -1;
    }

    vector<Variable> v;
    Variable v1;
    v1.digit = 7;
    v1.lowerBound = -10;
    v1.upperBound = 10;
    v.push_back(v1);
    v.push_back(v1);

    int t;
    vector<Cromosoma> Gt;

    int populationSize = atoi(argv[1]);
    int MaxGen = atoi(argv[2]);
    int tipoCruce = atoi(argv[3]);
    double pCruce = atof(argv[4]);
    double pMutacion = atof(argv[5]);

    srand(time(NULL));

    vector<double> bestX;
    vector<double> meanX;



    /** Inicio algoritmo genetico **/

    t = 0;
    inicializa(Gt, v, populationSize);
    evaluar(Gt, v);
    sort(Gt.begin(), Gt.end(), cmp);
    while(t < MaxGen){
        t = t + 1;
        Cromosoma gBest = Gt[Gt.size() - 1];
        //// INICIO Almacena los resultados para posterior analisis
        double meanFitness = 0;
        for(int i = 0; i < Gt.size(); i++){
            meanFitness += Gt[i].fitness;
        }
        bestX.push_back(gBest.fitness);
        meanX.push_back(meanFitness / Gt.size());
        //// FIN Almacena los resultados para posterior analisis
        seleccionar(Gt);
        operar(Gt, pCruce, pMutacion, tipoCruce);
        Gt.push_back(gBest);
        evaluar(Gt, v);
        sort(Gt.begin(), Gt.end(), cmp);
    }

    /** Fin algoritmo genetico **/


    /* Imprime los resultados */
    printf("gen  best     voff     von\n");
    for(int i = 0; i < (int)bestX.size(); i++){
        double voff = 0;
        for(int j = 0; j <= i; j++){
            voff += bestX[j];
        }
        voff = voff / (i + 1.0);

        double von = 0;
        for(int j = 0; j <= i; j++){
            von += meanX[j];
        }
        von = von / (i + 1.0);

        printf("%4d,%.6lf,%.6lf,%.6lf;\n", i, bestX[i], voff, von);

    }
    vector<double> x = decodificar(Gt[Gt.size() - 1].str, v);
    printf("Solucion: (%lf, %lf)\n", x[0], x[1]);



    return 0;
}




/******************/
/** La funcion f **/
/******************/
double f(vector<double> x){
    double x1 = x[0];
    double x2 = x[1];
    double aux1 = (sin(sqrt(x1 * x1 + x2 * x2))) * (sin(sqrt(x1 * x1 + x2 * x2))) - 0.5;
    double aux2 = (1.0 + 0.001 * (x1 * x1 + x2 * x2)) * (1.0 + 0.001 * (x1 * x1 + x2 * x2));
    return 0.5 - aux1 / aux2;
}




/****************/
/** Inicializa **/
/****************/
void inicializa(vector<Cromosoma> &Gt, vector<Variable> v, int populationSize){
    for(int k = 0; k < populationSize; k++){
        string str = "";
        for(int i = 0; i < v.size(); i++){
            double x = randomDouble(v[i].upperBound, v[i].lowerBound);
            str = str + codificar(x, v[i]);
        }

        Cromosoma g;
        g.str = str;

        Gt.push_back(g);
    }
}




/*************/
/** Evaluar **/
/*************/
void evaluar(vector<Cromosoma> &Gt, vector<Variable> v){
    for(int i = 0; i < Gt.size(); i++){
        vector<double> valores = decodificar(Gt[i].str, v);
        Gt[i].fitness = f(valores);
    }
}




/*****************/
/** Seleccionar **/
/*****************/

void seleccionar(vector<Cromosoma> &Gt){
    vector<Cromosoma> nuevoGt;

    double total = 0;
    for(int i = 0; i < Gt.size(); i++){
        total += Gt[i].fitness;
    }
    for(int k = 0; k < Gt.size() - 1; k++){
        double r = randomDouble(0, total);
        double sum = 0;
        int i = 0;
        while(sum < r && i < Gt.size()){
            sum += Gt[i].fitness;
            i++;
        }
        nuevoGt.push_back(Gt[i - 1]);
    }

    Gt = nuevoGt;
}




/************/
/** Operar **/
/************/

void operar(vector<Cromosoma> &Gt, double pCruce, double pMutacion, int tipoCruce){
    for(int i = 0; i < Gt.size() && i + 1 < Gt.size(); i+= 2){
        int pareja = rand() % Gt.size();
        swap(Gt[pareja], Gt[i + 1]);
        if(randomDouble(0, 1) < pCruce){
            switch(tipoCruce){
                case CRUCE_1PUNTO:
                    cruce1punto(Gt[i].str, Gt[i + 1].str);
                    break;
                case CRUCE_2PUNTOS:
                    cruce2puntos(Gt[i].str, Gt[i + 1].str);
                    break;
                case CRUCE_UNIFORME:
                    cruceUniforme(Gt[i].str, Gt[i + 1].str);
                    break;
            }
        }
   }
    for(int i = 0; i < Gt.size(); i++){
        if(randomDouble(0, 1) < pMutacion){
            mutacion1bit(Gt[i].str);
        }
    }
}




/****************/
/** Operadores **/
/****************/

void cruce1punto(string &str1, string &str2){
    int punto = rand() % str1.size();
    for(int i = punto; i < str2.length(); i++){
        swap(str1[i], str2[i]);
    }
}

void cruce2puntos(string &str1, string &str2){
    int punto1 = rand() % str1.size();
    int punto2 = rand() % str1.size();
    if(punto1 > punto2){
        swap(punto1, punto2);
    }
    for(int i = punto1; i < punto2; i++){
        swap(str1[i], str2[i]);
    }
}

void cruceUniforme(string &str1, string &str2){
    string str;
    for(int i = 0; i < str1.length(); i++){
        if(rand() % 2){
            swap(str1[i], str2[i]);
        }
    }
}

void mutacion1bit(string &str){
    int i = rand() % str.length();
    str[i] = (str[i] == '0') ? '1' : '0';
}




/******************/
/** Codificacion **/
/******************/

string codificar(double x, Variable v){
    double precision = pow(10.0, -v.digit);
    double delta = pow(2, ceil(log2((v.upperBound - v.lowerBound) / precision))) / (v.upperBound - v.lowerBound);
    int bits = (x - v.lowerBound) * delta;
    int iteration = ceil(log2((v.upperBound - v.lowerBound) * delta));
    string str;
    for(int i = 0; i < iteration; i++){
        str.push_back('0' + (bits & 1));
        bits = bits >> 1;
    }
    return str;
}

double decodificar(string str, Variable v){
    double precision = pow(10.0, -v.digit);
    double delta = pow(2, ceil(log2((v.upperBound - v.lowerBound) / precision))) / (v.upperBound - v.lowerBound);
    int bits = 0;
    for(int i = str.length() - 1; i >= 0; i--){
        bits = bits | (str[i] == '1' ? 1 : 0);
        bits = bits << 1;
    }
    return bits / delta + v.lowerBound;
}

string codificar(vector<double> x, vector<Variable> v){
    string str;
    for(int i = 0; i < v.size(); i++){
        str = str + codificar(x[i], v[i]);
    }
    return str;
}

vector<double> decodificar(string str, vector<Variable> v){
    vector<double> x;
    int i0 = 0;
    int i1 = 0;
    for(int i = 0; i < v.size(); i++){
        double precision = pow(10.0, -v[i].digit);
        double delta = pow(2, ceil(log2((v[i].upperBound - v[i].lowerBound) / precision))) / (v[i].upperBound - v[i].lowerBound);
        int length = ceil(log2((v[i].upperBound - v[i].lowerBound) * delta));
        string substring = str.substr(i0, length);
        x.push_back(decodificar(substring, v[i]));
        i0 = i0 + length;
    }
    return x;
}
