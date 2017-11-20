// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// %% QIEAR
// %% http://ieeexplore.ieee.org/document/5586193/
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// %% Developed by jecs89
// %% website: jecs89.site
// %% github: https://github.com/jecs89/EA
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// %% DON'T FORGET THE CREDITS OF THE CODE
// %% IF YOU ARE INTERESTED TO PERFORM EXPERIMENTS
// %% WRITE ME AND WE CAN WORK TOGETHER
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// A own modification was used in update step 
// because QIEAR could increase indefinitely.

#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <math.h>
#include <random>
#include "iomanip"
#include <sstream>
#include <algorithm>
#include <thread>

using namespace std;

const double PI = atan(1.0)*4;
const double my_inf = powf(2,64) - 1;

double global_opt; //10^-3

double max_value;

int iter_max;

//Twister
typedef std::mt19937 MyRNG;  // the Mersenne Twister with a popular choice of parameters
uint32_t seed_val;           // populate somehow

MyRNG rng;                   // e.g. keep one global instance (per thread)

void initialize(){
  rng.seed(seed_val);
}
//Twister

template<typename T>
void print( string title, vector< vector<T> >& Q_t ){
    cout << title << endl ;
    int spc = 15;
    for( int i = 0 ; i < Q_t.size(); ++i){
        for( int j = 0 ; j < Q_t[0].size(); ++j){
            cout << setw(spc) << Q_t[i][j] ;
        }
        cout << endl;
    }
}
template<typename T>
void print( string title, vector<T> & Q_t ){
    cout << title << endl ;
    int spc = 15;
    for( int i = 0 ; i < Q_t.size(); ++i){
        cout << setw(spc) << i << setw(spc) << Q_t[i] << endl;
    }
}
    
pair<double,int> get_max(vector<double>& F, int type){

	double n_max = F[0]; int idx = 0;

	for( int i = 1 ; i < F.size() ; ++i){
		if( type == 1 ) {
			if( n_max < F[i] ){
				n_max = F[i];
				idx = i;
			}
		}
		else if( type == -1 ){
			if( n_max > F[i] ){
				n_max = F[i];
				idx = i;
			}
		}
	}

	

	pair<double,int> ans;
	ans.first = n_max; ans.second = idx;

	return ans;
}


double standard_deviation(vector<double> data){
    double mean=0.0, sum_deviation=0.0;
    
    int i;
    for(i=0; i<data.size();++i){
        mean+=data[i];
    }
    mean=mean/data.size();

    for(i=0; i<data.size();++i){
   		sum_deviation += (data[i]-mean)*(data[i]-mean);
    }

    return sqrtf(sum_deviation/data.size());
}

void write_vector2csvfile(vector<double>& best_F, ofstream& file, vector<double>& p){
	string name_vector = "qc_";
    string s_tmp = static_cast<ostringstream*>( &(ostringstream() << p[0]) )->str();	name_vector = name_vector + s_tmp + "_";       
    s_tmp = static_cast<ostringstream*>( &(ostringstream() << p[1]) )->str();	name_vector = name_vector + s_tmp + "_";       
    s_tmp = static_cast<ostringstream*>( &(ostringstream() << p[4]) )->str();	name_vector = name_vector + s_tmp + "_";       
    s_tmp = static_cast<ostringstream*>( &(ostringstream() << p[6]) )->str();	name_vector = name_vector + s_tmp + "_";       
    s_tmp = static_cast<ostringstream*>( &(ostringstream() << p[8]) )->str();	name_vector = name_vector + s_tmp + "_";       
    s_tmp = static_cast<ostringstream*>( &(ostringstream() << p[9]) )->str();	name_vector = name_vector + s_tmp + "_";  
    s_tmp = static_cast<ostringstream*>( &(ostringstream() << p[10]) )->str();	name_vector = name_vector + s_tmp ;      

	file << name_vector ;
	for(int i = 0 ; i < best_F.size() ; ++i){
		file << "\t," << best_F[i];
	}
	file << "\n";
}

vector<vector<double>> domain_limits(30,vector<double>(2));

void set_domain(){
    for( int i = 0; i < domain_limits.size(); i++){
        domain_limits[i][0] = -30;
        domain_limits[i][1] = 30;

    }
}


void Q_generation( vector<double>& p, vector< vector<double> >& Q_t){
	vector< vector<double> > Q_actual ( p[0], vector<double>(p[1]*2) );

    //Same domain

    if( p[13] == 0 ){

        for( int d = 0; d < domain_limits.size(); d++){
            p[3] = domain_limits[d][1];
            p[2] = domain_limits[d][0];

            double pulso_ancho = ( p[3] - p[2] ) / 1;

            for(int i = 0; i < p[0]; ++i){
                for( int j = 0; j < p[1]*2; ++j){

                     if( j % 2 != 0){ //sigma
                         Q_actual[i][j] = pulso_ancho / 2.0 ;
                     }
                     else{ //mu
                         Q_actual[i][j] = - p[3] + pulso_ancho/2.0;
                     }           
                }
            }

        }

    }

    else if( p[13] == 1 ){
        for( int d = 0; d < domain_limits.size(); d++){
            p[3] = domain_limits[d][1];
            p[2] = domain_limits[d][0];
            
            // Domain splitted
        	double pulso_ancho = ( p[3] - p[2] ) / p[0];

        	for(int i = 0; i < p[0]; ++i){
        		if (i == 0){
        			Q_actual[i][0] = -p[3] + pulso_ancho/2 * (i+1); 
        		}
        		else{
        			Q_actual[i][0] = Q_actual[i-1][0] + pulso_ancho;
        		}

        		for( int j = 0; j < p[1]*2; ++j){
        			if( j % 2 != 0){
        				Q_actual[i][j] = pulso_ancho;
        			}
        			else{
        				Q_actual[i][j] = Q_actual[i][0];
        			}			
        		}
        	}

        }
    }

	Q_t = Q_actual;
}

void C_generation( vector< vector<double> >& C_t, vector<double>& p, vector< vector<double> >& Q_t){

	//default_random_engine rng( random_device{}() ); 
	// mt19937 rng(time(0));		
	uniform_real_distribution<double> dist( 0, p[0] ); 

	// default_random_engine rng2( random_device{}() ); 		
	// mt19937 rng2(time(0));
	uniform_real_distribution<double> dist2( 0, 1 ); 

    // Picking random q individuals
	// for(int i = 0; i < p[8]; ++i){
	// 	int idx = floor( dist(rng) );
	// 	indexes_t[i] = idx;
	// 	// cout << idx << endl;
	// 	for(int j = 0; j < p[1]; ++j){
	// 		C_t[i][j] = dist2(rng2) * Q_t[idx][2*j+1] + ( Q_t[idx][0] - Q_t[idx][2*j+1]/2.0 );
	// 	}
	// }

    // Picking q ~ classical
    int nq = int(p[8] / p[0]);

    int q = 0;

    for(int i = 0; i < p[8]; i+=nq ){
        // int idx = floor( dist(rng) );
        // indexes_t[i] = idx;
        // cout << idx << endl;
        for( int k = 0 ; k < nq ; ++k ){
            for(int j = 0; j < p[1]; ++j){
                C_t[i+k][j] = dist2(rng) * Q_t[q][2*j+1] + ( Q_t[q][0] - Q_t[q][2*j+1]/2.0 );
            }
        }
        q++;
    }

}


void Eval(vector<double>& F, vector< vector<double> >& C_t, vector<double>& p){

    int function = p[12];
    int shifted = p[11];
    
    double bias;

    if( function == 1 ){ //Ackley 0
        bias = ( shifted == 1 ) ? -20.0 : 0.0 ;    
        double a = 20, b = 0.2;
        for(int i = 0 ; i < C_t.size() ; ++i){
            double exp_1 = 0, exp_2 = 0;
            for(int j = 0 ; j < C_t[0].size() ; ++j){
                exp_1 += ( C_t[i][j] - bias ) * ( C_t[i][j] - bias );
                exp_2 += cos( 2.0*PI*( C_t[i][j] - bias ) );
            }    
        F[i] = -a * exp(-b*sqrtf(1.0/C_t[0].size() * exp_1)) - exp(1.0/C_t[0].size()*exp_2) + a + exp(1);
        }    
    }    
    else if( function == 2 ){ //Rastr 0
        bias = ( shifted == 1 ) ? -5.0 : 0.0 ;
        for(int i = 0 ; i < C_t.size() ; ++i){
            double sum = 0;
            double exp_1 = 0, exp_2 = 0;
            for(int j = 0 ; j < C_t[0].size() ; ++j){
                exp_1 = ( C_t[i][j] - bias )*( C_t[i][j] - bias );
                exp_2 = cos( 2.0*PI*( C_t[i][j] - bias ) );
                sum += exp_1 - 10*exp_2;
            }    
            F[i] = 10*C_t[0].size() + sum;
        }    
    }    
    else if( function == 3 ){ //REVISAR, para más variables //Rosenbrock 1
        bias = ( shifted == 1 ) ? -8.0 : 0.0 ;
        for(int i = 0 ; i < C_t.size() ; ++i){
            double sum = 0;
            double exp_1 = 0, exp_2 = 0;
            for(int j = 0 ; j < C_t[0].size() - 1; ++j){
                exp_1 = ( C_t[i][j+1] - bias ) - powf(( C_t[i][j] - bias ),2);
                exp_2 = powf(( C_t[i][j] - bias )-1.0,2);
                sum += 100*powf(exp_1,2) + exp_2;
                // cout << sum << endl;
            }    
            F[i] = sum;
        }    
    }    
    else if( function == 4 ){ //Scheffel 420.9687
        bias = ( shifted == 1 ) ? -420 : 0.0 ;
        for(int i = 0 ; i < C_t.size() ; ++i){
            double sum = 0;
            double exp_1 = 0;
            for(int j = 0 ; j < C_t[0].size() ; ++j){
                exp_1 = (C_t[i][j]-bias) * sin(sqrtf(abs(C_t[i][j]-bias)));
                sum += exp_1;
            }    
            F[i] = 418.9829*C_t[0].size() - sum;
        }    
    }    
    else if( function == 5 ){ //Sphere 0
        bias = ( shifted == 1 ) ? -50.0 : 0.0 ;
        for(int i = 0 ; i < C_t.size() ; ++i){
            double sum = 0;
            double exp_1 = 0;
            for(int j = 0 ; j < C_t[0].size() ; ++j){
                exp_1 = powf(( C_t[i][j] - bias ),2);
                sum += exp_1;
            }    
            F[i] = sum ;
        }    
    }    
}  

struct my_double {
    vector<double> x;
    double prob;
    bool operator<( const my_double& val ) const {
        return prob < val.prob;
    }
};

struct Xgreater{
    bool operator()( const my_double& lx, const my_double& rx ) const {
        return lx.prob < rx.prob;
    }
};

void qiear( vector<double>& best_X, double& best_val, vector<double>& p, vector< vector<double> > & save_data, double& save_iteration ){

	vector< vector<double> > Q_t ( p[0], vector<double>(p[1]*2) );	//Quantum population
	vector<double> f ( p[8] );		//Evaluation
    vector<double> f_t_1 ( p[8] );      //Evaluation of t-1 iteration

	//Classical population
	vector< vector<double> > C_t ( p[8], vector<double>(p[1]) );	
	vector< vector<double> > E_t ( p[8], vector<double>(p[1]) );

    //Temporal for resizing
    vector<double> res_f(p[8]) ;      //Evaluation
    vector<double> res_f_t_1(p[8]);      //Evaluation of t-1 iteration

    vector< vector<double> > res_C_t( p[8], vector<double>(p[1]) ) ;
    vector< vector<double> > res_E_t( p[8], vector<double>(p[1]) ) ;

	//First iteration
	Q_generation(p,Q_t); //Creating Q population

    double best_global = my_inf;

    int t = 0;

    int size_pop = p[8];

    vector<my_double> C_td;
    vector<my_double> E_td;

	while( t < iter_max ){
        // cout << "Iter: " << t << endl;
		
		C_generation( E_t, p, Q_t );

        // print( "C_generation:" , E_t );

		if( t == 0 ){
			C_t = E_t;
			
			Eval(f,C_t, p);

            f_t_1 = f; //copy

            max_value = get_max( f , 1 ).first;

            // print( "E_t:" , E_t );
            // print( "f:" , f );

            // print( "C_t:" , C_t );
            // print( "f t-1:" , f_t_1 );
		}
		else{
			//Recombination using DE
			uniform_real_distribution<double> distF( 0, 1 ); 	//To get new individual
			uniform_int_distribution<int> distC( 0, p[8]-1 ); 	//To get index of classical individuals
			uniform_real_distribution<double> dist_crossover( 0, 1 );	//Get prob of crossover
			
			double F = distF(rng), rate_crossover = 0, prob_population = 0; 
			int a , b , c;	//Coefficients
			vector<double> x_c1(p[1]);	//To save mutation
            vector<double> x_c2(p[1]);   //To save mutation

			for( int i = 0 ; i < p[8] ; ++i){
				rate_crossover = dist_crossover(rng);
				prob_population = dist_crossover(rng);
				do{
					a = distC(rng);
				}while( a == i );
				do{
					b = distC(rng);
				}while( b == i || b == a );

				for( int j = 0 ; j < p[1] ; ++j){
                    x_c1[j] = prob_population * E_t[a][j] + ( 1 - prob_population ) * C_t[b][j];
                    x_c2[j] = prob_population * C_t[b][j] + ( 1 - prob_population ) * E_t[a][j];
				}

				//Crossover
				if( rate_crossover <= p[10] ){
					E_t[a] = x_c1;
                    E_t[b] = x_c2;
				}
			}

            Eval( f_t_1 , C_t , p ); //Saving previous C_t

            // print( "C_t:" , C_t );
            // print( "f t-1:" , f_t_1 );

            Eval(f,E_t, p); //Evaluation of E_t
            C_t = E_t;

            // print( "E_t:" , E_t );
            // print( "f:" , f );
		}

        if( t > 0.75*p[6] ){
            //Resizing
            C_td.resize(p[8]);
            E_td.resize(p[8]);

            for( int i = 0 ; i < p[8] ; i++ ){
                C_td[i].x.resize(p[1]);
                E_td[i].x.resize(p[1]);

                for( int j = 0 ; j < p[1] ; j++ ){
                    C_td[i].x[j] = C_t[i][j];
                    C_td[i].prob = f_t_1[i];

                    E_td[i].x[j] = E_t[i][j];
                    E_td[i].prob = f[i];
                }
            }

            sort( C_td.begin(), C_td.end() , Xgreater() );
            // for( int i = 0 ; i < p[8] ; i++ ){
            //     for( int j = 0 ; j < p[1] ; j++ ){
            //         cout << setw(15) << C_td[i].x[j] ;
            //     }
            //     cout << endl;
            // }
            sort( E_td.begin(), E_td.end() , Xgreater() );
        }



        //Updating
        //Rule 1/5th
        uniform_real_distribution<double> distF( 0, 1 );
        double prop_mod = distF(rng);
        
        if( prop_mod < 0.5 ){

            int counter = 0;

            //Comparison between(C_t) and best_global of last pop
            for( int i = 0 ; i < p[8] ; ++i ){
                    if( f[i] < best_global ){
                        counter++;
                    }
            }

            double my_phi = double(counter) / ( p[8] )    ;

            double factor_mult = 0.0;

            if( my_phi < 0.2 ){
                factor_mult = p[7];

                size_pop = size_pop * 0.9;
            }
            else if( my_phi > 0.2 ){
                factor_mult = 1.0 / p[7];

                size_pop = size_pop * 1.1;
            }
            else if( my_phi == 0.2 ){
                factor_mult = 1.0;
            }

            for( int i = 0 ; i < Q_t.size() ; ++i ){
                for( int j = 1 ; j < Q_t[0].size() ; j+=2){

                    Q_t[i][j] = Q_t[i][j]*factor_mult;
                }
            }
        }

		
		if( t > 0 && (t+1)%(int)p[9] == 0 ){
            //Moving square pulse
            double lambda = 0.5;

            f_t_1 = f;

            for (int i = 0; i < p[0]; ++i){
                
                pair<double,int> ans =  get_max(f_t_1,-1) ;

                for (int j = 0, k = 0; j < p[1]*2; j+=2, ++k){
                    Q_t[i][j] = Q_t[i][j] + lambda * ( C_t[ ans.second ][k] - Q_t[i][j] ) ;
                }

                f_t_1[ ans.second ] = my_inf;
            }

        }

        Eval( f, C_t, p); //Evaluation of last C_t

        pair<double,int> ans =  get_max(f,-1) ;

        best_global = ( ans.first < best_global ) ? ans.first : best_global;

        // if( best_global < global_opt )
        //     break;

        if( t > 0.75*p[6] ){
            p[8] = ( size_pop < 10 ) ? 10 : size_pop;
            C_t.resize(p[8]);
            E_t.resize(p[8]);

            f.resize(p[8]);
            f_t_1.resize(p[8]);

            for( int i = 0 ; i < p[8] ; i++ ){
                C_t[i].resize( p[1] );
                E_t[i].resize( p[1] );
            }

            for( int i = 0 ; i < p[8] ; i++ ){
                for( int j = 0 ; j < p[1] ; j++ ){
                    C_t[i][j] = C_td[i].x[j];
                    E_t[i][j] = E_td[i].x[j];
                }
            }

            // Eval( f_t_1 , C_t , p );
            Eval( f , E_t , p );
        }


        // cout << size_pop << endl;
    
        t++;
	}

    save_iteration = t;

    Eval( f, C_t, p); //Evaluation of last C_t

    pair<double,int> ans =  get_max(f,-1) ;
		
	best_val = ans.first;
	best_X   = C_t[ans.second];
}


void write_vector2file(vector<double>& best_F, ofstream& file, vector<double>& p){
	string name_vector = "q";
    string s_tmp = static_cast<ostringstream*>( &(ostringstream() << p[0]) )->str();	name_vector = name_vector + s_tmp ;       
    // s_tmp = static_cast<ostringstream*>( &(ostringstream() << p[1]) )->str();	name_vector = name_vector + s_tmp + "_";       
    // s_tmp = static_cast<ostringstream*>( &(ostringstream() << p[4]) )->str();	name_vector = name_vector + s_tmp + "_";       
    // s_tmp = static_cast<ostringstream*>( &(ostringstream() << p[6]) )->str();	name_vector = name_vector + s_tmp + "_";       
    // s_tmp = static_cast<ostringstream*>( &(ostringstream() << p[8]) )->str();	name_vector = name_vector + s_tmp + "_";       
    // s_tmp = static_cast<ostringstream*>( &(ostringstream() << p[9]) )->str();	name_vector = name_vector + s_tmp + "_";  
    // s_tmp = static_cast<ostringstream*>( &(ostringstream() << p[10]) )->str();	name_vector = name_vector + s_tmp ;      

	file << name_vector << " = [";
	for(int i = 0 ; i < best_F.size() ; ++i){
		file << "\t" << best_F[i];
		
	}
	file << " ];\n";
}

template<typename T>
void print( string title, vector<T> & Q_t, int spc, ofstream & os ){
    os << title ;
    // os << "% ";
    // int spc = 10;
    for( int i = 0 ; i < Q_t.size(); ++i){
        // cout << setw(spc) << i << setw(spc) << Q_t[i] << endl;
        os << setw(spc) << Q_t[i] ;
    }    
    os << endl;
    
}

