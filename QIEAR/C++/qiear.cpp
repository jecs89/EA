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

#include "qiear.h"

using namespace std;

void set_params( vector<double>& p){
	p[0] = 5; p[4] = 2; p[6] = 75; p[7] = 0.82; p[8] = 100;  p[10] = 0.66;
	p[5] = 0.5*p[8]; p[9] = 2;

	//ackley
	if( p[12] == 1){
		p[2] = -32;  p[3] = 32;  p[5] = 0.5*p[8]; p[9] = 1;
	}
	//rastr
	else if( p[12] == 2){
		p[2] = -5; 	 p[3] = 5;   p[5] = 0.5*p[8]; p[9] = 1;
	}
	//rosenbrock
	else if( p[12] == 3){
		p[2] = -100; p[3] = 100; p[5] = 0.5*p[8]; p[9] = 10;
	}
	//scheffel
	else if( p[12] == 4){
		p[2] = -100; p[3] = 100; p[5] = 0.5*p[8]; p[9] = 10;
	}
	//sphere
	else if( p[12] == 5){
		p[2] = -100; p[3] = 100; p[5] = 0.5*p[8]; p[9] = 1;
	}

}

void euclidean_distance( vector<double> X, double& dist , vector<double>& p ){

	double value = 0;
	if( p[11] == 1 ){
		if( p[12] == 1 )
			value = -20.0;
		else if( p[12] == 2 )
			value = -5.0;
		else if( p[12] == 3 )
			value = -8.0;
		else if( p[12] == 4 )
			value = 420.9687;
		else if( p[12] == 5 )
			value = -50.0;
	}

	cout << value << endl;
	vector<double> real_opt( p[1] , value );

	double d = 0.0;

	for( int i = 0 ; i < p[1] ; ++i ){
		d += ( real_opt[i] - X[i] ) * ( real_opt[i] - X[i] ) ;
	}

	dist = sqrt(d);

}

int NE = 0;

void qiear_thread( vector<double>& p ){

	time_t timer = time(0); 

	ofstream tiempos("time.txt",ofstream::app);
	ofstream vbest   ("vbest.txt",ofstream::app);
	ofstream result ("results.m",ofstream::app);
	ofstream measures("measures.txt",ofstream::app);

	// p[0] = 5;  p[6] = 75; p[7] = 0.82; p[8] = 100;  p[10] = 0.66;
	// p[5] = 0.5*p[8]; p[9] = 1;


	vector<int> T(5); T[0] = NE; T[1] = 50; T[2] = 100; T[3] = 500; T[4] = 1000;

	vector<double> avgxiter ( p[6] );
	vector<double> final_best (T[0]);

	vector< vector <double> > save_data( p[6] , vector<double>(4) );
	//(0) best ever, (1) best iteration , (2) median, (3) worst

	vector<double> save_iteration( T[0] ); //save iteration


	string name_vector = "q";
	string s_tmp = name_vector;       

	

	for( int k = 0 ; k < 1; ++k){

		cout << "QIEA-R: Fitness vs iterations " << T[k] << endl;
		vector<double> best_val ( T[k] ); //Saving best of every experiment
		vector< vector<double> > best_X ( T[k] , vector<double> ( p[1] )); //Saving X of the best

		//Dynamic use of Threads
		vector<thread> ths;
		// Launching threads
		for( int i = 0 ; i < T[k] ; ++i){
			ths.emplace_back( qiear, ref(best_X[i]), ref(best_val[i]) , ref(p), ref(save_data), ref(save_iteration[i]) );
		}

		//Joining threads
		for ( auto t = ths.begin(); t != ths.end(); ++t){ //int i = 0; i < ths.size(); i++ ){
			t->join();// ths[i].join();	
		}

		print( "%params:" , p , 5, result );
		print( "%params:" , p , 5, vbest );
		print( "%params:" , p , 5, tiempos );
		print( "%params:" , p , 5, measures );

		//Getting the best and std_dev
		pair<double,int> pbest = get_max(best_val,-1);
		// std_dev << p[11] << "\n";

		double avg_best = 0;
		//Getting the average
		for( int i = 0 ; i < best_val.size() ; ++i ){
			avg_best += best_val[i];
		}
		avg_best = avg_best / best_val.size();

		/*vector<double> temp_best_val ( T[k] ); //temporal
		temp_best_val = best_val;
		sort(temp_best_val.begin(), temp_best_val.end());

		//writing the best, avg + %, median
		measures << name_vector << " = [ " << pbest.first << ", \t" << avg_best << " ± " << standard_deviation(best_val) << ", \t" << ( temp_best_val[(int)T[k]/2] + temp_best_val[(int)T[k]/2 - 1 ] )/ 2.0 << "];\n";
		*/
	    
	    //write vector of the best
		vbest << name_vector << " = [ ";
		for( int i = 0 ; i < p[1] ; ++i){
			// result << std::fixed; setprecision(4) <<
			vbest << setw(15) << best_X[pbest.second][i];
		}
		vbest << "];\n";

		//write the result of experiments
		result << name_vector << " = [ ";
		for( int i = 0 ; i < T[k] ; ++i ){
			result << "\t" << best_val[i];
		}
		result << "];\n";


		//plotting save_data
		// ofstream plot_data("plot_data.m",ofstream::app);

		// plot_data << "%figure;\ndata = [ ";
		// for( int i = 0 ; i < p[6] ; ++i ){
		// 	plot_data << setprecision(10) << save_data[ i ][ 0 ] << " " << setprecision(10)<< save_data[ i ][ 1 ] << " " << setprecision(10)<< save_data[ i ][ 2 ] << ";\n";
		// }
		// plot_data << "];\n";
		// // plot_data << "idx=1:" << p[6] << ";\n";
		// // plot_data << "plot( idx , data(:, 1) , '*', 'linewidth', 3 , idx , data(:, 2) , '+', 'linewidth', 3 );\n";
		
		// plot_data << "figure;\nhold on;\n" ;
		// plot_data << "for i=1:"<< p[6] <<"\nplot( i , log2(data(i, 1)) , 'b*', 'linewidth', 2 , i , log2(data(i, 2)) , 'r+', 'linewidth', 2 , i , log2(data(i, 3)) , 'go', 'linewidth', 2 );\n";
		// plot_data << "title( '" << p[12] << "-" << p[11] << "');\n";
		// plot_data << "pause(0.005);\nend\n";

		// plot_data.close();

		////
		//getting iterations
		double avg_t = 0;
		//Getting the average
		for( int i = 0 ; i < save_iteration.size() ; ++i ){
			avg_t += save_iteration[i];
		}
		avg_t = avg_t / save_iteration.size();

		measures << name_vector << " = [ " << avg_best << " / " << max_value << " , " << avg_t << "]; \n";

		/*
		sort(save_iteration.begin(), save_iteration.end());

		//writing the best, avg + %, median
		measures << name_vector << " = [ " << save_iteration[0] << ", \t" << avg_best << " ± " << standard_deviation(save_iteration) << ", \t" << ( save_iteration[(int)T[0]/2] + save_iteration[(int)T[0]/2 - 1 ] )/ 2.0 << "];\n";
		
		//distance
		vector<double> best_vector(p[1]);
		best_vector = best_X[pbest.second];

		double distance = 0;
		euclidean_distance( best_vector, distance, p );

		cout << "distance: " << distance << endl;
		*/
		////
		
	}

	
	measures.close();
	vbest.close();
	result.close();

	time_t timer2 = time(0);
	// tiempos << p[11] << "\n"; 
	tiempos << name_vector << " = [ " <<  difftime(timer2, timer) << " ";
	tiempos << "];\n";
	tiempos.close();

}

int main(int argc, char *argv[]){

	vector<double> p(14);

	cout << "QIEA-R Model" << "\t" ;
	
	//Reading params
	p[1]  = atoi(argv[2]); //Dimension
	p[11] = atoi(argv[3]); //shifted
	p[12] = atoi(argv[1]); //Function
	p[13] = atoi(argv[4]); //Tupe of Q
	// cout << p[12] << endl;

	set_params(p);

	p[6] = atoi(argv[5]); //# Iter

	NE = atoi(argv[6]); //# Iter	

	global_opt = atoi(argv[7]);
    iter_max = p[6];
	
	qiear_thread( p );




	// for( int a = 5 ; a <= 30 ; a+=5 ){
	// 	for( int b = 75 ; b <= 200 ; b+=25 ){
	// 		for( double c = 0.1 ; c <= 1 ; c+=0.1 ){
	// 			for( int d = 125 ; d <= 200 ; d+=25 ){
	// 				for( int e = 1 ; e <= 10 ; e+=2 ){
	// 					for( double f = 0.1 ; f <= 1.0 ; f+=0.1 ){

	// 						p[0] = a;  p[6] = b; p[7] = c; p[8] = d;  p[10] = f;
	// 						p[5] = 0.5*p[8]; p[9] = e;

	// 						qiear_thread( p );
	// 					}
	// 				}
	// 			}
	// 		}
	// 	}
	// }

	return 0;
}
	