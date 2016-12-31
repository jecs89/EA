%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% QIEAR
%% http://ieeexplore.ieee.org/document/5586193/
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Developed by jecs89
%% website: jecs89.site
%% github: https://github.com/jecs89/EA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DON'T FORGET THE CREDITS OF THE CODE
%% IF YOU ARE INTERESTED TO PERFORM EXPERIMENTS
%% WRITE ME AND WE CAN WORK TOGETHER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ bestV bestF ] = qiear()
	tic;
	figure; hold on;
	% printf( 'Quantum Algorithm with Real Codification\n')
	% t=cputime;
	p = params();

	Q_t = Q_Generation( p );
	myprint('Q_t',Q_t);

	% f_t = zeros(p(1),1);
	% f_t1 = zeros(p(1),1);

	%%Best
	% best_f = zeros(p(7),1);
	% best_values = zeros(p(7),p(2));
 
 	if( p(5) == 2) %% minimization
 		best_global = inf;
 	elseif( p(5) == 1) % maximization
 		best_global = -inf;
 	end

 	t = 1;

 	C_t = zeros( p(9) , p(2) );
 	E_t = zeros( p(9) , p(2) );

	while ( t <= p(7))

		% printf( '\nIteracion %i \n', t);
		E_t = C_Generation( p, Q_t );
		% E_t'

		if ( t == 1 )
			C_t = E_t;

			f_t = Eval( C_t , p );

			f_t1 = f_t;

		elseif ( t != 1 )

			x_c1 = zeros(1,p(2));
			x_c2 = zeros(1,p(2));

			for i=1:p(9)
				not_done = true;
  				while not_done
					a = floor(random( 'unif' , 0 , p(9) , [ 1 1 ] ) + 1 );
					not_done = (a == i);
				end

				not_done = true;
  				while not_done
					b = floor(random( 'unif' , 0 , p(9) , [ 1 1 ] ) + 1 );
					not_done = (b == i || b == a);
				end

				prob_population = rand;
				for j=1:p(2)
                    x_c1(j) = prob_population * E_t(a,j) + ( 1 - prob_population ) * C_t(b,j);
                    x_c2(j) = prob_population * C_t(a,j) + ( 1 - prob_population ) * C_t(b,j);
				end

				%%Crossover
				if( rand <= p(11) )
					E_t(a,:) = x_c1;
					E_t(b,:) = x_c2;
				end
			end
			% f_t1 = Eval( C_t , p );
			f_t   = Eval( E_t , p );
			C_t = E_t;
		end

		%C_t

		%Updatings waves
		prop_mod = rand;

		if( prop_mod > 0 )
			counter = 0;

			%%Comparison between(C_t) and best_global of last pop
			if( p(5) == 2 )
	            for i = 1 : p(9)
	                if( f_t(i) < best_global )
	                    counter++;
	                end
	                
	            end
	        elseif ( p(5) == 1 )
	        	for j = 1 : p(9)
	               	if( f_t(j) > best_global )
	               		counter++;
	               	end
	            end
	        end


            my_phi = counter / p(9);
            factor_mult = 0;

            if( my_phi < 0.2 )
                factor_mult = p(8);
            elseif( my_phi > 0.2 )
                factor_mult = 1.0 / p(8);
            elseif( my_phi == 0.2 )
                factor_mult = 1.0;
            end

            % factor_mult

            for i = 1 : p(1)
            	for j = 2:2:p(2)*2
                    Q_t(i,j) = Q_t(i,j) * factor_mult;
                end
            end

		end

		% Q_t
        %%Moving square pulse
		if( t > 1 && mod( t , p(12) ) == 0 ) 
            lambda = 0.5;
            f_t1 = f_t;

            if( p(5) == 2 )
            
	            for i=1:p(1)
	                
	                % [tmpval tmpidx] = min(f_t1);
	                % f_t1(tmpidx) = inf;

	                [val idx] = min(f_t1);
	                f_t1(idx) = inf;

	                % if( abs(tmpval - val ) < 3 )
	                % 	[val idx] = min(f_t1);
	                % 	f_t1(idx) = inf;
	                % end

	               	j = 1; k = 1;
	               	
	               	while( j <= p(2)*2 )
	                    Q_t(i,j) = Q_t(i,j) + lambda * ( C_t(idx,k) - Q_t(i,j) ) ;
	                    j = j+2;
	                    k++;
	                end
	            end
	        elseif ( p(5) == 1 )
	        	for i=1:p(1)
	                
	                [val idx] = max(f_t1);
	                f_t1(idx) = -inf;

	               	j = 1; k = 1;
	               	
	               	while( j <= p(2)*2 )
	                    Q_t(i,j) = Q_t(i,j) + lambda * ( C_t(idx,k) - Q_t(i,j) ) ;
	                    j = j+2;
	                    k++;
	                end
	            end
	        end
        end

		% myprint('Q_t',Q_t);

		% f_t = Eval( C_t , p );
	        	
        if( p(5) == 2 )
			[val idx] = min(f_t);
	    elseif ( p(5) == 1 )
	    	[val idx] = max(f_t);
	    end


	    if( p(5) == 2 )
			if( val < best_global )
				best_global = val;
			end
	    elseif ( p(5) == 1 )
	    	if( val > best_global )
				best_global = val;
			end
	    end

	    % best_global

		% [ C_t(idx,:) best_global ]

		plot(t,best_global);

		t++;
	end

	% myprint('C_t',C_t);

	f_t = Eval( C_t , p );
	if( p(5) == 2 )
		[bestF idx] = min(f_t);
	elseif( p(5) == 1 )
		[bestF idx] = max(f_t);
	end

	bestV = C_t(idx);

	myprint('Q_t',Q_t);

	toc;
endfunction



%%Parameters
function p = params()

	p = [];

	%%Parameters (1) Q_pop_size, (2) variables, (3) min_dom, (4) max_dom, (5) operation 1 min - 2 max, (6) n_worst, 
	%%(7) T, (8) phi, (9) C_pop_size, (10) function, p(11) update freq
	p(1) = 5; p(2) = 10; p(3) = -32; p(4) = 32; p(5) = 2; p(6) = 0.2*p(1); p(7) = 50; p(8) = 0.82; p(9) =500;
	p(10) = 4; p(11) = 0.66; p(12) = 1;

	if( p(10) == 1 )
		printf('ACkley\n');
	end

endfunction

%%Generation of Quantic Population
function Q_t = Q_Generation( p )

	Q_pop_size = p(1); variables = p(2); min_dom = p(3); max_dom = p(4); n_worst = p(6); T = p(7);

	Q_actual = zeros( Q_pop_size, variables*2 );

	pulso_ancho = ( max_dom - min_dom) / Q_pop_size;

	%%Filling matrix Q with sigma and mu
	for i=1:Q_pop_size
		if (i == 1)
				Q_actual(i,1) = min_dom + pulso_ancho/2 * i; 
			else
				Q_actual(i,1) = Q_actual(i-1,1) + pulso_ancho;
		end

		for j=2: size(Q_actual,2)
			if( mod(j,2) == 0)
				Q_actual(i,j) = pulso_ancho;
			else
				Q_actual(i,j) = Q_actual(i,1);
			end
		end
	end

	Q_t = Q_actual;

endfunction

%%Classic Generation
function C_t = C_Generation(p, Q_t)

	C_size = p(9); variables = p(2); Q_actual = Q_t;

	Xij = zeros( C_size, variables );

	%%Picking random q individual
	% idx = floor(random( 'unif' , 0 , 5 , [ 1 C_size ] ) + 1 );
	% for i=1:C_size
	% 	for j=1:variables
	% 		Xij(i,j) = rand * Q_actual(idx(i),j*2) + ( Q_actual(idx(i),1) - Q_actual(idx(i),j*2)/2.0 );
	% 	end
	% end

	%%Picking q ~ classical

	nq = p(9) / p(1) ;
	q = 1;

	i = 1;
	while( i < p(9) )
		for k=1:nq
			for j=1:p(2)
				Xij(i+k-1,j) = rand * Q_actual(q,j*2) + ( Q_actual(q,1) - Q_actual(q,j*2)/2.0 );
			end
		end
		q++;
		i = i + nq;
	end

	%% Check how to set less 0.5 floor and greater 0. ceil
	% for i=1:p(9)
	% 	if ( Xij(i,j) - )
	% 	Xij = floor( Xij );
	% end

	% printf( 'Xij: ');
	% Xij'

	C_t = Xij;

endfunction

function f = Eval( Xij , p )
	
	f = zeros(size(Xij,1),1);

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%>Function de Tarea AG
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% for i=1:size(Xij,1)
	% 	arg = 0;
	% 	for j=1:size(Xij,2)
	% 		arg += Xij(i,j)*Xij(i,j) ;
	% 	end
	% 	den = power(sin(sqrt(arg)), 2)- 0.5;
	% 	num = power((1+ 0.0001*arg),2);
	% 	f(i) = 0.5 - den/num;
	% end

	if( p(10) == 1 )
		%ackley
		for i=1:size(Xij,1)
			temp1 = 0;
			temp2 = 0;
			for j=1:size(Xij,2)
				temp1 += Xij(i,j) * Xij(i,j) ;
				temp2 += cos( 2*pi*Xij(i,j) ) ;
			end
			
			z = -20 * exp(-0.2*sqrt(1.0/p(2) * temp1)) - exp(1.0/p(2)*temp2) + 20 + exp(1);

			f(i) = z;
		end

	elseif( p(10) == 2 )
		%rastrigin
		x = Xij(:,1);
		y = Xij(:,2);
		
		temp1=x.^2-10*cos(2*pi.*x);
		temp2=y.^2-10*cos(2*pi.*y);
		z=temp1+temp2+20;
		f = z;

	elseif( p(10) == 3 )
		x = 0;
	elseif( p(10) == 4 )
		%schwefel
		x = Xij(:,1);
		y = Xij(:,2);
		z=(x-x.^2).^2+(x-1).^2+(x-y.^2).^2+(y-1).^2;

		f = z;

	elseif( p(10) == 5 )
		% sphere 
		x = Xij(:,1);
		y = Xij(:,2);
		z=x.^2+y.^2;
	
		f = z;
	end
endfunction


function myprint( msg , Xij )
	% format short g;
	printf(msg);
	printf('\n');
	for i=1:size(Xij,1)
		for j=1:size(Xij,2)
			printf( '%.3f\t', Xij(i,j) );
		end
		printf('\n');
	end
endfunction