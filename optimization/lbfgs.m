function [ fmin, xmin, Nit, cptm, err ] = lbfgs(x0,f,g,PlotFnc,Nmax)

%--------------------------------------------------------------------------
% L.Dede', 25 April 2007
% MOX - Politecnico di Milano
%--------------------------------------------------------------------------
%        [ fmin, xmin, Nit, cptm, err ] = lbfgs( f, g, x0, tol );
% 
% Find the minimum of the function f(x) \in R, with x \in R^N, given the 
% gradient vector g(x) \in R^N; the L-BFGS method is used. 
% The method is implemented as:
%            Liu, Nocedal. Math. Programming 45 (1989), 407-436
% and it uses the scaling techinque for the basic matrix proposed in:
%            Al-Baali. Numer. Algorithms 22 (1999) 99-112
% Linear search based on standard Armijo condition (in substitution of the 
% Wolfe-Powell conditions:
%           Wei, Li, Qi. App.Math.Comput. 179 (2006) 407-430
%
%
% Input:
%        f, g  =  strings for function and  gradient vector (use inline)
%        x0    =  initial guess (vector)
%        tol   =  stopping criterium tolerance ( ||g|| < tol )
% Output:
%        fmin  =  minimum value of f
%        xmin  =  vector x s.t. f(x)=fmin
%        Nit   =  number of iterations
%        cptm  =  cputime
%        err   =  errors (vector of ||g|| at each step k)
%--------------------------------------------------------------------------

% Parameters choice
%-------------------
tol = 1e-8;
M = 5;   % usually  \in (3,7) [Al-Baali]
rra = 0.5;   dda = 0.25; % parameters for Armijo cond. rra\in(0,1) dda\in(0,0.5)
Nja = 10;

% Initializations
%-----------------
N = length( x0 );
if ( size( x0, 2 ) > 1 )
    x0 = x0';
end

t = cputime;

k = 1;
Mh = min( k, M - 1 );

gk = sparse( N, 1 );
sk = sparse( N, 1 );
xk = sparse( N, 1 );
dgk = sparse( N, 1 );

gk( :, 1 ) = g( x0 );
sk( :, 1 ) = - gk;
xk( :, 1 ) = x0;
err( 1 ) = norm( gk( :, 1 ) );

% Iteration steps
%-----------------
while 1     
    disp('iteration')
    max(sk( :, k ))
    
%     % standard Armijo condition
    ja = 1;
    while ( ( ( f( xk( :, k ) + rra^ja * sk( :, k ) ) - f( xk( :, k ) ) ) > dda * rra^ja * gk( :, k )' * sk( :, k ) ) & ( ja < Nja ) )        
        ja = ja + 1;        
    end        
    alphak( k ) = rra^ja;
%     alphak( k ) = 1;
    
    % iterative step
    xk( :, k + 1 ) = xk( :, k ) + alphak( k ) * sk( :, k );
   
    % plot( 1: length( xk( :, 1 ) ), xk( :, k ), 'b', 1: length( xk( :, 1 ) ), xk( :, k+1 ), 'r' );
    % pause( 0.5 )    

    % stopping criterium
    gk( :, k + 1 ) = g( xk( :, k + 1 ) );
    err( k + 1 ) = norm( gk( :, k + 1 ) );    
    if ( err( k + 1 ) < tol || k >= Nmax)
        break
    end
   
    % matrix \times vector
    dgk( :, k ) = gk( :, k + 1 ) - gk( :, k );                
    %rhok( k ) = 1 / max( eps, alphak( k ) * dgk( :, k )' * sk( :, k ) );
    rhok( k ) = 1 / max( eps, abs(alphak( k ) * dgk( :, k )' * sk( :, k )) );
    if rhok( k ) > 1e10 || alphak( k ) * dgk( :, k )' * sk( :, k ) < 0
        aaa = 0;
    end
        
    % approximate inverse Hessian
       % 1
    vk = gk( :, k + 1 );
    
    for j = k : -1 : ( k - Mh + 1 )
        vk = vk - rhok( j ) * alphak( j ) * ( sk( :, j )' * vk ) * dgk( :, j );
    end
    
    nuk = max( alphak( k - Mh + 1 )' * ( sk( :, k - Mh + 1 )' * dgk( :, k - Mh + 1 ) ) / ( dgk( :, k - Mh + 1 )' * dgk( :, k - Mh + 1 ) ), ...
               alphak( k )' * ( sk( :, k )' * dgk( :, k ) ) / ( dgk( :, k )' * dgk( :, k ) ) );
    vk = nuk * vk;
    
    for j = ( k - Mh + 1 ) : k
        vk = vk - rhok( j ) * alphak( j ) * ( dgk( :, j )' * vk ) * sk( :, j );
    end
    
       % 1 : Mh-1
    wk = zeros( N, 1 );
    
    for j = 1 : ( Mh - 1 )
        
        wk( :, j ) = gk( :, k + 1 );                
        
        for i = k : -1 : ( k - Mh + j )
            wk( :, j ) = wk( :, j ) - rhok( i ) * alphak( i ) * ( sk( :, i )' * wk( :, j ) ) * dgk( :, i );            
        end
        
        wk( :, j ) = rhok( k - Mh + j )  * alphak( k - Mh + j )^2 * ...
                     ( sk( :, k - Mh + j )' * wk( :, j ) ) * sk( :, k - Mh + j );
        
        for i = ( k - Mh + j ) : k
            wk( :, j ) = wk( :, j ) - rhok( i ) * alphak( i ) * ( dgk( :, i )' * wk( :, j ) ) * sk( :, i );
        end
       
    end    
    
       % end
    sk( :, k + 1 ) = - vk - sum( wk, 2 ) ...
                     - rhok( k ) * alphak( k )^2 * ( sk( :, k )' * gk( :, k + 1 ) ) * sk( :, k );        
                     
    % update
    k = k + 1;    
    Mh = min( k, M - 1 );
       
    % remove not used steps
    xk( :, k - 1 ) = sparse( N, 1 );
    gk( :, k - 1) = sparse( N, 1 );
    if k > M + 1        
        dgk( :, k - M - 1 ) = sparse( N, 1 );
        sk( :, k - M - 1 ) = sparse( N, 1 );      
    end
    
    PlotFnc(xk( :, k));
end

% Output
%-------- 
cptm = cputime - t;
xmin = xk( :, end );
fmin = f( xmin );
Nit = k;

end