///////////////////////////////////////////////////////////
// Custom Gauss-Seidel Method
//
// Description: Solves a linear system of the form
// Ax = b through the Gauss-Seidel method. This 
// implementation uses a custom approach to invert
// matrices.
///////////////////////////////////////////////////////////
// Input:
//   A:         an n x n matrix (?)
//   b:         an n x 1 vector (?)
//   x0:        initial approximation for the solution 
//              of the system
//   error_tol: the error tolerance
//   max_iter:  maximum number of iterations
//   norm_type: prefered norm
///////////////////////////////////////////////////////////
// Output:
//   x:        approximate solution of the system
//   norm_dif: norm difference of the last two iterations
//   counter:  number of iterations executed
//   res_norm: residue norm
///////////////////////////////////////////////////////////

function [x, norm_dif, counter, res_norm] = solve_gauss_seidel_custom(A, b, x0, error_tol, max_iter, norm_type)
    // Determine dimensions
    [n]=size(A,1);

    // Partition A into triangular matrices (A = L + U)
    U = triu(A,1);
    L = tril(A);
    
    // First Iteration

    // Set up linear system
    x = zeros(n,1);
    system_vector = -U*x0 + b;

    // Solve through forward substitution
    for i=1:n
        // Sums products of known terms with elements of L
        sum_terms = sum( L(i,1:(i-1))*x(1:(i-1),1) );

        // Find xi through known elements
        x(i,:) = (system_vector(i,:) - sum_terms)/L(i,i);
        disp(i);
    end

    // Compute norm difference and initialize counter
    norm_dif = norm(x - x0, norm_type);
    counter = 1;

    // Other Iterations
    while (norm_dif > error_tol & counter < max_iter)
        // Set up linear system
        x_previous = x;
        system_vector = -U*x_previous + b;

        // Forward Substitution
        for i=1:n
            // Sums products of known terms with elements of L
            sum_terms = sum( L(i,1:(i-1))*x(1:(i-1),1) );

            // Find xi through known elements
            x(i,1) = (system_vector(i,1) - sum_terms)/L(i,i);
        end

        // Compute norm difference and update counter
        norm_dif = norm(x - x_previous, norm_type);
        counter = counter + 1;
    end

    // Compute residue norm
    res_norm = norm(b - A*x);
endfunction
