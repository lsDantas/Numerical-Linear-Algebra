///////////////////////////////////////////////////////////
// Gauss-Seidel Method
//
// Description: Solves a linear system of the form
// Ax = b through the Gauss-Seidel method. This 
// implementation uses the standard Scilab matrix 
// inversion method.
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

function [x, norm_dif, counter, res_norm] = solve_gauss_seidel(A, b, x0, error_tol, max_iter, norm_type)
    // Partition A into triangular matrices (A = L + U)
    U = triu(A,1);
    inv_L = inv(tril(A));

    // Define method matrices
    method_matrix = inv_L*U;
    constant_term = inv_L*b

    // First Iteration

    // Use method matrix to update x
    x = constant_term - method_matrix*x0;

    // Compute norm difference and initialize counter
    norm_dif = norm(x - x0, norm_type);
    counter = 1;

    // Other Iterations
    while (norm_dif > error_tol & counter < max_iter)
        // Use method matrix to update x
        x_previous = x;
        x = constant_term - method_matrix*x_previous;

        // Compute norm difference and update counter
        norm_dif = norm(x - x_previous, norm_type);
        counter = counter + 1;
    end

    // Compute residue norm
    res_norm = norm(b - A*x);
endfunction
