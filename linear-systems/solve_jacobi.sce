///////////////////////////////////////////////////////////
// Jacobi Method
//
// Description: Solves a linear system of the form
// Ax = b through the Jacobi method.
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

function [x, norm_dif, counter, res_norm] = solve_jacobi(A, b, x0, error_tol, max_iter, norm_type)
    // Partition A across diagonal (A = L + D + U)
    LU = A - diag(diag(A))
    inv_D = diag(1./diag(A))

    // Define method matrices
    method_matrix = inv_D*LU;
    constant_term = inv_D*b;
    
    // First Iteration
    x = constant_term - method_matrix*x0;
    norm_dif = norm(x - x0, norm_type);
    counter = 1;

    // Other Iterations
    while (norm_dif>error_tol & counter<max_iter)
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
