///////////////////////////////////////////////////////////
// Shifted Inverse Power Iteration
// 
// Description: Approximates closest eigenvalue to a
// particular value along with its corresponding 
// eigenvector. Uses the shifted inverse power iteration 
// method. 
///////////////////////////////////////////////////////////
// Input:
//   A:        a real n x n diagonalizable matrix 
//   x0:       initial approximation for eigenvector
//   epsilon:  precision used as halt criterion
//   alpha:    reference value for eigenvalue search
//   max_iter: maximum number of iterations
///////////////////////////////////////////////////////////
// Output:
//   lambda:      closest eigenvalue to alpha
//   x1:          eigenvector corresponding to
//                lambda
//   counter:     number of iterations executed
//   norm_error:  Infinity norm of difference
//                between subsequent iterations
///////////////////////////////////////////////////////////

function [lambda, x1, counter, norm_error] = power_iteration_shifted_inverse(A, x0, epsilon, alpha, max_iter)
    // Determine matrix dimensions and start counter
    [n]=size(A,1);
    counter = 1;

    // First iteration 
    
    // Obtain x1 through linear system
    x0 = x0/norm(x0,2);
    [x1, LU, P] = Gaussian_Elimination_4(A - alpha * eye(n,n), x0);
    x1 = x1/norm(x1,2);

    // Approximate lambda
    lambda = x1'*(A*x1);

    // Verify Precision
    norm_error = norm(x1 - x0, %inf);
    if norm_error < epsilon
        disp("Desired precision reached!");
        return;
    end

    // Store iteration and update counter
    x0 = x1;
    counter = counter + 1;

    // Other Iterations
    while counter < max_iter
        // Obtain x1 through linear system
        x1 = solve_with_LU(LU, x0, P);
        x1 = x1/norm(x1,2);

        // Approximate Lambda
        lambda = x1'*(A*x1);

        // Verify precision
        norm_error = norm(x1 - x0, %inf);
        if norm_error < epsilon
            disp("Desired precision reached!");
            return;
        end

        // Store iteration and update counter
        x0 = x1;
        counter = counter + 1;
    end
endfunction