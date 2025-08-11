lx = 50;
u_start = 0;
u_end = 3;
u_list = u_start+0.01:0.01:u_end;
y_ratio = [1 1 2];
x_ratio = [1 1 2];
spin = 2; %Chern insulator

%get the dimension of the subsystems
%on a lx*lx square lattice, subsystem C has 1 <= y <= lx, a < x <= lx
%subsystem A and B has 1 <= y <= b, 1 <= x <= a and b < y <= lx, 1 <= x <= a respectively
%where 1 < a, b < lx is some cutoff determined by x_ratio and y_ratio
d_Ax = lx*x_ratio(1)/x_ratio(3);
d_Bx = lx*x_ratio(1)/x_ratio(3);
d_Cx = lx*x_ratio(2)/x_ratio(3);
d_Ay = lx*y_ratio(1)/y_ratio(3);
d_By = lx*y_ratio(2)/y_ratio(3);
d_Cy = lx;
d_AB = (d_Ay + d_By) * d_Ax * spin;

%we construct the list of indices of sites that belong to subsystem
%A, B and C respectively, note that in the list of sites of size
%lx*ly, site (x, y) is in the position labeled by (y-1)*lx+x
A_sec = get_sec_index(1:d_Ax, 1:d_Ay, lx);
B_sec = get_sec_index(1:d_Ax, d_Ay+1: lx, lx);
C_sec = get_sec_index(d_Ax+1:lx, 1:lx, lx);
AB_sec = get_sec_index(1:d_Ax, 1:lx, lx);

%define lists to store data
kappa_list = zeros(1, length(u_list));
G_list = zeros(1, length(u_list));
R2_list = zeros(1, length(u_list));


for u_index = 1:length(u_list)
    u = u_list(u_index);
    disp(u);
        
    %get the correlation matrix <f_if_j^dagger> for the original state
    correlationObj = Chern_insulator_correlation_matrix_2d_APBC;
    D_full = correlationObj.Chern_insulator_correlation_matrix_func(u, lx);
        
    %if D_full does not look hermitian to matlab due to floating point
    %precision, we make it hermitian by manually setting D_ij =
    %conj(D_ji) so that the matrix of eigenvectors calculated by matlab 
    %would be unitary
    if not(ishermitian(D_full))
        if max(ctranspose(D_full)-D_full, [], "all") > 10^(-10)
            disp("error, D_full not hermitian");
        else
            D_full = make_hermitian(D_full);
        end
    end
    C_full = eye(lx*lx*spin) - transpose(D_full);
        
    D_orig = D_full(AB_sec, AB_sec); %correlation matrix restricted to subsystem AB
    C_orig = eye(d_AB) - transpose(D_orig);
        
    %the covariance matrix has the form (C-C^T)*1 + （1-C-C^T)*sigma_y
    I2 = eye(2);
    sigma_y = [0 -1j; 1j 0];
    Gamma_orig = kron(C_orig - transpose(C_orig), I2) + kron(eye(d_AB) - C_orig - transpose(C_orig), sigma_y);
        
    Gamma_real = Gamma_orig/1j;
    %schur decomposition in matlab for real matrices is the real schur
    %decomposition which returns O*E*transpose(O) = Gamma_real and 
    %O*transpose(O) = I
    %real schur decomposition coincides with block diagonalization for 
    %skew symmetric matrices
    [O, E] = schur(Gamma_real);
        
    %check O is orthogonal
    if norm(O*transpose(O) - eye(d_AB*2)) > 10^(-10)
        disp("O not orthogonal");
        disp(norm(O*transpose(O) - eye(d_AB*2)));
    end
        
    %construct the covariance matrix for rho^(2)
    %the construction exactly corresponds to the form of matrix on the note
    M_diag = zeros(2*d_AB, 2*d_AB);
    M_2_diag = zeros(2*d_AB, 2*d_AB);
    %this is to test that the Schur decomposition correctly gives the
    %eigenvalues, can comment out the corresponding lines in actual runs
    Gamma_test = zeros(2*d_AB, 2*d_AB);
    for i = 1:d_AB
        gamma_k_current = -E(2*i-1, 2*i);
        factor = 1+gamma_k_current^2;
        M_diag(2*i-1, 2*i) = -1j*2*gamma_k_current/factor;
        M_diag(2*i, 2*i-1) = 1j*2*gamma_k_current/factor;
        M_2_diag(2*i-1, 2*i-1) = -1j*(1-gamma_k_current^2)/factor;
        M_2_diag(2*i, 2*i) = -1j*(1-gamma_k_current^2)/factor;
            
        %a sanity check kept for (1, 2)-reflected entropy
        Gamma_test(2*i-1, 2*i) = -1j*gamma_k_current;
        Gamma_test(2*i, 2*i-1) = 1j*gamma_k_current;
 
    end
        
    if norm(O*Gamma_test*transpose(O) - Gamma_orig) > 10^(-10)
        disp("Incorrect block diagonalization");
        disp(norm(O*Gamma_test*transpose(O)-Gamma_orig));
    end
    %disp(max(transpose(O)*O - eye(d_AB*2), [], "all"));
    M = O*M_diag*transpose(O);
    M_2 = O*M_2_diag*transpose(O);

    %matrix defined to indicate the location of the submatrices
    I11 = [1 0; 0 0];
    I12 = [0 1; 0 0];
    I21 = [0 0; 1 0];
    I22 = [0 0; 0 1];
    Gamma_ABAB = kron(I11, M) + kron(I12, M_2) + kron(I21, -M_2) + kron(I22, -M);
        
    %get the submatrix of the covariance matrix corresponding to the
    %subsystem AA'
    %by the way of labeling the sites, all sites in subsystem A has
    %smaller label than sites in subsystem B
    AA_index = [1:d_Ax*d_Ay*2*spin d_AB*2+1:d_AB*2+d_Ax*d_Ay*2*spin];
    Gamma_AA = Gamma_ABAB(AA_index, AA_index);
        
    %calculate G(A:B:C), which can be expressed as 1/2 the (2, 2)-Renyi
    %reflected entropy plus the 2-Renyi entropy of rho_AB
    R2 = calculate_entropy_covariance(Gamma_AA, 1);
    G = R2/2 + calculate_entropy_correlation(C_orig, 1);
    disp(G);
    G_list(u_index) = G;
    R2_list(u_index) = R2;
        
    %kappa = G - 1/2(S_A^(2) + S_B^(2) + S_C^(2))
    C_A = C_full(A_sec, A_sec);
    C_B = C_full(B_sec, B_sec);
    C_C = C_full(C_sec, C_sec);
    sum_S_2 = calculate_entropy_correlation(C_A, 1/2) + calculate_entropy_correlation(C_B, 1/2) + calculate_entropy_correlation(C_C, 1/2);
    kappa = G - sum_S_2;
    disp(kappa);
    kappa_list(u_index) = kappa;  
end

%export the result to a csv file
file_name = sprintf("kappa_u_%d_%d_Chern_L_%d_APBC.csv", u_start, u_end, lx);
data_matrix = transpose([u_list; kappa_list; G_list; R2_list]);
writematrix(data_matrix, file_name);


function sec_index = get_sec_index(x_range, y_range, lx)
%x_range and y_range are defined such that sites (x_range(i), y_range(j))
%is in the subsystem for all 1 <= i <= len(x_range) and 1 <= j <= len(y_range)
    sec_index = zeros(1, length(x_range)*length(y_range));
    for i = 1:length(y_range)
        for j = 1:length(x_range)
            %Chern insulator, spin up and down
            sec_index(((i-1)*length(x_range)+j)*2-1) = ((y_range(i)-1)*lx+x_range(j))*2-1;
            sec_index(((i-1)*length(x_range)+j)*2) = ((y_range(i)-1)*lx+x_range(j))*2;
        end
    end
end

%this function returns a times the 2-renyi entropy for a given correlation matrix C
function renyi2_entropy = calculate_entropy_correlation(C, a)
    eigval = eig(C);
    renyi2_entropy = 0;
    for i = 1:length(eigval)
        renyi2_entropy = renyi2_entropy - a*log(eigval(i)^2 + (1-eigval(i))^2);
    end
end

%this function returns a times the 2-renyi entropy for a given covariance matrix Gamma
function renyi2_entropy = calculate_entropy_covariance(Gamma, a)
    eigval = eig(Gamma);
    %Due to floating point precision, Gamma may not be considered by
    %matlab as an exactly Hermitian matrix, so the eigenvalues returned
    %will contain an imaginary part a+0j, which messes up the sorting
    for i=1:length(eigval)
        %sanity check, make sure the eigenvalues are real
        if abs(imag(eigval(i)) - 0) > 10^(-8)
            disp(abs(imag(eigval(i)) - 0));
            disp("error, eigenvalues of covariance matrix not real");
         end
     end
     eigval_sorted = sort(real(eigval));
     %disp(eigval_sorted);
     renyi2_entropy = 0;
     for i = 1:length(eigval_sorted)/2
         epsilon = (1 - eigval_sorted(i))/2;
         renyi2_entropy = renyi2_entropy - a*log(epsilon^2 + (1-epsilon)^2); 
         %eigenvalues of the covariance matrix come in pairs, we only
         %use one of the pair of the eigenvalues
         if abs(-eigval_sorted(i) - eigval_sorted(length(eigval_sorted)+1-i)) > 10^(-10)
             disp("error, eigenvalues of covariance matrix not in pair");
         end
     end
end

%for a matrix M that is hermitian but doesn't "look" hermitian for matlab 
%due to floating point precision, this function returns a modified matrix M
%that is hermitian by changing all M_ij for i <= j to be equal to the
%complex conjugate of M_ji
function HermitianM = make_hermitian(M)
    HermitianM = zeros(size(M, 1), size(M, 2));
    for i = 1:size(M, 1)
        HermitianM(i, i) = real(M(i, i));
        for j = i+1:size(M, 1)
            HermitianM(j, i) = M(j, i);
            HermitianM(i, j) = conj(M(j, i));
        end
    end
end


