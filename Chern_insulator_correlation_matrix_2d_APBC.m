classdef Chern_insulator_correlation_matrix_2d_APBC
    methods
        function correlation_matrix = Chern_insulator_correlation_matrix_func(~, u, l_x)
            %we conside a 2D square lattice with dimension lx*lx
            %u is the parameter that tunes the model across different Chern
            %number regimes
            
            spin = 2;
            H = zeros(l_x*l_x*spin, l_x*l_x*spin);
            
            %hopping in the x-direction
            for y = 0:l_x-1
                for x = 1:l_x-1
                    %the tau_x term: -i/2*(f_r^\dagger \tau_x f_(r+a_x) - h.c.)
                    H((y*l_x+x)*spin-1, (y*l_x+x+1)*spin) = -1j/2;
                    H((y*l_x+x)*spin, (y*l_x+x+1)*spin-1) = -1j/2;
                    H((y*l_x+x+1)*spin-1, (y*l_x+x)*spin) = 1j/2;
                    H((y*l_x+x+1)*spin, (y*l_x+x)*spin-1) = 1j/2;
                    
                    %the tau_z term: 1/2*(f_r^\dagger \tau_z f_(r+a_x) + h.c.)
                    H((y*l_x+x)*spin-1, (y*l_x+x+1)*spin-1) = 1/2;
                    H((y*l_x+x)*spin, (y*l_x+x+1)*spin) = -1/2;
                    H((y*l_x+x+1)*spin-1, (y*l_x+x)*spin-1) = 1/2;
                    H((y*l_x+x+1)*spin, (y*l_x+x)*spin) = -1/2;
                end
                %anti-periodic boundary condition
                %the tau_x term: -i/2*(f_r^\dagger \tau_x f_(r+a_x) - h.c.)
                H((y*l_x+l_x)*spin-1, (y*l_x+1)*spin) = 1j/2;
                H((y*l_x+l_x)*spin, (y*l_x+1)*spin-1) = 1j/2;
                H((y*l_x+1)*spin-1, (y*l_x+l_x)*spin) = -1j/2;
                H((y*l_x+1)*spin, (y*l_x+l_x)*spin-1) = -1j/2;
                    
                %the tau_z term: 1/2*(f_r^\dagger \tau_z f_(r+a_x) + h.c.)
                H((y*l_x+l_x)*spin-1, (y*l_x+1)*spin-1) = -1/2;
                H((y*l_x+l_x)*spin, (y*l_x+1)*spin) = 1/2;
                H((y*l_x+1)*spin-1, (y*l_x+l_x)*spin-1) = -1/2;
                H((y*l_x+1)*spin, (y*l_x+l_x)*spin) = 1/2;
            end
            %hopping in the y-direction
            for x = 1:l_x
                for y = 0:l_x-2
                    %the tau_y term: -i/2*(f_r^\dagger \tau_y f_(r+a_y) - h.c.)
                    H((y*l_x+x)*spin-1, ((y+1)*l_x+x)*spin) = -1j*(-1j/2);
                    H((y*l_x+x)*spin, ((y+1)*l_x+x)*spin-1) = 1j*(-1j/2);
                    H(((y+1)*l_x+x)*spin-1, (y*l_x+x)*spin) = 1j*(-1j/2);
                    H(((y+1)*l_x+x)*spin, (y*l_x+x)*spin-1) = -1j*(-1j/2);
                    
                    %the tau_z term: 1/2*(f_r^\dagger \tau_z f_(r+a_y) + h.c.)
                    H((y*l_x+x)*spin-1, ((y+1)*l_x+x)*spin-1) = 1/2;
                    H((y*l_x+x)*spin, ((y+1)*l_x+x)*spin) = -1/2;
                    H(((y+1)*l_x+x)*spin-1, (y*l_x+x)*spin-1) = 1/2;
                    H(((y+1)*l_x+x)*spin, (y*l_x+x)*spin) = -1/2;
                end
                %anti-periodic boundary condition
                %the tau_y term: -i/2*(f_r^\dagger \tau_y f_(r+a_y) - h.c.)
                H(((l_x-1)*l_x+x)*spin-1, x*spin) = 1j*(-1j/2);
                H(((l_x-1)*l_x+x)*spin, x*spin-1) = -1j*(-1j/2);
                H(x*spin-1, ((l_x-1)*l_x+x)*spin) = -1j*(-1j/2);
                H(x*spin, ((l_x-1)*l_x+x)*spin-1) = 1j*(-1j/2);
                    
                %the tau_z term: 1/2*(f_r^\dagger \tau_z f_(r+a_y) + h.c.)
                H(((l_x-1)*l_x+x)*spin-1, x*spin-1) = -1/2;
                H(((l_x-1)*l_x+x)*spin, x*spin) = 1/2;
                H(x*spin-1, ((l_x-1)*l_x+x)*spin-1) = -1/2;
                H(x*spin, ((l_x-1)*l_x+x)*spin) = 1/2;
            end
            
            %the u-term
            for x = 1:l_x
                for y = 0:l_x-1
                    H((y*l_x+x)*spin-1, (y*l_x+x)*spin-1) = u;
                    H((y*l_x+x)*spin, (y*l_x+x)*spin) = -u;
                end
            end

            [Eigvec, Eigval] = eig(H); %obtain the eigenvectors and eigenvalues of H
                                        %the output is defined by H*Eigvec = Eigvec*Eigval
            %define the correlation matrix D_ij = <f_i f_j^\dagger> in the diagonal basis
            %the ground state only consists of states with energy <= 0
            D_diag = zeros(l_x*l_x*spin, l_x*l_x*spin);
            for i = 1:l_x*l_x*spin
                if Eigval(i, i) > 0 
                    D_diag(i, i) = 1;
                end 
                if abs(Eigval(i, i)) < 10^(-10)
                    disp("warning, zero eigenvalue");
                end
            end

            %if VHV^\dagger diagonalizes H, D = V^\dagger*D_diag*V
            correlation_matrix = Eigvec * D_diag * ctranspose(Eigvec);
        end
    end
end
