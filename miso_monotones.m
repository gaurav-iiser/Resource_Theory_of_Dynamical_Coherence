addpath(genpath('QETLAB-0.9'))
%{
%%%%%%%  Interconversion using MISC for qubit channels using the complete set of monotones as given in the paper 'Dynamical Resource Theory of Quantum Cohernce' by Gaurav Saxena, Eric Chitambar, and Gilad Gour.
This code is written by Gaurav Saxena.
This code is designed specifically for qubit channels.  %%%%%%%%%%
%}
%{
Pseudo Code:
    1. The problem is to maximize the inner product of the choi matrix of the
        DISO superchannel(J) with the tensor product of the choi matrix of a given
        channel and another quantum channel.
    2. Initialize different types of \Phi and \Lambda. Initially take
    \Lambda to be a replacement channel.
    3. The constraints are:
        a. J >= 0 (i.e. the choi matrix belongs to the cone of positive semi-definite
           matrices.)
        b,c. Impose the condition on J that makes it a superchannel.
        d. The condition that makes a superchannel a DISO.



Requirements:
1) Matlab
2) Qetlab
3) CVX

%}


%clear all
fileID = fopen('data.txt','w');

 %J_Lambda = RandomSuperoperator(2, 1, 0,  1)
    %%{
    plus_state = (1/sqrt(2))*[1; 1];
    plus_state_density_mat = plus_state * plus_state';
    J_Lambda = Tensor(eye(2), plus_state_density_mat ); %Replacement Channel (outputs a plus state)
    %J_Lambda = Tensor(eye(2), RandomDensityMatrix(2,1) ); %Replacement Channel (outputs a random state)
     %J_Lambda = MaxEntangled(2) * MaxEntangled(2)';
    %}
    
clear B
for count = 1:5
    
    
    %J_Phi = RandomSuperoperator(2, 1, 0,  1);
    %J_Phi = Tensor(eye(2),RandomDensityMatrix(2,1));%%%%% Choi matrix of the replacement channel with a random density matrix as output %%%%%%%%%%%
    
    %zero_zero = [1 0; 0 0];
    J_Phi = Tensor(eye(2), plus_state_density_mat );                  %%%%% Choi matrix of the replacement channel with output |+><+| %%%%%%%%%%% 
    %J_Phi = MaxEntangled(2) * MaxEntangled(2)';     %%%%%          Choi matrix of the identity channel              %%%%%%%%%%%
    %J_Phi = [ 1 0 0 1; 0 0 0 0; 0 0 0 0; 1 0 0 1];  %%%%%          Choi matrix of the identity channel   Unnormalized %%%%%%%%%%%
    %{
    J_Phi = zeros(4);
        random_num1 = rand;
        random_num2 = rand;
        J_Phi(1,1) = random_num1;
        J_Phi(2,2) = 1 - random_num1;
        J_Phi(3,3) = random_num2;
        J_Phi(4,4) = 1 - random_num2;
    %}

   

    J_Phi_Lambda = Tensor(J_Phi, J_Lambda);

    u_A1 = eye(dim) / trace( eye(dim) );
    %count = 0;

    cvx_begin sdp quiet
        dimension = dim;
        variable J(dimension^4,dimension^4) semidefinite
        maximize ( trace(J' * J_Phi_Lambda ) )
        subject to
            PartialTrace(J, [4],[dimension, dimension, dimension, dimension]) == Tensor( PartialTrace(J, [2,4],[dimension, dimension, dimension, dimension]) , u_A1)
            PartialTrace(J, [1,4],[dimension, dimension, dimension, dimension]) == eye(4)
            J == J'
            for f_1 = 1:dim
                for f_2 = 1:dim
                    f = [ f_1  f_2 ];
                    for i = 1:dim               %%for B_0
                        for j = 1:dim           %%for B_0
                           for u = 1:dim        %%for B_1
                                for v = 1:dim   %%for B_1
                                    if (i ~= j) | (u ~= v)
                                        trace(J' * Tensor(alpha(f), basis(i,j), basis(u,v) ) ) == 0
                                    end
                                end
                           end
                        end
                    end
                end
            end
    cvx_end
    A = sort(diag(J_Phi),'descend');
    B(count) = cvx_optval;
    fprintf(fileID,'%d \t %d \t %d \t %d \t : %d \n',A(1),A(2),A(3),A(4),  cvx_optval);

end
fclose(fileID);


%%%%%%%%%%%%%%%%%%%%%%%%%% 
%  Function Definitions  %
%%%%%%%%%%%%%%%%%%%%%%%%%% 

function pauli_basis = Pauli_basis
    I = eye(2);
    sigma_1 = [0, 1; 1  0];
    sigma_2 = [0 -1i; 1i  0];
    sigma_3 = [1  0; 0 -1];
    pauli_basis = I;
    pauli_basis(:,:,2) = sigma_1;
    pauli_basis(:,:,3) = sigma_2;
    pauli_basis(:,:,4) = sigma_3;
end

function a = alpha(f)
   d = dim;
   a = zeros(d^2, d^2);
   for k = 1: dim
      a = a + Tensor(basis(k,k),basis(f(k),f(k))); %function definition needs to go in
   end
end

function basis_elem = basis(i,j)
        basis_elem = Tensor(x(i), x(j)');
end

function basis_vector = x(a)
    basis_vector = zeros(dim,1);
    basis_vector(a) = 1;
end

function d = dim
    d = 2;
end
