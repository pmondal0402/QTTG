function [S] = get_entropy(rho) 
  % ----------------------------------------------
  % Returns log of density matrix rho 
  % epsilon is a small number
  % Note :
  %     Verified computing entropy for 3 e-3 spin
  % ----------------------------------------------

  eps = 10^-50 ; 
  [V, E] = eig(full(rho)) ;

  % If eigenvalue is less than 10^-6, replace it with eps
  % This is done for numerical accuracy
  for ii = 1:length(E)
    if E(ii, ii) <= 10^-8 
       E(ii, ii) = eps ;
    end
  end
  
  S = -trace(E * diag( log(diag(E))/log(2) ) ) ; 
  
