%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOCAL INTEGRAL F+
% 
% This function returns the nxn matrix (n is the number of nodes) that has
% to be multiplied by the density vector to get F+
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function matrix=Fplus(x,sigma)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% CHOOSE BETWEEN TRAPEZOIDAL OR SIMPSON RULES
%
choice=1; % Trapezoidal rule
%
%choice=2; % Simpson rule
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n=length(x);
deltax=x(2:end)-x(1:end-1); %diff(x)
matrix=zeros(n,n);

% TRAPEZOIDAL RULE

if choice==1
   
    for i=1:n-1       
        
    % FIRST: always x(n) is lower than x(i)+sigma
    % If node i enters here, the rest of the matrix is filled with a similar
    ... structure formed by 1/2 and 1
        
        if sigma+x(i)>=x(n)
              
            for j=i:n-1 % Going down in the rows from row i
            
                matrix(j,j)=deltax(j)/2; % First node in the diagonal
                matrix(j,n)=deltax(n-1)/2; % Last node in last column
                
                if j~=n-1 % In row n-1 there are only 2 elements
                    
                  % Taking all the colums
             
                   matrix(j,j+1:n-1)=(deltax(j:n-2)+deltax(j+1:n-1))/2;                 
                   
                end
            
            end
            
            break % Finish the for loop since the matrix has been filled
        
            
        % SECOND: x(n) is higher than x(i)+sigma
        % More delicate since the rows are not filled until the last column 
        % Each row is filled individually
        
        elseif sigma+x(i)<x(n)
            
            counter=0;
            
            for j=i+1:n % Going to the right in the column
                
                if x(i)+sigma<x(j) % When this loop is entered, j is the 
                    ...first node in which x(i)+sigma<x(j). Said in another
                    ...way, the first one outisde the integration domain.
                    
                % Three types of rows can occur:
                    
                    if counter==0 % Row with only two elements
                        
                        chi=sigma-(x(j-1)-x(i));
                        
                        matrix(i,i)=chi-chi^2/(2*deltax(j-1));
                        matrix(i,j)=chi^2/(2*deltax(j-1));
                        
                    elseif counter==1 % Row with only three elements
                        
                        chi=sigma-(x(j-1)-x(i));
                        
                        matrix(i,i)=deltax(i)/2;
                        matrix(i,j-1)=deltax(j-2)/2+chi-chi^2/(2*deltax(j-1));
                        matrix(i,j)=chi^2/(2*deltax(j-1));
                        
                    elseif counter>1 % Row with more than three elements
                        
                        chi=sigma-(x(j-1)-x(i));
                      
                        matrix(i,i)=deltax(i)/2;
                        matrix(i,i+1:j-2)=0.5*(deltax(i+1:j-2)+deltax(i:j-3));
                        matrix(i,j-1)=deltax(j-2)/2+chi-chi^2/(2*deltax(j-1));
                        matrix(i,j)=chi^2/(2*deltax(j-1));
                        
                    end
                    
                    break % The column has been finished
                    
                else % This counter distinguishes the type of row
                    counter=counter+1;
                end
                
            end
            
        end
        
    end
                            
                        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% SIMPSON'S RULE

elseif choice==2
    for i=1:n-1
        
        % FIRST: always x(n) is lower than x(i)+sigma
        % If node i enters here, the rest of the matrix is filled with a similar
        % structure
        
        if sigma+x(i)>=x(n)
            
            for j=i:n-1 % Going down in the rows from row i
                
                numbernodes=n-j+1;
                
                if mod(numbernodes,2)==1 % Odd number of nodes: exact Simpson
                    
                    matrix(j,j)=1/6*((deltax(j)+deltax(j+1))*(2*deltax(j)-...
                        deltax(j+1))/deltax(j)); % First node in the diagonal
                    matrix(j,n)=1/6*((deltax(n-1)+deltax(n-2))*(2*deltax(n-1)-...
                        deltax(n-2))/deltax(n-1)); % Last node in last column
                    
                    % The second node is put here since it does not enter
                    % in the loop
                    matrix(j,j+1)=1/6*((deltax(j)+deltax(j+1))^3/deltax(j)/deltax(j+1));
                    
                    % The rest of the nodes are computed in a loop
                    
                    % This loop is entered only if there are not 3 nodes in
                    % the row (row n-2)
                    
                    if j~=n-2 % Row j==n-1 will never have odd number of nodes...
                        % and will not enter the previous if
                        
                        for k=1:(numbernodes-1)/2-1
                            
                            % Firstly, from the 3 points, the ones of the
                            % extremes (sum of two parts)
                            matrix(j,j+2*k)=1/6*((deltax(j+2*k-2)+deltax(j+2*k-1))...
                                *(2*deltax(j+2*k-1)-deltax(j+2*k-2))/deltax(j+2*k-1))...
                                +1/6*((deltax(j+2*k)+deltax(j+2*k+1))*(2*deltax(j+2*k)...
                                -deltax(j+2*k+1))/deltax(j+2*k));
                            
                            % Secondly, from the 3 points, the ones in the
                            % middle (only one part)
                            matrix(j,j+2*k+1)=1/6*((deltax(j+2*k+1)+deltax(j+2*k))^3/...
                                deltax(j+2*k+1)/deltax(j+2*k));
                            
                        end
                        
                    end
                    
                elseif mod(numbernodes,2)==0 % Even number of nodes: exact Simpson
                    % until column n-1, and trapezoidal between columns n-1 and n
                    
                    % We add the part between n-1 and n with trapezoidal
                    % rule
                    
                    matrix(j,n-1)=deltax(n-1)/2;
                    matrix(j,n)=deltax(n-1)/2;
                    
                    % Now we proceed if we are not in the n-1 row (it has
                    % only two elements)
                    
                    if numbernodes~=2
                        
                        matrix(j,j)=1/6*((deltax(j)+deltax(j+1))*(2*deltax(j)-...
                            deltax(j+1))/deltax(j)); % First node in the diagonal
                        matrix(j,n-1)=matrix(j,n-1)+ 1/6*((deltax(n-2)+deltax(n-3))*...
                            (2*deltax(n-2)-deltax(n-3))/deltax(n-2)); % Node in n-1 column
                        
                        
                        % The rest is computed in a similar way as before
                        
                        % The second node is put here since it does not enter
                        % in the loop
                        matrix(j,j+1)=1/6*((deltax(j)+deltax(j+1))^3/deltax(j)/deltax(j+1));
                        
                        % The rest of the nodes are computed in a loop
                        
                        % This loop is entered only if there are not 4 nodes in
                        % the row (row n-3 has 4 elements)
                        
                        if j~=n-3
                            
                            for k=1:(numbernodes-2)/2-1
                                
                                % Firstly, from the 3 points, the ones of the
                                % extremes (sum of two parts)
                                matrix(j,j+2*k)=1/6*((deltax(j+2*k-2)+deltax(j+2*k-1))...
                                    *(2*deltax(j+2*k-1)-deltax(j+2*k-2))/deltax(j+2*k-1))...
                                    +1/6*((deltax(j+2*k)+deltax(j+2*k+1))*(2*deltax(j+2*k)...
                                    -deltax(j+2*k+1))/deltax(j+2*k));
                                
                                % Secondly, from the 3 points, the ones in the
                                % middle (only one part)
                                matrix(j,j+2*k+1)=1/6*((deltax(j+2*k+1)+deltax(j+2*k))^3/...
                                    deltax(j+2*k+1)/deltax(j+2*k));
                                
                            end
                            
                        end
                        
                    end
                    
                end
                
            end
            
            break % Finish the for loop since the matrix has been filled
            
            
            % SECOND: x(n) is higher than x(i)+sigma
            % More delicate since the rows are not filled until the last column
            % Each row is filled individually
            
        elseif sigma+x(i)<x(n)
            
            counter=0;
            
            for j=i+1:n % Going to the right in the column
                
                if x(i)+sigma<=x(j) % When this loop is entered, j is the
                    ...first node in which x(i)+sigma<=x(j). Said in another
                        ...way, the first one outisde the integration domain.
                        
                    % Three types of rows can occur:
                    
                    if counter==0 % Row with only two elements
                        % Trapezoidal rule is applied
                        
                        chi=sigma-(x(j-1)-x(i));
                        
                        matrix(i,i)=chi-chi^2/(2*deltax(j-1));
                        matrix(i,j)=chi^2/(2*deltax(j-1));
                        
                    elseif counter==1 % Row with only three elements
                        % Simpson rule, and the value at x=chi is obtained
                        % by averaging j-1 and j
                        
                        chi=sigma-(x(j-1)-x(i));
                        
                        matrix(i,i)=1/6*((deltax(i)+chi)*(2*deltax(i)-...
                            chi)/deltax(i));
                        matrix(i,j-1)=1/6*((chi+deltax(j-2))^3/deltax(j-2)/chi)+...
                            1/6*(1-chi/deltax(j-1))*((chi+deltax(j-2))*(2*chi-...
                            deltax(j-2))/chi);
                        matrix(i,j)=1/6*(chi/deltax(j-1))*((chi+deltax(j-2))*(2*chi-...
                            deltax(j-2))/chi);
                        
                    elseif counter>1 % Row with more than three elements
                        
                        chi=sigma-(x(j-1)-x(i));
                        
                        % It is necessary again to see whether there are
                        % odd or even nodes
                        
                        % In that count the node at chi is included
                        
                        numbernodeschi=j-i+1;
                        
                        
                        if mod(numbernodeschi,2)==1 % Odd number of nodes: exact Simpson
                            
                            matrix(i,i)=1/6*((deltax(i)+deltax(i+1))*(2*deltax(i)-...
                                deltax(i+1))/deltax(i)); % First node in the diagonal
                            matrix(i,i+1)=1/6*((deltax(i)+deltax(i+1))^3/deltax(i)/deltax(i+1)); % Second node
                            
                            matrix(i,j-2)=1/6*((deltax(j-2)+chi)*(2*deltax(j-2)-...
                                chi)/deltax(j-2))+1/6*((deltax(j-3)+deltax(j-4))*(2*deltax(j-3)-...
                                deltax(j-4))/deltax(j-3)); % Antepenultimate node
                            matrix(i,j-1)=1/6*((deltax(j-2)+chi)^3/deltax(j-2)/chi)+...
                                1/6*(1-chi/deltax(j-1))*((chi+deltax(j-2))*(2*chi-...
                                deltax(j-2))/chi); % Penultimate node
                            matrix(i,j)=1/6*(chi/deltax(j-1))*((chi+deltax(j-2))*(2*chi-...
                                deltax(j-2))/chi); % Last node
                            
                            
                            % The rest of the nodes are computed in a loop
                            
                            
                            % This loop is entered if there are more than 5
                            % nodeschi
                            
                            if numbernodeschi>5
                                
                                for k=1:(numbernodeschi-1)/2-2 % Up to the last three nodes
                                    
                                    % Firstly, from the 3 points, the ones of the
                                    % extremes (sum of two parts)
                                    matrix(i,i+2*k)=1/6*((deltax(i+2*k-2)+deltax(i+2*k-1))...
                                        *(2*deltax(i+2*k-1)-deltax(i+2*k-2))/deltax(i+2*k-1))...
                                        +1/6*((deltax(i+2*k)+deltax(i+2*k+1))*(2*deltax(i+2*k)...
                                        -deltax(i+2*k+1))/deltax(i+2*k));
                                    
                                    % Secondly, from the 3 points, the ones in the
                                    % middle (only one part)
                                    matrix(i,i+2*k+1)=1/6*((deltax(i+2*k+1)+deltax(i+2*k))^3/...
                                        deltax(i+2*k+1)/deltax(i+2*k));
                                    
                                end
                                
                            end
                            
                            
                        elseif mod(numbernodeschi,2)==0 % Even number of nodes: exact Simpson
                            % until column j-1, and trapezoidal between columns j-1 and j
                            
                            % We add the part between j-1 and j with trapezoidal
                            % rule
                            
                            
                            matrix(i,i)=1/6*((deltax(i)+deltax(i+1))*(2*deltax(i)-...
                                deltax(i+1))/deltax(i)); % First node
                            matrix(i,i+1)=1/6*((deltax(i)+deltax(i+1))^3/deltax(i)/deltax(i+1)); % Second node
                            
                            matrix(i,j-1)=+chi-chi^2/(2*deltax(j-1))+1/6*((deltax(j-2)+deltax(j-3))*...
                                (2*deltax(j-2)-deltax(j-3))/deltax(j-2));% Penultimate node
                            matrix(i,j)=chi^2/(2*deltax(j-1)); % Last node
                            
                        
                            
                            % The rest of the nodes are computed in a loop
                            
                            % This loop is entered only if there are more
                            % than 4 nodeschi
                            
                            if numbernodeschi>4
                                
                                for k=1:(numbernodeschi-2)/2-1
                                    
                                    % Firstly, from the 3 points, the ones of the
                                    % extremes (sum of two parts)
                                    matrix(i,i+2*k)=1/6*((deltax(i+2*k-2)+deltax(i+2*k-1))...
                                        *(2*deltax(i+2*k-1)-deltax(i+2*k-2))/deltax(i+2*k-1))...
                                        +1/6*((deltax(i+2*k)+deltax(i+2*k+1))*(2*deltax(i+2*k)...
                                        -deltax(i+2*k+1))/deltax(i+2*k));
                                    
                                    % Secondly, from the 3 points, the ones in the
                                    % middle (only one part)
                                    matrix(i,i+2*k+1)=1/6*((deltax(i+2*k+1)+deltax(i+2*k))^3/...
                                        deltax(i+2*k+1)/deltax(i+2*k));
                                    
                                end
                                
                            end
                            
                            
                            
                        end
                        
                    end
                    
                    break % The column has been finished
                    
                    
                    
                else % This counter distinguishes the type of row
                    counter=counter+1;
                end
                

                
            end
            
        end
        
    end


end


end