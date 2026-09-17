function [r2_old, bhat_selezionato, Hselezionato, indici_selezionati] = stepwise(regressoricel, dati, parametro)
    % FORWARD STEPWISE REGRESSION
    % Iteratively selects the best regressors based on the adjusted R-squared metric.
    %
    % Inputs:
    %   regressoricel - Cell array containing candidate regressor matrices
    %   dati          - Dataset (used to determine the number of observations)
    %   parametro     - Dependent variable/parameter to model
    %
    % Outputs:
    %   r2_old              - Final adjusted R-squared value of the model
    %   bhat_selezionato    - Estimated regression coefficients for the selected model
    %   Hselezionato        - Matrix of selected regressors
    %   indici_selezionati  - Indices of the selected regressors from the cell array

    Hselezionato = [];
    indici_selezionati = [];
    k = 1;
    r2_new = 0;
    r2_old = 0;
    
    % Continue adding regressors as long as the adjusted R-squared improves
    while r2_new >= r2_old
        r2_old = r2_new; 
        
        for i = 1:length(regressoricel)
            % Check if the current regressor has not been selected yet
            if isempty(find(indici_selezionati == i, 1)) 
                % Build a candidate regressor matrix by adding the current option
                Hcandidato = [Hselezionato regressoricel{i}]; 
                
                % Evaluate the model using the external regression function
                [r2adj(i), Ftest(i), bhat, ~] = regressione(length(dati), size(Hcandidato, 2), Hcandidato, parametro);
            end
        end
        
        % Find the maximum adjusted R-squared among the evaluated candidates
        [r2_new, indicemax] = max(r2adj);
        
        % If the new adjusted R-squared improves, update the model structure
        if r2_new >= r2_old
            % Add the best performing regressor to the selected matrix
            Hselezionato = [Hselezionato regressoricel{indicemax}];
            
            % Update the record of selected indices
            indici_selezionati = [indici_selezionati indicemax];
            
            % Re-estimate coefficients and fitted values for the updated model
            [~, ~, bhat, yhat] = regressione(length(dati), size(Hselezionato, 2), Hselezionato, parametro);
            bhat_selezionato = bhat;
            yhat_selezionato = yhat;
        end
        
        k = k + 1;
        clear r2adj;
    end
end