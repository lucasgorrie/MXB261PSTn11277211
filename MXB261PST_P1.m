%% n11277211, P1 %%

% Rules as described in task specification: %
% (1) Particles can move S/W/E %
% (2) Western and Eastern occupied cells cannot be moved to %
% (3a) If at any point the South cell is occupied and a particle is to move %
% south, its journey has ended %
% (3b) If a particle has reached the bottom row, its journey has ended %
% (4) Cyclic W/E boundaries %

% PMF cases as described in task specification: %
% (i)   s =      w =       e = 1/3  %
% (ii)  s = 2/3, w = 1/6,  e = 1/6  %
% (iii) s = 3/5, w = 3/10, e = 1/10 %
% (iv)  s = 3/5, w = 1/10, e = 3/10 %

% Compute walk function %
% Inputs: %
% N ∈ Z^+ = number of particles % 
% P ∈ [1, 'rand'] = particle starting position, P = 1 := X_0 = 50, P = 'rand' := X_0 ~ ⌈U(1, 99)⌉ %
% PMF of dimensions [1 X 3] array of probabilities in the order of [ S, W, E ] %
% Outputs: %
% U ∈ R^2 of dimensions [N X 3], where the third column describes whether a particle has stopped moving %
% At the point of return, this third column will always be all ones. %
% I would like to parameterise the dimensions of the sim space but I'm not %
% sure if that is strictly allowed.                                        %
function U = walk(N, P, PMF)

    % Initialize particles %
    U = zeros(N, 3);
    U(:, 1) = 99;  % Top row

    walkFinished = false;

    % Initial columns %
    if isstring(P)
        if strcmp(P, "rand")
            for i = 1:N
                c = ceil(rand()*99); % X_0 ~ ⌈U(1, 99)⌉
                U(i, 2) = c;
            end
        end
    elseif P == 1
        U(:, 2) = 50;
    end

    CMF = cumsum(PMF); % CMF
    
    % Steps %
    while ~walkFinished
        for i = 1:N
    
            % Skip this particle if it is deemed to have stopped moving. %
            % This will need to be checked again later, but checking it  %
            % now will save a good amount of computation time            %
            if U(i, 3) == 1
                continue 
            end

            ParticleIsSouth = false;
            ParticleIsWest = false;
            ParticleIsEast = false;

            % List of indices of all particles at the same X %
            ParticlesOfSameX = find(U(:, 2) == U(i, 2));
            % List of indices of all particles at the same Y %
            ParticlesOfSameY = find(U(:, 1) == U(i, 1));
            % List of indices of all particles at X - 1 %
            ParticlesAtPrevX = find(U(:, 2) == U(i, 2) - 1);
            % List of indices of all particles at X + 1 %
            ParticlesAtNextX = find(U(:, 2) == U(i, 2) + 1);
            % List of indices of all particles at Y - 1 %
            ParticlesAtPrevY = find(U(:, 1) == U(i, 1) - 1);

            % >0 particles directly southward if set intersection does not return null set %
            if ~isempty(intersect(ParticlesOfSameX, ParticlesAtPrevY))
                ParticleIsSouth = true;
            end

            % Check for particles directly eastward and westward via set intersection %
            if ~isempty(intersect(ParticlesOfSameY, ParticlesAtPrevX))
                ParticleIsWest = true;
            end
            if ~isempty(intersect(ParticlesOfSameY, ParticlesAtNextX))
                ParticleIsEast = true;
            end

            % Direction to move %
            MoveIsValid = false;
            while ~MoveIsValid

                r = rand(); % r ~ U(0, 1)
                dir = find((r <= CMF), 1);

                if dir == 1
                    MoveIsValid = true;
                end

                % Enforce rule (2) %
                if dir == 2 && ~ParticleIsWest
                    MoveIsValid = true;
                end
                if dir == 3 && ~ParticleIsEast
                    MoveIsValid = true;
                end

            end
    
            % Enforce rule (3a) %
            if (dir == 1) && ParticleIsSouth
                U(i, 3) = 1;
            end
    
            % Check a second time if this particle should be skipped %
            if U(i, 3) == 1
                continue 
            end
    
            % Step %
            % South %
            if dir == 1
                U(i, 1) = U(i, 1) - 1;
            end

            % West %
            if dir == 2
                
                % Enforce rule (4) %
                if U(i, 2) == 1
                    U(i, 2) = 99;
                else
                    U(i, 2) = U(i, 2) - 1;
                end

            end

            % East %
            if dir == 3

                % Enforce rule (4) %
                if U(i, 2) == 99
                    U(i, 2) = 1;
                else
                    U(i, 2) = U(i, 2) + 1;
                end

            end

            % Check if particle has reached the bottom %
            % In other words, enforce rule (3b)        %
            if U(i, 1) == 1
                U(i, 3) = 1;
            end

        end

        % Check if all particles have finished their journey %
        if all(U(:, 3))
            walkFinished = true;
        end

    end
end

% PMF cases %
C1 = [1/3 1/3 1/3];
C2 = [2/3 1/6 1/6];
C3 = [3/5 3/10 1/10];
C4 = [3/5 1/10 3/10];

% Desired arguments %
P1 = 1;
P2 = "rand";
N1 = 100;
N2 = 200;

% Cases - this generally takes 5-10 seconds on my machine %
U1a = walk(N1, P1, C1);
U1b = walk(N1, P1, C2);
U1c = walk(N1, P1, C3);
U1d = walk(N1, P1, C4);

U2a = walk(N2, P1, C1);
U2b = walk(N2, P1, C2);
U2c = walk(N2, P1, C3);
U2d = walk(N2, P1, C4);

U3a = walk(N1, P2, C1);
U3b = walk(N1, P2, C2);
U3c = walk(N1, P2, C3);
U3d = walk(N1, P2, C4);

U4a = walk(N2, P2, C1);
U4b = walk(N2, P2, C2);
U4c = walk(N2, P2, C3);
U4d = walk(N2, P2, C4);

% Visualisations - histograms of *height* as dependent variable %
% P1 %
figure;
t = tiledlayout(2, 2);  % For 2x2 facet setup
title(t, sprintf("N = %d, P = %d", N1, P1))

nexttile
histogram(U1a(:, 1), 99);
title("case (i)")
xlabel("Y")
ylabel("N")

nexttile
histogram(U1b(:, 1), 99);
title("case (ii)")
xlabel("Y")
ylabel("N")

nexttile
histogram(U1c(:, 1), 99);
title("case (iii)")
xlabel("Y")
ylabel("N")

nexttile
histogram(U1d(:, 1), 99);
title("case (iv)")
xlabel("Y")
ylabel("N")

% P2 %
figure
t = tiledlayout(2, 2);  % For 2x2 facet setup
title(t, sprintf("N = %d, P = %d", N2, P1))

nexttile
histogram(U2a(:, 1), 99);
title("case (i)")
xlabel("Y")
ylabel("N")

nexttile
histogram(U2b(:, 1), 99);
title("case (ii)")
xlabel("Y")
ylabel("N")

nexttile
histogram(U2c(:, 1), 99);
title("case (iii)")
xlabel("Y")
ylabel("N")

nexttile
histogram(U2d(:, 1), 99);
title("case (iv)")
xlabel("Y")
ylabel("N")

% P3 %
figure
t = tiledlayout(2, 2);  % For 2x2 facet setup
title(t, sprintf("N = %d, P = %s", N1, P2))

nexttile
histogram(U3a(:, 1), 99);
title("case (i)")
xlabel("Y")
ylabel("N")

nexttile
histogram(U3b(:, 1), 99);
title("case (ii)")
xlabel("Y")
ylabel("N")

nexttile
histogram(U3c(:, 1), 99);
title("case (iii)")
xlabel("Y")
ylabel("N")

nexttile
histogram(U3d(:, 1), 99);
title("case (iv)")
xlabel("Y")
ylabel("N")

% P4 %
figure
t = tiledlayout(2, 2);  % For 2x2 facet setup
title(t, sprintf("N = %d, P = %s", N2, P2))

nexttile
histogram(U4a(:, 1), 99);
title("case (i)")
xlabel("Y")
ylabel("N")

nexttile
histogram(U4b(:, 1), 99);
title("case (ii)")
xlabel("Y")
ylabel("N")

nexttile
histogram(U4c(:, 1), 99);
title("case (iii)")
xlabel("Y")
ylabel("N")

nexttile
histogram(U4d(:, 1), 99);
title("case (iv)")
xlabel("Y")
ylabel("N")