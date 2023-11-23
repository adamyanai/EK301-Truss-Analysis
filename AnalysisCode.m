function trussanalysis(C, Sx, Sy, X, Y, L, appliedLoad)
    [joints, members] = size(C);

    % Initialize moments
    xMoment = zeros(joints, members);
    yMoment = zeros(joints, members);
    
    lengths = zeros(members);

    % Calculate moments
    for i = 1:members
        [jointOne, jointTwo] = findJoints(C, i);

        xMoment(jointOne, i) = (X(jointTwo) - X(jointOne)) / calculateDistance(X(jointOne), X(jointTwo), Y(jointOne), Y(jointTwo));
        xMoment(jointTwo, i) = -1 * xMoment(jointOne, i);

        yMoment(jointOne, i) = (Y(jointTwo) - Y(jointOne)) / calculateDistance(X(jointOne), X(jointTwo), Y(jointOne), Y(jointTwo));
        yMoment(jointTwo, i) = -1 * yMoment(jointOne, i);

        lengths(i) = calculateDistance(X(jointOne), X(jointTwo), Y(jointOne), Y(jointTwo));
    end

    % Combine moments, Sx, and Sy into matrix A
    A = [xMoment Sx; yMoment Sy];

    % Solve for forces using matrix inversion
    forces = inv(A) * L;

    % Compute member lengths
    memberLengths = sqrt(L(1:joints).^2 + L(joints+1:end).^2);

    % Extract member forces and compute reactions
    memberForces = forces(1:members);
    memberReactions = memberForces / memberLengths;
    
    % Display output information
    disp('\% EK301, Section A5, Group #: Kelly, Gaby, Adam 11/16/23xx.');
    fprintf('\nLoad: %.2f oz\n\n', appliedLoad);
    
    % Display member forces
    fprintf('Member forces in oz\n');
    for i = 1:members
        fprintf('m%d: %.3f (%s)\n', i, abs(memberForces(i)), (memberForces(i) > 0.1) * 'T' + (memberForces(i) < -0.1) * 'C');
    end
    
    % Display reaction forces
    fprintf('\nReaction forces in oz:\n');
    fprintf('Sx1: %.2f\n', forces(end-2));
    fprintf('Sy1: %.2f\n', forces(end-1));
    fprintf('Sy2: %.2f\n', forces(end));

    totalLength = sum(lengths);

    % Compute truss cost
    trussCost = calculateTrussCost(totalLength);
    fprintf('\nCost of truss: $%.2f\n', trussCost);
    
    % Compute theoretical max load/cost ratio
    theoreticalLoadCostRatio = appliedLoad / trussCost;
    fprintf('Theoretical max load/cost ratio in oz/$: %.4f\n\n', theoreticalLoadCostRatio);
    
    % Code for Wf calculation
    for i = 1 : members
        newlengths(i) = lengths(i, 1);
    end
        for i = 1 : members
            forcesnew(i) = forces(i);
        end
        Rm = forcesnew ./ appliedLoad;
        Pcrit = 3654.533*newlengths.^(-2.119);
        zeros1 = 0;
        j = 0;
        for i = 1 : members
            if Pcrit(i) == min(Pcrit)
                j = j + 1;
                zeros1(j) = i;
            end
        end
        k = 1;
        Rm = Rm';
        for i = 1 : members
            if (-0.1 < Rm(i) && Rm(i) < 0.1)
                Rm(i) = 1;
            end
        end
        Rm1 = abs(Rm);
        minimum = 0;
        minloc = 0;
        y = 1;
        for i = 1 : members
            if Rm1(i) == min(Rm1)
                minimum(y) = Rm1(i);
                minloc(y) = i;
                y = y +1;
            end
        end
        fprintf('The following member(s) will buckle first: %d\n', minloc)
        fprintf('\n')
        Wf = zeros(size(minimum));
        for i = 1 : members
            if ismember(i, minloc)
                Wf(k) = -min(Pcrit) ./ Rm1(i);
                k = k +1;
            end
        end

    % Displaying Wf
    fprintf('The maximum load (at which those members fail) is: %.2f\n', Wf(1))

    % Helper function to find joints connected by a member
    function [jointOne, jointTwo] = findJoints(connections, memberIndex)
        connectedJoints = find(connections(:, memberIndex) == 1);
        jointOne = connectedJoints(1);
        jointTwo = connectedJoints(2);
    end

    % Helper function to calculate distance between two points
    function distance = calculateDistance(x1, x2, y1, y2)
        distance = sqrt((x2 - x1)^2 + (y2 - y1)^2);
    end

    % Helper function to calculate truss cost
    function cost = calculateTrussCost(lengths)
        totalLength = sum(lengths);
        costPerJoint = 10;
        totalJoints = sum(joints);
        cost = totalLength + totalJoints * costPerJoint;
    end
end
