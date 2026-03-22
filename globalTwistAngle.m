% =========================================================================
% Global Twist Angle Computation for a Cyclic Peptide
% =========================================================================
%
% This code is intended to implement the global twist-angle formulation
% for cyclic peptides based on projected normal vectors and accumulated
% signed rotational angles.
%
% Input files:
%   - coordinates_C.csv   : coordinates of C atoms
%   - coordinates_CA.csv  : coordinates of C-alpha atoms
%   - coordinates_N.csv   : coordinates of N atoms
%
% Each CSV file is assumed to contain four columns:
%   Name, X, Y, Z
%
% Output:
%   - totalAngles : accumulated twist angle for each 14-residue block
%   - theta       : local signed twist angles from the final processed block
%
% Authors: Oh, M., Oh, C., and Tabassum, E.
% =========================================================================
clc; clear; close;
format long

% Read data (no headers in csv)
C = readtable('coordinates_C.csv', 'ReadVariableNames', false);
Ca = readtable('coordinates_CA.csv', 'ReadVariableNames', false);
N = readtable('coordinates_N.csv', 'ReadVariableNames', false);

% Give column names
C.Properties.VariableNames  = {'Name','X','Y','Z'};
Ca.Properties.VariableNames = {'Name','X','Y','Z'};
N.Properties.VariableNames = {'Name','X','Y','Z'};

% Extract coordinates
C_xyz  = C{:, {'X','Y','Z'}};
Ca_xyz = Ca{:, {'X','Y','Z'}};
N_xyz  = N{:, {'X','Y','Z'}};

% Set block size
blockSize = 14;

% Find number of coordinate rows
nRows = size(C_xyz, 1);

% Preallocate vector to store total twist angle for each block
totalAngles = nan(nRows/blockSize, 1);

% Loop over data block by block
idx = 1;
for j = 1:blockSize:nRows
    % Find indices for current block
    rows = j:min(j + blockSize - 1, nRows);

    % Extract coordinates for current block
    C_block = C_xyz(rows, :);
    Ca_block = Ca_xyz(rows, :);
    N_block = N_xyz(rows, :);
    
    % Compute normal vectors
    CCa = Ca_block - C_block;
    NCa = Ca_block - N_block;
    normalVectors = cross(CCa, NCa, 2);

    % Normalize normal vectors
    tol = 1e-12;
    den = vecnorm(normalVectors, 2, 2);
    den(den < tol) = tol;
    normalVectors = normalVectors ./ den; 

    % Construct backbone vectors
    CaCopied_block = Ca_block([2:end, 1], :);
    backboneVectors = CaCopied_block - Ca_block;

    % Normalize backbone vectors
    ben = vecnorm(backboneVectors, 2, 2);
    ben(ben < tol) = tol;
    backboneVectors = backboneVectors ./ ben;
    
    % Preallocate vector to store local twist angles
    num = size(normalVectors, 1)-1;
    theta = nan(num, 1);

    % Set large parameter used in tanh-based orientation adjustment
    lambda = 1e5;

    % Compute local signed twist angle between consecutive projected normal vectors
    for i = 1:num
        % Adjust orientation of next normal vector so that its direction is consistent with current one
        ni = normalVectors(i, :);
        nii = normalVectors(i+1, :);
        Nii = tanh(lambda*dot(ni, nii))*nii;
        
        % Normalize adjusted normal vector
        ten = vecnorm(Nii, 2);
        if ten < tol
            ten = tol;
        end
        Nii = Nii / ten;
    
        % Project both normal vectors onto plane perpendicular to backbone direction
        ui = backboneVectors(i, :);
        ai = ni - dot(ni, ui)*ui;
        bi = Nii - dot(Nii, ui)*ui;
        
        % Normalize projected vectors
        na = vecnorm(ai, 2);
        if na < tol
            na = tol;
        end
        ai = ai / na;
        nb = vecnorm(bi, 2);
        if nb < tol
            nb = tol;
        end
        bi = bi / nb;
        
        % Compute signed angle between ai and bi
        nume = dot(ui, cross(ai, bi));
        deno = dot(ai, bi);
        if abs(nume) < tol && abs(deno) < tol
            theta(i) = 0;
        else
            theta(i) = atan2(nume, deno);
        end
    end
    
    % Sum all local twist angles in current block to obtain total accumulated twist angle for that block
    totalAngles(idx) = sum(theta);

    % Move to next block
    idx = idx + 1;
end