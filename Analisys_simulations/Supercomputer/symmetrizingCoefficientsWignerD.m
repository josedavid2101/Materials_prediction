%-------------------------------------------------------------------------%
% Filename:  symmetrizingCoefficientsWignerD.m
% Author:    Oliver Johnson
% Date:      5/24/2023
%
% symmetrizingCoefficientsWignerD computes the "A" coefficients required to
% express the fully (crystal and sample) symmetrized Wigner-D functions in
% terms of the unsymmetrized Wigner-D functions. Note This is specific to
% the definition of the Wigner-D functions used in MTEX 5.9.0.
%
% Inputs:
%   S - The order(s) of the Wigner-D functions. If S is a scalar, the 
%       coefficients for just that order will be returned. If S is a 
%       vector, the coefficients for all orders contained in S will be
%       returned.
%   CS - a crystalSymmetry object defining the crystal symmetry point-group
%        with respect to which symmetrization will be performed.
%   SS - a specimenSymmetry object defining the sample symmetry point-group
%        with respect to which symmetrization will be performed.
%
% Outputs:
%   A_C - A numel(S) cell array where A_C{i} contains the 
%         M(S(i))-by-(2*S(i)+1) matrix of left-symmetrizing (crystal) 
%         coefficients for the Wigner-D functions of order S(i).
%   A_S - A numel(S) cell array where A_S{i} contains the 
%         (2*S(i)+1)-by-N(S(i)) matrix of right-symmetrizing (sample) 
%         coefficients for the Wigner-D functions of order S(i).
%
% Given a matrix of unsymmetrized Wigner-D functions of order s (evaluated 
% at some point), the correpsonding symmetrized Wigner-D matrix is obtained 
% by
%
%   D_sym = A_C{s+1}*D*A_S{s+1};
%
% Given a matrix of coefficients for a harmonic expansion in the
% unsymmetrized Wigner-D funcitons of order s, the matrix for the 
% corresponding coefficients of a harmonic expansion in the symmetrized 
% Wigner-D functions of the same order is obtained by
%
%   fhat_sym = conj(A_C{s+1})*fhat*conj(A_S{s+1});
% 
% The reverse transformation is 
%
%   fhat = (A_C{s+1}.')*fhat_sym*(A_S{s+1}.');
%
% NOTE: For use with MTEX 5.9.0
%-------------------------------------------------------------------------%
function [A_C,A_S] = symmetrizingCoefficientsWignerD(S,CS,SS)

%% Extract the group generators

R_C = getGroupGenerators(CS);
R_S = getGroupGenerators(SS);

%% Compute symmetrizing coefficients

i = 1;
for s = S(:).'

    % get identity matrix for order s
    I_s = eye(2*s+1);

    % Build B matrices from Wigner-D matrices and identity
    B_C = cell(0,1);
    for j = 1:prod(size(R_C))
        B_C{j} = sqrt(2*s+1)*Wigner_D(s,R_C(j))-sqrt(2*s+1)*I_s; % note that D_s = sqrt(2*s+1)*Wigner_D(s,R_C(j)), so this is B_C{j} = D_s(R_C)-sqrt(2*s+1)*I_s
    end
    B_CT = cellfun(@(X) permute(X,[2,1]),B_C,'UniformOutput',false); % transpose each submatrix before concatenation
    B_CT = cat(1,B_CT{:});

    B_S = cell(0,1);
    for j = 1:prod(size(R_S))
        B_S{j} = sqrt(2*s+1)*Wigner_D(s,R_S(j))-sqrt(2*s+1)*I_s; % note that D_s = sqrt(2*s+1)*Wigner_D(s,R_S(j)), so this is B_S{j} = D_s(R_S)-sqrt(2*s+1)*I_s
    end
    B_S = cat(1,B_S{:});

    % solve system of equations
    A_C{i} = null(B_CT).'; % note the untranspose so A_C{i} forms a row-basis
    A_S{i} = null(B_S); % note there is NO transpose so A_S{i} forms a column-basis

    i = i+1;
end

end

function R = getGroupGenerators(sym)

switch sym.id
    case 1 % C1 (triclinic)
        R = rotation.byAxisAngle(vector3d(0,0,1),0*degree);
    case 16 % D2h (orthorhombic)
        R = [-rotation.byAxisAngle(vector3d(0,0,1),0*degree);...
              rotation.byAxisAngle(vector3d(0,0,1),180*degree);...
              rotation.byAxisAngle(vector3d(0,1,0),180*degree)];
    case 45 % Oh
        R = [-rotation.byAxisAngle(vector3d(0,0,1),0*degree);...
              rotation.byAxisAngle(vector3d(0,0,1),180*degree);...
              rotation.byAxisAngle(vector3d(0,1,0),180*degree);...
              rotation.byAxisAngle(vector3d(1,1,1),120*degree);...
              rotation.byAxisAngle(vector3d(1,1,0),180*degree)];
    otherwise
        error(['Symmetry ',sym.pointGroup,' has not yet been implemented.'])
end

end