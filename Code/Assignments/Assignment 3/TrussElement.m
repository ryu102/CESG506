function [Fe, Ke] = TrussElement(Xi, Xj, ui, uj, EA)

    % unmdeformed
    len0V = Xj - Xi;
    len0 = sqrt(len0V'*len0V);

    % deformed
    lenV = Xj + uj - Xi - ui;
    len = sqrt(lenV'*lenV);
    nvec = lenV/len;

    % find strain
    lambda = len/len0;
    strain = log(lambda);

    % internal force vector
    f = EA * strain;
    Fe = f * nvec;

    % tangent stiffness matrix
    NtensorN = nvec*nvec';
    Ke = EA/len * NtensorN + f/len * (eye(length(ui)) - NtensorN);

end
