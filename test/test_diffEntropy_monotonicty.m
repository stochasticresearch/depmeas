%**************************************************************************
%* 
%* Copyright (C) 2017  Kiran Karra <kiran.karra@gmail.com>
%*
%* This program is free software: you can redistribute it and/or modify
%* it under the terms of the GNU General Public License as published by
%* the Free Software Foundation, either version 3 of the License, or
%* (at your option) any later version.
%*
%* This program is distributed in the hope that it will be useful,
%* but WITHOUT ANY WARRANTY; without even the implied warranty of
%* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%* GNU General Public License for more details.
%*
%* You should have received a copy of the GNU General Public License
%* along with this program.  If not, see <http://www.gnu.org/licenses/>.
%*
%**************************************************************************

%% Experiments to show the relationship between
% Spearman's Rho, Differential Entropy, Kendall's Tau & tau_{KL}
clear;
clc;

M = 10000;

sRhoVec = 0.01:.01:0.99;
kTauVec = 0.01:.01:0.99;

entropyVecSRho = zeros(4,length(sRhoVec));
entropyVecKTau = zeros(4,length(kTauVec));

% compute sRho vs entropy
for ii=1:length(sRhoVec)
    sRho = sRhoVec(ii);
    kTau = kTauVec(ii);
    
    copulaType = 'Gaussian'; storeIdx = 1;
    rho = copulaparam(copulaType,sRho,'type','Spearman');
    tau = copulaparam(copulaType,kTau,'type','Kendall');
    % simulate samples
    Xsrho = copularnd(copulaType,rho,M);
    Xktau = copularnd(copulaType,tau,M);
    % compute entropy
    entropyVecSRho(storeIdx,ii) = sum(log(copulapdf(copulaType,Xsrho,rho)))/M;
    entropyVecKTau(storeIdx,ii) = sum(log(copulapdf(copulaType,Xktau,rho)))/M;
    
    copulaType = 'Frank'; storeIdx = 2;
    rho = copulaparam(copulaType,sRho,'type','Spearman');
    tau = copulaparam(copulaType,kTau,'type','Kendall');
    % simulate samples
    Xsrho = copularnd(copulaType,rho,M);
    Xktau = copularnd(copulaType,tau,M);
    % compute entropy
    entropyVecSRho(storeIdx,ii) = sum(log(copulapdf(copulaType,Xsrho,rho)))/M;
    entropyVecKTau(storeIdx,ii) = sum(log(copulapdf(copulaType,Xktau,rho)))/M;
    
    copulaType = 'Gumbel'; storeIdx = 3;
    rho = copulaparam(copulaType,sRho,'type','Spearman');
    tau = copulaparam(copulaType,kTau,'type','Kendall');
    % simulate samples
    Xsrho = copularnd(copulaType,rho,M);
    Xktau = copularnd(copulaType,tau,M);
    % compute entropy
    entropyVecSRho(storeIdx,ii) = sum(log(copulapdf(copulaType,Xsrho,rho)))/M;
    entropyVecKTau(storeIdx,ii) = sum(log(copulapdf(copulaType,Xktau,rho)))/M;
    
    copulaType = 'Clayton'; storeIdx = 4;
    rho = copulaparam(copulaType,sRho,'type','Spearman');
    tau = copulaparam(copulaType,kTau,'type','Kendall');
    % simulate samples
    Xsrho = copularnd(copulaType,rho,M);
    Xktau = copularnd(copulaType,tau,M);
    % compute entropy
    entropyVecSRho(storeIdx,ii) = sum(log(copulapdf(copulaType,Xsrho,rho)))/M;
    entropyVecKTau(storeIdx,ii) = sum(log(copulapdf(copulaType,Xktau,rho)))/M;
   
end

figure;
% plot
subplot(2,2,1);
plot(sRhoVec,entropyVecSRho(1,:),kTauVec,entropyVecKTau(1,:));
title('Gaussian Copula');
legend('\rho','\tau'); grid on;
xlabel('\kappa'); ylabel('-H')

subplot(2,2,2);
plot(sRhoVec,entropyVecSRho(2,:),kTauVec,entropyVecKTau(2,:));
title('Frank Copula');
legend('\rho','\tau'); grid on;
xlabel('\kappa'); ylabel('-H')

subplot(2,2,3);
plot(sRhoVec,entropyVecSRho(3,:),kTauVec,entropyVecKTau(3,:));
title('Gumbel Copula');
legend('\rho','\tau'); grid on;
xlabel('\kappa'); ylabel('-H')

subplot(2,2,4);
plot(sRhoVec,entropyVecSRho(4,:),kTauVec,entropyVecKTau(4,:));
title('Clayton Copula');
legend('\rho','\tau'); grid on;
xlabel('\kappa'); ylabel('-H')

% Plot all of the different copula's on one plot to reproduce Fig.2 in SMS
% paper
figure;
plot(sRhoVec,entropyVecSRho(1,:), ...
     sRhoVec,entropyVecSRho(2,:), ...
     sRhoVec,entropyVecSRho(3,:), ...
     sRhoVec,entropyVecSRho(4,:) );
legend('Gaussian','Frank','Gumbel','Clayton')
xlabel('\rho');
ylabel('E(-H)')
grid on;

% plot just the priors
p_prior = zeros(3,length(sRhoVec));

for ii=1:length(sRhoVec)
    rho_s = sRhoVec(ii);
    p_prior(1,ii) = hcbn_prior('Gaussian',rho_s);
    p_prior(2,ii) = hcbn_prior('Gumbel',rho_s);
    p_prior(3,ii) = hcbn_prior('Clayton',rho_s);
end

subplot(1,2,1);
plot(sRhoVec,p_prior);
grid on;
xlabel('rho_s');
ylabel('Prior');

% Compute & Plot the Posteriors
p_posterior = zeros(3,length(sRhoVec));
p_posterior(1,:) = entropyVecSRho(1,:) + log(p_prior(1,:));
p_posterior(2,:) = entropyVecSRho(3,:) + log(p_prior(2,:));
p_posterior(3,:) = entropyVecSRho(4,:) + log(p_prior(3,:));
subplot(1,2,2);
plot(sRhoVec,p_posterior);
grid on;
xlabel('rho_s');
ylabel('Posterior');
legend('Gaussian', 'Gumbel', 'Clayton');
