function trackdata = cluster_contour_tracker(videodata,template_struct,numIterations,coilIntensityCorrectionOn,coilSensitivityMatFileName,newtonMethodOn,contourCleanUpOn,plotOn,notifyOn, frames, parallelOn, workers)

N=size(videodata,1);

if isempty(frames)
    frames = 1:size(videodata,3);
end;

%if notifyOn, 
fprintf('spanContourTracker - processing\n'); 
%end

kMatrixFull2DFT = repmat(((-N/2):(N/2-1))/N,N,1); 
kMatrixFull2DFT = kMatrixFull2DFT + 1i*(kMatrixFull2DFT.');
wMatrixFull2DFT = ones(N,N)/(N^2);

kMatrixFull2DFT_ = kMatrixFull2DFT(:);
wMatrixFull2DFT_ = wMatrixFull2DFT(:);

if coilIntensityCorrectionOn 
    load(coilSensitivityMatFileName,'magnitudeCoilSensMap'); 
else
    magnitudeCoilSensMap=ones(N);
end

if parallelOn
    videodata_reordered = videodata(:,:,frames);
    p=parpool('TorqueProfile1',workers);
    addAttachedFiles(p,{'gridkb.m','calckbkernel.m','kb.m','ift.m','interpft2.m'});
    parfor indframe=1:length(frames)
        
        %f = frames(indframe);
        
        %fprintf('Starting frame %i of %i\n',f,length(frames));
        
        img = videodata_reordered(:,:,indframe);
        
        amax=max(max(img));
        amin=min(min(img));
        
        imObserved=flip(mat2gray(img,[amin amax]),1);
        
        %do coil intensity correction
        if coilIntensityCorrectionOn
            imObservedCorrected = imObserved ./ magnitudeCoilSensMap;
        else
            imObservedCorrected = imObserved;
        end
        
        %go back into spatial frequency domain
        fObserved = fftshift(fft2(fftshift(imObservedCorrected)));
        %fDataType = 'full2DFT';
        
        %choose template
        
        trackingScoreEvolutionArray=1000*ones(1,size(template_struct,2));
        
        for i=1:size(template_struct,2);
            
            model=template_struct(i).model;
            
            for s=1:size(model.segment,2)
                model.segment{s}.v = model.segment{s}.v * N;
            end
            
            %add other items to the model data structure
            model.wMatrixSqrt = sqrt(wMatrixFull2DFT_);
            model.mBar        = fObserved(:) .* model.wMatrixSqrt; % TODO: THIS IS JUST DIVISION BY N
            model.kMatrix     = kMatrixFull2DFT_;
            model.N           = N;
            model = updatePsiScore(model);
            model = updateMeanScore(model);
            model = updateGradient(model);
            
            %run the segmentation for the current image
            [~, trackingScoreEvolution] = processOneImage(model,[10 20 19 1],newtonMethodOn,contourCleanUpOn,plotOn,notifyOn);
            trackingScoreEvolutionArray(i)=trackingScoreEvolution(1,end);
            
        end;
        
        [~,indBest] = min(trackingScoreEvolutionArray);
        
        %fprintf('Frame: %i, Selected template: %i\n',f,indBest);
        
        model=template_struct(indBest).model;
        
        for s=1:size(model.segment,2)
            model.segment{s}.v = model.segment{s}.v * N;
        end
        
        
        %add other items to the model data structure
        model.wMatrixSqrt = sqrt(wMatrixFull2DFT(:));
        model.mBar        = fObserved(:) .* model.wMatrixSqrt; % TODO: THIS IS JUST DIVISION BY N
        model.kMatrix     = kMatrixFull2DFT(:);
        model.N           = N;
       % model.templateNo  = indBest;
        model = updatePsiScore(model);
        model = updateMeanScore(model);
        model = updateGradient(model);
        
        %run the segmentation for the current image
        [modelEvolution, trackingScoreEvolution] = processOneImage(model,numIterations,newtonMethodOn,contourCleanUpOn,plotOn,notifyOn);
        
        
        
        %save the tracking score and the model evolution in the data structure
        %         frame{indframe}.modelEvolution          = modelEvolution;
        %         frame{indframe}.image                   = imObserved;
        %         frame{indframe}.imageCorrected          = imObservedCorrected;
        %         frame{indframe}.frameNo = f;
        %save the result
        
        trackdata{indframe}.contours = modelEvolution{5};
        trackdata{indframe}.frameNo  = frames(indframe);
        trackdata{indframe}.template = indBest;
        
        %fprintf('Frame: %i, done\n',f);
            
    end
                
    %delete(gcp('nocreate'));
    delete(p);
        
else
    
    for indframe=1:length(frames)
        
        f = frames(indframe);
        fprintf('Frame: %i\n',f);
        
        img = videodata(:,:,f);
        
        amax = max(max(img));
        amin = min(min(img));
        
        imObserved=flip(mat2gray(img,[amin amax]));
        
        %do coil intensity correction
        if coilIntensityCorrectionOn
            imObservedCorrected = imObserved ./ magnitudeCoilSensMap;
        else
            imObservedCorrected = imObserved;
        end
        
        %go back into spatial frequency domain
        fObserved = fftshift(fft2(fftshift(imObservedCorrected)));
        fDataType = 'full2DFT';
        
        %load the model template from the templateMatFileName
        
        
        %choose template
        
        trackingScoreEvolutionArray=1000*ones(1,size(template_struct,2));
        
        for i=1:size(template_struct,2);
            
            model=template_struct(i).model;
            
            for s=1:size(model.segment,2)
                model.segment{s}.v = model.segment{s}.v * N;
            end
            
            
            %add other items to the model data structure
            model.wMatrixSqrt = sqrt(wMatrixFull2DFT(:));
            model.mBar        = fObserved(:) .* model.wMatrixSqrt; % TODO: THIS IS JUST DIVISION BY N
            model.kMatrix     = kMatrixFull2DFT(:);
            model.N           = N;
            model = updatePsiScore(model);
            model = updateMeanScore(model);
            model = updateGradient(model);
            
            %run the segmentation for the current image
            [~, trackingScoreEvolution] = processOneImage(model,[10 20 19 1],newtonMethodOn,contourCleanUpOn,plotOn,notifyOn);
            trackingScoreEvolutionArray(i)=trackingScoreEvolution(1,end);
            
        end;
        
        [~,indBest] = min(trackingScoreEvolutionArray);
        
        
        model=template_struct(indBest).model;
        
        for s=1:size(model.segment,2)
            model.segment{s}.v = model.segment{s}.v * N;
        end
        
        
        %add other items to the model data structure
        model.wMatrixSqrt = sqrt(wMatrixFull2DFT(:));
        model.mBar        = fObserved(:) .* model.wMatrixSqrt; % TODO: THIS IS JUST DIVISION BY N
        model.kMatrix     = kMatrixFull2DFT(:);
        model.N           = N;
        model = updatePsiScore(model);
        model = updateMeanScore(model);
        model = updateGradient(model);
        
        %run the segmentation for the current image
        [modelEvolution, trackingScoreEvolution] = processOneImage(model,numIterations,newtonMethodOn,contourCleanUpOn,plotOn,notifyOn);
        
        
        %save the tracking score and the model evolution in the data structure
        %frame{indframe}.modelEvolution         = modelEvolution;
        %frame{indframe}.frameNo = f;
        trackdata{indframe}.contours = modelEvolution{5};
        trackdata{indframe}.frameNo  = f;
    
    end
    
end

%data.frame=frame;
%save(trackingMatFileName,'data');
clear trackingScoreEvolution;
return

function [modelEvolution, trackingScoreEvolution] = processOneImage(model,numIterations,newtonMethodOn,contourCleanUpOn,plotOn,notifyOn)
modelEvolution{1}      = model;
trackingScoreEvolution = [inf ; 0 ; 1 ; datenum(now)];
%**********
% stage 1-3
%**********
for optLevel=1:3
    iter = 0;
    epsilonMultiplier = 1;
    while iter<numIterations(optLevel)
        if notifyOn disp(sprintf('\nLevel %d, Iteration %d, epsilonMultiplier %1.3f',optLevel,iter,epsilonMultiplier)); end
        modelNew = gradientDescentStep(model,optLevel,epsilonMultiplier,plotOn);
        if (modelNew.score < trackingScoreEvolution(1,end))
            model = modelNew;
            if contourCleanUpOn
                model = cleanUpIntersections(model);
            end
            trackingScoreEvolution = [trackingScoreEvolution [modelNew.score                ; optLevel ; epsilonMultiplier ; datenum(now)]];
            epsilonMultiplier = 1;
        else
            trackingScoreEvolution = [trackingScoreEvolution [trackingScoreEvolution(1,end) ; optLevel ; epsilonMultiplier ; datenum(now)]];
            epsilonMultiplier = epsilonMultiplier/2;
            if epsilonMultiplier < 0.0005; break; end;
        end
        if plotOn
            figure(1); subplot(1,2,2);
            plot(trackingScoreEvolution(1,:));
            v=axis;
            axis([1 sum(numIterations) 0 1.05*v(4)]);
        end
        iter = iter + 1;
    end
    modelEvolution{optLevel+1} = model;
end
%**********
% stage 4
%**********
model = modelEvolution{4};
optLevel = 4;
if newtonMethodOn
    %use the optimization toolbox
    if notifyOn
        options = optimset('Diagnostics','off','Display','iter','LargeScale','on','PrecondBandWidth',Inf,'GradObj','on','Hessian','on','MaxIter',1);
    else
        options = optimset('Diagnostics','off','Display','off' ,'LargeScale','on','PrecondBandWidth',Inf,'GradObj','on','Hessian','on','MaxIter',1);
    end
    iter = 0;
    while iter<numIterations(optLevel)
        if notifyOn disp(sprintf('\nLevel %d, Iteration %d, using Newton Method',optLevel,iter)); end
        p = zeros(2*(size(model.segment{1}.v,1) + size(model.segment{2}.v,1) + size(model.segment{3}.v,1)),1);
        [p model.score] = fminunc(@(p) objectiveFunctionStage4(p,model),p,options );
        %put the result of the toolbox function into the model data structure +
        %update the model, incl. means and score
        for s=1:(size(model.segment,2)-1)
            Ns = length(model.segment{s}.v);
            model.segment{s}.v = model.segment{s}.v + [p(1:2:2*Ns) p(2:2:2*Ns)];
            p(1:2*Ns) = [];
        end
        if contourCleanUpOn
            model = cleanUpIntersections(model);
        end
        trackingScoreEvolution = [trackingScoreEvolution [model.score ; optLevel ; epsilonMultiplier ; datenum(now)]];
        if plotOn
            figure(1); subplot(1,2,2)
            plot(trackingScoreEvolution(1,:));
            v=axis;
            axis([1 sum(numIterations) 0 1.05*v(4)]);
        end
        iter = iter + 1;
    end
else
    %use Gradient descent
    epsilonMultiplier = 1;
    iter = 0;
    bestmodel = model;
    while iter<numIterations(optLevel)
        if notifyOn disp(sprintf('\nLevel %d, Iteration %d, epsilonMultiplier %1.3f',optLevel,iter,epsilonMultiplier)); end
        model = gradientDescentStep(model,optLevel,epsilonMultiplier,plotOn);
        if contourCleanUpOn
            model = cleanUpIntersections(model);
        end
                if (model.score < trackingScoreEvolution(1,end))
            bestmodel = model;
        end;
        trackingScoreEvolution = [trackingScoreEvolution [model.score                ; optLevel ; epsilonMultiplier ; datenum(now)]];
        %             epsilonMultiplier = 1;
        %         else
        %             trackingScoreEvolution = [trackingScoreEvolution [trackingScoreEvolution(1,end) ; optLevel ; epsilonMultiplier ; datenum(now)]];
        %             epsilonMultiplier = epsilonMultiplier/2;
        %             if epsilonMultiplier < 0.0005; break; end;
        %         end
        %         model = gradientDescentStep(model,optLevel,epsilonMultiplier,plotOn);
        %         if contourCleanUpOn
        %             model = cleanUpIntersections(model);
        %         end
        %         trackingScoreEvolution = [trackingScoreEvolution [model.score ; optLevel ; epsilonMultiplier ; datenum(now)]];
        if plotOn
            figure(1); subplot(1,2,2)
            plot(trackingScoreEvolution(1,:));
            v=axis;
            axis([1 sum(numIterations) 0 1.05*v(4)]);
        end
        iter = iter + 1;
    end
end
model = bestmodel;
model = updatePsiScore(model);
model = updateMeanScore(model);
trackingScoreEvolution = [trackingScoreEvolution [model.score ; optLevel ; epsilonMultiplier ; datenum(now)]];
modelEvolution{5} = model;

%********************************************************************
%clean up the returned data structure to eliminate excessive memory
%requirements; retain only the vertex vectors and the means
%********************************************************************
modelEvolutionTemp = [];
for m=1:size(modelEvolution,2)
    for s=1:size(modelEvolution{m}.segment,2)
        modelEvolutionTemp{m}.segment{s}.v  = modelEvolution{m}.segment{s}.v;
        modelEvolutionTemp{m}.segment{s}.i  = modelEvolution{m}.segment{s}.i;
        modelEvolutionTemp{m}.segment{s}.mu = modelEvolution{m}.segment{s}.mu;
    end
end
modelEvolution = modelEvolutionTemp;
return

function [J gradJ HessianJ] = objectiveFunctionStage4(p,model)
%extract the vertices from the parameter vector using the size information
%from the global variable
for s=1:(size(model.segment,2)-1)
    Ns = length(model.segment{s}.v);
    model.segment{s}.v = model.segment{s}.v + [p(1:2:2*Ns) p(2:2:2*Ns)];
    p(1:2*Ns) = [];
end
% update the geomtry and the score of the model
% estimate the gradient and the Hessian 
model = updatePsiScore(model);
model = updateMeanScore(model);
model = updateGradientHessian(model);
J        = model.score;
gradJ    = model.grad;
HessianJ = model.Hessian;
return

function model = updateGradientHessian(model);
%get the non-zero column of the derivative of Psi for all segments and
%concatenate in to one vector; then weight with sqrt(w)-matrix
numSegments                              = size(model.segment,2);
numParametersPerSegment                  = [];
parameterNumberToSegmentNumberConversion = [];
dPsi_dp_nonZeroColumn                    = [];
for s=1:(numSegments-1)
    %determine the number of parameters in the current segment and make a
    %conversion LUT for parameter-number - to - segment-number
    numParametersPerSegment(s) = 2*size(model.segment{s}.v,1);
    parameterNumberToSegmentNumberConversion = ...
        [parameterNumberToSegmentNumberConversion s*ones(1,numParametersPerSegment(s))];
    %get the first derivative of the shape function of the current segment
    [dfPoly_dxS dfPoly_dyS] = spanPoly2FourierDerivative(model.segment{s}.v(:,1),model.segment{s}.v(:,2),real(model.kMatrix(:)),imag(model.kMatrix(:)));
    dPsi_dp_nonZeroColumn  = [dPsi_dp_nonZeroColumn ...
        reshape([squeeze(dfPoly_dxS) ; squeeze(dfPoly_dyS)],[length(model.wMatrixSqrt(:)),numParametersPerSegment(s)])];
end
numParameters = sum(numParametersPerSegment);
dPsi_dp_nonZeroColumn = repmat(model.wMatrixSqrt(:),[1,numParameters]).*dPsi_dp_nonZeroColumn;
%compute dA_dp and dinvA_dp for all parameters
dPsiH_dp_nonZeroColmns_times_Psi  = dPsi_dp_nonZeroColumn'*model.Psi;
invA_dA_dp_stackedInColumn = [];
dinvA_dp_stackedInColumn   = [];
for ii=1:size(dPsi_dp_nonZeroColumn,2)
    dA_dp = zeros(numSegments);
    dA_dp(parameterNumberToSegmentNumberConversion(ii),1:numSegments) = dPsiH_dp_nonZeroColmns_times_Psi(ii,1:numSegments);
    dA_dp = dA_dp + dA_dp';
    dinvA_dp =  -model.invPsiPsi * dA_dp * model.invPsiPsi;
    dinvA_dp_stackedInColumn   = [dinvA_dp_stackedInColumn   ; dinvA_dp             ];
    invA_dA_dp_stackedInColumn = [invA_dA_dp_stackedInColumn ; model.invPsiPsi*dA_dp];
end
%re-compute the means and other intermediate constants
PsiHtimesMbar                               = model.Psi'*model.mBar;
mu                                          = model.invPsiPsi*PsiHtimesMbar;
mBarminusPsimu                              = model.mBar - model.Psi * mu;
mBarminusPsiMuH_times_Psi                   = mBarminusPsimu' * model.Psi;
dPsiH_dpi_dPsi_dpj                          = dPsi_dp_nonZeroColumn'*dPsi_dp_nonZeroColumn;
invA_dA_dp_dinvA_dp_stackedInColumn         = invA_dA_dp_stackedInColumn * dinvA_dp_stackedInColumn';
dPsiH_dp_nonZeroColumn_times_mBar           = dPsi_dp_nonZeroColumn'*model.mBar;
mBarminusPsiMuH_times_dPsi_dp_nonZeroColumn = mBarminusPsimu' * dPsi_dp_nonZeroColumn;
% recompute the gradient of the mean and the gradient of the objective
% function
model.gradMu = reshape(dinvA_dp_stackedInColumn*PsiHtimesMbar,numSegments,numParameters) + ...
    model.invPsiPsi(:,parameterNumberToSegmentNumberConversion).*repmat(dPsiH_dp_nonZeroColumn_times_mBar.',[numSegments,1]);
model.grad = -2*real(mBarminusPsiMuH_times_Psi*model.gradMu + mBarminusPsiMuH_times_dPsi_dp_nonZeroColumn.*mu(parameterNumberToSegmentNumberConversion).');
%compute the first order terms
HessianFirstOrderTerm = [];
for ii=1:numParameters
    for jj=ii:numParameters
        %contribution from second derivative of Mu
            %Term1
            temp1_Term1 = invA_dA_dp_dinvA_dp_stackedInColumn(((ii-1)*numSegments + 1) :(ii*numSegments),((jj-1)*numSegments + 1) :(jj*numSegments));
            HessianMuFirstOrderTerm1 = -(temp1_Term1 + temp1_Term1')*PsiHtimesMbar;
            %Term2
            HessianMuFirstOrderTerm2 =...
                -dinvA_dp_stackedInColumn(((ii-1)*numSegments+1) : (ii*numSegments) , parameterNumberToSegmentNumberConversion(jj)) * dPsiH_dp_nonZeroColumn_times_mBar(jj)...
                +dinvA_dp_stackedInColumn(((jj-1)*numSegments+1) : (jj*numSegments) , parameterNumberToSegmentNumberConversion(ii)) * dPsiH_dp_nonZeroColumn_times_mBar(ii);
            %Term3
            temp1_Term3 = zeros(numSegments);
            temp1_Term3(parameterNumberToSegmentNumberConversion(ii),parameterNumberToSegmentNumberConversion(jj)) = dPsiH_dpi_dPsi_dpj(ii,jj);
            temp1_Term3 = temp1_Term3 + temp1_Term3';
            HessianMuFirstOrderTerm3 = -model.invPsiPsi * temp1_Term3 * mu;
            %Put them together
            HessianFirstOrderTerms_1_2_3(ii,jj) = mBarminusPsiMuH_times_Psi * (HessianMuFirstOrderTerm1 + HessianMuFirstOrderTerm2 + HessianMuFirstOrderTerm3);
        %contribution from gradient-mu terms
            HessianFirstOrderTerm4(ii,jj) = ...
                mBarminusPsiMuH_times_dPsi_dp_nonZeroColumn(ii) * model.gradMu(parameterNumberToSegmentNumberConversion(ii),jj) + ...
                mBarminusPsiMuH_times_dPsi_dp_nonZeroColumn(jj) * model.gradMu(parameterNumberToSegmentNumberConversion(jj),ii);
    end
end
HessianFirstOrderTerms_1_2_3_4 = HessianFirstOrderTerms_1_2_3 + HessianFirstOrderTerm4;
HessianFirstOrderTerms_1_2_3_4 = HessianFirstOrderTerms_1_2_3_4 + HessianFirstOrderTerms_1_2_3_4' - diag(diag(HessianFirstOrderTerms_1_2_3_4));
% the second outside term
temp1_Term5 = repmat(mu(parameterNumberToSegmentNumberConversion).',[length(model.wMatrixSqrt(:)),1]) .* dPsi_dp_nonZeroColumn;
temp1_Term5 = temp1_Term5 + model.Psi*model.gradMu;
HessianFirstOrderTerm5 = temp1_Term5.' * conj(temp1_Term5);
% now put everything together    
HessianFirstOrderTerm = HessianFirstOrderTerms_1_2_3_4 - HessianFirstOrderTerm5;
%now compute the second order terms of the Hessian of Mu 
wMatrixSqrtPsi                            = repmat(model.wMatrixSqrt(:),[1,numSegments]) .* model.Psi;
wMatrixSqrtmBar                           = model.wMatrixSqrt(:) .* model.mBar;
mBarminusPsiMuH_times_Psi_times_invPsiPsi = mBarminusPsiMuH_times_Psi * model.invPsiPsi;
HessianSecondOrderTerm = [];
for s=1:(numSegments-1)
    %calculate the Hessian of the shape function in sparse mode
    Hs = spanPoly2FourierHessianSparse(model.segment{s}.v(:,1),model.segment{s}.v(:,2),...
        real(model.kMatrix),imag(model.kMatrix));
    %calculate second order terms of the Hessian of the mean in sparse mode
    mBarminusPsiMuH_times_wMatrixSqrtMu = (conj(mu(s)*model.wMatrixSqrt).* mBarminusPsimu)';
    HsTemp = zeros(size(Hs,1),size(Hs,2));
    for r=1:size(Hs,1)
        for c=1:size(Hs,2)
            %the contribution of Psi
            secondOrderTermOfPsi = mBarminusPsiMuH_times_wMatrixSqrtMu * squeeze(Hs(r,c,:));
            %the contribution of mu
            temp1_Term4_and_5      = zeros(numSegments);
            temp1_Term4_and_5(s,:) = squeeze(Hs(r,c,:))'*wMatrixSqrtPsi;
            temp1_Term4_and_5      = -(temp1_Term4_and_5 + temp1_Term4_and_5')*mu;
            temp1_Term4_and_5(s,1) = temp1_Term4_and_5(s,1) + squeeze(Hs(r,c,:))'*wMatrixSqrtmBar;
            secondOrderTermOfMu    = mBarminusPsiMuH_times_Psi_times_invPsiPsi * temp1_Term4_and_5;
            %sum the second order contributions
            HsTemp(r,c) = secondOrderTermOfPsi + secondOrderTermOfMu;
        end
    end
    Htemp= zeros(2*length(model.segment{s}.v));
    Htemp(1:2:end,1) = .5* HsTemp(:,1); %d2S_dx2_vector(:);
    Htemp(2:2:end,1) = .5* HsTemp(:,5); %d2S_dy2_vector(:);
    Htemp(1:2:end,2) =     HsTemp(:,2); %d2S_dxdy_vector(:);
    Htemp(2:2:end,2) =     HsTemp(:,6); %d2S_dydxplus1_vector(:);
    Htemp(1:2:end,3) =     HsTemp(:,3); %d2S_dxdxplus1_vector(:);
    Htemp(2:2:end,3) =     HsTemp(:,7); %d2S_dydyplus1_vector(:);
    Htemp(1:2:end,4) =     HsTemp(:,4); %d2S_dxdyplus1_vector(:);
    for r=1:2*length(model.segment{s}.v)
        Htemp(r,:) = circshift(Htemp(r,:),[0,r-1]);
    end
    Htemp = Htemp + Htemp.'; % this is the final result!
    %append to Hessian
    HessianSecondOrderTerm(...
        size(HessianSecondOrderTerm,1)+1:size(HessianSecondOrderTerm,1)+size(Htemp,1),...
        size(HessianSecondOrderTerm,2)+1:size(HessianSecondOrderTerm,2)+size(Htemp,2),:) = ...
        Htemp;
end
model.Hessian = -2*real(HessianFirstOrderTerm + HessianSecondOrderTerm);
return

function model = updateGradient(model);
mBarminusPsimuwMatrixSqrt = model.mBarminusPsimu .* model.wMatrixSqrt;
model.gradMu = [];
for s=1:(size(model.segment,2)-1)
    %gradient of objective function J
    [dfPoly_dxS dfPoly_dyS] = spanPoly2FourierDerivative(model.segment{s}.v(:,1),model.segment{s}.v(:,2),real(model.kMatrix(:)),imag(model.kMatrix(:)));
    dfPoly_dxS = squeeze(dfPoly_dxS);
    dfPoly_dyS = squeeze(dfPoly_dyS);
    temp1 = mBarminusPsimuwMatrixSqrt' * [dfPoly_dxS  dfPoly_dyS];

    %gradient of the shape function only
    deltavPartOneS = model.segment{s}.mu * temp1;
    
    %gradient of the mean vector only
    gradMuPartOneS =  model.invPsiPsi(:,s) * conj(temp1);
    gradMuPartTwoS = -(model.pinvPsi.* repmat(model.wMatrixSqrt(:),[1,size(model.segment,2)]).') * [dfPoly_dxS  dfPoly_dyS] * model.segment{s}.mu;
    gradMuS = gradMuPartOneS + gradMuPartTwoS;
    deltavPartTwoS = model.mBarminusPsimu'* model.Psi * gradMuS;
    
    %put together the gradient of J
    model.segment{s}.deltav = -2*real(reshape(deltavPartOneS + deltavPartTwoS,[length(deltavPartOneS)/2,2]));
    
    %reshuffle gradMuS because [p1 p2 p3 p4 p5 p6 ...] = [x1 y1 x2 y2 x3 y3
    %...], and store it separately
    gradMuS = reshape([gradMuS(:,1:end/2) ; gradMuS(:,end/2+1:end)],size(gradMuS));
    model.gradMu = [model.gradMu gradMuS];
end
return

function model = updateMeanScore(model);
mu = model.pinvPsi*(model.mBar);
for s=1:(size(model.segment,2))
    model.segment{s}.mu = mu(s);
end
model.mBarminusPsimu = model.mBar - model.Psi*mu;
model.score = model.mBarminusPsimu'*model.mBarminusPsimu;
return

function model = updatePsiScore(model);
model.Psi = [];
for s=1:(size(model.segment,2)-1)
%     [fPolyR, fPolyI] = spanPoly2Fourier_mex(... % TODO: wmatrixsqrt is basically divsion by N
%                 model.segment{s}.v(:,1),...
%                 model.segment{s}.v(:,2),...
%                 real(model.kMatrix(:,1)),imag(model.kMatrix(:,1)));
%            fPoly=model.wMatrixSqrt.*reshape((fPolyR+1i*fPolyI),size(real(model.kMatrix(:,1))));
    fPoly = model.wMatrixSqrt.*spanPoly2Fourier(...
                model.segment{s}.v(:,1),...
                model.segment{s}.v(:,2),...
                real(model.kMatrix),imag(model.kMatrix));    
model.Psi = [model.Psi fPoly];
end
fBackground     = (model.wMatrixSqrt).*(model.N^2).*sinc(real(model.kMatrix)*model.N).*sinc(imag(model.kMatrix)*model.N);
model.Psi       = [model.Psi fBackground];
model.PsiPsi    = model.Psi'*model.Psi;
model.invPsiPsi = inv(model.PsiPsi);
model.pinvPsi   = model.invPsiPsi*model.Psi';
for s=1:(size(model.segment,2))
    mu(s,1) = model.segment{s}.mu;
end
model.mBarminusPsimu = model.mBar - model.Psi * mu;
model.score = model.mBarminusPsimu'*model.mBarminusPsimu;
return

function model = gradientDescentStep(model,optLevel,epsilonMultiplier,plotOn)
%step 1: estimate the gradient of the score (objective function) with
%respect to all contours' vertex coordinates (parameters)
model = updateGradient(model);

%step 2: GVF w /boosting
segment = model.segment;
if optLevel==1
    %disp('L=1');
    vAllSegmentsCombined    = [];
    lAllSegmentsCombined    = [];
    gradAllSegmentsCombined = [];
    for s=1:(size(segment,2)-1)
        v = segment{s}.v; x = v(:,1); y = v(:,2);
        l = .5*(sqrt((circshift(x,-1) - x).^2 + (circshift(y,-1) - y).^2) + sqrt((circshift(x,1) - x).^2 + (circshift(y,1) - y).^2));
        gradSegment = segment{s}.deltav;
        
        vAllSegmentsCombined    = [vAllSegmentsCombined    ; v];
        lAllSegmentsCombined    = [lAllSegmentsCombined    ; l];
        gradAllSegmentsCombined = [gradAllSegmentsCombined ; gradSegment];
    end
    vBar = lAllSegmentsCombined'*vAllSegmentsCombined/sum(lAllSegmentsCombined);
    [gradT gradR gradS gradC Fx Fy Mz Fr] = ...
            projectGradient(gradAllSegmentsCombined,vAllSegmentsCombined,vBar,lAllSegmentsCombined);

    timeStep  = epsilonMultiplier * .0005;
    lambdaT   = 1;
    lambdaR   = 1;
    lambdaS   = 0;

    A = (1+timeStep*lambdaS*(-Fr))*[cos(timeStep*lambdaR*(-Mz)) -sin(timeStep*lambdaR*(-Mz));sin(timeStep*lambdaR*(-Mz)) cos(timeStep*lambdaR*(-Mz))];
    B = timeStep*lambdaT*[-Fx -Fy];

    for s=1:(size(segment,2)-1)
        for iii=1:size(segment{s}.v,1)
            segment{s}.v(iii,1:2) = (segment{s}.v(iii,1:2)-vBar)*A + B + vBar;
        end
    end
elseif ((optLevel==2) | (optLevel==3))
    for s=1:(size(segment,2)-1)
        sectionsId              = segment{s}.i;
        boostingParameterMatrix = segment{s}.p;
        v                       = segment{s}.v;
        x = v(:,1);
        y = v(:,2);
        l = .5*(sqrt((circshift(x,-1) - x).^2 + (circshift(y,-1) - y).^2) + sqrt((circshift(x,1) - x).^2 + (circshift(y,1) - y).^2));
        if optLevel==2
            %disp('L=2');
            gradSection = segment{s}.deltav;
            lSection    = l;
            vSection    = v;
            vSectionBar = lSection'*vSection/sum(lSection);
            [gradSectionT gradSectionR gradSectionS gradSectionC Fx Fy Mz Fr] = projectGradient(gradSection,vSection,vSectionBar,lSection);

            timeStep  = epsilonMultiplier * .001;
            lambdaT   = 1;
            lambdaR   = 1;
            lambdaS   = 0;

            A = (1+timeStep*lambdaS*(-Fr))*[cos(timeStep*lambdaR*(-Mz)) -sin(timeStep*lambdaR*(-Mz));sin(timeStep*lambdaR*(-Mz)) cos(timeStep*lambdaR*(-Mz))];
            B = timeStep*lambdaT*[-Fx -Fy];

            vSectionTemp = [];
            for iii=1:size(vSection,1)
                    vSectionTemp(iii,1:2) = (vSection(iii,1:2)-vSectionBar)*A + B + vSectionBar;
            end        
            segment{s}.v = vSectionTemp;
        elseif optLevel==3
            %disp('L=3');
            vTemp = zeros(size(segment{s}.v));
            for sId=1:max(sectionsId)
                gradSection   = segment{s}.deltav(sectionsId == sId,:);
                lSection      = l(sectionsId == sId,:);
                vSection      = segment{s}.v(sectionsId == sId,:);
                vSectionBar   = lSection'*vSection/sum(lSection);

                if ((optLevel==3) & (s==3) & (sId==2))
                    %disp('VELUM!!!');                
                    %use Palate and Nasal Cavity end points to find rotation center
                    vPalate      = segment{s}.v(sectionsId == 1,:);
                    vNasalcavity = segment{s}.v(sectionsId == 3,:);
                    vSectionBar  = 1/4*(vPalate(end,:) + vNasalcavity(1,:) + vSection(1,:) + vSection(end,:));
                    %figure(10);clf;
                    %figure(1); subplot(121); hold off;
                    %plot(vSection(:,1),vSection(:,2),'kx-');hold on;
                    %plot(vPalate(:,1),vPalate(:,2),'rx-');hold on;
                    %plot(vNasalcavity(:,1),vNasalcavity(:,2),'gx-');hold on;
                    %plot(vSectionBar(1),vSectionBar(2),'ko');
                    %pause;
                end

                if ((optLevel==3) & (s==1) & (sId==1))
                    %disp('EPIGLOTTIS!!!');
                    %use Neck and Tongue end points to find rotation center
                    vNeck   = segment{s}.v(sectionsId == 6,:);
                    vTongue = segment{s}.v(sectionsId == 2,:);
                    vSectionBar  = 1/4*(vNeck(end,:) + vTongue(1,:) + vSection(1,:) + vSection(end,:));
                    %figure(11);clf;
                    %figure(1); subplot(122); hold off;
                    %plot(vSection(:,1),vSection(:,2),'kx-');hold on;
                    %plot(vNeck(:,1),vNeck(:,2),'rx-');hold on;
                    %plot(vTongue(:,1),vTongue(:,2),'gx-');hold on;
                    %plot(vSectionBar(1),vSectionBar(2),'ko');
                    %pause;
                end

                [gradSectionT gradSectionR gradSectionS gradSectionC Fx Fy Mz Fr] = projectGradient(gradSection,vSection,vSectionBar,lSection);

                timeStep = epsilonMultiplier * 0.001;
                lambdaT  = boostingParameterMatrix(sId,1);
                lambdaR  = boostingParameterMatrix(sId,2);
                lambdaS  = boostingParameterMatrix(sId,3);

                A = (1+timeStep*lambdaS*(-Fr))*[cos(timeStep*lambdaR*(-Mz)) -sin(timeStep*lambdaR*(-Mz));sin(timeStep*lambdaR*(-Mz)) cos(timeStep*lambdaR*(-Mz))];
                B = timeStep*lambdaT*[-Fx -Fy];
                vSectionTemp = [];
                for iii=1:size(vSection,1)
                    vSectionTemp(iii,1:2) = (vSection(iii,1:2)-vSectionBar)*A + vSectionBar + B;
                end        
                vTemp(sectionsId == sId,:) = vSectionTemp;
            end
            % blending at the connections between the segments
            vDifference = vTemp - segment{s}.v;
            vD1 = vDifference .* repmat((sectionsId ~= circshift(sectionsId,[-1 0])),1,2);
            vD2 = vDifference .* repmat((sectionsId ~= circshift(sectionsId,[ 1 0])),1,2);
            vD3 = (vD1 + circshift(vD2,[-1 0]))/2;
            vD4 = vD3 + circshift(vD3,[ 1 0]);
            vDifference = vDifference.*(~repmat(((sectionsId ~= circshift(sectionsId,[-1 0])) ~= (sectionsId ~= circshift(sectionsId,[ 1 0]))),1,2));
            vDifference = vDifference + vD4;
            %contour data update
            segment{s}.v = segment{s}.v + vDifference;
        end
    end
elseif optLevel==4
    % ASTERIOS do not deform hard palate in step 4
    %disp('L=4');
      timeStep = epsilonMultiplier * .05;
%     for s=[1 2]
%         segment{s}.v = segment{s}.v - timeStep * segment{s}.deltav;
%     end;
%     s=3;
%     sectionID=segment{s}.i;
%     for secID=2:6;
%         segment{s}.v(sectionID==secID,:) = segment{s}.v(sectionID==secID,:) - timeStep * segment{s}.deltav(sectionID==secID,:);
%     end;
    for s=1:(size(segment,2)-1)
        segment{s}.v = segment{s}.v - timeStep * segment{s}.deltav;
    end
else
    disp('L=not valid');
end
model.segment  = segment;

%step 3: MMSE estimate the parameters for all inside regions and the outside region
model = updatePsiScore(model);
model = updateMeanScore(model);
%disp(sprintf('Score after = %4.5f',model.score));

%step 4: plot some stuff
if plotOn 
    plotAll(1,model.mBar ./ model.wMatrixSqrt,'gridded',model.N,model.kMatrix,model.wMatrixSqrt.^2,model.segment);
end
return

function [gradT gradR gradS gradC Fx Fy Mz Fr] = projectGradient(grad,v,vBar,massVector)
% gradT - gradient projected on translational movement 
% gradR - gradient projected on rotational movement around vBar
% gradS - gradient projected in radial stretch centered at vBar
% gradC - remaining gradient, i.e. complement to sum of all projections
% Fx Fy - mean forces in x and y direction
% Mz    - z component of the mean torque
% Fr    - radial mean stretch force 

vMinusvBar      = [v(:,1)-vBar(1) v(:,2)-vBar(2)];
%Translational component
% gradT           = repmat(massVector'*grad / sum(massVector),[size(grad,1),1]);
% gradT           = repmat(ones(size(massVector))'*grad / sum(massVector),[size(grad,1),1]);
Fx = sum(grad(:,1),1);
Fy = sum(grad(:,2),1);
gradT = [Fx*ones(size(grad,1),1) Fy*ones(size(grad,1),1)];

%Rotational component
% integralNum     = ones(size(massVector))'*cross([grad zeros(size(grad,1),1)], [vMinusvBar zeros(size(vMinusvBar,1),1)]);
% integralNum     = massVector'*cross([grad zeros(size(grad,1),1)], [vMinusvBar zeros(size(vMinusvBar,1),1)]);
% integralDenom   = sum(massVector'* sum(vMinusvBar.^2,2),1);
% angularMomentum = integralNum / integralDenom;
% gradR           = cross([vMinusvBar zeros(size(vMinusvBar,1),1)],repmat(angularMomentum,size(vMinusvBar,1),1));
M  = cross([grad zeros(size(grad,1),1)], [vMinusvBar zeros(size(vMinusvBar,1),1)]);
M  = sum(M,1);
Mz = M(3);
gradR = cross([vMinusvBar zeros(size(vMinusvBar,1),1)],repmat(M,size(vMinusvBar,1),1));
gradR = gradR(:,1:2);

%Scaling component
% integralNum     = massVector'*sum(grad.*vMinusvBar,2);
% integralNum     = ones(size(massVector))'*sum(grad.*vMinusvBar,2);
% integralDenom   = sum(massVector'* sum(vMinusvBar.^2,2),1);
% overallStretch  = integralNum / integralDenom;
% gradS           = overallStretch * vMinusvBar;
Fr    = sum(sum(grad.*vMinusvBar,2),1);
gradS = Fr * vMinusvBar;

gradC = grad - gradT - gradR - gradS;
return

function im = reconImage(f,fDataType,N,kMatrix,wMatrix)
if strcmp(fDataType,'gridded')
    fgridded = gridkb(kMatrix, f, wMatrix, 2*N, 3, 2);
    im = ift(fgridded);
    im = im( N/2+1 : 3*N/2,N/2+1 : 3*N/2 );
    [yDeapod,xDeapod] = meshgrid(-N/2 : N/2-1, -N/2 : N/2-1);
    im = im ./ (sinc(.5*xDeapod/N).^2 .* sinc(.5*yDeapod/N).^2 );
    im = 2 * im * (N^2);
elseif strcmp(fDataType,'full2DFT')
    im = fftshift(ifft2(fftshift(f)));
elseif strcmp(fDataType,'half2DFT')
    disp('half2DFT recon has not been implemented yet.');
    pause;
    return;
end
return

function plotAll(figureNumber,fObserved,fDataType,N,kMatrix,wMatrix,segment)
figure(1); subplot(1,2,1)
imObserved = reconImage(fObserved,fDataType,N,kMatrix,wMatrix);
[lx,ly]=size(imObserved);
%imObservedUpsampled = interpft2(imObserved,5*lx,5*ly);
imObservedUpsampled = imObserved;
%figure(1); subplot(241);
colormap(gray);
iid=image([(-lx/2):(lx/2-1)], [(-ly/2):(ly/2-1)],abs(real(flipud(imObservedUpsampled))*50)); hold on;
for s=1:(size(segment,2)-1)
    sectionsId = segment{s}.i;
    v          = segment{s}.v;
    colors = ['r' 'g' 'b' 'y' 'c' 'm' 'k'];
    for sId=1:max(sectionsId)
        %plot( v(sectionsId==sId,1),-v(sectionsId==sId,2),[colors(sId) 'x-']); hold on;
        plot( v(sectionsId==sId,1),-v(sectionsId==sId,2),[colors(sId)],  'LineWidth', 2); hold on;
    end
end
axis([-.5 .5 -.5 .5]*N);
axis equal; axis off;
drawnow;
return

function model = cleanUpIntersections(model)
%overlap removal
for s=1:(size(model.segment,2)-1)
    s1=s;                                 %this is the first boundary to be checked
    s2=mod((s1+1),(size(model.segment,2)-1))+1; %this is the second boundary
    v1 = model.segment{s1}.v;
    v2 = model.segment{s2}.v;

    pointOfV2IsInsideV1 = inpolygon(v2(:,1),v2(:,2),v1(:,1),v1(:,2));
    transitions = filter([1 -1],1, [pointOfV2IsInsideV1 ; pointOfV2IsInsideV1(1)]);
    id1 = find(transitions==1);
    for ss=1:length(id1)
        shiftAmount = id1(ss)-2;
        v2Shifted          = circshift(v2         ,[-shiftAmount 0]);
        transitionsShifted = circshift(transitions,[-shiftAmount 0]);
        id2 = find(transitionsShifted==-1);
        vOutside1 = v2Shifted(1  ,:);
        vOutside2 = v2Shifted(id2(1),:);
        vInside   = v2Shifted(2:id2(1)-1,:);
        newX = linspace(vOutside1(1),vOutside2(1),size(vInside,1)+2);
        newY = linspace(vOutside1(2),vOutside2(2),size(vInside,1)+2);
        v2Shifted(2:id2(1)-1,1) = newX(2:end-1);
        v2Shifted(2:id2(1)-1,2) = newY(2:end-1);
        v2 = circshift(v2Shifted,[shiftAmount 0]);
    end    
    
    pointOfV1IsInsideV2 = inpolygon(v1(:,1),v1(:,2),v2(:,1),v2(:,2));
    transitions = filter([1 -1],1, [pointOfV1IsInsideV2 ; pointOfV1IsInsideV2(1)]);
    id1 = find(transitions==1);
    for ss=1:length(id1)
        shiftAmount = id1(ss)-2;
        v1Shifted          = circshift(v1         ,[-shiftAmount 0]);
        transitionsShifted = circshift(transitions,[-shiftAmount 0]);
        id2 = find(transitionsShifted==-1);
        vOutside1 = v1Shifted(1  ,:);
        vOutside2 = v1Shifted(id2(1),:);
        vInside   = v1Shifted(2:id2(1)-1,:);
        newX = linspace(vOutside1(1),vOutside2(1),size(vInside,1)+2);
        newY = linspace(vOutside1(2),vOutside2(2),size(vInside,1)+2);
        v1Shifted(2:id2(1)-1,1) = newX(2:end-1);
        v1Shifted(2:id2(1)-1,2) = newY(2:end-1);
        v1 = circshift(v1Shifted,[shiftAmount 0]);
    end    
    
    model.segment{s1}.v =  v1;
    model.segment{s2}.v =  v2;
end

%special case: lower lip vs. tongue (ASTERIOS)
    v1 = model.segment{1}.v(model.segment{1}.i == 2,:);
    v2 = model.segment{1}.v(model.segment{1}.i == 4,:);

    pointOfV2IsInsideV1 = inpolygon(v2(:,1),v2(:,2),v1(:,1),v1(:,2));
    transitions = filter([1 -1],1, [pointOfV2IsInsideV1 ; pointOfV2IsInsideV1(1)]);
    id1 = find(transitions==1);
    for ss=1:length(id1)
        shiftAmount = id1(ss)-2;
        v2Shifted          = circshift(v2         ,[-shiftAmount 0]);
        transitionsShifted = circshift(transitions,[-shiftAmount 0]);
        id2 = find(transitionsShifted==-1);
        if numel(id2)>0
            vOutside1 = v2Shifted(1  ,:);
            vOutside2 = v2Shifted(id2(1),:);
            vInside   = v2Shifted(2:id2(1)-1,:);
            newX = linspace(vOutside1(1),vOutside2(1),size(vInside,1)+2);
            newY = linspace(vOutside1(2),vOutside2(2),size(vInside,1)+2);
            v2Shifted(2:id2(1)-1,1) = newX(2:end-1);
            v2Shifted(2:id2(1)-1,2) = newY(2:end-1);
            v2 = circshift(v2Shifted,[shiftAmount 0]);
        end;
    end    
    
    pointOfV1IsInsideV2 = inpolygon(v1(:,1),v1(:,2),v2(:,1),v2(:,2));
    transitions = filter([1 -1],1, [pointOfV1IsInsideV2 ; pointOfV1IsInsideV2(1)]);
    id1 = find(transitions==1);
    for ss=1:length(id1)
        shiftAmount = id1(ss)-2;
        v1Shifted          = circshift(v1         ,[-shiftAmount 0]);
        transitionsShifted = circshift(transitions,[-shiftAmount 0]);
        id2 = find(transitionsShifted==-1);
        vOutside1 = v1Shifted(1  ,:);
        vOutside2 = v1Shifted(id2(1),:);
        vInside   = v1Shifted(2:id2(1)-1,:);
        newX = linspace(vOutside1(1),vOutside2(1),size(vInside,1)+2);
        newY = linspace(vOutside1(2),vOutside2(2),size(vInside,1)+2);
        v1Shifted(2:id2(1)-1,1) = newX(2:end-1);
        v1Shifted(2:id2(1)-1,2) = newY(2:end-1);
        v1 = circshift(v1Shifted,[shiftAmount 0]);
    end    
    
    model.segment{1}.v(model.segment{1}.i == 2,:) =  v1;
    model.segment{1}.v(model.segment{1}.i == 4,:) =  v2;


%contour spike removal
angleThreshold = 22;
for s=1:(size(model.segment,2)-1)
    v = model.segment{s}.v;
    for vId = 1:size(v,1)
        v = circshift(v,[ vId 0]);
        vTriang = [v(end,:) ; v(1,:) ; v(2,:)];
%         triangArea = abs(polyarea(vTriang(:,1),vTriang(:,2))); %%%
%         triangPeri = sum(sqrt(sum((diff([vTriang ; vTriang(1,:)],1,1).^2),2)),1);
%         triangCirc = 4*pi*triangArea/(triangPeri^2);
%         circularityThreshold = 0.0001;
%         if triangCirc < circularityThreshold
%             v(1,:) = (v(end,:) + v(2,:))/2;
%             v(end,:) = (v(end-1,:) + v(1,:))/2;
%             v( 2 ,:) = (v(  3  ,:) + v(1,:))/2;
%         end %%%
        v1 = v(end,:) - v(1,:);
        l1 = sqrt(v1*v1');
        v2 =  v(2,:)  - v(1,:);
        l2 = sqrt(v2*v2');
        enclosedAngle = acosd((v1*v2')/(l1*l2));
        if (enclosedAngle < angleThreshold) | ((360 - enclosedAngle) < angleThreshold)
            v(1,:) = (v(end,:) + v(2,:))/2;
            v(end,:) = (v(end-1,:) + v(1,:))/2;
            v( 2 ,:) = (v(  3  ,:) + v(1,:))/2;
        end
        v = circshift(v,[-vId 0]);
    end
    model.segment{s}.v = v;
end

% for s=1:(size(model.segment,2)-1)
%     
%     for ii=1:max(model.segment{s}.i)
%         
%         v = model.segment{s}.v(model.segment{s}.i == ii,:);
%         nPoints=size(v,1);
%         v1=InterpolateContourPoints2D(v,nPoints);
%         %[xs,ys]=smooth_contours(v1(:,1),v1(:,2),ceil(size(v1,1)/2));
%         %model.segment{s}.v(model.segment{s}.i == ii,:)= [xs,ys];
%         model.segment{s}.v(model.segment{s}.i == ii,:) = v1;
%         
%     end;
%     
% end;

return

function S = spanPoly2Fourier(x,y,kx,ky)
%function S = spanPoly2Fourier(x,y,kx,ky)
% Computes the 2D-DFT of the polygonal shape defined by the corner points 
% (x,y) at the k-space location defined by the kx and ky matrices.
%
% kx,ky : [ -.5             , .5             )
% x     : [ -1/(2*delta_kx) , 1/(2*delta_kx) ) 
% y     : [ -1/(2*delta_ky) , 1/(2*delta_ky) ) 
%
% The algorithm is from:
% K.McInturff and P.S.Simon, "The Fourier Transform of Linearly Varying
% Functions with Polygonal Support", IEEE Transactions on Antennas and
% Propagation, 39(9), 1991.
%
% Erik Bresch, USC SPAN Group, 2007

gammaN        = [real(x(:)') ; real(y(:)')];
gammaNplusOne = circshift(gammaN,[0,-1]);
alphaN        = gammaNplusOne - gammaN;
betaN         = gammaNplusOne + gammaN;
K             = [real(kx(:)) real(ky(:))];

KxAndKyAreZero = (K(:,1)==0) & (K(:,2)==0);

Kbar = (sum(K.^2,2) + KxAndKyAreZero);
Kbar = [real(ky(:))./Kbar -real(kx(:))./Kbar];

nu=size(gammaN,1);


S = sum(((Kbar * alphaN) .* exp(-1i*pi*K*betaN) .* sinc(K*alphaN)),2)/(2*1i*pi);  % TODO Could I remove sinc ??? (besselj)
S(KxAndKyAreZero==1) = polyarea(x,y);
S = reshape(S,size(kx));
return

function [dS_dx, dS_dy] = spanPoly2FourierDerivative(x,y,kx,ky)
%function [dS_dx dS_dy] = spanPoly2FourierDerivative(x,y,kx,ky)
% Computes the derivative of the 2D-DFT of the polygonal shape defined 
% by the corner points (x,y) at the k-space location defined by the kx and
% ky matrices for dx_n and dy_n.
%
% dS_dx(:,:,i) = dS/dxi for the (kx,ky) locations specified in "kx","ky".
% dS_dy(:,:,i) = dS/dyi for the (kx,ky) locations specified in "kx","ky".
%
% kx,ky : [ -.5             , .5             )
% x     : [ -1/(2*delta_kx) , 1/(2*delta_kx) ) 
% y     : [ -1/(2*delta_ky) , 1/(2*delta_ky) ) 
%
% The algorithm uses the equations for the Fourier transform of polygons
% from:
% K.McInturff and P.S.Simon, "The Fourier Transform of Linearly Varying
% Functions with Polygonal Support", IEEE Transactions on Antennas and
% Propagation, 39(9), 1991.
%
% Erik Bresch, USC SPAN Group, 2007

gammai        = [real(x(:)') ; real(y(:)')];
gammaiplusOne = circshift(gammai,[0,-1]);
alphai        = gammaiplusOne - gammai;
betai         = gammaiplusOne + gammai;
K             = [real(kx(:)) real(ky(:))];

numberOfVertices     = size(gammai,2);
numberOfKspacepoints = size(K,1);

KxAndKyAreZero  = (K(:,1)==0)&(K(:,2)==0);
KxAndKyZeroRows = repmat(KxAndKyAreZero,1,numberOfVertices);

Kbar = j*2*pi*(sum(K.^2,2) + KxAndKyAreZero);
Kbar = [real(ky(:))./Kbar -real(kx(:))./Kbar];

exp_minusjpi_K_betai             = exp(-j*pi*K*betai);
Kbar_alphai_exp_minusjpi_K_betai = Kbar*alphai .* exp_minusjpi_K_betai;
sinc_K_alphai                    = sinc(K*alphai);
sincDerivative_K_alphai          = sincDerivative(K*alphai);

exp_minusjpi_K_betai_sinc_K_alphai                       = exp_minusjpi_K_betai .* sinc_K_alphai;
Kbar_alphai_exp_minusjpi_K_betai_sinc_K_alphai           = Kbar_alphai_exp_minusjpi_K_betai .* sinc_K_alphai;
Kbar_alphai_exp_minusjpi_K_betai_sincDerivative_K_alphai = Kbar_alphai_exp_minusjpi_K_betai .* sincDerivative_K_alphai;

c_i = repmat(-Kbar(:,1),1,numberOfVertices)     .* exp_minusjpi_K_betai_sinc_K_alphai;
d_i = repmat((-j*pi*K(:,1)),1,numberOfVertices) .* Kbar_alphai_exp_minusjpi_K_betai_sinc_K_alphai;
e_i = repmat(-K(:,1),1,numberOfVertices)        .* Kbar_alphai_exp_minusjpi_K_betai_sincDerivative_K_alphai;

f_i = repmat(-Kbar(:,2),1,numberOfVertices)     .* exp_minusjpi_K_betai_sinc_K_alphai;
g_i = repmat((-j*pi*K(:,2)),1,numberOfVertices) .* Kbar_alphai_exp_minusjpi_K_betai_sinc_K_alphai;
h_i = repmat(-K(:,2),1,numberOfVertices)        .* Kbar_alphai_exp_minusjpi_K_betai_sincDerivative_K_alphai;

dS_dx = (c_i - circshift(c_i,[0,1]) + d_i + circshift(d_i,[0,1]) + e_i - circshift(e_i,[0,1]));
dS_dy = (f_i - circshift(f_i,[0,1]) + g_i + circshift(g_i,[0,1]) + h_i - circshift(h_i,[0,1]));

dS_dx_AtOrigin = .5*(circshift( real(y(:)'),[0 -1]) - circshift( real(y(:)'),[0 +1] ) );
dS_dy_AtOrigin = .5*(circshift( real(x(:)'),[0 +1]) - circshift( real(x(:)'),[0 -1] ) );

dS_dx = (~KxAndKyZeroRows).*dS_dx  + KxAndKyZeroRows .* repmat(dS_dx_AtOrigin,numberOfKspacepoints,1);
dS_dy = (~KxAndKyZeroRows).*dS_dy  + KxAndKyZeroRows .* repmat(dS_dy_AtOrigin,numberOfKspacepoints,1);

dS_dx = reshape(dS_dx,[size(kx) numberOfVertices]);
dS_dy = reshape(dS_dy,[size(ky) numberOfVertices]);


return

function dsinc_dx = sincDerivative(x)
xEqZero = (x==0);
x(xEqZero) = 1;
dsinc_dx = (((pi*x).*cos(pi*x) - sin(pi*x))./(pi* (x.^2))) .* (1-xEqZero);
return

function K=InterpolateContourPoints2D(P,nPoints)
% This function resamples a few points describing a countour , to a smooth
% contour of uniform sampled points.
%
% K=InterpolateContourPoints(P,nPoints)
%
% input,
%  P : Inpute Contour, size N x 2  (with N>=4)
%  nPoints : Number of Contour points as output
% 
% output,
%  K : Uniform sampled Contour points, size nPoints x 2
%
% Function is written by D.Kroon University of Twente (July 2010)
%
% example, 
%  % Show an image
%   figure, imshow(imread('moon.tif'));
%  % Select some points with the mouse
%   [y,x] = getpts;
%  % Make an array with the clicked coordinates
%   P=[x(:) y(:)];
%  % Interpolate inbetween the points
%   Pnew=InterpolateContourPoints2D(P,100)
%  % Show the result
%   hold on; plot(P(:,2),P(:,1),'b*');
%   plot(Pnew(:,2),Pnew(:,1),'r.');
%
% Function is written by D.Kroon University of Twente (July 2010)

% Interpolate points inbetween
%O(:,1)=interp([P(end-3:end,1);P(:,1);P(:,1);P(1:4,1)],10);
%O(:,2)=interp([P(end-3:end,2);P(:,2);P(:,2);P(1:4,2)],10);
% O(:,1)=interp(P(:,1),2);
% O(:,2)=interp(P(:,2),2);
%O=O(41:end-39,:); 

O=P;

% Calculate distance between points
dis=[0;cumsum(sqrt(sum((O(2:end,:)-O(1:end-1,:)).^2,2)))];

% Resample to make uniform points
K(:,1) = interp1(dis,O(:,1),linspace(0,dis(end),nPoints));
K(:,2) = interp1(dis,O(:,2),linspace(0,dis(end),nPoints));
%K=K(round(end/4):round(end/4)+nPoints-1,:);

 
 return