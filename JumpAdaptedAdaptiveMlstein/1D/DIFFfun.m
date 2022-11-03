function W_terms =DIFFfun(sigma, input, dt, W)
% if additive == 1
%     W_terms = sigma*W;
% else  % additive == 0
%     diffusion = sigma*input;
%     D_diffusion = sigma;
%     W_terms = diffusion*W+0.5*D_diffusion*diffusion*(W^2-dt);
%     %     diffusion= sigma*(1-input^2);
%     %     D_diffusion =  -2*sigma*input;
% end

% W_terms=sigma*W;
% 
% diffusion = sigma;
% W_terms=diffusion*W;

% diffusion = sigma*input;
% W_terms=diffusion*W+0.5*sigma*diffusion*(W^2-dt);
% 
diffusion = sigma*(1-input^2);
D_diffusion =  -2*sigma*input;
W_terms=diffusion*W+0.5*D_diffusion*diffusion*(W^2-dt);
return


