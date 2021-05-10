clear

do_save = 0;
do_plot = 0;
do_table = 1;

figure_size_alpha = [100 100 250 220];
% figure_size_f = [100 100 320 280];
figure_size_f = [100 100 1100 280];

i_mod = 1;
mods{i_mod} = struct('dir','test_int_N1_hlayF3_dof13_2019-01-28_10-29-25','name','free1'); i_mod = i_mod+1;
mods{i_mod} = struct('dir','test_int_N1_hlayF3_dof13_2019-02-19_18-38-29','name','free2'); i_mod = i_mod+1;
mods{i_mod} = struct('dir','test_int_N1_hlayF3_dof13_2019-01-27_11-21-37','name','free3'); i_mod = i_mod+1;
mods{i_mod} = struct('dir','test_int_N1_hlayF3_dof12_2019-02-01_11-43-58','name','strong'); i_mod = i_mod+1;
mods{i_mod} = struct('dir','test_int_N1_hlayF3_dof13_2019-02-20_18-54-43','name','weak'); i_mod = i_mod+1;

%%

problem = problem_get('testcase','testcase_1var_exp.ini');
mod_HF = problem.get_model(problem);

for i_mod = 1:length(mods)
    mod_learned = read_model_fromfile(problem,mods{i_mod}.dir);
    
    if do_plot || do_save
        alphafig = figure('units','pixels','position',figure_size_alpha);
        fig_alpha_plot(mod_learned)
        set(alphafig,'PaperPositionMode','auto')
        if do_save
            pause()
            print(alphafig, sprintf('fig/testcase_1var_exp_alphaplot_%s', mods{i_mod}.name),'-depsc','-painters');
    %         aa = get(gcf,'position')
        end
        pause(1e-16)
    end
    if do_table
        table = array2table([mod_learned.alpha_learned', mod_learned.alpha_original'],'VariableNames',{'alpha_learned','alpha_original'});
        writetable(table,sprintf('fig/testcase_1var_exp_alphaplot_%s.dat', mods{i_mod}.name),'Delimiter','tab')
    end
    
    if isfield(mod_learned,'alpha_to_alpha')
            nptx = 10;
            npta = 10;
            x_min = mod_HF.problem.y_min;
            x_max = mod_HF.problem.y_max;
            a_min = mod_HF.alpha_min;
            a_max = mod_HF.alpha_max;
            xx = linspace(x_min,x_max,nptx);
            aa = linspace(a_min,a_max,npta);
            [XX,AA] = meshgrid(xx,aa);
            FF_HF = zeros(npta,nptx);
            FF_NN = zeros(npta,nptx);
            for ix = 1:nptx
                for ia = 1:npta
                    FF_HF(ia,ix) = mod_HF.f_alpha(xx(ix),[],aa(ia));
                    FF_NN(ia,ix) = mod_learned.f_alpha(xx(ix),[],mod_learned.alpha_to_alpha(aa(ia)));
                end
            end

            Fmin = min(min(FF_HF(:)),min(FF_NN(:)));
            Fmax = max(max(FF_HF(:)),max(FF_NN(:)));

        if do_table
            pgfplots_write3Ddatafile(AA,XX,FF_HF,sprintf('fig/testcase_1var_exp_modelplot_HF_%s.dat', mods{i_mod}.name));
            pgfplots_write3Ddatafile(AA,XX,FF_NN,sprintf('fig/testcase_1var_exp_modelplot_NN_%s.dat', mods{i_mod}.name));
            pgfplots_write3Ddatafile(AA,XX,abs(FF_HF-FF_NN),sprintf('fig/testcase_1var_exp_modelplot_EE_%s.dat', mods{i_mod}.name));
%             array = [reshape(AA,[],1) reshape(XX,[],1) reshape(FF_HF,[],1)];
%             writematrix(array,sprintf('fig/testcase_1var_exp_modelplot_HF_%s.dat', mods{i_mod}.name),'Delimiter','tab')
%             array = [reshape(AA,[],1) reshape(XX,[],1) reshape(FF_NN,[],1)];
%             writematrix(array,sprintf('fig/testcase_1var_exp_modelplot_NN_%s.dat', mods{i_mod}.name),'Delimiter','tab')
%             array = [reshape(AA,[],1) reshape(XX,[],1) reshape(abs(FF_HF-FF_NN),[],1)];
%             writematrix(array,sprintf('fig/testcase_1var_exp_modelplot_EE_%s.dat', mods{i_mod}.name),'Delimiter','tab')
%             table = array2table(table,'VariableNames',{'state','alpha','rhs_HF', 'rhs_NN'});
%             writetable(table,sprintf('fig/testcase_1var_exp_modelplot_%s.dat', mods{i_mod}.name),'Delimiter','tab')
        end
        if do_plot || do_save
            modelfig = figure('units','pixel','position',figure_size_f);

    %         figure('units','pixel','position',figure_size_f);
            subplot(1,3,1)
            surf(XX,AA,FF_HF)
            title('HF model')
            axis([x_min x_max a_min a_max Fmin Fmax])
            xlabel('\theta'); ylabel('\alpha'); zlabel('$g(\theta,\alpha)$', 'interpreter','latex')

            subplot(1,3,2)
    %         figure('units','pixel','position',figure_size_f);
            surf(XX,AA,FF_NN)
            title('Learned model')
            axis([x_min x_max a_min a_max Fmin Fmax])
            xlabel('\theta'); ylabel('\alpha'); zlabel('$\hat{g}(\theta,\Gamma\{g\}(\alpha))$', 'interpreter','latex')

            subplot(1,3,3)
    %         figure('units','pixel','position',figure_size_f);
            surf(XX,AA,abs(FF_HF-FF_NN))
            title('Error')
            xlabel('\theta'); ylabel('\alpha'); zlabel('$|g(\theta,\alpha) - \hat{g}(\theta,\Gamma\{g\}(\alpha)) |$', 'interpreter','latex')

            if do_save
                pause()
                print(modelfig, sprintf('fig/testcase_1var_exp_modelplot_%s', mods{i_mod}.name),'-depsc','-painters');
            end
            pause(1e-16)
        end
    end
end
