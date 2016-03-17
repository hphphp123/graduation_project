function [parameters_noise,parameters_enhance] = initialise_parameters(ns_ps)
len_val = length(ns_ps);%���ò���
alpha_d_val=0.85;
alpha_s_val=0.8;
U_val=2;
V_val=8;
Bmin_val=1.66;
gamma0_val=4.6;
gamma1_val=3;
psi0_val=1.67;
alpha_val=0.92;
beta_val=3;
b_val=hanning(3);
B_val=sum(b_val);
b_val=b_val/B_val;
eps_min_val=10^(-15/10);
Sf_val=zeros(len_val,1);
Sf_tild_val=zeros(len_val,1);
Sf_val(1) = ns_ps(1);
for f=2:len_val-1
    Sf_val(f)=sum(b_val.*[ns_ps(f-1);ns_ps(f);ns_ps(f+1)]); 
end
Sf_val(len_val)=ns_ps(len_val);
Sf_tild_val = zeros(len_val,1);

len2a=len_val;
psi_val=ones(len2a,1);

parameters_noise =  struct('n',2,'len',len_val,'noise_ps',Sf_val,'noise_tild',Sf_val,'gamma',ones(len_val,1),'Sf',Sf_val,...
            'Smin',Sf_val,'S',Sf_val,'S_tild',Sf_val,'GH1',ones(len_val,1),'Smin_tild',Sf_val,'Smin_sw',Sf_val,'Smin_sw_tild',Sf_val,...
            'stored_min',max(Sf_val)*ones(len_val,U_val),'stored_min_tild',max(Sf_val)*ones(len_val,U_val),'u1',1,'u2',1,'j',1,...
            'alpha_d',alpha_d_val,'alpha_s',alpha_s_val,'U',U_val,'V',V_val,'Bmin',Bmin_val,'gamma0',gamma0_val,'gamma1',gamma1_val,'psi0',psi0_val,'alpha',alpha_val,'beta',beta_val,...
            'b',b_val,'Sf_tild',Sf_tild_val,'eps_min',eps_min_val);
parameters_enhance = struct('len',len_val,'Gmin',10^(-10/10),'noise_ps',Sf_val,'eps_min',eps_min_val,'beta',0.7,'psi',psi_val,...
                    'psi_min',0.1,'psi_max',0.316,'psi_pmin',1,'psi_pmax',10,'psi_fr',1000,'psi_peak',0,'qmax',0.95,'sigout',ones(len_val,1),...
                   'alpha',alpha_val,'GH1',ones(len_val,1),'gamma',ones(len_val,1),'eps',ones(len_val,1));

  