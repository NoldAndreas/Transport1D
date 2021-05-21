import matplotlib.pyplot as plt
import numpy as np
import json
from scipy.integrate import solve_bvp
import pandas as pd

class Transport:
    def __init__(self, paramsFilename,L,supply_vs_demand):
        #***************************
        #Parameters from input file
        #***************************

        with open(paramsFilename) as f:
          self.params_input = json.load(f)

        self.params_input['L']                = L;
        self.params_input['supply_vs_demand'] = supply_vs_demand;

    def __initialization(self):
        params   = self.params_input;

        if(params['name'] == "somatic, passive protein transport"):
            params['translation'] = 'somatic';
            params['transport_p'] = 'passive';
        elif(params['name'] == "somatic, active protein transport"):
            params['translation'] = 'somatic';
            params['transport_p'] = 'active';
        elif(params['name'] == "dendritic, passive mRNA transport"):
            params['translation'] = 'dendritic';
            params['transport_p'] = 'passive';
            params['transport_m'] = 'passive';
        elif(params['name'] == "dendritic, active mRNA transport"):
            params['translation'] = 'dendritic';
            params['transport_p'] = 'passive';
            params['transport_m'] = 'active';


        L                = params['L'];
        supply_vs_demand = params['supply_vs_demand'];

        gamma_tl = params['gamma_tl'];#Translation rate (proteins per second per mRNA)

        # Active transport parameters of mRNA
        beta_plus_m  = params['beta_plus_m'];  #Switching rate from anterograde to resting (per second)
        beta_minus_m = params['beta_minus_m']; #Switching rate from retrograde to resting (per second)
        alpha_m      = params['alpha_m'];      #Switching rate from resting to anterograde (equal to to retrograde) (per second)

        # Active transport parameters of protein
        beta_plus_p  = params['beta_plus_p'];  #Switching rate from anterograde to resting (per second)
        beta_minus_p = params['beta_minus_p']; #Switching rate from retrograde to resting (per second)
        alpha_p      = params['alpha_p'];      #Switching rate from resting to anterograde (equal to to retrograde) (per second)

        #Translation
        translation_type  = params['translation']; #{mRNA or protein}

        #************************
        # Initialization
        #************************
        lambda_m  = np.log(2)/params['halflife_m']; #Degradation ratio for mRNA
        lambda_p  = np.log(2)/params['halflife_p']; #Degradation ratio for Protein

        eta_max_p  = params['eta_max'];
        eta_max_m  = eta_max_p/gamma_tl*lambda_m;

        #Compute effective velocities of mRNA transport
        #Set passive/active transport:
        if(params['transport_m'] == 'active'):
            gamma_m = 1/beta_plus_m + 1/beta_minus_m + 1/alpha_m;
            v_m     = params['v_m']/gamma_m*(1/beta_plus_m - 1/beta_minus_m);
            D_m     = (params['v_m'] - v_m)**2/(gamma_m*beta_plus_m**2) + (params['v_m'] + v_m)**2/(gamma_m*beta_minus_m**2);
        elif(params['transport_m'] == 'passive'):
            v_m     = 0;
            D_m     = params['D_m'];
        else:
            raise ValueError("transport_m must be either active or passive")

        #Compute effective velocities of protein transport
        if(params['transport_p'] == 'active'):
            gamma_p = 1/beta_plus_p + 1/beta_minus_p + 1/alpha_p;
            v_p     = params['v_p']/gamma_p*(1/beta_plus_p - 1/beta_minus_p);
            D_p     = (params['v_p'] - v_p)**2/(gamma_p*beta_plus_p**2) + (params['v_p'] + v_p)**2/(gamma_p*beta_minus_p**2);
        elif(params['transport_p'] == 'passive'):
            v_p     = 0;
            D_p     = params['D_p'];
        else:
            raise ValueError("transport_p must be either active or passive")

        #Set influxes
        if(translation_type=='dendritic'):
            J_p = 0;
            J_m = supply_vs_demand*L*eta_max_m;
        elif(translation_type=='somatic'):
            J_m      = 0;
            J_p      = supply_vs_demand*L*eta_max_p;
        else:
            print('No supply source defined');


        #Save
        self.params = {'D_m':D_m,'D_p':D_p,'v_m':v_m,'v_p':v_p,\
                       'lambda_m':lambda_m,'lambda_p':lambda_p,\
                       'J_p':J_p,'J_m':J_m,\
                       'eta_max_m':eta_max_m,'eta_max_p':eta_max_p,\
                       'gamma_tl':gamma_tl,\
                       'eta_0':params['eta_0'],\
                       'L':L};

    def GetSolution(self):
        if(hasattr(self,'res_p')):
            x     = np.linspace(0,self.params['L'],1000);
            p     = self.res_p.sol(x)[0];
            m     = self.res_m.sol(x)[0];
        else:
            x = [0,0];
            p = [0,0];
            m = [0,0];
        return pd.DataFrame({'x':x,'protein':p,'mRNA':m});

    def GetParametersAndResults(self):
        p_all = self.params_input;
        p_all['particles_in_active_transport_perLocalizedProtein']           = self.particles_in_active_transport_perLocalizedProtein;
        p_all['ratio_synapses_supplied']                 = self.ratio_synapses_supplied;
        p_all['ExcessmRNA_transcribed_perSynapticProtein']     = self.ExcessmRNA_transcribed_perSynapticProtein;
        p_all['ExcessProteins_translated_perSynapticProtein']  = self.ExcessProteins_translated_perSynapticProtein;
        return p_all;

    def __computeCostItems(self):

        x     = np.linspace(0,self.params['L'],1000);
        p     = self.res_p.sol(x)[0];
        m     = self.res_m.sol(x)[0];

        eta_max_p = self.params['eta_max_p'];
        eta_max_m = self.params['eta_max_m'];
        eta_0     = self.params['eta_0'];
        gamma_tl  = self.params['gamma_tl'];
        lambda_m  = self.params['lambda_m'];
        lambda_p  = self.params['lambda_p'];

        def u_m(m):
            return np.tanh(m/eta_max_m*eta_0)*eta_max_m;

        def u_p(p):
            return np.tanh(p/eta_max_p*eta_0)*eta_max_p;

        proteinsToSynapses            = np.trapz(u_p(p),x=x);
        ratio_synapses_supplied       = proteinsToSynapses/eta_max_p/self.params['L'];
        particles_in_active_transport = 0;
        if(self.params_input['transport_p'] == 'active'):
            proteins_all = np.trapz(p)
            particles_in_active_transport += proteins_all;
        if(self.params_input['transport_m'] == 'active'):
            mRNA_all = np.trapz(m);
            particles_in_active_transport += mRNA_all;

        m0                                          = np.trapz(u_m(m)/lambda_m,x=x);
        self.particles_in_active_transport_perLocalizedProtein  = particles_in_active_transport/(proteinsToSynapses/lambda_p);
        self.ratio_synapses_supplied                = ratio_synapses_supplied;
        self.ExcessProteins_translated_perSynapticProtein = np.max([0,(self.params['J_p'] + gamma_tl*m0 )/proteinsToSynapses - 1]);
        self.ExcessmRNA_transcribed_perSynapticProtein    = (self.params['J_m'] + self.params['J_p']*lambda_m/gamma_tl)/(proteinsToSynapses*lambda_m/gamma_tl)- 1;

        return ratio_synapses_supplied,particles_in_active_transport;

    def SolveTransportEqs(self,N):

        self.__initialization();

        D_m, D_p = self.params['D_m'], self.params['D_p']
        v_m, v_p = self.params['v_m'], self.params['v_p']
        J_m, J_p = self.params['J_m'], self.params['J_p']
        lambda_m, lambda_p   = self.params['lambda_m'], self.params['lambda_p']
        eta_max_m, eta_max_p = self.params['eta_max_m'], self.params['eta_max_p']

        gamma_tl = self.params['gamma_tl'];
        eta_0    = self.params['eta_0'];
        L        = self.params['L'];
        #locals().update(self.params);

        # uptake of protein
        # u_p(p) = 1*tanh(p/epsilon_uptake)
        #
        # uptake of mRNA into translating state: (tuned to satisfy protein demand)
        # u_m(m) = lambda/gamma_tl*tanh(m/epsilon_uptake)

        def u_m(m):
            return np.tanh(m/eta_max_m*eta_0)*eta_max_m;

        def u_p(p):
            return np.tanh(p/eta_max_p*eta_0)*eta_max_p;

        #Setting up equations / matrices
        #
        #Solve for mRNA
        # 0   = D_eff_m * d^2m/dx^2 - v_eff_m*dm/dx - lambda_m*m - u_m(m)
        # m_0 = u_m(m)/lambda_m;
        #
        #Solve for protein
        # 0   = D_eff_p * d^2p/dx^2 - v_eff_p*dp/dx - lambda_p*p - u_p(p) + gamma_tl*m_0
        # p_0 = u_p(p)/lambda_p;

        #y[0] = m
        #y[1] = dm/dx
        def fun_m(x, y):
            return np.vstack((y[1],-1/D_m*(-v_m*y[1]-lambda_m*y[0]-u_m(y[0]) )));

        #y[0] = p
        #y[1] = dp/dx
        def fun_p(x, y):
            m0 = u_m(res_m.sol(x)[0])/lambda_m;
            return np.vstack((y[1],-1/D_p*(-v_p*y[1]-lambda_p*y[0]-u_p(y[0])+gamma_tl*m0)));

        #****************************
        #Boundary conditions:
        #****************************
        #
        # at x=0: J_m = - D_eff_m*dm/dx + v_eff_m*m (influx)
        #
        # at x=L: 0   = - D_eff_m*dm/dx + v_eff_m*m (influx)(no-outflux)
        def bc_m(ya, yb):
            return np.array([-D_m*ya[1]+v_m*ya[0]-J_m, -D_m*yb[1]+v_m*yb[0]]);

        def bc_p(ya, yb):
            return np.array([-D_p*ya[1]+v_p*ya[0]-J_p, -D_p*yb[1]+v_p*yb[0]]);

        x     = np.linspace(0,L,N)
        y_0   = np.zeros((2, x.size))

        res_m = solve_bvp(fun_m, bc_m, x, y_0,max_nodes=1e4);
        if(res_m.status != 0):
            print(res_m.message);
            raise RuntimeError("solve_bvp for p did not converge");
        else:
            print("Max residuals m:"+str(np.max(res_m.rms_residuals)));

        res_p = solve_bvp(fun_p, bc_p, x, y_0,max_nodes=1e4);
        if(res_p.status != 0):
            print(res_p.message);
            raise RuntimeError("solve_bvp for p did not converge");
        else:
            print("Max residuals p:"+str(np.max(res_p.rms_residuals)));

        self.res_m = res_m;
        self.res_p = res_p;

        self.__computeCostItems();

        return res_m,res_p;
