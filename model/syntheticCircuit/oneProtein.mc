%% ---------------------------------------------------------
%% Definition of the class OneProtein.
%% ---------------------------------------------------------

use ./model/wildType/dynamic.mc;

class OneProtein

  extends WildType;

  %% ---------------------------------------------------------
  %% Parameters. 
  %% ---------------------------------------------------------

  parameter lp_A(
    value = 195,
    units = 'aa',
    comment = 'Length of protein A.',
    nameTex = 'l_p^A',
    reference = '**'
    );
  
  parameter dm_A(
    value = 0.16,
    units = 'min^{-1}',
    comment = 'Mean degradation rate of protein A mRNA.',
    nameTex = 'd_m^A',
    reference = '**'
    );
  
  parameter ku_A(
    value = 117.2305,
    range = [6 135],
    units = 'min^{-1}',
    comment = 'Dissotiation rate RBS-ribosome for protein A mRNA.',
    nameTex = 'k_u^A',
    valueTex = '[6 135]',
    reference = '\cite{deSmit2003,Kierzek2001}'
    );
  
  parameter kb_A(
    value = 5.1754,
    range = [3 15],
    units = 'cell \cdot min^{-1} \cdot molec^{-1}',
    comment = 'Association rate RBS-ribosome for protein A mRNA.',
    nameTex = 'k_b^A',
    valueTex = '[3 15]',
    reference = '\cite{deSmit2003,Kierzek2001,SALIS201119}'
    );
  
  parameter N_A(
    value = 1,
    units = 'adim',
    comment = 'Number of copies of gene $A$.',
    nameTex = 'N_A',
    valueTex = '[1 70]',
    reference = '**'
    );
  
  parameter omega_A(
    value = 10, %5.0519,
    units = 'molec \cdot min^{-1} \cdot cell^{-1}',
    comment = 'Average transcription rate for protein $A$.',
    nameTex = '\omega_A',
    valueTex = '[0 5]',
    reference = '**'
    );
  
  %% ---------------------------------------------------------
  %% Variables. 
  %% ---------------------------------------------------------

  variable KkC0_A(
    value = kb_A/(ku_A+ke),
    units = 'cell \cdot molec^{-1}',
    comment = 'Effective RBS affinity of protein A mRNA.',
    isPlot = false,
    nameTex = 'K^A_{C_0}(s)',
    equationTex = 'K^A_{C_0}(s) = k_b^A /(k_u^A + k_e(s))'
    );
  
  variable Emk_A(
    value = 0.62*lp_A/l_e,
    units = 'adim',
    comment = 'Ribosomes density related term for protein A mRNA.',
    isPlot = false,
    nameTex = 'E_m^A',
    equationTex = 'E_m^A = 0.62 l_p^A/l_e'
    );
  
  variable J_A(
    value = Emk_A*omega_A/(dm_A/KkC0_A+mu*r),
    units = 'adim',
    comment = 'Average J value of one protein A gene.',
    title = '$J_A$',
    nameTex = 'J_A(\mu,r,s)',
    equationTex = 'J_A(\mu,r,s) = E_m^A \omega_A / (d_m^A/K_{C_0}^A(s) + \mu r)'
    );
  
  variable m_A(
    start = 0,
    units = 'fg \cdot cell^{-1}',
    comment = 'Total mass of protein $A$ in the cell.',
    title = 'm$_A$',
    nameTex = 'm_A',
    equationTex = '\dot m_A &= \left[m_p(\mu) \frac{N_A J_A(\mu,r,s)}{J_{sum}(\mu,r,s)} - m_A \right] \mu'
  );
  
  %% ---------------------------------------------------------
  %% Equations. 
  %% ---------------------------------------------------------
  
  equation der_m_A == (m_p*N_A*J_A/J_sum - m_A)*mu;
  
  % Override these equations taking into account protein A.
  m_sum_eq.eqn = 'm_sum == m_r + m_nr + m_A';
  J_sum_eq.eqn = 'J_sum == N_r*J_r + N_nr*J_nr + N_A*J_A';

end OneProtein;

%% ---------------------------------------------------------
%% Stand alone model.
%% ---------------------------------------------------------

parameter s(value = 3.6, isTex = false);
parameter m_p(value = 450, isTex = false);

OneProtein wt;

connect(wt__s, s);
connect(wt__m_p, m_p);
