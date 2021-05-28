%% ---------------------------------------------------------
%% Definition of the class WildType.
%% ---------------------------------------------------------

class WildType

  import ./model/wildType/param.mc;

  %% ---------------------------------------------------------
  %% Input variables.
  %% ---------------------------------------------------------

  variable s(
    comment = 'Concentration of substrate in the biorreactor.',
    units = 'g \cdot L^{-1}',
    title = 'Substrate concentration $s$',
    isTex = false
    );

  variable m_p(
    units = 'fg \cdot cell^{-1}',
    comment = 'Total protein mass of the cell.',
    title = 'Cell mass',
    nameTex = 'm_p',
    isTex = false
    );
  
  %% ---------------------------------------------------------
  %% General variables. 
  %% ---------------------------------------------------------

  variable nu(
    value = nu_max*s/(s+K_s),
    units = 'aa \cdot min^{-1}',
    comment = 'Effective translation rate per ribosome.',
    title = '$\nu(s)$',
    nameTex = '\nu(s)',
    equationTex = '\nu(s) = \nu_{max} s / (s+K_s)'
    );
 
  variable mu(
    start = 0.01,
    units = 'min^{-1}',
    comment = 'Specific cell growth rate.',
    title='Specifi cell growth rate.',
    ylabel='[min$^{-1}$]',
    ylim = [0 0.03],
    nameTex = '\mu',
    equationTex = '\mu &= \frac{m_{aa}}{m_p(\mu)} \nu(s) \phi_b(\mu,r,s) \phi_t r_t(m_r)'
    );
 
  variable ke(
    value = nu/l_e,
    units = 'min^{-1}',
    comment = 'Translation initation rate.',
    isPlot = false,
    nameTex = 'k_e(s)',
    equationTex = 'k_e(s) = \nu(s) / l_e'
    );
  
  variable m_sum(
    value = m_r + m_nr,
    units = 'fg \cdot cell^{-1}',
    comment = 'Total protein mass calculated from the mass of individual proteins in the cell',
    title = 'm$_{sum}$',
    plotIn = 'cellModel__m_p',
    isPlot = false,
    isTex = false
    );
  
  variable J_sum(
    value = N_r*J_r + N_nr*J_nr,
    units = 'adim',
    comment = 'Total sum of all the J in the cell.',
    title = '$J_{sum}$',
    nameTex = 'J_{sum}(\mu,r,s)',
    equationTex = 'J_{sum}(\mu,r,s) = \sum_{i=r,nr,A} N_i J_i(\mu,r,s)'
    );
  
  %% ---------------------------------------------------------
  %% Ribosomal proteins.
  %% ---------------------------------------------------------
  
  variable KkC0_r(
    value = kb_r/(ku_r+ke),
    units = 'cell \cdot molec^{-1}',
    comment = 'Effective RBS affinity of ribosomal mRNA.',
    isPlot = false,
    nameTex = 'K^r_{C_0}(s)',
    equationTex = 'K^r_{C_0}(s) = k_b^r/(k_u^r + k_e(s))'
    );
  
  variable Emk_r(
    value = 0.62*lp_r/l_e,
    units = 'adim',
    comment = 'Ribosomes density related term for ribsomal mRNA.',
    isPlot = false,
    nameTex = 'E_{m}^r',
    equationTex = 'E_m^r = 0.62 l_p^r / l_e'
    );
  
  variable J_r(
    value = Emk_r*omega_r/(dm_r/KkC0_r+mu*r),
    units = 'adim',
    comment = 'Average J value of one ribosomal \textit{E.coli.} gene.',
    title = '$J_r$',
    nameTex = 'J_r(\mu,r,s)',
    equationTex = 'J_r(\mu,r,s) = E_m^r \omega_r / (d_m^r/K^r_{C_0}(s) + \mu r)'
    );
 
  variable r_t(
    value = m_r/ribosomeWeight,
    units = 'molec \cdot cell^{-1}',
    comment = 'Number of mature and inmmature ribosomes.',
    title = 'r$_t$',
    isPlot = false,
    nameTex = 'r_t(m_r)',
    equationTex = 'r_t(m_r) = m_r/r_w'
    );
  
  variable phi_b(
    value = J_sum/(1+J_sum),
    units = 'adim',
    comment = 'Fraction of translating ribosomes of $\phi_t r_t(m_r)$.',
    title = '$\phi_b$',
    isPlot = false,
    nameTex = '\phi_b(\mu,r,s)',
    equationTex = '\phi_b(\mu,r,s) = J_{sum}(\mu,r,s)/(1 + J_{sum}(\mu,r,s))'
    );
  
  variable r(
    start = 350,
    units = 'molec \cdot cell^{-1}',
    comment = 'Free mature ribosomes in the cell.',
    title = 'r',
    nameTex = 'r',
    equationTex = 'r &= \frac{\phi_t r_t(m_r)}{1+J_{sum}(\mu,r,s)}'
    );

  variable m_r(
    start = 10,
    units = 'fg \cdot cell^{-1}',
    comment = 'Total mass of ribosomal proteins in the cell.',
    title = 'm$_r$',
    nameTex = 'm_r',
    equationTex = '\dot m_r &= \left[m_p(\mu) \frac{N_r J_r(\mu,r,s)}{J_{sum}(\mu,r,s)} - m_r \right] \mu'
  );
  
  %% ---------------------------------------------------------
  %% Non ribosomal proteins.
  %% ---------------------------------------------------------
  
  variable KkC0_nr(
    value = kb_nr/(ku_nr+ke),
    units = 'cell \cdot molec^{-1}',
    comment = 'Effective RBS affinity of non-ribosomal mRNA.',
    isPlot = false,
    nameTex = 'K^{nr}_{C_0}(s)',
    equationTex = 'K^{nr}_{C_0} = k_b^{nr}/(k_u^{nr} + k_e(s))'
    );
  
  variable Emk_nr(
    value = 0.62*lp_nr/l_e,
    units = 'adim',
    comment = 'Ribosomes density related term for non-ribsomal mRNA.',
    isPlot = false,
    nameTex = 'E_{m}^{nr}',
    equationTex = 'E_{m}^{nr} = 0.62 l_p^{nr} / l_e'
    );
  
  variable J_nr(
    value = Emk_nr*omega_nr/(dm_nr/KkC0_nr+mu*r),
    units = 'adim',
    comment = 'Average J value of one non-ribosomal \textit{E.coli.} gene.',
    title = '$J_{nr}$',
    nameTex = 'J_{nr}(\mu,r,s)',
    equationTex = 'J_{nr}(\mu,r,s) = E_m^{nr} \omega_{nr} / (d_m^{nr}/K_{C_0}^{nr}(s) + \mu r)'
    );
  
  variable m_nr(
    start = 10,
    units = 'fg \cdot cell^{-1}',
    comment = 'Total mass of non ribosomal proteins in the cell.',
    title = 'm$_{nr}$',
    nameTex = 'm_{nr}',
    equationTex = '\dot m_{nr} &= \left[m_p(\mu) \frac{N_{nr} J_{nr}(\mu,r,s)}{J_{sum}(\mu,r,s)} - m_{nr}\right] \mu'
  );
  
  %% ---------------------------------------------------------
  %% Equations.
  %% ---------------------------------------------------------
  
  equation r == phi_t*r_t/(1+J_sum);
  equation mu == m_aa/m_p*nu*phi_b*phi_t*r_t;
  equation der_m_r == (m_p*N_r*J_r/J_sum - m_r)*mu;
  equation der_m_nr == (m_p*N_nr*J_nr/J_sum - m_nr)*mu;

end WildType;

%% ---------------------------------------------------------
%% Stand alone model.
%% ---------------------------------------------------------

parameter s(value = 3.6, isTex = false);
parameter m_p(value = 450, isTex = false);

WildType wt;

connect(wt__s, s);
connect(wt__m_p, m_p);