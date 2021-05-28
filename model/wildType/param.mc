%% ---------------------------------------------------------
%% General parameters.
%% ---------------------------------------------------------

parameter K_s( 
    comment = 'Half activation threshold of growth rate.',
    units = 'g \cdot L^{-1}',
    value = 0.1802,
    nameTex = 'K_s',
    reference = '\cite{Zhuang2013}'
    % reference = 'Dynamic strain scanning optimization: An efficient strain design strategy for balanced yield, titer, and productivity. DySScO strategy for strain design'
    );

parameter nu_max(
    value = 1260,
    units = 'aa \cdot min^{-1}',
    comment = 'Maximum effective translation rate per ribosome.',
    nameTex = '\nu_{max}',
    reference = '\cite{Milo2009}'
    );
 
parameter m_aa(
  value = 182.6e-9,
  units = 'fg \cdot aa^{-1}',
  comment = 'Average aminoacid mass.',
  nameTex = 'm_{aa}',
  reference = '\cite{Sundararaj2004}',
  valueTex = '$182.6 \cdot 10^{-9}$'
  );

parameter l_e(
  value = 25,
  units = 'aa',
  comment = 'Ribosome occupancy length.',
  nameTex = 'l_e',
  reference = 'estimated \cite{Fernandes2017,Eriksen2017,Picard2013,Siwiak2013}'

  );

parameter phi_t(
  value = 0.7796, % 0.8
  units = 'adim',
  comment = 'Fraction of mature available ribosomes relative to the total.',
  nameTex = '\phi_t',
  reference = '\cite{doi:10.1111/febs.13258,Bremer:2008}*'
  );

%% ---------------------------------------------------------
%% RIBOSOMAL GENES PARAMETERS.
%% ---------------------------------------------------------

parameter lp_r(
  value = 195,
  units = 'aa',
  comment = 'Mean length of ribosomal proteins.',
  nameTex = 'l_p^{r}',
  reference = 'calculated from \cite{Hausser2019}'
  );

parameter dm_r(
  value = 0.16,
  units = 'min^{-1}',
  comment = 'Mean degradation rate of ribosomal mRNA.',
  nameTex = 'd_m^{r}',
  reference = 'calculated from \cite{Hausser2019}'
  );

parameter ku_r(
  value = 135, % 132.68
  range = [6 135],
  units = 'min^{-1}',
  comment = 'Dissotiation rate RBS-ribosome for ribosomal mRNA.',
  nameTex = 'k_u^{r}',
  reference = '\cite{nobel2020resources}*'
  );

parameter kb_r(
  value = 8.8530, % 7.82
  range = [3 15],
  units = 'cell \cdot min^{-1} \cdot molec^{-1}',
  comment = 'Association rate RBS-ribosome for ribosomal mRNA.',
  nameTex = 'k_b^r',
  reference = '\cite{nobel2020resources}*'
  );

parameter N_r(
  value = 55,
  units = 'adim',
  comment = 'Number of proteins that make up a ribosome.',
  nameTex = 'N_r',
  reference = '\cite{nobel2020resources}'
  );

parameter omega_r(
  value = 4.8658,
  units = 'molec \cdot min^{-1} \cdot cell^{-1}',
  comment = 'Average transcription rate for ribosomal proteins.',
  nameTex = '\omega_r',
  reference = '\cite{nobel2020resources}*'
  );

parameter ribosomeWeight(
  %value = m_r/(N_r*lp_r*m_aa), % Esto da 0.002 fg cuando con los datos de Bremer debería ser 0.0045 fg.
  value = 0.0045,
  units = 'fg',
  comment = 'Weight of a ribosome.',
  nameTex = 'r_w',
  reference = '\cite{Bremer:2008}'
  );

%% ---------------------------------------------------------
%% NON-RIBOSOMAL GENES PARAMETERS.
%% ---------------------------------------------------------

parameter lp_nr(
  value = 333,
  units = 'aa',
  comment = 'Mean length of non-ribosomal proteins.',
  nameTex = 'l_p^{nr}',
  reference = 'calculated from \cite{Hausser2019}'
  );

parameter dm_nr(
  value = 0.2,
  units = 'min^{-1}',
  comment = 'Mean degradation rate of non-ribosomal mRNA.',
  nameTex = 'd_m^{nr}',
  reference = 'calculated from \cite{Hausser2019}'
  );

parameter ku_nr(
  value = 6.1297, % 6
  range = [6 135],
  units = 'min^{-1}',
  comment = 'Dissotiation rate RBS-ribosome for non-ribosomal mRNA.',
  nameTex = 'k_u^{nr}',
  reference = '\cite{nobel2020resources}*'
  );

parameter kb_nr(
  value = 14.9971, % 15
  range = [3 15],
  units = 'cell \cdot min^{-1} \cdot molec^{-1}',
  comment = 'Association rate RBS-ribosome for non-ribosomal mRNA.',
  nameTex = 'k_b^{nr}',
  reference = '\cite{nobel2020resources}*'
  );

parameter N_nr(
  value = 1735,
  units = 'adim',
  comment = 'Number of non ribosomal proteins expressed at one time.',
  nameTex = 'N_{nr}',
  reference = '\cite{nobel2020resources}'
  );

parameter omega_nr(
  value = 0.03,
  units = 'molec \cdot min^{-1} \cdot cell^{-1}',
  comment = 'Average transcription rate for non ribosomal proteins.',
  nameTex = '\omega_{nr}',
  reference = '\cite{nobel2020resources}*'
  );
