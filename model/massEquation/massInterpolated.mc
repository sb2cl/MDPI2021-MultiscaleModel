%% ---------------------------------------------------------
%% Definition of the class MassInterpolated.
%% ---------------------------------------------------------

class MassInterpolated

  %% ---------------------------------------------------------
  %% Input variables.
  %% ---------------------------------------------------------

  variable mu(
    units = '1/min',
    comment = 'Growth rate.',
    isPlot = false,
    isTex = false
    );
 
  %% ---------------------------------------------------------
  %% Parameters and variables of the model.
  %% ---------------------------------------------------------

  %Parameter mp_p1(value = 18.41,isTex = false);
  %Parameter mp_p2(value = 138.9,isTex = false);
  %Parameter mp_p3(value = 241.3,isTex = false);
  %Parameter mp_mean(value = 0.01756,isTex = false);
  %Parameter mp_std(value = 0.008775,isTex = false);
  %
  %Variable mu_norm(
  %  value = (mu-mp_mean)/mp_std,
  %  isPlot = false,
  %  isTex = false
  %  );

  parameter c_1(
    value = 239089,
    units = 'fg \cdot cell^{-1} \cdot min^2',
    comment = 'First coefficient of mass equation.',
    nameTex = 'c_1',
    reference = '\cite{nobel2020resources}*'
    );

  parameter c_2(
    value = 7432,
    units = 'fg \cdot cell^{-1} \cdot min',
    comment = 'Second coefficient of mass equation.',
    nameTex = 'c_2',
    reference = '\cite{nobel2020resources}*'
    );

  parameter c_3(
    value = 37.06,
    units = 'fg \cdot cell^{-1}',
    comment = 'Third coefficient of mass equation.',
    nameTex = 'c_3',
    reference = '\cite{nobel2020resources}*'
    );

  variable m_p(
    units = 'fg \cdot cell^{-1}',
    comment = 'Total protein mass of the cell.',
    isPlot = false,
    nameTex = 'm_p(\mu)',
    %equationTex = 'm_p = 18.41 ((\mu - 0.01756)/0.008775)^2 + 138.9 (\mu - 0.01756)/0.008775 + 241.3'
    equationTex = 'm_p(\mu) = c_1 \mu^2 + c_2 \mu + c_3',
    isSubstitution = true
    );
 
  %% ---------------------------------------------------------
  %% Equations.
  %% ---------------------------------------------------------
  
  %Equation (m_p == mp_p1*mu_norm*mu_norm + mp_p2*mu_norm + mp_p3, isSubstitution = true);
  equation (m_p == c_1*mu*mu + c_2*mu + c_3, isSubstitution = true);

end MassInterpolated;
