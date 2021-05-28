classdef EquationClass < ModelPartClass
  %% EQUATIONCLASS This class defines equations for OneModel.
  %

  properties
    % [char] Char array equation.
    eqn
    % [sym] Simbolic equation.
    eqnSym
    % bool Is the equation just a substitution of variables?
    isSubstitution
  end % properties

  properties (Dependent)
    % [sym] Left part of the equation.
    left
    % [sym] Right part of the equation.
    right 
    % [char] Right part of the equation.
    rightStr
    % [sym] Variables of the equation.
    vars
    % {[char]} Variables of the equation.
    varsStr
    % [sym] Derivatives variables in the equation.
    ders
    % {[char]} Derivatives variables in the equation.
    dersStr
    % [bool] Is the equation algebraic?
    isAlgebraic
  end % properties (Dependent)

  methods 

    function [obj] = EquationClass(mc, name)
      %% Constructor of EquationClass.
      %
      % param: mc   OneModel object.
      %        name Name of the equation. 

      obj = obj@ModelPartClass(mc);

      % Is the equation is defined with an empty name?
      if isempty(name)
        % Generate a name for it.
        name = sprintf('eq_%d',length(mc.eqns));
      end

      % Define the name taking into account the namespace.
      name = [obj.namespace name];

      obj.name = name;
      obj.isSubstitution = false;

    end % EquationClass

    function [out] =  getFreeVars(obj,knownVars)
      %% GETFREEVARS Return variables that are free in the equation.
      %
      % param: knownVars [sym] Array with known vars.
      %
      % return: out [sym] Free vars in the equation.

      vars = obj.vars;

      out = sym([]);

      for i = 1:length(vars)
        isFree = true;
        for j = 1:length(knownVars)
          if vars(i) == knownVars(j)
            isFree = false;
            break;
          end
        end
        if isFree
          out(end+1) = vars(i);
        end
      end
      
    end % getFreeVars

    function [] = set.eqn(obj,eqn)
      %% SET.EQN Set interface for eqn propierty.
      %
      % param: eqn
      %
      % return: void

      if ~ischar(eqn) 
        error('eqn must be a char array.');
      end

      % Check if '= 'is used instead of '=='.
      expression = '(?<!=)=(?!=)';
      [tokens,matches] = regexp(eqn,expression,'tokens','match');
      if ~isempty(matches)
          error('''='' symbol is used for variable assignation, use ''=='' instead in Equations.');
      end
      
      % Check if there is one and only one '=='.
      expression = '==';
      [tokens,matches] = regexp(eqn,expression,'tokens','match');
      if isempty(matches)
          error('there must be one ''=='' in the Equation.');
      elseif length(matches) > 1
          error('there must be only one ''=='' in the Equation.');
      end

      % Take into account the namespace.
      vars = StrSymbolic.symvar(eqn);
      
      for i = 1:length(vars)
          % Check if the var does not have a namespace defined by the user.
          if ~contains(vars{i},'__')
              % Add the namespace.
              varsNamespace{i} = [obj.namespace vars{i}];
          else
              varsNamespace{i} = vars{i};
          end
      end
      
      eqn = StrSymbolic.subs(eqn,vars,varsNamespace);

      obj.eqn = eqn;
      obj.eqnSym = str2sym(eqn);
      
    end % set.eqn

    function [] = set.isSubstitution(obj,isSubstitution)
      %% SET.ISSUBSTITUTION Set interface for isSubstitution propierty.
      %
      % param: isSubstitution
      %
      % return: void
      
      if ~islogical(isSubstitution)
        error('isSubstitution must be logical.');
      end

      obj.isSubstitution = isSubstitution;
      
    end % set.isSubstitution

    function [out] =  get.left(obj)
      %% GET.LEFT Get left equal part of the equation.
      %
      % return: out Simbolic left part of the equation.

      ind_equal = strfind(obj.eqn,'=');

      out = obj.eqn(1:ind_equal(1)-1);

      out = str2sym(out).';
    end % get.left

    function [out] =  get.right(obj)
      %% GET.RIGHT Get right equal part of the equation.
      %
      % return: out Simbolic right part of the equation.

      ind_equal = strfind(obj.eqn,'=');

      out = obj.eqn(ind_equal(2)+1:end);

      out = str2sym(out).';
    end % get.right

    function [out] =  get.rightStr(obj)
      %% GET.RIGHT Get right equal part of the equation.
      %
      % return: out [char] String right part of the equation.

      ind_equal = strfind(obj.eqn,'=');

      out = obj.eqn(ind_equal(2)+1:end);

    end % get.right


    function [out] = get.vars(obj)
      %% GET.VARS Variables of the equation.
      %
      % return: out [sym] Variables of the equation.

      out = symvar(obj.eqnSym);
      
    end % get.vars

    function [out] = get.varsStr(obj)
      %% GET.VARSSTR Variables of the equation.
      %
      % return: out {[char]} Variables of the equation.

      out = StrSymbolic.symvar(obj.eqn);
      
    end % get.vars

    function [out] = get.ders(obj)
      %% GET.DERS Get the variables that have a derivative in the equation.
      %
      % return: out [sym] Derivative variables.

      vars = obj.vars;

      out = sym([]);

      for i = 1:length(vars)
        [matches] = regexp(char(vars(i)), 'der_', 'tokens');

        if ~isempty(matches)
          out(end+1) = sym(erase(char(vars(i)), 'der_'));
        end
      end

    end % get.ders

    function [out] = get.dersStr(obj)
      %% GET.DERS Get the variables that have a derivative in the equation.
      %
      % return: out {[char]} Derivative variables.

      vars = obj.vars;

      out = {};

      for i = 1:length(vars)
        [matches] = regexp(char(vars(i)), 'der_', 'tokens');

        if ~isempty(matches)
          out{end+1} = erase(char(vars(i)), 'der_');
        end
      end

    end % get.ders

    function [out] = get.isAlgebraic(obj)
      %% GET.ISALGEBRAIC Is the equation algebraic?
      %
      % return: out bool isAlgebraic

      vars = obj.vars;

      expression = 'der_(\w*)';

      for i = 1:length(vars)
        [tokens,matches] = regexp(char(vars(i)),expression,'tokens','match');
        if ~isempty(tokens)
          out = false;
          return
        end
      end
      
      out = true;

    end % get.isAlgebraic

  end % methods

end % classdef
