(* Created with the Wolfram Language for Students - Personal Use Only : www.wolfram.com *)
{A -> (#2^(-2) + As[#1, #2]*#2^2 + (2*\[Xi][#1])/#2 + \[Xi][#1]^2 + 
    Log[#2^(-1)]*((-2*\[ScriptCapitalB]^2*#2^2)/3 + 
      (4*\[ScriptCapitalB]^2*#2^3*\[Xi][#1])/3 - 2*\[ScriptCapitalB]^2*#2^4*
       \[Xi][#1]^2) - 2*Derivative[1][\[Xi]][#1] & ), 
 B -> ((-2*Log[#1])/3 - (2*#2)/(3*#1) + Bs[#1, #2]*#2^4 + 
    (#2^2*(1 + 2*#1*\[Xi][#1]))/(3*#1^2) - 
    (2*#2^3*(1 + 3*#1*\[Xi][#1] + 3*#1^2*\[Xi][#1]^2))/(9*#1^3) + 
    Log[#2^(-1)]*((\[ScriptCapitalB]^2*#2^4)/3 - 
      (4*#2^5*(\[ScriptCapitalB]^2 + 15*\[ScriptCapitalB]^2*#1*\[Xi][#1]))/
       (45*#1) + (2*#2^6*(2*\[ScriptCapitalB]^2 + 10*\[ScriptCapitalB]^2*#1*
          \[Xi][#1] + 75*\[ScriptCapitalB]^2*#1^2*\[Xi][#1]^2))/
       (45*#1^2)) & ), dplus[1][B] -> 
  (-1/(3*#1) + #2/(3*#1^2) + (#2^2*(-1 - #1*\[Xi][#1]))/(3*#1^3) + 
    Log[#2^(-1)]*((-2*\[ScriptCapitalB]^2*#2^3)/3 + 2*\[ScriptCapitalB]^2*
       #2^4*\[Xi][#1] + #2^5*((2*\[ScriptCapitalB]^2)/(45*#1^2) - 
        4*\[ScriptCapitalB]^2*\[Xi][#1]^2)) + #2^3*dplus[1][Bs[#1, #2]] & ), 
 S -> (#1^(1/3)/#2 - #2/(9*#1^(5/3)) + #2^3*Ss[#1, #2] + 
    (1 + 3*#1*\[Xi][#1])/(3*#1^(2/3)) + (#2^2*(5 + 9*#1*\[Xi][#1]))/
     (81*#1^(8/3)) + Log[#2^(-1)]*((2*\[ScriptCapitalB]^2*#2^4)/
       (45*#1^(2/3)) - (4*\[ScriptCapitalB]^2*#2^5*(1 + 6*#1*\[Xi][#1]))/
       (135*#1^(5/3))) & ), dplus[1][S] -> 
  (#1^(1/3)/(2*#2^2) + (10*#2)/(81*#1^(8/3)) + (1 + 3*#1*\[Xi][#1])/
     (3*#1^(2/3)*#2) + (-1 + 2*#1*\[Xi][#1] + 3*#1^2*\[Xi][#1]^2)/
     (6*#1^(5/3)) + Log[#2^(-1)]*(-(\[ScriptCapitalB]^2*#1^(1/3)*#2^2)/3 + 
      (2*#2^3*(-2*\[ScriptCapitalB]^2 + 15*\[ScriptCapitalB]^2*#1*\[Xi][#1]))/
       (45*#1^(2/3)) + (#2^4*(\[ScriptCapitalB]^2 + 36*\[ScriptCapitalB]^2*#1*
          \[Xi][#1] - 135*\[ScriptCapitalB]^2*#1^2*\[Xi][#1]^2))/
       (135*#1^(5/3))) + #2^2*dplus[1][Ss[#1, #2]] & )}
