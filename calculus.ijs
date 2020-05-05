NB. Calculus: derivative, partial derivative, secant slope
NB. deriv, pderiv, sslope
NB. Interface is like d., D., D:   sample sentences at end of file

cocurrent 'jcalculus'

NB. Derivative operator, like d.
deriv =: 2 : 0
assert. '' -: $n   NB. n must be an atom
order =. <. n
assert. n = order
assert. 3 = 4!:0 <'u'   NB. Must be in user's namespace
ufix =. u. f.  NB. u in user's namespace
resstg =. 5!:5 <'ufix'   NB. User's verb in character form
while. order do.   NB. If the order is not 1 or _1, keep increasing the order.
NB. Find the derivative of the AR.
  try.
    resstg =. 0:`derivstg`intstg@.(*order) resstg
  catch.
    break.  NB. keep resstg/order from the last successful calculation
  end.
  order =. order - *order   NB. Decrement # orders remaining
end.
NB. resstg has the string form and order is the number of derivs unfinished.  Return appropriately
if. order < 0 do. 13!:8 (3)  end.  NB. domain error if we don't know the integral
NB. Parse the verb to create a verb
NB. Return the verb itself if it's the right order, otherwise error.  We could do secant approx instead of error
if. order=0 do. resstg vnofu else. 13!:8 (3) end.
)

NB. Routine to handle unevaluable derivatives by approximation, or D:
NB. x is amount to change eval point (0 is replaced by 1e_7)
NB. y is eval point
NB. Dyad only
derivsecant =: 2 : 0
NB. Get function to use: u. or higher-order secant 
if. n = 1 do. func =. u.@] else. func =. u. derivsecant_jcalculus_ (n-1) end.
NB. x must be an atom or conform to shape of y
if. 0=#@$x do. x =. ($y)$x end.  NB. replicate atom
assert. x -:&$ y  NB. shapes must agree
x =. x + 1e_7 * 0 = x  NB. replace 0 by epsilon
newy =. y +"(#@$y) (,~$x) $ (#~  1 j. #) ,x  NB. array of moved points, each an array with 1 nonzero
f0 =. x func y  NB. the function at the initial point
((x func"(#@$y) newy) -"(#@$f0) f0) % x   NB. evaluate function at moved points, calc slope.  x used only for higher orders
)

NB. x and y are linear reps of functions
NB. Result is linear rep of their product
NB. Someday we may combine polynomials etc
ftymes =: 4 : 0
'((',x,')*(',y,'))'
)
fmp =: 4 : 0  NB. same but matrix  product
'((',x,') +/ . * (',y,'))'
)

NB. x and y are linear reps of functions
NB. Result is linear rep of their sum
NB. Someday we may combine polynomials etc
fplus =: 4 : 0
'((',x,')+(',y,'))'
)

NB. x and y are linear reps of functions
NB. Result is linear rep of their difference
NB. Someday we may combine polynomials etc
fminus =: 4 : 0
'((',x,')-(',y,'))'
)

NB. x and y are linear reps of functions
NB. Result is linear rep of their quotient
NB. Someday we may combine polynomials etc
fdiv =: 4 : 0
'((',x,')%(',y,'))'
)

NB. Operand y of AR x, as an AR
opar =: 4 : 0"_ 0
y { 1 {:: x
)

NB. true if operand y of x is a noun
opisnoun =: 4 : 0"_ 0
(,'0') -: (1;y;0) {:: x
)

NB. Operand y of AR x, as a string
opstr =: 4 : 0"_ 0
op =. (x opar y) 5!:0
5!:5 <'op'
)

NB. AR of u (noun or verb)
arofu =: 1 : '5!:1 <''u'''

NB. convert AR u to verb/noun
vnofaru =: 5!:0

NB. Convert string u to verb/noun
vnofu =: 1 : 0
". 'Xvcv98df9d =. ' , u  NB. This sets noun result
if. 0 ~: 4!:0 <'Xvcv98df9d' do. Xvcv98df9d f. end.  NB. This sets verb result
)

NB. Convert verb/noun u to string.  Must fix first in case u is a name
strofu =: f. 1 : '5!:5 <''u'''

NB. x and y are string forms of verb
NB. Result is x@y in string form
atops =: '(',[,')@(',],')'"_

NB. y is f;g;h, string forms
NB. Result is fork in string form
forks =: '(',(0&{::),') (',(1&{::),') (',(2&{::),')'"_

NB. x is string form, y is numeric rank
NB. Result is x"y in string form
ranks =: '(',[,')"',":@]

NB. x is string form, y is numeric power
NB. Result is x^:y in string form
powops =: '(',[,')^:',":@]

arofstringu =: vnofu arofu

NB. symbolic derivative, in string form
NB. y is string form of a verb
NB. Result is string form of the derivative verb
NB. If we don't know the derivative, we raise an error
NB. We pass around strings rather than ARs for convenience of coding & table maintenance
derivstg =: 3 : 0
yar =. y arofstringu
NB. Top-down recursion through the AR
if. 1 (>: L.) yar do.
  NB. Unboxed contents.  It's a primitive or a name
  NB. Look it up in the derivative table
  if. (#primvb) > primno =. primvb i. yar do. primno {:: primderiv return. end.
  NB. If undifferentiable primitive, fail
  13!:8 (3)
end.
NB. Not a primitive or named verb.  Must be a 2-box list.  The first box indicates the (possibly invisible)
NB. modifier; process it if we know it
yar =. > yar  NB. discard outer boxing
select. {.yar
case. <,'0' do.  NB. noun - should not occur
case. <,'2' do.  NB. hook - not supported
case. <,'3' do.  NB. fork
  'farg garg harg' =. yar&opstr&.> 0 1 2
  if. (<'[:') -: yar opar 0 do. derivstg garg atops harg return. end.   NB. [: g h
  if. yar opisnoun 0 do. derivstg forks (farg ranks 0);garg;harg return. end. NB. n v v
  df =. derivstg farg [ dh =. derivstg harg
  select. garg
  case. <,'+' do. df fplus dh return.
  case. <,'-' do. df fminus dh return.
  case. <,'*' do. (df ftymes harg) fplus (farg ftymes dh) return.  NB. product rule
  case. <,',' do. '(',df,') , (',dh,')' return.
  case. <,'%' do. ((df ftymes harg) fminus (farg ftymes dh)) fdiv '*:' atops harg return.  NB. quotient rule
  case. <,'^' do. NB. no support for f ^ h
  end.

case. <,'&' do.  NB. & - first check for bonded constant
  if. yar opisnoun 0 do.   NB. m&v  (v must be a verb)
    nounarg =. (yar opar 0) vnofaru  NB. m
    verbarg =. (yar opar 1)  NB. v as an AR
    if. verbarg -: <'p.' do.   NB. Handle p. as a special case
      if. 0 (< L.) nounarg do. nounarg =. p. nounarg end.  NB. convert multiplier or multinomial form to coeffs
      (": (* #\)@}. nounarg),'&p.' return.   NB. take the derivative and return
    end.
    if. *@#@$ nounarg do. 13!:8 (3) end. NB. except for p. only atomic m is recognized
    NB. Look it up in the derivative table, execute it on m, return string result
    if. (#primmandv) > primno =. primmandv i. verbarg do.
      (nounarg 1 : (primno {:: primmandvderiv)) strofu return.
    end.
    NB. If undifferentiable primitive, fail
    13!:8 (3)
  elseif. yar opisnoun 1 do.   NB. u&n  (u must be a verb)
    nounarg =. (yar opar 1) vnofaru  NB. n
    verbarg =. (yar opar 0)  NB. u as an AR
    if. *@#@$ nounarg do. 13!:8 (3) end. NB. except for p. only atomic m is recognized
    NB. Look it up in the derivative table, execute it on m, return string result
    if. (#primuandn) > primno =. primuandn i. verbarg do.
      (nounarg 1 : (primno {:: primuandnderiv)) strofu return.
    end.
    NB. If undifferentiable primitive, fail
    13!:8 (3)
  end.
  NB. fallthrough cases are u&v

fcase. ;:'&:@@:' do.  NB. chain rule, including fallthrough from &
  if. yar opisnoun 1 do. 13!:8 (3) end.  NB. v must be a verb
  uop =. yar opstr 0   NB. the verb as a string
  vop =. yar opstr 1   NB. the verb as a verb
  (derivstg vop) ftymes (derivstg uop) atops vop return.

case. <'^:' do.  NB. power
  if. yar opisnoun 0 do. 13!:8 (3) end.  NB. u must be a verb
  if. -. yar opisnoun 1 do. 13!:8 (3) end.  NB. v must be a noun
  nop =. (yar opar 1) vnofaru  NB. extract the noun
  if. '' -.@-: $nop do. 13!:8 (3) end.  NB. must be an atom
  if. nop ~: <.nop do. 13!:8 (3) end. NB. must be integral
  uop =. (yar opar 0) vnofaru   NB. the verb as a verb
  if. nop < 0 do.  NB. If power is negative, replace the function with its inverse
    uop =. ((< (1;0)&{:: yar) 5!:0) b. _1  NB. inverse as a string
    nop =. -nop   NB. # applications always positive
  else.  NB. power is positive, extract verb as a string
    uop =. uop strofu  NB. the verb as a string
  end.
  select. nop
  case. 1 do. derivstg uop
  case. 0 do. '1"0'
  case. do. derivstg (uop powops <:nop) atops uop  NB. Large powers become u^:(<:n)@u
  end.
  return.

case. <,'"' do.  NB. rank
  if. -. yar opisnoun 1 do. 13!:8 (3) end.  NB. v must be a noun
  nop =. (yar opar 1) vnofaru   NB. extract the noun
  uop =. (yar opar 0) vnofaru  NB. v, either verb or noun
  if. yar opisnoun 0 do.
    if. nop = #@$ uop do. '0"0' return. end.  NB. value"0 or equivalent
  else. (derivstg uop strofu) ranks nop return.
  end.

case. <,'~' do.  NB. reflexive - only a few verbs supported
  NB. Look it up in the derivative table
  if. (#primrefvb) > primno =. primrefvb i. 1 {:: yar do. primno {:: primrefderiv return. end.
end.
13!:8 (3)   NB. Unknown or unprocessed modifier, fail
)

NB. Canned table of derivatives for the primitive verbs
'primvb primderiv' =: <"1 |: (({. ; deb@:}.)~ i.&' ');._2 (0 : 0)
-  _1"0
-. _1:0
<: 1"0
>: 1"0
[ 1"0
] 1"0
+: 2"0
*: +:
-: 1r2"0
o. 1p1"0
% -@%@*:
%: -:@%@%:
^ ^
^. %
j. 0j1"0
r. j.@r.
0: 0"0
1: 0"0
2: 0"0
3: 0"0
4: 0"0
5: 0"0
6: 0"0
7: 0"0
8: 0"0
9: 0"0
_1: 0"0
_2: 0"0
_3: 0"0
_4: 0"0
_5: 0"0
_6: 0"0
_7: 0"0
_8: 0"0
_9: 0"0
_: 0"0
__: 0"0
)

NB. Canned table of derivatives for m&v
NB. The derivative will be interpreted as an adverb, where m is the m from m&v
NB. That adverb might itself execute another verb, eg for o.
'primmandv primmandvderiv' =: <"1 |: (({. ; deb@:}.)~ i.&' ');._2 (0 : 0)
+ 1"0
* m"0
- _1"0
% (-m)&%@*:
%: (%m)&*@(^&(<:%m))
^. (%^.m)&%
^ (^.m)&*@(m&^)
o. m circlederiv
)

NB. Canned table of derivatives for u&n
NB. The derivative will be interpreted as an adverb, where m is the n from u&n
'primuandn primuandnderiv' =: <"1 |: (({. ; deb@:}.)~ i.&' ');._2 (0 : 0)
+ 1"0
* m"0
- 1"0
% (%m)"0
^. (-^.m)&%@(* *:@^.)
^ m&*@(^&(m-1))
)

NB. Canned table of derivatives for reflexives
'primrefvb primrefderiv' =: <"1 |: (({. ; deb@:}.)~ i.&' ');._2 (0 : 0)
+ 2"0
* +:
- 0"0
^. 0"0
% 0"0
^ (^~ * >:@^.)"0
)

NB. Derivatives of m&o., with verb results
circlederiv =: 1 : 0
select. m
case. 0 do. - % 0&o.
case. 1 do. 2&o.
case. 2 do. -@(1&o.)
case. 3 do. %@*:@(2&o.)
case. 5 do. 6&o.
case. 6 do. 5&o.
case. 7 do. %@*:@(6&o.)
case. do. 13!:8 (3) NB. if unknown type, fail
end.
)


NB. symbolic integration, in string form
NB. y is string form of a verb
NB. Result is string form of the integral verb (without the constant)
NB. If we don't know the integral, we raise an error
NB. We pass around strings rather than ARs for convenience of coding & table maintenance
intstg =: 3 : 0
yar =. y arofstringu
NB. If y is a polynomial, integrate it.  This will handle all combinations producing a polynomial
if. # ypoly =. topoly yar do.  (0j17 ": 0 , (% #\) ypoly),'&p.' return. end.
NB. Top-down recursion through the AR
if. 1 (>: L.) yar do.
  NB. Unboxed contents.  It's a primitive or a name
  NB. Look it up in the derivative table
  if. (#intprimvb) > primno =. intprimvb i. yar do. primno {:: primint return. end.
  NB. If undifferentiable primitive, fail
  13!:8 (3)
end.
NB. Not a primitive or named verb.  Must be a 2-box list.  The first box indicates the (possibly invisible)
NB. modifier; process it if we know it
yar =. > yar  NB. discard outer boxing
select. {.yar
case. <,'0' do.  NB. noun - should not occur
case. <,'2' do.  NB. hhok - not supported
case. <,'3' do.  NB. fork
  'farg garg harg' =. yar&opstr&.> 0 1 2
  if. (<'[:') -: yar opar 0 do. 13!:8 (3) end.   NB. [: g h
  if. yar opisnoun 0 do. 13!:8 (3) return. end. NB. n v v
  pf =. topoly farg [ ph =. topoly harg
  select. garg
  case. <,'+' do. if. pf *.&(*&#) ph do. intstg pf pplus ph else. forks (intstg pf);garg;(intstg ph) end. return.
  case. <,'-' do. if. pf *.&(*&#) ph do. intstg pf pminus ph else. forks (intstg pf);garg;(intstg ph) end. return.
  case. <,'*' do. if. pf *.&(*&#) ph do. intstg pf ptymes ph end. return.
  end.

case. <,'&' do.  NB. & - first check for bonded constant
  if. yar opisnoun 0 do.   NB. m&v  (v must be a verb)
    nounarg =. (yar opar 0) vnofaru  NB. m
    verbarg =. (yar opar 1)  NB. v as an AR
    if. *@#@$ nounarg do. 13!:8 (3) end. NB. except for p. handled above, only atomic m is recognized
    NB. Look it up in the integral table, execute it on m, return string result
    if. (#primintmandvvb) > primno =. primintmandvvb i. verbarg do.
      (nounarg 1 : (primno {:: primintmandv)) strofu return.
    end.
    NB. If undifferentiable primitive, fail
    13!:8 (3)
  elseif. yar opisnoun 1 do.   NB. u&n  (u must be a verb)
    nounarg =. (yar opar 1) vnofaru  NB. n
    verbarg =. (yar opar 0)  NB. u as an AR
    if. *@#@$ nounarg do. 13!:8 (3) end. NB. except for p. only atomic m is recognized
    NB. Look it up in the derivative table, execute it on m, return string result
    if. (#primintuandnvb) > primno =. primintuandnvb i. verbarg do.
      (nounarg 1 : (primno {:: primintuandn)) strofu return.
    end.
    NB. If undifferentiable primitive, fail
    13!:8 (3)
  end.

  NB. fallthrough cases are u&v
fcase. ;:'&:@@:' do.  NB. f@g (verb only)
  if. # upoly =. topoly yar opar 0 do.  NB. if u is a polynomial
    if. (-: 1 {.~ -@#) upoly do.
      NB. Handle ^&m for the cases we know
      if. 1 = #upoly do. ']' return. end.  NB. 1@g, just like 1
      vop =. yar opstr 1   NB. v as a string
      if. 2 = #upoly do. intstg vop return. end.
      if. (#primintpatopvb) > primno =. primintpatopvb i. <vop do.
        ((<: #upoly) 1 : (primno {:: primintpatopint)) strofu return.
      end.
    end.
    NB. If undifferentiable primitive, fail
    13!:8 (3)
  end.

case. <'^:' do.  NB. power
  if. yar opisnoun 0 do. 13!:8 (3) end.  NB. u must be a verb
  if. -. yar opisnoun 1 do. 13!:8 (3) end.  NB. v must be a noun
  nop =. (yar opar 1) vnofaru  NB. extract the noun
  if. '' -.@-: $nop do. 13!:8 (3) end.  NB. must be an atom
  if. nop ~: <.nop do. 13!:8 (3) end. NB. must be integral
  uop =. (yar opar 0) vnofaru   NB. the verb as a verb
  if. nop < 0 do.  NB. If power is negative, replace the function with its inverse
    uop =. ((< (1;0)&{:: yar) 5!:0) b. _1  NB. inverse as a string
    nop =. -nop   NB. # applications always positive
  else.  NB. power is positive, extract verb as a string
    uop =. uop strofu  NB. the verb as a string
  end.
  select. nop
  case. 1 do. intstg uop
  case. 0 do. '0 0 0.5&p.'
  case. do. intstg (uop powops <:nop) atops uop  NB. Large powers become u^:(<:n)@u
  end.
  return.

case. <,'"' do.  NB. rank
  if. -. yar opisnoun 1 do. 13!:8 (3) end.  NB. v must be a noun
  nop =. (yar opar 1) vnofaru   NB. extract the noun
  uop =. (yar opar 0) vnofaru  NB. v, either verb or noun
  if. yar opisnoun 0 do.
    if. nop = #@$ uop do. (":uop),'&*' return. end.  NB. value"0 or equivalent
  else. (intstg uop strofu) ranks nop return.
  end.

case. <,'~' do.  NB. reflexive - only a few verbs supported
  NB. Look it up in the integral table
  if. (#intrefvb) > primno =. intrefvb i. 1 {:: yar do. primno {:: refint return. end.
end.
13!:8 (3)   NB. Unknown or unprocessed modifier, fail
)


NB. Canned table of integrals for the primitive verbs
'intprimvb primint' =: <"1 |: (({. ; deb@:}.)~ i.&' ');._2 (0 : 0)
-  0 0 _0.5&p.
-. 0 1 _0.5&p.
<: 0 _1 0.5&p.
>: 0 1 0.5&p.
[ 0 0 0.5&p.
] 0 0 0.5&p.
+: *:
*: 0 0 0 0.3333333333333333333&p.
-: 0 0 0.25&p.
o. 0 0 0.5p1&p.
% ^.
%: %: * 0 0.666666666666666666&p.
^ ^
^. (]*^.) - ]
j. 0 0 0j0.5&p.
r. -@j.@r.
0: 0&*
1: 1&*
2: 2&*
3: 3&*
4: 4&*
5: 5&*
6: 6&*
7: 7&*
8: 8&*
9: 9&*
_1: _1&*
_2: _2&*
_3: _3&*
_4: _4&*
_5: _5&*
_6: _6&*
_7: _7&*
_8: _8&*
_9: _9&*
_: _&*
__: __&*
)


NB. Canned table of integrals for ^&n@:v
NB. The derivative will be interpreted as an adverb, where m is the n from ^&n
'primintpatopvb primintpatopint' =: <"1 |: (({. ; deb@:}.)~ i.&' ');._2 (0 : 0)
^. (] * ^&m@^.) - m&* @(^&(m-1)@^. deriv _1)
1&o. %&(-m )@(^&(m-1)@(1&o.) * 2&o.) + ((m-1)%m)&*@(^&(m-2)@(1&o.) deriv _1)
2&o. %&m    @(^&(m-1)@(2&o.) * 1&o.) + ((m-1)%m)&*@(^&(m-2)@(2&o.) deriv _1)
3&o. %&(m-1)@(^&(m-1)@(3&o.)       ) -              ^&(m-2)@(3&o.) deriv _1 
7&o. %&(1-m)@(^&(m-1)@(7&o.)       ) +              ^&(m-2)@(7&o.) deriv _1 
)

NB. Canned table of integrals for m&v
NB. The integral will be interpreted as an adverb, where m is the m from m&v
NB. That adverb might itself execute another verb, eg for o.
'primintmandvvb primintmandv' =: <"1 |: (({. ; deb@:}.)~ i.&' ');._2 (0 : 0)
% m&*@^.
^ %&(^. m)@(m&^)
o. m circleint
)

NB. Canned table of integrals for u&n
NB. The integral will be interpreted as an adverb, where m is the n from u&n
'primintuandnvb primintuandn' =: <"1 |: (({. ; deb@:}.)~ i.&' ');._2 (0 : 0)
 %&(m+1)@(^&(m+1))`^.@.(m-:_1)
)

NB. Canned table of integrals for reflexive verbs
'intrefvb refint' =: <"1 |: (({. ; deb@:}.)~ i.&' ');._2 (0 : 0)
+ *:
* 0 0 0 0.33333333333333331&p.
- 0:
^. 0:
% ]
)

NB. Integrals of m&o., with verb results
circleint =: 1 : 0
select. m
case. 1 do. -@:(2&o.)
case. 2 do. 1&o.
case. 3 do. -@^.@(2&o.)
case. 5 do. 6&o.
case. 6 do. 5&o.
case. 7 do. ^.@(6&o.)
case. do. 13!:8 (3) NB. if unknown type, fail
end.
)


NB. operations on polynomials
NB. Generally an empty list means 'not a polynomial'
pplus =: +/@,:
pminus =: -/@,:
ptymes =: +//.@:(*/)
pexp =: 4 : 0  NB. y must have only 1 term
if. y ~: <. y do. $0 return. end.
if. y < 0 do. $0 return. end.
if. y = 0 do. ,1 return. end.
x0 =. x
for_lsb. }. #: {. y do. x =. x0 ptymes^:lsb x ptymes x end.
)

NB. y is a string or an AR, result is the polynomial
topoly =: 3 : 0
if. 0 (>: L.) yar =. y do. yar =. y arofstringu end.  NB. AR of y
NB. Top-down recursion through the AR
if. 1 (>: L.) yar do.
  NB. Unboxed contents.  It's a primitive or a name
  NB. Look it up in the derivative table
  if. (#polyprimvb) > primno =. polyprimvb i. yar do. primno {:: primpoly return. end.
  NB. If undifferentiable primitive, fail
  '' return.
end.
NB. Not a primitive or named verb.  Must be a 2-box list.  The first box indicates the (possibly invisible)
NB. modifier; process it if we know it
yar =. > yar  NB. discard outer boxing
select. {.yar
case. <,'0' do.  NB. noun
  if. '' -: $nounarg =. 1{::yar do. ,nounarg return. end.  NB. If an atom, treat as polynomial constant
  '' return.  NB. error otherwise
case. <,'2' do.  NB. hook - not supported
case. <,'3' do.  NB. fork
  'farg garg harg' =. yar&opstr&.> 0 1 2
  if. (<'[:') -: yar opar 0 do. topoly garg atops harg return. end.   NB. [: g h
  f =. topoly farg [ h =. topoly harg
  if. 0 = f *.&# h do. '' return. end.  NB. f and h must be polynomials
  select. garg
  case. <,'+' do. f pplus h return.
  case. <,'-' do. f pminus h return.
  case. <,'*' do. f ptymes h return.
  case. <,'^' do. if. 1=#f do. f pexp h return. end.
  end.

case. <,'&' do.
  if. yar opisnoun 0 do.   NB. m&v  (v must be a verb)
    nounarg =. (yar opar 0) vnofaru  NB. m
    verbarg =. (yar opar 1)  NB. v as an AR
    if. verbarg -: <'p.' do.   NB. Handle p. as a special case
      if. 0 (< L.) nounarg do.   NB.  roots or multinomial
        if. 1 < #nounarg do. p. nounarg return. end.  NB. roots: convert to coeffs and return
        if. (2 = #@$ > nounarg) *. (2 = {>@$ > nounarg) do. p. nounarg return. end.  NB. multinomial convert to coeffs and return
        $0 return.  NB. Bad shape, not polynomial
      end.
      (,nounarg) return.   NB. coeffs, return them
    end.
    if. *@#@$ nounarg do. 13!:8 (3) end. NB. except for p. only atomic m is recognized
    NB. Look it up in the polynomial table, execute it on m, return string result
    if. (#primpmandvvb) > primno =. primpmandvvb i. verbarg do.
      (nounarg 1 : (primno {:: primpmandv))  return.
    end.
    $0 return.
    
  elseif. yar opisnoun 1 do.   NB. u&n  (u must be a verb)
    nounarg =. (yar opar 1) vnofaru  NB. n
    verbarg =. (yar opar 0)  NB. u as an AR
    if. *@#@$ nounarg do. 13!:8 (3) end. NB. except for p. only atomic m is recognized
    NB. Look it up in the polynomial table, execute it on m, return string result
    if. (#primpuandnvb) > primno =. primpuandnvb i. verbarg do.
      (nounarg 1 : (primno {:: primpuandn))  return.
    end.
    NB. ^&n is a polynomial if n is nonnegative integer
    if. (verbarg -: <,'^') *. (nounarg=<.nounarg) *. (nounarg>:0) do. (->:nounarg) {. 1 return. end.
    NB. If undifferentiable primitive, fail
    $0 return.
  end.

  NB. fall through to...
case. ;:'&:@@:' do.
  'uarg varg' =. yar&opstr&.> 0 1
  f =. topoly uarg [ h =. topoly varg
  if. 0 = f *.&# h do. '' return. end.  NB. f and h must be polynomials
  +/ f * h pexp"_ 0 i. #f return.

case. <'^:' do.  NB. power
  'uarg varg' =. yar&opstr&.> 0 1
  f =. topoly uarg [ h =. topoly varg
  if. 0 = f *.&# h do. '' return. end.  NB. f and h must be polynomials
  if. h ~: <. h do. $0 return. end.
  if. h < 0 do. $0 return. end.  NB. h must be positive
  select. h
  case. 1 do. f
  case. 0 do. 0 1
  case. do. topoly (uarg powops <:h) atops uarg  NB. Large powers become u^:(<:n)@u
  end.
  return.

case. <,'"' do.  NB. rank
  'uarg varg' =. yar&opstr&.> 0 1
  f =. topoly uarg [ h =. topoly varg
  if. 0 = f *.&# h do. '' return. end.  NB. f and h must be polynomials
  f return.  NB. If an atom, treat as polynomial constant

case. <,'~' do.  NB. reflexive - these cases are rare & we ignore them
end.
$0   NB. Unknown or unprocessed modifier, fail
)


NB. Canned table of derivatives for the primitive verbs
'polyprimvb primpoly' =: <"1 |: (({. ; ".@}.)~ i.&' ');._2 (0 : 0)
-  0 _1
-. 1 _1
<: _1 1
>: 1 1
[ 0 1
] 0 1
+: 0 2
*: 0 0 1
-: 0 0.5
o. 0 1p1
j. 0 0j1
0: ,0
1: ,1
2: ,2
3: ,3
4: ,4
5: ,5
6: ,6
7: ,7
8: ,8
9: ,9
_1: ,_1
_2: ,_2
_3: ,_3
_4: ,_4
_5: ,_5
_6: ,_6
_7: ,_7
_8: ,_8
_9: ,_9
_: ,_
__: ,__
)

NB. Canned table of polynomials for m&v
NB. The polynomial will be interpreted as an adverb, where m is the m from m&v
NB. That adverb might itself execute another verb, eg for o.
'primpmandvvb primpmandv' =: <"1 |: (({. ; deb@:}.)~ i.&' ');._2 (0 : 0)
+ m,1
* 0,m
- m,_1
)

NB. Canned table of polynomials for u&n
NB. The polynomial will be interpreted as an adverb, where m is the n from u&n
'primpuandnvb primpuandn' =: <"1 |: (({. ; deb@:}.)~ i.&' ');._2 (0 : 0)
+ m,1
* 0,m
- (-m),1
% 0,(%m)
)


NB. Partial Derivative operator, like D.
pderiv =: 2 : 0
assert. '' -: $n   NB. n must be an atom
order =. <. n
assert. n = order
assert. 3 = 4!:0 <'u'   NB. Must be in user's namespace
if. order < 0 do. 13!:8 (3)  end.  NB. integrals not allowed
ufix =. u. f.  NB. u in user's namespace
resstg =. 5!:5 <'ufix'   NB. User's verb in character form
while. order do.   NB. If the order is not 1 or _1, keep increasing the order.
NB. Find the derivative of the AR.  If it fails at any point, revert to secant approximation for the remnant
  try.
    resstg =. pderivstg resstg
  catch.
    break.  NB. keep resstg/order from the last successful calculation
  end.
  order =. <: order   NB. Decrement # orders remaining
end.
NB. resstg has the string form and order is the number of derivs unfinished.  Return appropriately
NB. Parse the verb to create a verb
NB. Return the verb itself if it's the right order, otherwise a call to the approximator
if. order=0 do. resstg vnofu else. 0&(resstg vnofu derivsecant_jcalculus_ order) end.
)

NB. symbolic partial  derivative, in string form
NB. y is string form of a verb
NB. Result is string form of the derivative verb
NB. If we don't know the derivative, we raise an error
NB. We pass around strings rather than ARs for convenience of coding & table maintenance
pderivstg =: 3 : 0
yar =. y arofstringu
NB. First, see if we can take the total derivative
try.
  d =. derivstg y  NB. This fails if there is no derivative
  NB.  if not rnk  0, append the multivariate appendix
  if. 0 ~: {. (d vnofu) b. 0 do. d =. '(* =/~@(i.@$))@(',d,')' end.
  d  return.
catch.
end.
NB. Try the forms that have no variables, on the original string verb
if. (#primpvb) > primno =. primpvb i. yar do. primno {:: primpderiv return. end.
if. 1 (>: L.) yar do. 13!:8 (3) end.   NB. If undifferentiable primitive, fail

NB. Not a primitive or named verb.  Must be a 2-box list.  The first box indicates the (possibly invisible)
NB. modifier; process it if we know it
yar =. > yar  NB. discard outer boxing
select. {.yar
case. <,'0' do.  NB. noun - should not occur
case. <,'2' do.  NB. hook - not supported
case. <,'3' do.  NB. fork
  'farg garg harg' =. yar&opstr&.> 0 1 2
  df =. pderivstg farg [ dh =. pderivstg harg
  select. garg
  case. <,'+' do. df fplus dh return.
  case. <,'-' do. df fminus dh return.
  case. <,'*' do. (df ftymes harg) fplus (farg ftymes dh) return.  NB. product rule
  case. <,'%' do. ((df ftymes harg) fminus (farg ftymes dh)) fdiv '*:' atops harg return.  NB. quotient rule
  end.

case. <,'&' do.  NB. &
  NB. The only thing we handle here are structural permutations where we can create a boolean
  NB. matrix indicating what is connected to what
  if. yar opisnoun 0 do.   NB. m&v  (v must be a verb)
    nounarg =. (yar opar 0) vnofaru  NB. m
    verbarg =. (yar opar 1)  NB. v as an AR
    if. verbarg e. ;:'|. |: { A. C.' do.   NB. perm type
      '(=/ ',y,')@(i.@$)' return.
    end.
    13!:8 (3)
  end.
  NB. fallthrough cases are u&v

fcase. ;:'&:@@:' do.  NB. chain rule, including fallthrough from &
  if. yar opisnoun 1 do. 13!:8 (3) end.  NB. v must be a verb
  uop =. yar opstr 0   NB. the verb as a string
  vop =. yar opstr 1   NB. the verb as a verb
  (pderivstg vop) fmp (pderivstg uop) atops vop return.

case. <,'"' do.  NB. rank
  if. -. yar opisnoun 1 do. 13!:8 (3) end.  NB. v must be a noun
  nop =. (yar opar 1) vnofaru   NB. extract the noun
  uop =. (yar opstr 0) vnofaru  NB. extract the verb u as a string
  if. -. yar opisnoun 0 do. (pderivstg uop),'"',":nop  end.  NB. u is a verb
  '$&0@0"',":nop return.

case. <,'~' do.  NB. reflexive - these cases are rare & we ignore them
NB. +/ omitted
end.
13!:8 (3)   NB. Unknown or unprocessed modifier, fail
)

NB. Canned table of derivatives for the primitive verbs
'primpvb primpderiv' =: <"1 |: (({. ; deb@:}.)~ i.&' ');._2 (0 : 0)
|. (|.=/])@(i.@$)
|: (|:=/])@(i.@$)
+/ ({. =/ */@}.@$ | ])@(i.@$)

)

NB. Secant slope, like D:
sslope =: 2 : 0
assert. _1 = 4!:0 <'m'  NB. u is verb
assert. 0 = 4!:0 <'n'  NB. v is a noun
assert. 0 = #$n  NB. n is an atom
assert. n > 0   NB. n > 0
(u. derivsecant_jcalculus_ n)"(u. f. b. 0)
)

0 : 0  NB. Sample sentences
*: deriv 1
*: deriv 2
*: pderiv 1
*:"_ pderiv 1
0.0005 +/@:(1 2&*) sslope 1 i. 2 3
)