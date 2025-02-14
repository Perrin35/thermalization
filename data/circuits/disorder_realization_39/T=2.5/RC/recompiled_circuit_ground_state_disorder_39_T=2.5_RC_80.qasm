OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.429739) q[0];
sx q[0];
rz(-1.9965594) q[0];
sx q[0];
rz(2.1249007) q[0];
rz(-3.9859803) q[1];
sx q[1];
rz(3.9581668) q[1];
sx q[1];
rz(8.858455) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96903548) q[0];
sx q[0];
rz(-0.96562562) q[0];
sx q[0];
rz(3.023554) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2015622) q[2];
sx q[2];
rz(-0.21182952) q[2];
sx q[2];
rz(0.43583696) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.8992246) q[1];
sx q[1];
rz(-1.4903233) q[1];
sx q[1];
rz(2.7182275) q[1];
rz(-pi) q[2];
rz(1.3641961) q[3];
sx q[3];
rz(-2.1461764) q[3];
sx q[3];
rz(2.8386338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.5113968) q[2];
sx q[2];
rz(-0.87367311) q[2];
sx q[2];
rz(1.1279747) q[2];
rz(-1.5026211) q[3];
sx q[3];
rz(-2.2962544) q[3];
sx q[3];
rz(-0.42475167) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86654919) q[0];
sx q[0];
rz(-2.6847222) q[0];
sx q[0];
rz(0.04091111) q[0];
rz(-0.54967898) q[1];
sx q[1];
rz(-2.7626541) q[1];
sx q[1];
rz(3.0612225) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5588829) q[0];
sx q[0];
rz(-2.2930995) q[0];
sx q[0];
rz(1.5764019) q[0];
rz(0.92795579) q[2];
sx q[2];
rz(-1.11737) q[2];
sx q[2];
rz(-2.3969899) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.077420635) q[1];
sx q[1];
rz(-1.6257833) q[1];
sx q[1];
rz(1.085678) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8714058) q[3];
sx q[3];
rz(-0.98434292) q[3];
sx q[3];
rz(1.0932066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.4139159) q[2];
sx q[2];
rz(-2.2985986) q[2];
sx q[2];
rz(0.28398871) q[2];
rz(1.5208288) q[3];
sx q[3];
rz(-1.5512543) q[3];
sx q[3];
rz(2.3782597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3182227) q[0];
sx q[0];
rz(-2.9871873) q[0];
sx q[0];
rz(0.41686091) q[0];
rz(2.2082224) q[1];
sx q[1];
rz(-2.2552172) q[1];
sx q[1];
rz(2.0534168) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6898854) q[0];
sx q[0];
rz(-2.5390115) q[0];
sx q[0];
rz(-1.129928) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6112664) q[2];
sx q[2];
rz(-1.010561) q[2];
sx q[2];
rz(0.76228415) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.5610366) q[1];
sx q[1];
rz(-1.2966178) q[1];
sx q[1];
rz(-0.31117532) q[1];
x q[2];
rz(1.635664) q[3];
sx q[3];
rz(-1.8577575) q[3];
sx q[3];
rz(0.6395517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.1209968) q[2];
sx q[2];
rz(-0.72046295) q[2];
sx q[2];
rz(0.32568112) q[2];
rz(2.7459512) q[3];
sx q[3];
rz(-0.49042693) q[3];
sx q[3];
rz(-2.4237848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2242103) q[0];
sx q[0];
rz(-2.2069187) q[0];
sx q[0];
rz(0.14064661) q[0];
rz(-2.7463101) q[1];
sx q[1];
rz(-1.544516) q[1];
sx q[1];
rz(-0.43221727) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0790063) q[0];
sx q[0];
rz(-0.64029988) q[0];
sx q[0];
rz(-2.2878134) q[0];
rz(-pi) q[1];
rz(-1.725196) q[2];
sx q[2];
rz(-2.6302166) q[2];
sx q[2];
rz(-1.7768154) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.4443214) q[1];
sx q[1];
rz(-1.2978076) q[1];
sx q[1];
rz(-1.4876025) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1915901) q[3];
sx q[3];
rz(-2.0696313) q[3];
sx q[3];
rz(2.102559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.29493368) q[2];
sx q[2];
rz(-2.0294971) q[2];
sx q[2];
rz(-2.6726932) q[2];
rz(2.9851959) q[3];
sx q[3];
rz(-0.96894914) q[3];
sx q[3];
rz(1.8739353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0676607) q[0];
sx q[0];
rz(-1.1336552) q[0];
sx q[0];
rz(2.3611948) q[0];
rz(-2.228915) q[1];
sx q[1];
rz(-2.3527805) q[1];
sx q[1];
rz(1.3891634) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3816504) q[0];
sx q[0];
rz(-1.3033322) q[0];
sx q[0];
rz(-0.36422361) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.19540968) q[2];
sx q[2];
rz(-0.59048803) q[2];
sx q[2];
rz(-0.54412847) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.011574419) q[1];
sx q[1];
rz(-0.60638529) q[1];
sx q[1];
rz(-2.9467086) q[1];
x q[2];
rz(1.0205998) q[3];
sx q[3];
rz(-2.3781099) q[3];
sx q[3];
rz(-0.17624763) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.7954365) q[2];
sx q[2];
rz(-0.6554335) q[2];
sx q[2];
rz(-0.1795086) q[2];
rz(1.4795715) q[3];
sx q[3];
rz(-0.69412762) q[3];
sx q[3];
rz(3.1365385) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96106225) q[0];
sx q[0];
rz(-1.5164277) q[0];
sx q[0];
rz(-1.4638715) q[0];
rz(-3.1278817) q[1];
sx q[1];
rz(-1.6746215) q[1];
sx q[1];
rz(-2.1965068) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7031539) q[0];
sx q[0];
rz(-0.8301211) q[0];
sx q[0];
rz(1.0891455) q[0];
rz(1.361998) q[2];
sx q[2];
rz(-1.8755302) q[2];
sx q[2];
rz(1.3688251) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.2234563) q[1];
sx q[1];
rz(-2.0386749) q[1];
sx q[1];
rz(1.0997869) q[1];
x q[2];
rz(2.1626445) q[3];
sx q[3];
rz(-2.0267539) q[3];
sx q[3];
rz(0.78612529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.98662394) q[2];
sx q[2];
rz(-1.0633609) q[2];
sx q[2];
rz(-1.0490136) q[2];
rz(1.5341885) q[3];
sx q[3];
rz(-2.1395855) q[3];
sx q[3];
rz(1.775942) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5904215) q[0];
sx q[0];
rz(-0.48749247) q[0];
sx q[0];
rz(-2.3792939) q[0];
rz(-0.47362348) q[1];
sx q[1];
rz(-1.7601687) q[1];
sx q[1];
rz(-0.90829888) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1518734) q[0];
sx q[0];
rz(-2.3494968) q[0];
sx q[0];
rz(-0.0072974722) q[0];
rz(-pi) q[1];
rz(2.9680266) q[2];
sx q[2];
rz(-1.3527292) q[2];
sx q[2];
rz(1.1898927) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6544853) q[1];
sx q[1];
rz(-1.5853549) q[1];
sx q[1];
rz(1.3338939) q[1];
x q[2];
rz(1.1486111) q[3];
sx q[3];
rz(-3.1109538) q[3];
sx q[3];
rz(-0.8493166) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.121754) q[2];
sx q[2];
rz(-2.8826931) q[2];
sx q[2];
rz(-2.1799977) q[2];
rz(-0.52729765) q[3];
sx q[3];
rz(-1.7617825) q[3];
sx q[3];
rz(2.5808047) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2335662) q[0];
sx q[0];
rz(-2.5416424) q[0];
sx q[0];
rz(2.7954234) q[0];
rz(1.2973805) q[1];
sx q[1];
rz(-0.66645122) q[1];
sx q[1];
rz(-0.60874879) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93701836) q[0];
sx q[0];
rz(-1.6547583) q[0];
sx q[0];
rz(2.2473826) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4791711) q[2];
sx q[2];
rz(-2.8067744) q[2];
sx q[2];
rz(-1.3864558) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.1947555) q[1];
sx q[1];
rz(-1.5018592) q[1];
sx q[1];
rz(0.26534715) q[1];
x q[2];
rz(-2.5777588) q[3];
sx q[3];
rz(-2.2632416) q[3];
sx q[3];
rz(3.0094224) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.79494563) q[2];
sx q[2];
rz(-0.55507675) q[2];
sx q[2];
rz(-0.65183276) q[2];
rz(0.73976222) q[3];
sx q[3];
rz(-1.0659822) q[3];
sx q[3];
rz(-0.8968001) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5480492) q[0];
sx q[0];
rz(-0.47320047) q[0];
sx q[0];
rz(2.5439673) q[0];
rz(-1.301544) q[1];
sx q[1];
rz(-1.5664682) q[1];
sx q[1];
rz(0.32599932) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7571288) q[0];
sx q[0];
rz(-1.8600896) q[0];
sx q[0];
rz(-2.3209653) q[0];
rz(-pi) q[1];
rz(2.4818899) q[2];
sx q[2];
rz(-2.7355425) q[2];
sx q[2];
rz(-2.7999807) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1797267) q[1];
sx q[1];
rz(-0.71007198) q[1];
sx q[1];
rz(-1.1674787) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3365485) q[3];
sx q[3];
rz(-2.9784103) q[3];
sx q[3];
rz(1.9022579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.7035383) q[2];
sx q[2];
rz(-1.2039098) q[2];
sx q[2];
rz(-0.34029141) q[2];
rz(1.641364) q[3];
sx q[3];
rz(-0.54268018) q[3];
sx q[3];
rz(-0.59804183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21094766) q[0];
sx q[0];
rz(-1.1177381) q[0];
sx q[0];
rz(-3.0979544) q[0];
rz(-2.4234407) q[1];
sx q[1];
rz(-1.3190045) q[1];
sx q[1];
rz(-1.7125548) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4627422) q[0];
sx q[0];
rz(-1.382811) q[0];
sx q[0];
rz(-1.7440656) q[0];
rz(-pi) q[1];
rz(-1.4347836) q[2];
sx q[2];
rz(-2.2439661) q[2];
sx q[2];
rz(1.7409131) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.82270993) q[1];
sx q[1];
rz(-1.4850933) q[1];
sx q[1];
rz(-2.6515168) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9472856) q[3];
sx q[3];
rz(-1.2293136) q[3];
sx q[3];
rz(-1.7870513) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.67937294) q[2];
sx q[2];
rz(-0.93928176) q[2];
sx q[2];
rz(-0.77793724) q[2];
rz(2.098162) q[3];
sx q[3];
rz(-1.611004) q[3];
sx q[3];
rz(1.9061609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0935681) q[0];
sx q[0];
rz(-2.1401736) q[0];
sx q[0];
rz(-2.370136) q[0];
rz(-1.7780766) q[1];
sx q[1];
rz(-1.9693146) q[1];
sx q[1];
rz(1.6065425) q[1];
rz(-2.4189998) q[2];
sx q[2];
rz(-1.9836139) q[2];
sx q[2];
rz(0.2000533) q[2];
rz(-1.0092953) q[3];
sx q[3];
rz(-2.3875368) q[3];
sx q[3];
rz(-0.33215678) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
