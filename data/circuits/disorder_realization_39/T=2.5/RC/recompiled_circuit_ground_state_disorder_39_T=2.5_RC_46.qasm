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
rz(-1.0166919) q[0];
rz(-0.84438762) q[1];
sx q[1];
rz(-0.8165741) q[1];
sx q[1];
rz(-2.5752697) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96903548) q[0];
sx q[0];
rz(-0.96562562) q[0];
sx q[0];
rz(-3.023554) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3728605) q[2];
sx q[2];
rz(-1.6467484) q[2];
sx q[2];
rz(1.4966485) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.8493718) q[1];
sx q[1];
rz(-1.148889) q[1];
sx q[1];
rz(-1.4825691) q[1];
x q[2];
rz(-1.7773966) q[3];
sx q[3];
rz(-2.1461764) q[3];
sx q[3];
rz(2.8386338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.6301959) q[2];
sx q[2];
rz(-2.2679195) q[2];
sx q[2];
rz(-1.1279747) q[2];
rz(1.6389716) q[3];
sx q[3];
rz(-0.84533826) q[3];
sx q[3];
rz(0.42475167) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86654919) q[0];
sx q[0];
rz(-2.6847222) q[0];
sx q[0];
rz(-0.04091111) q[0];
rz(0.54967898) q[1];
sx q[1];
rz(-0.37893852) q[1];
sx q[1];
rz(3.0612225) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98438063) q[0];
sx q[0];
rz(-1.5750021) q[0];
sx q[0];
rz(-2.4192817) q[0];
rz(-0.54687107) q[2];
sx q[2];
rz(-1.0016707) q[2];
sx q[2];
rz(-2.6324181) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3894905) q[1];
sx q[1];
rz(-2.6536149) q[1];
sx q[1];
rz(1.6882903) q[1];
rz(-pi) q[2];
rz(2.8714058) q[3];
sx q[3];
rz(-0.98434292) q[3];
sx q[3];
rz(-2.0483861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.4139159) q[2];
sx q[2];
rz(-2.2985986) q[2];
sx q[2];
rz(-2.8576039) q[2];
rz(1.5208288) q[3];
sx q[3];
rz(-1.5903383) q[3];
sx q[3];
rz(0.76333299) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82336998) q[0];
sx q[0];
rz(-0.1544054) q[0];
sx q[0];
rz(-2.7247317) q[0];
rz(2.2082224) q[1];
sx q[1];
rz(-0.88637543) q[1];
sx q[1];
rz(1.0881759) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8899208) q[0];
sx q[0];
rz(-1.8150738) q[0];
sx q[0];
rz(1.0142465) q[0];
rz(-pi) q[1];
rz(0.56060411) q[2];
sx q[2];
rz(-1.5365155) q[2];
sx q[2];
rz(-2.3115668) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.8510906) q[1];
sx q[1];
rz(-2.7298285) q[1];
sx q[1];
rz(2.3985833) q[1];
rz(-0.21622323) q[3];
sx q[3];
rz(-0.29400405) q[3];
sx q[3];
rz(2.2764429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0205959) q[2];
sx q[2];
rz(-2.4211297) q[2];
sx q[2];
rz(-2.8159115) q[2];
rz(0.39564141) q[3];
sx q[3];
rz(-2.6511657) q[3];
sx q[3];
rz(-2.4237848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91738236) q[0];
sx q[0];
rz(-0.93467394) q[0];
sx q[0];
rz(3.000946) q[0];
rz(-2.7463101) q[1];
sx q[1];
rz(-1.544516) q[1];
sx q[1];
rz(-0.43221727) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88971602) q[0];
sx q[0];
rz(-1.1036627) q[0];
sx q[0];
rz(2.6863195) q[0];
x q[1];
rz(1.4163966) q[2];
sx q[2];
rz(-0.51137608) q[2];
sx q[2];
rz(1.7768154) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.85104698) q[1];
sx q[1];
rz(-1.49069) q[1];
sx q[1];
rz(2.8677031) q[1];
rz(-pi) q[2];
rz(2.3234576) q[3];
sx q[3];
rz(-2.3664118) q[3];
sx q[3];
rz(-1.1213746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.846659) q[2];
sx q[2];
rz(-2.0294971) q[2];
sx q[2];
rz(-0.46889949) q[2];
rz(-2.9851959) q[3];
sx q[3];
rz(-2.1726435) q[3];
sx q[3];
rz(1.8739353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.073931996) q[0];
sx q[0];
rz(-1.1336552) q[0];
sx q[0];
rz(0.78039783) q[0];
rz(-0.91267768) q[1];
sx q[1];
rz(-2.3527805) q[1];
sx q[1];
rz(1.7524293) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4311543) q[0];
sx q[0];
rz(-1.9214993) q[0];
sx q[0];
rz(-1.8560657) q[0];
rz(-pi) q[1];
rz(-1.4413799) q[2];
sx q[2];
rz(-2.1485818) q[2];
sx q[2];
rz(-2.8313864) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.720019) q[1];
sx q[1];
rz(-1.6813845) q[1];
sx q[1];
rz(-2.5441267) q[1];
x q[2];
rz(-0.88646226) q[3];
sx q[3];
rz(-1.9406978) q[3];
sx q[3];
rz(-0.97742888) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.3461561) q[2];
sx q[2];
rz(-2.4861591) q[2];
sx q[2];
rz(0.1795086) q[2];
rz(1.6620212) q[3];
sx q[3];
rz(-0.69412762) q[3];
sx q[3];
rz(0.0050541335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1805304) q[0];
sx q[0];
rz(-1.5164277) q[0];
sx q[0];
rz(-1.4638715) q[0];
rz(-0.013710984) q[1];
sx q[1];
rz(-1.4669712) q[1];
sx q[1];
rz(-2.1965068) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4714519) q[0];
sx q[0];
rz(-1.9197122) q[0];
sx q[0];
rz(-0.80100153) q[0];
rz(0.58269549) q[2];
sx q[2];
rz(-2.7740363) q[2];
sx q[2];
rz(0.75395298) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.764115) q[1];
sx q[1];
rz(-0.65122737) q[1];
sx q[1];
rz(2.4098544) q[1];
x q[2];
rz(2.2919214) q[3];
sx q[3];
rz(-0.73007089) q[3];
sx q[3];
rz(1.3645614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.1549687) q[2];
sx q[2];
rz(-1.0633609) q[2];
sx q[2];
rz(-2.092579) q[2];
rz(1.6074041) q[3];
sx q[3];
rz(-2.1395855) q[3];
sx q[3];
rz(-1.775942) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5904215) q[0];
sx q[0];
rz(-2.6541002) q[0];
sx q[0];
rz(0.76229873) q[0];
rz(0.47362348) q[1];
sx q[1];
rz(-1.381424) q[1];
sx q[1];
rz(-0.90829888) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57595162) q[0];
sx q[0];
rz(-1.5759908) q[0];
sx q[0];
rz(-0.7920825) q[0];
x q[1];
rz(-2.9680266) q[2];
sx q[2];
rz(-1.3527292) q[2];
sx q[2];
rz(1.9517) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.087203793) q[1];
sx q[1];
rz(-1.3339195) q[1];
sx q[1];
rz(-0.014976784) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5987464) q[3];
sx q[3];
rz(-1.5833491) q[3];
sx q[3];
rz(-0.2994699) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.0198387) q[2];
sx q[2];
rz(-2.8826931) q[2];
sx q[2];
rz(-2.1799977) q[2];
rz(0.52729765) q[3];
sx q[3];
rz(-1.7617825) q[3];
sx q[3];
rz(-2.5808047) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2335662) q[0];
sx q[0];
rz(-2.5416424) q[0];
sx q[0];
rz(-2.7954234) q[0];
rz(1.2973805) q[1];
sx q[1];
rz(-2.4751414) q[1];
sx q[1];
rz(0.60874879) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52972165) q[0];
sx q[0];
rz(-2.4606315) q[0];
sx q[0];
rz(1.7044071) q[0];
rz(-0.66242156) q[2];
sx q[2];
rz(-2.8067744) q[2];
sx q[2];
rz(-1.3864558) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.37571463) q[1];
sx q[1];
rz(-0.27395136) q[1];
sx q[1];
rz(0.25744827) q[1];
x q[2];
rz(0.56383384) q[3];
sx q[3];
rz(-0.878351) q[3];
sx q[3];
rz(0.13217029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.346647) q[2];
sx q[2];
rz(-0.55507675) q[2];
sx q[2];
rz(0.65183276) q[2];
rz(2.4018304) q[3];
sx q[3];
rz(-2.0756105) q[3];
sx q[3];
rz(-0.8968001) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5480492) q[0];
sx q[0];
rz(-0.47320047) q[0];
sx q[0];
rz(-2.5439673) q[0];
rz(1.301544) q[1];
sx q[1];
rz(-1.5664682) q[1];
sx q[1];
rz(2.8155933) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38446389) q[0];
sx q[0];
rz(-1.8600896) q[0];
sx q[0];
rz(-2.3209653) q[0];
x q[1];
rz(2.8140961) q[2];
sx q[2];
rz(-1.8153037) q[2];
sx q[2];
rz(-1.8482908) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.1797267) q[1];
sx q[1];
rz(-0.71007198) q[1];
sx q[1];
rz(1.9741139) q[1];
x q[2];
rz(-0.80504412) q[3];
sx q[3];
rz(-0.16318233) q[3];
sx q[3];
rz(-1.9022579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.43805435) q[2];
sx q[2];
rz(-1.2039098) q[2];
sx q[2];
rz(0.34029141) q[2];
rz(-1.641364) q[3];
sx q[3];
rz(-2.5989125) q[3];
sx q[3];
rz(2.5435508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21094766) q[0];
sx q[0];
rz(-1.1177381) q[0];
sx q[0];
rz(3.0979544) q[0];
rz(2.4234407) q[1];
sx q[1];
rz(-1.8225881) q[1];
sx q[1];
rz(-1.7125548) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71009174) q[0];
sx q[0];
rz(-2.8866308) q[0];
sx q[0];
rz(-2.405317) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4638955) q[2];
sx q[2];
rz(-1.6770098) q[2];
sx q[2];
rz(-2.8863557) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.70245269) q[1];
sx q[1];
rz(-1.0826775) q[1];
sx q[1];
rz(1.6678651) q[1];
rz(-pi) q[2];
rz(-2.7765482) q[3];
sx q[3];
rz(-1.2170346) q[3];
sx q[3];
rz(0.34788528) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.67937294) q[2];
sx q[2];
rz(-0.93928176) q[2];
sx q[2];
rz(-0.77793724) q[2];
rz(-2.098162) q[3];
sx q[3];
rz(-1.5305887) q[3];
sx q[3];
rz(1.9061609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.048024561) q[0];
sx q[0];
rz(-2.1401736) q[0];
sx q[0];
rz(-2.370136) q[0];
rz(1.363516) q[1];
sx q[1];
rz(-1.9693146) q[1];
sx q[1];
rz(1.6065425) q[1];
rz(0.58495782) q[2];
sx q[2];
rz(-2.3282606) q[2];
sx q[2];
rz(-1.797779) q[2];
rz(2.2424768) q[3];
sx q[3];
rz(-1.9439144) q[3];
sx q[3];
rz(1.6685566) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
