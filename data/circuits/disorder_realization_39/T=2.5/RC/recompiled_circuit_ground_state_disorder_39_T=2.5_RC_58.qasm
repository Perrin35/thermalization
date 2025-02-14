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
rz(-0.84438762) q[1];
sx q[1];
rz(-0.8165741) q[1];
sx q[1];
rz(-2.5752697) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53439683) q[0];
sx q[0];
rz(-1.473794) q[0];
sx q[0];
rz(-2.1792381) q[0];
rz(-pi) q[1];
rz(1.9400305) q[2];
sx q[2];
rz(-0.21182952) q[2];
sx q[2];
rz(-2.7057557) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.8493718) q[1];
sx q[1];
rz(-1.9927036) q[1];
sx q[1];
rz(1.4825691) q[1];
rz(-pi) q[2];
rz(2.5563566) q[3];
sx q[3];
rz(-1.7437616) q[3];
sx q[3];
rz(-1.1542785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.5113968) q[2];
sx q[2];
rz(-0.87367311) q[2];
sx q[2];
rz(2.013618) q[2];
rz(1.5026211) q[3];
sx q[3];
rz(-0.84533826) q[3];
sx q[3];
rz(2.716841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2750435) q[0];
sx q[0];
rz(-0.4568704) q[0];
sx q[0];
rz(-3.1006815) q[0];
rz(0.54967898) q[1];
sx q[1];
rz(-2.7626541) q[1];
sx q[1];
rz(0.080370195) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59118862) q[0];
sx q[0];
rz(-0.72232095) q[0];
sx q[0];
rz(-0.0063615464) q[0];
rz(2.2136369) q[2];
sx q[2];
rz(-2.0242226) q[2];
sx q[2];
rz(0.7446028) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7521022) q[1];
sx q[1];
rz(-0.48797777) q[1];
sx q[1];
rz(1.4533024) q[1];
rz(2.1744022) q[3];
sx q[3];
rz(-1.7949824) q[3];
sx q[3];
rz(0.32550117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.72767672) q[2];
sx q[2];
rz(-2.2985986) q[2];
sx q[2];
rz(0.28398871) q[2];
rz(-1.5208288) q[3];
sx q[3];
rz(-1.5903383) q[3];
sx q[3];
rz(2.3782597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82336998) q[0];
sx q[0];
rz(-0.1544054) q[0];
sx q[0];
rz(0.41686091) q[0];
rz(-0.93337026) q[1];
sx q[1];
rz(-2.2552172) q[1];
sx q[1];
rz(-1.0881759) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4517072) q[0];
sx q[0];
rz(-2.5390115) q[0];
sx q[0];
rz(2.0116647) q[0];
x q[1];
rz(3.0771824) q[2];
sx q[2];
rz(-2.5800534) q[2];
sx q[2];
rz(-0.68622996) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0644898) q[1];
sx q[1];
rz(-1.2716115) q[1];
sx q[1];
rz(1.283518) q[1];
rz(-pi) q[2];
rz(1.635664) q[3];
sx q[3];
rz(-1.8577575) q[3];
sx q[3];
rz(-2.502041) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0205959) q[2];
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
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91738236) q[0];
sx q[0];
rz(-0.93467394) q[0];
sx q[0];
rz(-3.000946) q[0];
rz(-0.39528254) q[1];
sx q[1];
rz(-1.5970767) q[1];
sx q[1];
rz(-0.43221727) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2435139) q[0];
sx q[0];
rz(-1.9742516) q[0];
sx q[0];
rz(2.0824758) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4163966) q[2];
sx q[2];
rz(-2.6302166) q[2];
sx q[2];
rz(1.7768154) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.6972713) q[1];
sx q[1];
rz(-1.2978076) q[1];
sx q[1];
rz(-1.4876025) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.95000259) q[3];
sx q[3];
rz(-2.0696313) q[3];
sx q[3];
rz(1.0390337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.846659) q[2];
sx q[2];
rz(-1.1120956) q[2];
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
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.073931996) q[0];
sx q[0];
rz(-1.1336552) q[0];
sx q[0];
rz(2.3611948) q[0];
rz(-0.91267768) q[1];
sx q[1];
rz(-0.78881216) q[1];
sx q[1];
rz(1.3891634) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7245623) q[0];
sx q[0];
rz(-0.44831845) q[0];
sx q[0];
rz(2.4858767) q[0];
x q[1];
rz(-2.5599586) q[2];
sx q[2];
rz(-1.6791145) q[2];
sx q[2];
rz(-1.9519639) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.9174172) q[1];
sx q[1];
rz(-0.97748102) q[1];
sx q[1];
rz(1.437287) q[1];
rz(-2.2551304) q[3];
sx q[3];
rz(-1.9406978) q[3];
sx q[3];
rz(-2.1641638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.3461561) q[2];
sx q[2];
rz(-0.6554335) q[2];
sx q[2];
rz(-2.9620841) q[2];
rz(1.4795715) q[3];
sx q[3];
rz(-0.69412762) q[3];
sx q[3];
rz(-0.0050541335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96106225) q[0];
sx q[0];
rz(-1.625165) q[0];
sx q[0];
rz(1.4638715) q[0];
rz(-0.013710984) q[1];
sx q[1];
rz(-1.4669712) q[1];
sx q[1];
rz(0.94508583) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7031539) q[0];
sx q[0];
rz(-2.3114716) q[0];
sx q[0];
rz(-2.0524471) q[0];
x q[1];
rz(-2.5588972) q[2];
sx q[2];
rz(-0.36755633) q[2];
sx q[2];
rz(-0.75395298) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.764115) q[1];
sx q[1];
rz(-2.4903653) q[1];
sx q[1];
rz(-2.4098544) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2919214) q[3];
sx q[3];
rz(-2.4115218) q[3];
sx q[3];
rz(1.7770313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.98662394) q[2];
sx q[2];
rz(-2.0782317) q[2];
sx q[2];
rz(-1.0490136) q[2];
rz(1.5341885) q[3];
sx q[3];
rz(-1.0020071) q[3];
sx q[3];
rz(-1.775942) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55117115) q[0];
sx q[0];
rz(-0.48749247) q[0];
sx q[0];
rz(0.76229873) q[0];
rz(2.6679692) q[1];
sx q[1];
rz(-1.7601687) q[1];
sx q[1];
rz(-0.90829888) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.98971924) q[0];
sx q[0];
rz(-2.3494968) q[0];
sx q[0];
rz(3.1342952) q[0];
x q[1];
rz(-1.7920831) q[2];
sx q[2];
rz(-1.7402116) q[2];
sx q[2];
rz(2.7986023) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.0543889) q[1];
sx q[1];
rz(-1.3339195) q[1];
sx q[1];
rz(0.014976784) q[1];
rz(-1.5987464) q[3];
sx q[3];
rz(-1.5582435) q[3];
sx q[3];
rz(-0.2994699) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.0198387) q[2];
sx q[2];
rz(-2.8826931) q[2];
sx q[2];
rz(-2.1799977) q[2];
rz(-2.614295) q[3];
sx q[3];
rz(-1.7617825) q[3];
sx q[3];
rz(-2.5808047) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2335662) q[0];
sx q[0];
rz(-0.59995025) q[0];
sx q[0];
rz(-0.34616923) q[0];
rz(1.2973805) q[1];
sx q[1];
rz(-2.4751414) q[1];
sx q[1];
rz(0.60874879) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70102126) q[0];
sx q[0];
rz(-0.89703575) q[0];
sx q[0];
rz(-0.10752039) q[0];
rz(1.7815963) q[2];
sx q[2];
rz(-1.3087052) q[2];
sx q[2];
rz(2.0767625) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.94683719) q[1];
sx q[1];
rz(-1.6397335) q[1];
sx q[1];
rz(0.26534715) q[1];
x q[2];
rz(2.5777588) q[3];
sx q[3];
rz(-2.2632416) q[3];
sx q[3];
rz(0.13217029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.346647) q[2];
sx q[2];
rz(-2.5865159) q[2];
sx q[2];
rz(-2.4897599) q[2];
rz(-0.73976222) q[3];
sx q[3];
rz(-2.0756105) q[3];
sx q[3];
rz(-0.8968001) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5480492) q[0];
sx q[0];
rz(-2.6683922) q[0];
sx q[0];
rz(-2.5439673) q[0];
rz(-1.301544) q[1];
sx q[1];
rz(-1.5751244) q[1];
sx q[1];
rz(-0.32599932) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38446389) q[0];
sx q[0];
rz(-1.281503) q[0];
sx q[0];
rz(-0.82062736) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4818899) q[2];
sx q[2];
rz(-0.40605011) q[2];
sx q[2];
rz(0.34161196) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.9218775) q[1];
sx q[1];
rz(-1.829521) q[1];
sx q[1];
rz(-2.2398276) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6889309) q[3];
sx q[3];
rz(-1.6836327) q[3];
sx q[3];
rz(-1.090534) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7035383) q[2];
sx q[2];
rz(-1.9376829) q[2];
sx q[2];
rz(2.8013012) q[2];
rz(1.5002286) q[3];
sx q[3];
rz(-0.54268018) q[3];
sx q[3];
rz(-2.5435508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21094766) q[0];
sx q[0];
rz(-2.0238545) q[0];
sx q[0];
rz(0.043638226) q[0];
rz(2.4234407) q[1];
sx q[1];
rz(-1.3190045) q[1];
sx q[1];
rz(1.7125548) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71009174) q[0];
sx q[0];
rz(-0.25496182) q[0];
sx q[0];
rz(2.405317) q[0];
rz(-pi) q[1];
rz(2.9731644) q[2];
sx q[2];
rz(-2.456924) q[2];
sx q[2];
rz(1.9569966) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.43914) q[1];
sx q[1];
rz(-1.0826775) q[1];
sx q[1];
rz(1.6678651) q[1];
rz(1.9472856) q[3];
sx q[3];
rz(-1.2293136) q[3];
sx q[3];
rz(1.3545413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.67937294) q[2];
sx q[2];
rz(-0.93928176) q[2];
sx q[2];
rz(-2.3636554) q[2];
rz(2.098162) q[3];
sx q[3];
rz(-1.611004) q[3];
sx q[3];
rz(-1.2354318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.048024561) q[0];
sx q[0];
rz(-1.0014191) q[0];
sx q[0];
rz(0.77145664) q[0];
rz(1.363516) q[1];
sx q[1];
rz(-1.9693146) q[1];
sx q[1];
rz(1.6065425) q[1];
rz(2.099299) q[2];
sx q[2];
rz(-0.92005554) q[2];
sx q[2];
rz(-1.0307606) q[2];
rz(-2.1322973) q[3];
sx q[3];
rz(-0.75405585) q[3];
sx q[3];
rz(2.8094359) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
