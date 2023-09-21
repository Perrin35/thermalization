OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.70513201) q[0];
sx q[0];
rz(-2.5897265) q[0];
sx q[0];
rz(3.119757) q[0];
rz(2.7472189) q[1];
sx q[1];
rz(-1.4596649) q[1];
sx q[1];
rz(-0.2149166) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2673187) q[0];
sx q[0];
rz(-1.312717) q[0];
sx q[0];
rz(-1.8549071) q[0];
rz(-0.46918842) q[2];
sx q[2];
rz(-1.4326296) q[2];
sx q[2];
rz(-2.5772622) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.79214) q[1];
sx q[1];
rz(-2.1556902) q[1];
sx q[1];
rz(-2.5930415) q[1];
rz(-pi) q[2];
rz(-1.0825726) q[3];
sx q[3];
rz(-2.2256652) q[3];
sx q[3];
rz(2.0522576) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.73137838) q[2];
sx q[2];
rz(-1.6822858) q[2];
sx q[2];
rz(2.5773876) q[2];
rz(-1.7764067) q[3];
sx q[3];
rz(-2.6919638) q[3];
sx q[3];
rz(1.8723429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0441701) q[0];
sx q[0];
rz(-1.928227) q[0];
sx q[0];
rz(2.2136097) q[0];
rz(1.9762951) q[1];
sx q[1];
rz(-1.6033283) q[1];
sx q[1];
rz(-2.2448418) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1363163) q[0];
sx q[0];
rz(-1.3747842) q[0];
sx q[0];
rz(-2.2851903) q[0];
rz(-2.9727544) q[2];
sx q[2];
rz(-1.4988006) q[2];
sx q[2];
rz(0.26693401) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.6078728) q[1];
sx q[1];
rz(-2.3405511) q[1];
sx q[1];
rz(-0.17287066) q[1];
rz(-pi) q[2];
rz(-1.1284626) q[3];
sx q[3];
rz(-2.3395174) q[3];
sx q[3];
rz(-2.341552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8759878) q[2];
sx q[2];
rz(-0.53498712) q[2];
sx q[2];
rz(-1.0401475) q[2];
rz(-1.4552207) q[3];
sx q[3];
rz(-1.2929595) q[3];
sx q[3];
rz(0.35475981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74137694) q[0];
sx q[0];
rz(-2.564036) q[0];
sx q[0];
rz(2.1133912) q[0];
rz(1.0785412) q[1];
sx q[1];
rz(-2.5787347) q[1];
sx q[1];
rz(-2.7064586) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1895771) q[0];
sx q[0];
rz(-0.36882419) q[0];
sx q[0];
rz(2.8394305) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.52559678) q[2];
sx q[2];
rz(-0.28738775) q[2];
sx q[2];
rz(1.4488066) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.0355465) q[1];
sx q[1];
rz(-1.4075081) q[1];
sx q[1];
rz(2.0590904) q[1];
rz(-pi) q[2];
rz(1.2262672) q[3];
sx q[3];
rz(-0.83762729) q[3];
sx q[3];
rz(-1.4602349) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.6083287) q[2];
sx q[2];
rz(-1.3003131) q[2];
sx q[2];
rz(0.30291525) q[2];
rz(-1.3251925) q[3];
sx q[3];
rz(-1.15851) q[3];
sx q[3];
rz(0.091025092) q[3];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19642297) q[0];
sx q[0];
rz(-1.7315995) q[0];
sx q[0];
rz(-0.91745013) q[0];
rz(-2.4687185) q[1];
sx q[1];
rz(-1.0854951) q[1];
sx q[1];
rz(2.8767169) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6661975) q[0];
sx q[0];
rz(-0.65984939) q[0];
sx q[0];
rz(0.049113627) q[0];
rz(-pi) q[1];
rz(2.8633966) q[2];
sx q[2];
rz(-1.858466) q[2];
sx q[2];
rz(0.96565914) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.781573) q[1];
sx q[1];
rz(-2.2847166) q[1];
sx q[1];
rz(0.54846958) q[1];
rz(-pi) q[2];
x q[2];
rz(2.360965) q[3];
sx q[3];
rz(-0.98140162) q[3];
sx q[3];
rz(-1.9115703) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.778487) q[2];
sx q[2];
rz(-2.6532756) q[2];
sx q[2];
rz(1.5650361) q[2];
rz(-2.1145084) q[3];
sx q[3];
rz(-0.74147195) q[3];
sx q[3];
rz(-1.1013793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.251579) q[0];
sx q[0];
rz(-0.13680923) q[0];
sx q[0];
rz(0.47873163) q[0];
rz(2.1084673) q[1];
sx q[1];
rz(-0.9712351) q[1];
sx q[1];
rz(-2.1889401) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6081776) q[0];
sx q[0];
rz(-1.0262283) q[0];
sx q[0];
rz(1.4334701) q[0];
x q[1];
rz(-2.7258337) q[2];
sx q[2];
rz(-1.2649049) q[2];
sx q[2];
rz(1.3410459) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.488727) q[1];
sx q[1];
rz(-1.1796724) q[1];
sx q[1];
rz(-2.8847242) q[1];
rz(-1.5289375) q[3];
sx q[3];
rz(-2.0484945) q[3];
sx q[3];
rz(2.4822513) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.4218563) q[2];
sx q[2];
rz(-0.36281261) q[2];
sx q[2];
rz(-2.6110113) q[2];
rz(1.7355708) q[3];
sx q[3];
rz(-2.0134182) q[3];
sx q[3];
rz(0.83166844) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56753165) q[0];
sx q[0];
rz(-1.4968137) q[0];
sx q[0];
rz(1.6249599) q[0];
rz(1.3051055) q[1];
sx q[1];
rz(-1.790698) q[1];
sx q[1];
rz(0.17257246) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0046376) q[0];
sx q[0];
rz(-1.1543659) q[0];
sx q[0];
rz(-0.01491551) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0748323) q[2];
sx q[2];
rz(-0.50903382) q[2];
sx q[2];
rz(-0.80392716) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3249358) q[1];
sx q[1];
rz(-1.4656193) q[1];
sx q[1];
rz(-1.9865958) q[1];
rz(0.63038007) q[3];
sx q[3];
rz(-1.7582338) q[3];
sx q[3];
rz(-2.558625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.3036348) q[2];
sx q[2];
rz(-1.6875608) q[2];
sx q[2];
rz(1.1266358) q[2];
rz(-2.3593694) q[3];
sx q[3];
rz(-1.2354847) q[3];
sx q[3];
rz(-1.8036028) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28850266) q[0];
sx q[0];
rz(-2.8350916) q[0];
sx q[0];
rz(2.4801168) q[0];
rz(2.181197) q[1];
sx q[1];
rz(-1.7405225) q[1];
sx q[1];
rz(2.3849934) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9264383) q[0];
sx q[0];
rz(-1.4190136) q[0];
sx q[0];
rz(0.093869165) q[0];
x q[1];
rz(-1.7724178) q[2];
sx q[2];
rz(-1.1668418) q[2];
sx q[2];
rz(2.902365) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7797864) q[1];
sx q[1];
rz(-0.27946073) q[1];
sx q[1];
rz(-0.96868412) q[1];
x q[2];
rz(0.19161253) q[3];
sx q[3];
rz(-1.629047) q[3];
sx q[3];
rz(0.91764698) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.0044272) q[2];
sx q[2];
rz(-2.9512773) q[2];
sx q[2];
rz(-2.9471617) q[2];
rz(0.91313177) q[3];
sx q[3];
rz(-1.3875995) q[3];
sx q[3];
rz(2.156179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59654355) q[0];
sx q[0];
rz(-0.61674917) q[0];
sx q[0];
rz(0.066666691) q[0];
rz(2.8170259) q[1];
sx q[1];
rz(-1.5044183) q[1];
sx q[1];
rz(0.98888046) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38935223) q[0];
sx q[0];
rz(-1.1712495) q[0];
sx q[0];
rz(2.7826392) q[0];
x q[1];
rz(-3.0481911) q[2];
sx q[2];
rz(-1.5867234) q[2];
sx q[2];
rz(1.2707368) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.0162184) q[1];
sx q[1];
rz(-1.6375293) q[1];
sx q[1];
rz(1.1557505) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.057468) q[3];
sx q[3];
rz(-1.1063965) q[3];
sx q[3];
rz(-1.7341136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.6797592) q[2];
sx q[2];
rz(-0.89670783) q[2];
sx q[2];
rz(-2.7339593) q[2];
rz(-2.3729825) q[3];
sx q[3];
rz(-1.8278443) q[3];
sx q[3];
rz(2.2176946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.75893629) q[0];
sx q[0];
rz(-1.8396682) q[0];
sx q[0];
rz(0.60920238) q[0];
rz(-0.095104782) q[1];
sx q[1];
rz(-1.8895878) q[1];
sx q[1];
rz(2.2682155) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3082335) q[0];
sx q[0];
rz(-1.5033659) q[0];
sx q[0];
rz(1.7504577) q[0];
rz(-2.9940669) q[2];
sx q[2];
rz(-1.4359056) q[2];
sx q[2];
rz(2.877176) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.1904618) q[1];
sx q[1];
rz(-1.8776263) q[1];
sx q[1];
rz(0.61913403) q[1];
rz(-pi) q[2];
x q[2];
rz(0.43975131) q[3];
sx q[3];
rz(-2.1034749) q[3];
sx q[3];
rz(-0.38910481) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0218899) q[2];
sx q[2];
rz(-0.78803524) q[2];
sx q[2];
rz(2.4592887) q[2];
rz(0.37426379) q[3];
sx q[3];
rz(-1.5687317) q[3];
sx q[3];
rz(-0.16690978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4436214) q[0];
sx q[0];
rz(-1.3598096) q[0];
sx q[0];
rz(-0.95296729) q[0];
rz(-2.3151746) q[1];
sx q[1];
rz(-0.73917878) q[1];
sx q[1];
rz(1.3964765) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5293225) q[0];
sx q[0];
rz(-1.4850052) q[0];
sx q[0];
rz(-0.32353185) q[0];
rz(-pi) q[1];
rz(1.868532) q[2];
sx q[2];
rz(-2.2072788) q[2];
sx q[2];
rz(2.0993078) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.1435946) q[1];
sx q[1];
rz(-1.8815787) q[1];
sx q[1];
rz(2.8423611) q[1];
rz(-2.9328437) q[3];
sx q[3];
rz(-0.26780805) q[3];
sx q[3];
rz(2.1684614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.46618) q[2];
sx q[2];
rz(-2.7853577) q[2];
sx q[2];
rz(0.15979016) q[2];
rz(2.8397078) q[3];
sx q[3];
rz(-0.92697898) q[3];
sx q[3];
rz(0.2872428) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0970584) q[0];
sx q[0];
rz(-0.67561588) q[0];
sx q[0];
rz(-1.5560879) q[0];
rz(-3.0083169) q[1];
sx q[1];
rz(-1.6242846) q[1];
sx q[1];
rz(-0.12856738) q[1];
rz(-2.2014387) q[2];
sx q[2];
rz(-1.4174145) q[2];
sx q[2];
rz(-2.6819475) q[2];
rz(-0.53363579) q[3];
sx q[3];
rz(-1.5279557) q[3];
sx q[3];
rz(0.91844311) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
