OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.8013826) q[0];
sx q[0];
rz(-0.27685452) q[0];
sx q[0];
rz(-1.0387596) q[0];
rz(-0.34173319) q[1];
sx q[1];
rz(3.9781456) q[1];
sx q[1];
rz(9.8415924) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3562389) q[0];
sx q[0];
rz(-2.2921341) q[0];
sx q[0];
rz(-1.712404) q[0];
rz(2.7567467) q[2];
sx q[2];
rz(-0.59078465) q[2];
sx q[2];
rz(0.42921517) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3089247) q[1];
sx q[1];
rz(-2.1410258) q[1];
sx q[1];
rz(-0.69242386) q[1];
rz(-pi) q[2];
x q[2];
rz(0.19617041) q[3];
sx q[3];
rz(-0.16465287) q[3];
sx q[3];
rz(0.97033721) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.9716399) q[2];
sx q[2];
rz(-1.6728741) q[2];
sx q[2];
rz(1.2539585) q[2];
rz(1.5597255) q[3];
sx q[3];
rz(-2.4032148) q[3];
sx q[3];
rz(2.019465) q[3];
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
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1935683) q[0];
sx q[0];
rz(-1.3688315) q[0];
sx q[0];
rz(-2.3728306) q[0];
rz(-1.7747152) q[1];
sx q[1];
rz(-1.3221075) q[1];
sx q[1];
rz(2.1034525) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65377849) q[0];
sx q[0];
rz(-3.1285435) q[0];
sx q[0];
rz(-2.292146) q[0];
rz(-pi) q[1];
rz(0.409276) q[2];
sx q[2];
rz(-2.3774494) q[2];
sx q[2];
rz(-2.0913569) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.86377108) q[1];
sx q[1];
rz(-1.9367332) q[1];
sx q[1];
rz(-2.4863431) q[1];
rz(-pi) q[2];
rz(-1.6065381) q[3];
sx q[3];
rz(-0.48887353) q[3];
sx q[3];
rz(-1.3066178) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.64608964) q[2];
sx q[2];
rz(-1.8547408) q[2];
sx q[2];
rz(2.3213279) q[2];
rz(-0.25343728) q[3];
sx q[3];
rz(-0.35496747) q[3];
sx q[3];
rz(-2.7804815) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(-2.8609817) q[0];
sx q[0];
rz(-0.51787037) q[0];
sx q[0];
rz(0.88731998) q[0];
rz(-2.6112828) q[1];
sx q[1];
rz(-0.92875004) q[1];
sx q[1];
rz(2.0022154) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4654008) q[0];
sx q[0];
rz(-1.7125107) q[0];
sx q[0];
rz(-0.34872524) q[0];
rz(-2.1043037) q[2];
sx q[2];
rz(-0.59743687) q[2];
sx q[2];
rz(-3.1410599) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.92147747) q[1];
sx q[1];
rz(-0.24610983) q[1];
sx q[1];
rz(-3.0819478) q[1];
rz(-pi) q[2];
rz(0.51872571) q[3];
sx q[3];
rz(-1.5253807) q[3];
sx q[3];
rz(0.96744999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.016971074) q[2];
sx q[2];
rz(-1.4554224) q[2];
sx q[2];
rz(2.7023081) q[2];
rz(0.47075054) q[3];
sx q[3];
rz(-3.0622523) q[3];
sx q[3];
rz(1.9973756) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73862326) q[0];
sx q[0];
rz(-2.2082177) q[0];
sx q[0];
rz(-0.11548197) q[0];
rz(0.11784095) q[1];
sx q[1];
rz(-1.203048) q[1];
sx q[1];
rz(-1.8353362) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7457434) q[0];
sx q[0];
rz(-1.9877399) q[0];
sx q[0];
rz(0.31503079) q[0];
rz(-0.66008644) q[2];
sx q[2];
rz(-1.8750817) q[2];
sx q[2];
rz(0.10709281) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1237643) q[1];
sx q[1];
rz(-0.99010795) q[1];
sx q[1];
rz(0.75947096) q[1];
x q[2];
rz(-1.4545069) q[3];
sx q[3];
rz(-2.5657885) q[3];
sx q[3];
rz(-0.35279122) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8370886) q[2];
sx q[2];
rz(-1.8279165) q[2];
sx q[2];
rz(-0.15466776) q[2];
rz(1.6020487) q[3];
sx q[3];
rz(-1.3402901) q[3];
sx q[3];
rz(-0.60607564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-2.5681169) q[0];
sx q[0];
rz(-2.0982168) q[0];
sx q[0];
rz(-2.306275) q[0];
rz(-0.71594816) q[1];
sx q[1];
rz(-0.42778152) q[1];
sx q[1];
rz(3.0221525) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1560005) q[0];
sx q[0];
rz(-1.553926) q[0];
sx q[0];
rz(1.6266247) q[0];
x q[1];
rz(-1.4177464) q[2];
sx q[2];
rz(-2.0156392) q[2];
sx q[2];
rz(-0.21013021) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.2970256) q[1];
sx q[1];
rz(-0.78999619) q[1];
sx q[1];
rz(-0.67761919) q[1];
x q[2];
rz(2.4923513) q[3];
sx q[3];
rz(-0.63707817) q[3];
sx q[3];
rz(-1.0604881) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.9014088) q[2];
sx q[2];
rz(-2.8808424) q[2];
sx q[2];
rz(0.060997941) q[2];
rz(-0.048246233) q[3];
sx q[3];
rz(-1.2333074) q[3];
sx q[3];
rz(-0.72859305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.401684) q[0];
sx q[0];
rz(-0.70867276) q[0];
sx q[0];
rz(2.0106864) q[0];
rz(2.6629958) q[1];
sx q[1];
rz(-2.5391948) q[1];
sx q[1];
rz(1.7344249) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81086195) q[0];
sx q[0];
rz(-1.6243582) q[0];
sx q[0];
rz(-1.5723482) q[0];
rz(-pi) q[1];
rz(-2.8596609) q[2];
sx q[2];
rz(-0.86840668) q[2];
sx q[2];
rz(-1.3739746) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.9764912) q[1];
sx q[1];
rz(-2.8860984) q[1];
sx q[1];
rz(-1.1063904) q[1];
x q[2];
rz(0.21295548) q[3];
sx q[3];
rz(-1.3939438) q[3];
sx q[3];
rz(0.18421728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.56167928) q[2];
sx q[2];
rz(-2.734197) q[2];
sx q[2];
rz(-1.178406) q[2];
rz(2.0922349) q[3];
sx q[3];
rz(-2.6047843) q[3];
sx q[3];
rz(1.2843081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2105763) q[0];
sx q[0];
rz(-1.1897621) q[0];
sx q[0];
rz(-1.3354906) q[0];
rz(1.8203075) q[1];
sx q[1];
rz(-2.0704465) q[1];
sx q[1];
rz(0.30805045) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46522147) q[0];
sx q[0];
rz(-1.9587059) q[0];
sx q[0];
rz(1.3419536) q[0];
rz(-0.63703434) q[2];
sx q[2];
rz(-0.61372988) q[2];
sx q[2];
rz(1.5368652) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.25702204) q[1];
sx q[1];
rz(-1.0858469) q[1];
sx q[1];
rz(-2.229175) q[1];
rz(-pi) q[2];
rz(-3.0940476) q[3];
sx q[3];
rz(-1.9564087) q[3];
sx q[3];
rz(1.4908718) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.41219741) q[2];
sx q[2];
rz(-0.9407548) q[2];
sx q[2];
rz(-0.13128734) q[2];
rz(0.73733759) q[3];
sx q[3];
rz(-1.8359343) q[3];
sx q[3];
rz(2.9837515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3103631) q[0];
sx q[0];
rz(-1.0897626) q[0];
sx q[0];
rz(0.53037733) q[0];
rz(1.4250379) q[1];
sx q[1];
rz(-1.1385463) q[1];
sx q[1];
rz(-2.4200965) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.90008) q[0];
sx q[0];
rz(-0.78966138) q[0];
sx q[0];
rz(-0.65171839) q[0];
rz(-pi) q[1];
rz(-2.6725298) q[2];
sx q[2];
rz(-1.7028168) q[2];
sx q[2];
rz(-1.309883) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2635195) q[1];
sx q[1];
rz(-1.7612447) q[1];
sx q[1];
rz(1.3506804) q[1];
rz(0.11105342) q[3];
sx q[3];
rz(-0.78563443) q[3];
sx q[3];
rz(1.812404) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.7044652) q[2];
sx q[2];
rz(-2.7361054) q[2];
sx q[2];
rz(-1.4777769) q[2];
rz(-1.1635228) q[3];
sx q[3];
rz(-1.082837) q[3];
sx q[3];
rz(-1.6400953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8116233) q[0];
sx q[0];
rz(-1.4412619) q[0];
sx q[0];
rz(-0.60923088) q[0];
rz(1.5902663) q[1];
sx q[1];
rz(-0.32569277) q[1];
sx q[1];
rz(-1.4564266) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1066172) q[0];
sx q[0];
rz(-1.2143232) q[0];
sx q[0];
rz(-2.8750505) q[0];
x q[1];
rz(-0.48522075) q[2];
sx q[2];
rz(-1.5169797) q[2];
sx q[2];
rz(-0.14419989) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.2472154) q[1];
sx q[1];
rz(-0.4900107) q[1];
sx q[1];
rz(1.9202597) q[1];
rz(-pi) q[2];
x q[2];
rz(0.84734579) q[3];
sx q[3];
rz(-1.0718126) q[3];
sx q[3];
rz(-0.33667281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.41813254) q[2];
sx q[2];
rz(-0.89833608) q[2];
sx q[2];
rz(-2.2612803) q[2];
rz(-0.20600016) q[3];
sx q[3];
rz(-2.7166631) q[3];
sx q[3];
rz(2.6515085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3708165) q[0];
sx q[0];
rz(-2.4801065) q[0];
sx q[0];
rz(-3.1296375) q[0];
rz(-1.6318343) q[1];
sx q[1];
rz(-2.6963574) q[1];
sx q[1];
rz(-0.59648046) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76145455) q[0];
sx q[0];
rz(-1.5517457) q[0];
sx q[0];
rz(-2.6721902) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8305186) q[2];
sx q[2];
rz(-1.2190281) q[2];
sx q[2];
rz(0.042366926) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.0582907) q[1];
sx q[1];
rz(-1.6495678) q[1];
sx q[1];
rz(2.0473256) q[1];
rz(-pi) q[2];
rz(0.29774547) q[3];
sx q[3];
rz(-1.7216873) q[3];
sx q[3];
rz(-2.6835364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.97800469) q[2];
sx q[2];
rz(-0.96077335) q[2];
sx q[2];
rz(-0.19700024) q[2];
rz(-1.152285) q[3];
sx q[3];
rz(-0.14324337) q[3];
sx q[3];
rz(0.44767374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
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
rz(-0.86354179) q[0];
sx q[0];
rz(-2.1320237) q[0];
sx q[0];
rz(-0.19436819) q[0];
rz(2.9327783) q[1];
sx q[1];
rz(-1.5369692) q[1];
sx q[1];
rz(2.1280638) q[1];
rz(-3.0473982) q[2];
sx q[2];
rz(-2.258156) q[2];
sx q[2];
rz(-1.5428262) q[2];
rz(0.87416762) q[3];
sx q[3];
rz(-0.71577358) q[3];
sx q[3];
rz(-1.7795455) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
