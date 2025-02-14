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
rz(-0.67796081) q[0];
sx q[0];
rz(-0.27422658) q[0];
sx q[0];
rz(1.8507313) q[0];
rz(-0.21386799) q[1];
sx q[1];
rz(-2.8254421) q[1];
sx q[1];
rz(-1.6419799) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8052657) q[0];
sx q[0];
rz(-1.299843) q[0];
sx q[0];
rz(-1.8150369) q[0];
x q[1];
rz(-0.77969061) q[2];
sx q[2];
rz(-1.3609972) q[2];
sx q[2];
rz(-2.9203134) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5747524) q[1];
sx q[1];
rz(-0.87927239) q[1];
sx q[1];
rz(-1.994446) q[1];
rz(-pi) q[2];
rz(-0.94120195) q[3];
sx q[3];
rz(-1.1735907) q[3];
sx q[3];
rz(1.3324228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.9700254) q[2];
sx q[2];
rz(-2.1000523) q[2];
sx q[2];
rz(2.2402666) q[2];
rz(2.6775635) q[3];
sx q[3];
rz(-1.8601067) q[3];
sx q[3];
rz(3.071781) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0932015) q[0];
sx q[0];
rz(-0.036660107) q[0];
sx q[0];
rz(-1.8638336) q[0];
rz(-3.0473862) q[1];
sx q[1];
rz(-2.6779046) q[1];
sx q[1];
rz(1.6069848) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8746895) q[0];
sx q[0];
rz(-2.1703566) q[0];
sx q[0];
rz(0.039681704) q[0];
rz(-0.056939967) q[2];
sx q[2];
rz(-1.7527662) q[2];
sx q[2];
rz(-1.6873311) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.31904768) q[1];
sx q[1];
rz(-1.8314059) q[1];
sx q[1];
rz(2.6316363) q[1];
rz(-pi) q[2];
rz(-1.1538131) q[3];
sx q[3];
rz(-1.1975884) q[3];
sx q[3];
rz(2.5905379) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.719912) q[2];
sx q[2];
rz(-1.9399425) q[2];
sx q[2];
rz(-0.16033944) q[2];
rz(-0.73873377) q[3];
sx q[3];
rz(-2.2952047) q[3];
sx q[3];
rz(-1.4275449) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5147603) q[0];
sx q[0];
rz(-1.4381831) q[0];
sx q[0];
rz(-2.6318188) q[0];
rz(-1.431142) q[1];
sx q[1];
rz(-1.9606083) q[1];
sx q[1];
rz(-2.3950155) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1603161) q[0];
sx q[0];
rz(-1.0183987) q[0];
sx q[0];
rz(-3.1180326) q[0];
rz(-0.66303894) q[2];
sx q[2];
rz(-2.2929077) q[2];
sx q[2];
rz(0.3438562) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.4164734) q[1];
sx q[1];
rz(-1.2095569) q[1];
sx q[1];
rz(2.3010985) q[1];
rz(-pi) q[2];
x q[2];
rz(0.88991092) q[3];
sx q[3];
rz(-1.2632252) q[3];
sx q[3];
rz(-2.6079082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.15472445) q[2];
sx q[2];
rz(-2.8747323) q[2];
sx q[2];
rz(1.9179087) q[2];
rz(1.0736505) q[3];
sx q[3];
rz(-1.4370388) q[3];
sx q[3];
rz(-0.8980208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8680442) q[0];
sx q[0];
rz(-0.88669625) q[0];
sx q[0];
rz(0.37937382) q[0];
rz(2.9647478) q[1];
sx q[1];
rz(-1.6601446) q[1];
sx q[1];
rz(-2.3462229) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.39033) q[0];
sx q[0];
rz(-2.7613233) q[0];
sx q[0];
rz(1.3185391) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3248939) q[2];
sx q[2];
rz(-1.3743322) q[2];
sx q[2];
rz(0.24345556) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.6713555) q[1];
sx q[1];
rz(-2.0908326) q[1];
sx q[1];
rz(0.29754559) q[1];
x q[2];
rz(0.42511149) q[3];
sx q[3];
rz(-0.13775857) q[3];
sx q[3];
rz(0.83168304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.4312326) q[2];
sx q[2];
rz(-1.7962339) q[2];
sx q[2];
rz(-1.360652) q[2];
rz(2.1121173) q[3];
sx q[3];
rz(-1.6607213) q[3];
sx q[3];
rz(1.3220538) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65115702) q[0];
sx q[0];
rz(-1.5999726) q[0];
sx q[0];
rz(1.5486451) q[0];
rz(-2.7410638) q[1];
sx q[1];
rz(-1.488204) q[1];
sx q[1];
rz(-1.4097479) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9841138) q[0];
sx q[0];
rz(-1.7703919) q[0];
sx q[0];
rz(0.73961135) q[0];
rz(2.7422041) q[2];
sx q[2];
rz(-2.7542973) q[2];
sx q[2];
rz(-0.072810955) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.21239195) q[1];
sx q[1];
rz(-2.1521795) q[1];
sx q[1];
rz(0.53668944) q[1];
rz(3.0733969) q[3];
sx q[3];
rz(-2.1917097) q[3];
sx q[3];
rz(-0.62486006) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.8283525) q[2];
sx q[2];
rz(-1.7007549) q[2];
sx q[2];
rz(-2.7102615) q[2];
rz(-3.0865772) q[3];
sx q[3];
rz(-0.31156817) q[3];
sx q[3];
rz(0.55606786) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7677652) q[0];
sx q[0];
rz(-2.8887833) q[0];
sx q[0];
rz(-2.1164236) q[0];
rz(-2.688736) q[1];
sx q[1];
rz(-0.57944524) q[1];
sx q[1];
rz(-0.90389171) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65176455) q[0];
sx q[0];
rz(-1.2625041) q[0];
sx q[0];
rz(2.2405008) q[0];
x q[1];
rz(-2.3329817) q[2];
sx q[2];
rz(-2.3249194) q[2];
sx q[2];
rz(-2.6220235) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.55906193) q[1];
sx q[1];
rz(-2.2901223) q[1];
sx q[1];
rz(1.2079617) q[1];
rz(-pi) q[2];
x q[2];
rz(0.76670209) q[3];
sx q[3];
rz(-0.23582102) q[3];
sx q[3];
rz(-0.23877777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.30770939) q[2];
sx q[2];
rz(-1.4557975) q[2];
sx q[2];
rz(-2.9252388) q[2];
rz(0.89573914) q[3];
sx q[3];
rz(-0.68550617) q[3];
sx q[3];
rz(-1.5007796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1322587) q[0];
sx q[0];
rz(-2.0348771) q[0];
sx q[0];
rz(-2.4427781) q[0];
rz(1.1837333) q[1];
sx q[1];
rz(-1.4212757) q[1];
sx q[1];
rz(-2.2898477) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4119371) q[0];
sx q[0];
rz(-1.0348399) q[0];
sx q[0];
rz(-0.044117731) q[0];
rz(1.5450655) q[2];
sx q[2];
rz(-0.99359578) q[2];
sx q[2];
rz(2.4943697) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.3892365) q[1];
sx q[1];
rz(-2.3515278) q[1];
sx q[1];
rz(0.52100928) q[1];
rz(-pi) q[2];
rz(2.5301039) q[3];
sx q[3];
rz(-2.4019255) q[3];
sx q[3];
rz(0.31114331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.8375887) q[2];
sx q[2];
rz(-1.8303266) q[2];
sx q[2];
rz(2.6336929) q[2];
rz(-1.0287644) q[3];
sx q[3];
rz(-2.2977836) q[3];
sx q[3];
rz(-0.60104162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69865882) q[0];
sx q[0];
rz(-1.3154987) q[0];
sx q[0];
rz(0.0096631924) q[0];
rz(1.6161605) q[1];
sx q[1];
rz(-1.7639672) q[1];
sx q[1];
rz(-0.38472167) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54668173) q[0];
sx q[0];
rz(-1.4230886) q[0];
sx q[0];
rz(-1.3315931) q[0];
x q[1];
rz(1.184395) q[2];
sx q[2];
rz(-1.7586054) q[2];
sx q[2];
rz(1.9822497) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.2710451) q[1];
sx q[1];
rz(-1.7779852) q[1];
sx q[1];
rz(1.8533587) q[1];
rz(-pi) q[2];
rz(-2.1423903) q[3];
sx q[3];
rz(-0.88547844) q[3];
sx q[3];
rz(2.9504321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.5821417) q[2];
sx q[2];
rz(-2.1820575) q[2];
sx q[2];
rz(0.01586308) q[2];
rz(1.9150241) q[3];
sx q[3];
rz(-0.80569402) q[3];
sx q[3];
rz(0.62172833) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3847619) q[0];
sx q[0];
rz(-0.66487304) q[0];
sx q[0];
rz(2.050198) q[0];
rz(-0.81820828) q[1];
sx q[1];
rz(-1.3312157) q[1];
sx q[1];
rz(-2.4726726) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92907897) q[0];
sx q[0];
rz(-0.90796472) q[0];
sx q[0];
rz(2.3431091) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7162076) q[2];
sx q[2];
rz(-2.2375319) q[2];
sx q[2];
rz(-0.83620706) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3260169) q[1];
sx q[1];
rz(-1.6417802) q[1];
sx q[1];
rz(1.0128338) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0505597) q[3];
sx q[3];
rz(-0.63891131) q[3];
sx q[3];
rz(2.0990438) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.5074671) q[2];
sx q[2];
rz(-0.52906817) q[2];
sx q[2];
rz(-0.96366209) q[2];
rz(-1.4108747) q[3];
sx q[3];
rz(-1.4286634) q[3];
sx q[3];
rz(-1.3655519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.016634781) q[0];
sx q[0];
rz(-0.5683012) q[0];
sx q[0];
rz(-1.4547263) q[0];
rz(2.0330873) q[1];
sx q[1];
rz(-1.5364372) q[1];
sx q[1];
rz(2.813521) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5660368) q[0];
sx q[0];
rz(-1.6854316) q[0];
sx q[0];
rz(0.057997142) q[0];
rz(1.6871618) q[2];
sx q[2];
rz(-1.7504217) q[2];
sx q[2];
rz(-2.3723974) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.1197443) q[1];
sx q[1];
rz(-1.9520063) q[1];
sx q[1];
rz(-0.66701835) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0193908) q[3];
sx q[3];
rz(-0.20272045) q[3];
sx q[3];
rz(2.3644476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.6803117) q[2];
sx q[2];
rz(-1.842247) q[2];
sx q[2];
rz(0.4168365) q[2];
rz(2.4428115) q[3];
sx q[3];
rz(-2.4650033) q[3];
sx q[3];
rz(-0.39066395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9417435) q[0];
sx q[0];
rz(-1.9244292) q[0];
sx q[0];
rz(1.695965) q[0];
rz(-0.70882123) q[1];
sx q[1];
rz(-2.2292021) q[1];
sx q[1];
rz(3.0880047) q[1];
rz(2.9403953) q[2];
sx q[2];
rz(-0.94379776) q[2];
sx q[2];
rz(-1.1371053) q[2];
rz(2.3117456) q[3];
sx q[3];
rz(-1.1329069) q[3];
sx q[3];
rz(-2.2923242) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
