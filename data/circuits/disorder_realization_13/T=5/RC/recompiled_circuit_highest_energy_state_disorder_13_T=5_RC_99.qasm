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
rz(-2.9450671) q[0];
sx q[0];
rz(-2.3152469) q[0];
sx q[0];
rz(2.2235121) q[0];
rz(-0.11153587) q[1];
sx q[1];
rz(-1.3476975) q[1];
sx q[1];
rz(-1.5741875) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8188326) q[0];
sx q[0];
rz(-0.010060223) q[0];
sx q[0];
rz(-1.9181817) q[0];
rz(-pi) q[1];
rz(1.1996881) q[2];
sx q[2];
rz(-2.1016309) q[2];
sx q[2];
rz(0.67000721) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.35167978) q[1];
sx q[1];
rz(-1.6761629) q[1];
sx q[1];
rz(-1.8743452) q[1];
rz(-pi) q[2];
rz(0.2229177) q[3];
sx q[3];
rz(-0.50361982) q[3];
sx q[3];
rz(0.38583392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.7142746) q[2];
sx q[2];
rz(-1.6552507) q[2];
sx q[2];
rz(-1.8753258) q[2];
rz(2.0990939) q[3];
sx q[3];
rz(-1.5407341) q[3];
sx q[3];
rz(-1.3955759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0365486) q[0];
sx q[0];
rz(-0.62750134) q[0];
sx q[0];
rz(3.0446766) q[0];
rz(2.3333683) q[1];
sx q[1];
rz(-0.40882912) q[1];
sx q[1];
rz(-1.854863) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5135749) q[0];
sx q[0];
rz(-2.6903408) q[0];
sx q[0];
rz(-1.6523163) q[0];
rz(-2.497235) q[2];
sx q[2];
rz(-0.8335118) q[2];
sx q[2];
rz(2.4887511) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.8342487) q[1];
sx q[1];
rz(-1.1888224) q[1];
sx q[1];
rz(2.3267507) q[1];
rz(-pi) q[2];
rz(-1.2169514) q[3];
sx q[3];
rz(-1.2520188) q[3];
sx q[3];
rz(0.48976605) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.6404932) q[2];
sx q[2];
rz(-2.2576136) q[2];
sx q[2];
rz(-2.7596149) q[2];
rz(0.52465087) q[3];
sx q[3];
rz(-1.7347521) q[3];
sx q[3];
rz(-2.9978571) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76661888) q[0];
sx q[0];
rz(-1.5856278) q[0];
sx q[0];
rz(-2.4554456) q[0];
rz(1.067591) q[1];
sx q[1];
rz(-0.99718863) q[1];
sx q[1];
rz(3.0608665) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2093182) q[0];
sx q[0];
rz(-0.22695146) q[0];
sx q[0];
rz(-1.7448241) q[0];
rz(2.3600134) q[2];
sx q[2];
rz(-0.8197166) q[2];
sx q[2];
rz(0.21910659) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.9517652) q[1];
sx q[1];
rz(-1.8009342) q[1];
sx q[1];
rz(1.5946424) q[1];
rz(-pi) q[2];
rz(-1.0992381) q[3];
sx q[3];
rz(-0.84065372) q[3];
sx q[3];
rz(-1.7748347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.8778235) q[2];
sx q[2];
rz(-0.70085415) q[2];
sx q[2];
rz(0.43886718) q[2];
rz(-1.2980488) q[3];
sx q[3];
rz(-0.38709199) q[3];
sx q[3];
rz(-0.41518655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91442672) q[0];
sx q[0];
rz(-0.54410797) q[0];
sx q[0];
rz(-2.4650204) q[0];
rz(-3.0034972) q[1];
sx q[1];
rz(-2.0510249) q[1];
sx q[1];
rz(0.20060435) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4355136) q[0];
sx q[0];
rz(-1.9301751) q[0];
sx q[0];
rz(0.91889221) q[0];
x q[1];
rz(-1.3820962) q[2];
sx q[2];
rz(-2.1917135) q[2];
sx q[2];
rz(3.0953654) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.8183867) q[1];
sx q[1];
rz(-0.96508677) q[1];
sx q[1];
rz(-1.731864) q[1];
x q[2];
rz(1.4273321) q[3];
sx q[3];
rz(-2.1273566) q[3];
sx q[3];
rz(-1.0866764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.6382008) q[2];
sx q[2];
rz(-1.8663422) q[2];
sx q[2];
rz(-1.4208043) q[2];
rz(-3.0270789) q[3];
sx q[3];
rz(-1.4736466) q[3];
sx q[3];
rz(2.1421471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3223709) q[0];
sx q[0];
rz(-0.79077661) q[0];
sx q[0];
rz(0.97212273) q[0];
rz(0.32360336) q[1];
sx q[1];
rz(-1.4515896) q[1];
sx q[1];
rz(1.9591029) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2401449) q[0];
sx q[0];
rz(-1.6109544) q[0];
sx q[0];
rz(1.3754649) q[0];
x q[1];
rz(1.576831) q[2];
sx q[2];
rz(-2.1409935) q[2];
sx q[2];
rz(1.6369199) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.8132434) q[1];
sx q[1];
rz(-2.9651838) q[1];
sx q[1];
rz(-1.117068) q[1];
rz(-1.0504405) q[3];
sx q[3];
rz(-2.0791868) q[3];
sx q[3];
rz(0.45724487) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.16367308) q[2];
sx q[2];
rz(-1.7719496) q[2];
sx q[2];
rz(0.53691205) q[2];
rz(2.4162857) q[3];
sx q[3];
rz(-3.091843) q[3];
sx q[3];
rz(-2.3427826) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86730114) q[0];
sx q[0];
rz(-2.0766356) q[0];
sx q[0];
rz(-1.6987479) q[0];
rz(2.9310215) q[1];
sx q[1];
rz(-1.6981533) q[1];
sx q[1];
rz(-0.64705667) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0083425) q[0];
sx q[0];
rz(-1.1971336) q[0];
sx q[0];
rz(1.295412) q[0];
rz(0.17433106) q[2];
sx q[2];
rz(-0.56340471) q[2];
sx q[2];
rz(2.8772815) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.9155026) q[1];
sx q[1];
rz(-1.9880297) q[1];
sx q[1];
rz(-1.9482062) q[1];
x q[2];
rz(1.8451377) q[3];
sx q[3];
rz(-1.8963433) q[3];
sx q[3];
rz(-0.59301585) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.1232542) q[2];
sx q[2];
rz(-2.047057) q[2];
sx q[2];
rz(-1.1798165) q[2];
rz(1.2849464) q[3];
sx q[3];
rz(-2.732087) q[3];
sx q[3];
rz(-0.063260945) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0354075) q[0];
sx q[0];
rz(-1.4466865) q[0];
sx q[0];
rz(2.2440198) q[0];
rz(1.3342185) q[1];
sx q[1];
rz(-1.370627) q[1];
sx q[1];
rz(2.1468377) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3883427) q[0];
sx q[0];
rz(-2.3412447) q[0];
sx q[0];
rz(2.0092756) q[0];
rz(-pi) q[1];
rz(2.1734235) q[2];
sx q[2];
rz(-1.954284) q[2];
sx q[2];
rz(0.17959514) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.89284507) q[1];
sx q[1];
rz(-2.020128) q[1];
sx q[1];
rz(2.8161418) q[1];
rz(-pi) q[2];
rz(-3.1174421) q[3];
sx q[3];
rz(-1.8514048) q[3];
sx q[3];
rz(-0.58089549) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.083954088) q[2];
sx q[2];
rz(-1.8123764) q[2];
sx q[2];
rz(0.29339054) q[2];
rz(-0.053226274) q[3];
sx q[3];
rz(-0.86881995) q[3];
sx q[3];
rz(-1.4330385) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5695802) q[0];
sx q[0];
rz(-1.7520289) q[0];
sx q[0];
rz(0.30297512) q[0];
rz(1.6527893) q[1];
sx q[1];
rz(-1.8257414) q[1];
sx q[1];
rz(2.4724919) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8527864) q[0];
sx q[0];
rz(-2.0159205) q[0];
sx q[0];
rz(1.2034125) q[0];
rz(-pi) q[1];
rz(-0.24208787) q[2];
sx q[2];
rz(-1.228516) q[2];
sx q[2];
rz(0.9108327) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.2337906) q[1];
sx q[1];
rz(-2.0747342) q[1];
sx q[1];
rz(-0.48598098) q[1];
x q[2];
rz(2.7176526) q[3];
sx q[3];
rz(-1.3518847) q[3];
sx q[3];
rz(-0.11737897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.9992493) q[2];
sx q[2];
rz(-1.3331022) q[2];
sx q[2];
rz(2.2744501) q[2];
rz(-2.0182746) q[3];
sx q[3];
rz(-0.85223782) q[3];
sx q[3];
rz(1.999202) q[3];
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
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2808696) q[0];
sx q[0];
rz(-0.46361247) q[0];
sx q[0];
rz(-3.0238357) q[0];
rz(0.87751687) q[1];
sx q[1];
rz(-2.120647) q[1];
sx q[1];
rz(2.1868475) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68761008) q[0];
sx q[0];
rz(-1.9508525) q[0];
sx q[0];
rz(-2.4847772) q[0];
rz(-2.0659201) q[2];
sx q[2];
rz(-1.9946163) q[2];
sx q[2];
rz(2.0002805) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.6462209) q[1];
sx q[1];
rz(-2.2471161) q[1];
sx q[1];
rz(-0.27797525) q[1];
x q[2];
rz(-0.8189965) q[3];
sx q[3];
rz(-1.0291417) q[3];
sx q[3];
rz(-2.0217232) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.8942326) q[2];
sx q[2];
rz(-2.8017513) q[2];
sx q[2];
rz(-1.6575238) q[2];
rz(-1.6190716) q[3];
sx q[3];
rz(-1.3831474) q[3];
sx q[3];
rz(-1.1744261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5840983) q[0];
sx q[0];
rz(-1.431594) q[0];
sx q[0];
rz(0.75310055) q[0];
rz(1.9901216) q[1];
sx q[1];
rz(-0.52422062) q[1];
sx q[1];
rz(1.6023191) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1881018) q[0];
sx q[0];
rz(-1.3878146) q[0];
sx q[0];
rz(-1.1142) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.98277715) q[2];
sx q[2];
rz(-2.3522107) q[2];
sx q[2];
rz(-3.0510117) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6021612) q[1];
sx q[1];
rz(-1.185158) q[1];
sx q[1];
rz(-1.3630609) q[1];
x q[2];
rz(2.8026641) q[3];
sx q[3];
rz(-0.68343168) q[3];
sx q[3];
rz(-0.11250699) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.3573542) q[2];
sx q[2];
rz(-2.9604762) q[2];
sx q[2];
rz(2.6757346) q[2];
rz(-1.0957796) q[3];
sx q[3];
rz(-2.1491137) q[3];
sx q[3];
rz(2.5994658) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.096238484) q[0];
sx q[0];
rz(-1.5662554) q[0];
sx q[0];
rz(1.5731496) q[0];
rz(-2.0737598) q[1];
sx q[1];
rz(-0.44034958) q[1];
sx q[1];
rz(-0.82028295) q[1];
rz(-1.3374511) q[2];
sx q[2];
rz(-1.3139616) q[2];
sx q[2];
rz(-1.4001605) q[2];
rz(2.3798306) q[3];
sx q[3];
rz(-0.69275012) q[3];
sx q[3];
rz(2.6150273) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
