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
rz(0.077724783) q[0];
sx q[0];
rz(-2.1169777) q[0];
sx q[0];
rz(1.8416564) q[0];
rz(2.6990702) q[1];
sx q[1];
rz(4.2808851) q[1];
sx q[1];
rz(12.130375) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3740957) q[0];
sx q[0];
rz(-2.7347235) q[0];
sx q[0];
rz(-1.8396729) q[0];
rz(-pi) q[1];
rz(2.9098815) q[2];
sx q[2];
rz(-0.75318906) q[2];
sx q[2];
rz(1.1540268) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.2548542) q[1];
sx q[1];
rz(-1.5553932) q[1];
sx q[1];
rz(-2.7387357) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0150839) q[3];
sx q[3];
rz(-0.96130575) q[3];
sx q[3];
rz(1.786527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.062833) q[2];
sx q[2];
rz(-1.6597513) q[2];
sx q[2];
rz(3.1256342) q[2];
rz(1.6503085) q[3];
sx q[3];
rz(-1.242638) q[3];
sx q[3];
rz(-2.7119467) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66351873) q[0];
sx q[0];
rz(-1.0915382) q[0];
sx q[0];
rz(1.7953405) q[0];
rz(-1.9740055) q[1];
sx q[1];
rz(-1.0941894) q[1];
sx q[1];
rz(1.8928554) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25062171) q[0];
sx q[0];
rz(-2.3912794) q[0];
sx q[0];
rz(-0.98558275) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9850626) q[2];
sx q[2];
rz(-1.4505149) q[2];
sx q[2];
rz(2.2379654) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.75227164) q[1];
sx q[1];
rz(-1.5157425) q[1];
sx q[1];
rz(1.9717713) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1571659) q[3];
sx q[3];
rz(-1.274144) q[3];
sx q[3];
rz(0.082060952) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.0060129082) q[2];
sx q[2];
rz(-0.69301444) q[2];
sx q[2];
rz(-2.6683624) q[2];
rz(3.0114975) q[3];
sx q[3];
rz(-1.3869163) q[3];
sx q[3];
rz(-0.010802833) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3267645) q[0];
sx q[0];
rz(-2.9846551) q[0];
sx q[0];
rz(0.21406315) q[0];
rz(-1.7474984) q[1];
sx q[1];
rz(-2.1653445) q[1];
sx q[1];
rz(-3.0551547) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28537108) q[0];
sx q[0];
rz(-0.32111327) q[0];
sx q[0];
rz(-1.1620528) q[0];
rz(-pi) q[1];
rz(0.8930703) q[2];
sx q[2];
rz(-2.0166335) q[2];
sx q[2];
rz(-2.0634212) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.1972034) q[1];
sx q[1];
rz(-1.7447339) q[1];
sx q[1];
rz(2.7157213) q[1];
rz(-pi) q[2];
x q[2];
rz(0.91032116) q[3];
sx q[3];
rz(-2.3757138) q[3];
sx q[3];
rz(-0.41492763) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.67132407) q[2];
sx q[2];
rz(-2.2094122) q[2];
sx q[2];
rz(2.0051125) q[2];
rz(1.350435) q[3];
sx q[3];
rz(-1.8108436) q[3];
sx q[3];
rz(2.7854846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7478624) q[0];
sx q[0];
rz(-2.6342454) q[0];
sx q[0];
rz(-0.90079975) q[0];
rz(-1.9081217) q[1];
sx q[1];
rz(-2.2869459) q[1];
sx q[1];
rz(2.5305117) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38310941) q[0];
sx q[0];
rz(-2.0462667) q[0];
sx q[0];
rz(-0.0080541797) q[0];
rz(-pi) q[1];
x q[1];
rz(2.99154) q[2];
sx q[2];
rz(-0.23153472) q[2];
sx q[2];
rz(-0.48047149) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.722695) q[1];
sx q[1];
rz(-1.9729821) q[1];
sx q[1];
rz(0.45940347) q[1];
x q[2];
rz(2.8318367) q[3];
sx q[3];
rz(-2.7742226) q[3];
sx q[3];
rz(-1.0760372) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.63371199) q[2];
sx q[2];
rz(-0.66358006) q[2];
sx q[2];
rz(-3.0207685) q[2];
rz(3.1145596) q[3];
sx q[3];
rz(-2.9398672) q[3];
sx q[3];
rz(2.2656608) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19288572) q[0];
sx q[0];
rz(-2.3159733) q[0];
sx q[0];
rz(-1.8008308) q[0];
rz(1.5277398) q[1];
sx q[1];
rz(-1.8283045) q[1];
sx q[1];
rz(0.54940474) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4260146) q[0];
sx q[0];
rz(-2.3142646) q[0];
sx q[0];
rz(1.6756776) q[0];
rz(-3.1099186) q[2];
sx q[2];
rz(-1.0397382) q[2];
sx q[2];
rz(-0.90716328) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.0548965) q[1];
sx q[1];
rz(-0.48829309) q[1];
sx q[1];
rz(1.6773305) q[1];
rz(0.033614393) q[3];
sx q[3];
rz(-1.1096802) q[3];
sx q[3];
rz(2.5534843) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.01650979) q[2];
sx q[2];
rz(-1.5089401) q[2];
sx q[2];
rz(-2.5588918) q[2];
rz(-1.6216283) q[3];
sx q[3];
rz(-0.71526066) q[3];
sx q[3];
rz(-1.8250072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80578605) q[0];
sx q[0];
rz(-1.0842706) q[0];
sx q[0];
rz(1.8871319) q[0];
rz(-1.1955903) q[1];
sx q[1];
rz(-1.7343727) q[1];
sx q[1];
rz(1.5171299) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82255581) q[0];
sx q[0];
rz(-1.5040888) q[0];
sx q[0];
rz(0.011724756) q[0];
x q[1];
rz(2.1468494) q[2];
sx q[2];
rz(-2.2855504) q[2];
sx q[2];
rz(-1.7896259) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.9395917) q[1];
sx q[1];
rz(-0.67145214) q[1];
sx q[1];
rz(0.47452773) q[1];
rz(-pi) q[2];
x q[2];
rz(0.51503564) q[3];
sx q[3];
rz(-1.5608369) q[3];
sx q[3];
rz(1.8731093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.17455165) q[2];
sx q[2];
rz(-1.5341362) q[2];
sx q[2];
rz(3.0957481) q[2];
rz(0.52715078) q[3];
sx q[3];
rz(-1.0322626) q[3];
sx q[3];
rz(-2.3003858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4884278) q[0];
sx q[0];
rz(-0.66075745) q[0];
sx q[0];
rz(-0.38594693) q[0];
rz(-1.7653607) q[1];
sx q[1];
rz(-1.0877345) q[1];
sx q[1];
rz(-1.2765346) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7273743) q[0];
sx q[0];
rz(-0.91394768) q[0];
sx q[0];
rz(-1.0947544) q[0];
x q[1];
rz(-1.6930503) q[2];
sx q[2];
rz(-2.7664693) q[2];
sx q[2];
rz(1.5135672) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.701653) q[1];
sx q[1];
rz(-1.5787485) q[1];
sx q[1];
rz(-1.6064744) q[1];
rz(-pi) q[2];
rz(-2.0790398) q[3];
sx q[3];
rz(-0.073598737) q[3];
sx q[3];
rz(1.910833) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.71629268) q[2];
sx q[2];
rz(-1.425681) q[2];
sx q[2];
rz(-1.3670134) q[2];
rz(-0.084271757) q[3];
sx q[3];
rz(-2.027498) q[3];
sx q[3];
rz(-3.058694) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0609584) q[0];
sx q[0];
rz(-3.0301889) q[0];
sx q[0];
rz(1.8633307) q[0];
rz(-3.0807965) q[1];
sx q[1];
rz(-2.1863329) q[1];
sx q[1];
rz(-1.0999058) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5583145) q[0];
sx q[0];
rz(-0.41039) q[0];
sx q[0];
rz(1.9463842) q[0];
rz(-pi) q[1];
rz(2.3387735) q[2];
sx q[2];
rz(-2.2051797) q[2];
sx q[2];
rz(-0.54976058) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.5775313) q[1];
sx q[1];
rz(-1.2857591) q[1];
sx q[1];
rz(-2.7708745) q[1];
x q[2];
rz(2.527929) q[3];
sx q[3];
rz(-1.1432768) q[3];
sx q[3];
rz(0.18117426) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.78201571) q[2];
sx q[2];
rz(-0.64543739) q[2];
sx q[2];
rz(2.6673356) q[2];
rz(-0.54840243) q[3];
sx q[3];
rz(-0.67265284) q[3];
sx q[3];
rz(-1.3801581) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7062374) q[0];
sx q[0];
rz(-2.7147003) q[0];
sx q[0];
rz(-1.2109582) q[0];
rz(-1.2779166) q[1];
sx q[1];
rz(-1.5958818) q[1];
sx q[1];
rz(0.011215297) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1007014) q[0];
sx q[0];
rz(-0.56849231) q[0];
sx q[0];
rz(-0.16162737) q[0];
rz(-pi) q[1];
rz(1.340788) q[2];
sx q[2];
rz(-1.9267779) q[2];
sx q[2];
rz(-1.9614416) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.6026759) q[1];
sx q[1];
rz(-2.3178812) q[1];
sx q[1];
rz(-3.0680821) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2809926) q[3];
sx q[3];
rz(-1.6609238) q[3];
sx q[3];
rz(0.33099701) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8212905) q[2];
sx q[2];
rz(-1.9052637) q[2];
sx q[2];
rz(-2.383929) q[2];
rz(2.441794) q[3];
sx q[3];
rz(-0.57707969) q[3];
sx q[3];
rz(2.2363766) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59117544) q[0];
sx q[0];
rz(-1.2400405) q[0];
sx q[0];
rz(-0.33429876) q[0];
rz(-2.7475157) q[1];
sx q[1];
rz(-0.95029345) q[1];
sx q[1];
rz(1.2324415) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.619433) q[0];
sx q[0];
rz(-1.30511) q[0];
sx q[0];
rz(2.9152318) q[0];
rz(-2.443529) q[2];
sx q[2];
rz(-1.4163829) q[2];
sx q[2];
rz(-2.7464641) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(3.0403172) q[1];
sx q[1];
rz(-2.7849232) q[1];
sx q[1];
rz(-0.89445313) q[1];
rz(-0.27421342) q[3];
sx q[3];
rz(-1.0875487) q[3];
sx q[3];
rz(-0.26457126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.95840994) q[2];
sx q[2];
rz(-0.56768688) q[2];
sx q[2];
rz(3.0779823) q[2];
rz(-2.375864) q[3];
sx q[3];
rz(-2.016341) q[3];
sx q[3];
rz(-0.061802797) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80339377) q[0];
sx q[0];
rz(-0.82427187) q[0];
sx q[0];
rz(-1.6765539) q[0];
rz(1.6784531) q[1];
sx q[1];
rz(-1.8747765) q[1];
sx q[1];
rz(2.3864559) q[1];
rz(1.7378308) q[2];
sx q[2];
rz(-0.93334953) q[2];
sx q[2];
rz(2.2095528) q[2];
rz(-0.83972558) q[3];
sx q[3];
rz(-1.6942377) q[3];
sx q[3];
rz(-0.89577976) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
