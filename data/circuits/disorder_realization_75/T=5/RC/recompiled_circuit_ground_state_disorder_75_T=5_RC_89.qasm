OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.6899941) q[0];
sx q[0];
rz(-2.8328083) q[0];
sx q[0];
rz(0.2398332) q[0];
rz(2.580515) q[1];
sx q[1];
rz(-0.74220389) q[1];
sx q[1];
rz(1.5129369) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6282677) q[0];
sx q[0];
rz(-1.4428992) q[0];
sx q[0];
rz(1.2012175) q[0];
rz(-pi) q[1];
rz(2.0432908) q[2];
sx q[2];
rz(-1.3368946) q[2];
sx q[2];
rz(0.59532524) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.4862772) q[1];
sx q[1];
rz(-1.4239171) q[1];
sx q[1];
rz(-0.49238654) q[1];
rz(-pi) q[2];
rz(-1.9001107) q[3];
sx q[3];
rz(-2.9639177) q[3];
sx q[3];
rz(-0.8575646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.1538887) q[2];
sx q[2];
rz(-0.99606267) q[2];
sx q[2];
rz(1.055701) q[2];
rz(-0.40134564) q[3];
sx q[3];
rz(-1.515712) q[3];
sx q[3];
rz(-1.5041941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6317247) q[0];
sx q[0];
rz(-2.8724176) q[0];
sx q[0];
rz(-1.7357695) q[0];
rz(-0.74613219) q[1];
sx q[1];
rz(-1.1765307) q[1];
sx q[1];
rz(3.0768118) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9822183) q[0];
sx q[0];
rz(-1.1752593) q[0];
sx q[0];
rz(-0.48898029) q[0];
x q[1];
rz(0.6798052) q[2];
sx q[2];
rz(-1.7035988) q[2];
sx q[2];
rz(-0.50253579) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.8967817) q[1];
sx q[1];
rz(-0.1644539) q[1];
sx q[1];
rz(-1.0599778) q[1];
rz(3.1094743) q[3];
sx q[3];
rz(-1.7815551) q[3];
sx q[3];
rz(-3.075066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.88227162) q[2];
sx q[2];
rz(-0.52053014) q[2];
sx q[2];
rz(-0.041291324) q[2];
rz(1.1193554) q[3];
sx q[3];
rz(-1.7740403) q[3];
sx q[3];
rz(1.1411427) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(1.816788) q[0];
sx q[0];
rz(-2.7140706) q[0];
sx q[0];
rz(0.69931716) q[0];
rz(-2.518867) q[1];
sx q[1];
rz(-2.3284349) q[1];
sx q[1];
rz(-2.8831388) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9759443) q[0];
sx q[0];
rz(-1.5830402) q[0];
sx q[0];
rz(-2.5077639) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0270789) q[2];
sx q[2];
rz(-2.2199759) q[2];
sx q[2];
rz(0.69930062) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.041641673) q[1];
sx q[1];
rz(-1.1219684) q[1];
sx q[1];
rz(-1.5931409) q[1];
rz(-3.1051226) q[3];
sx q[3];
rz(-1.1410645) q[3];
sx q[3];
rz(3.1231073) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.2059325) q[2];
sx q[2];
rz(-1.460133) q[2];
sx q[2];
rz(-2.4350731) q[2];
rz(-0.62890729) q[3];
sx q[3];
rz(-0.8258515) q[3];
sx q[3];
rz(-2.0188873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.127447) q[0];
sx q[0];
rz(-1.4268459) q[0];
sx q[0];
rz(2.1674147) q[0];
rz(1.6575419) q[1];
sx q[1];
rz(-2.0482792) q[1];
sx q[1];
rz(1.0008224) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38353048) q[0];
sx q[0];
rz(-0.62072004) q[0];
sx q[0];
rz(1.8412983) q[0];
x q[1];
rz(2.3685969) q[2];
sx q[2];
rz(-2.2064035) q[2];
sx q[2];
rz(-2.5781812) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.9467869) q[1];
sx q[1];
rz(-1.1071702) q[1];
sx q[1];
rz(-0.63321106) q[1];
x q[2];
rz(0.14814143) q[3];
sx q[3];
rz(-0.32399789) q[3];
sx q[3];
rz(-1.2202386) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.68690825) q[2];
sx q[2];
rz(-1.3409706) q[2];
sx q[2];
rz(-2.3243813) q[2];
rz(0.59598437) q[3];
sx q[3];
rz(-1.7242566) q[3];
sx q[3];
rz(1.3982841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
rz(-0.78438321) q[0];
sx q[0];
rz(-1.7359808) q[0];
sx q[0];
rz(1.0719517) q[0];
rz(1.1373854) q[1];
sx q[1];
rz(-1.5147361) q[1];
sx q[1];
rz(0.17328182) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.539285) q[0];
sx q[0];
rz(-2.0208911) q[0];
sx q[0];
rz(2.0099239) q[0];
rz(-pi) q[1];
rz(2.2229175) q[2];
sx q[2];
rz(-2.5498418) q[2];
sx q[2];
rz(0.39338912) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.73663354) q[1];
sx q[1];
rz(-1.4990605) q[1];
sx q[1];
rz(-0.81291764) q[1];
rz(-pi) q[2];
rz(1.6638443) q[3];
sx q[3];
rz(-1.4085839) q[3];
sx q[3];
rz(-1.4402267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.5002354) q[2];
sx q[2];
rz(-2.6884029) q[2];
sx q[2];
rz(2.267061) q[2];
rz(-2.4747961) q[3];
sx q[3];
rz(-1.6801445) q[3];
sx q[3];
rz(2.3165406) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85103971) q[0];
sx q[0];
rz(-0.14851004) q[0];
sx q[0];
rz(-0.17459757) q[0];
rz(-2.5945276) q[1];
sx q[1];
rz(-0.51536307) q[1];
sx q[1];
rz(1.863716) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96748501) q[0];
sx q[0];
rz(-0.293403) q[0];
sx q[0];
rz(1.3551329) q[0];
rz(-2.3894887) q[2];
sx q[2];
rz(-1.9255203) q[2];
sx q[2];
rz(0.51774509) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.3681952) q[1];
sx q[1];
rz(-0.85933472) q[1];
sx q[1];
rz(1.0991286) q[1];
rz(0.90323351) q[3];
sx q[3];
rz(-2.0732911) q[3];
sx q[3];
rz(-2.5046405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.7586907) q[2];
sx q[2];
rz(-0.26472696) q[2];
sx q[2];
rz(-0.36941377) q[2];
rz(-0.82768011) q[3];
sx q[3];
rz(-1.8281432) q[3];
sx q[3];
rz(0.15414342) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5675885) q[0];
sx q[0];
rz(-2.2249157) q[0];
sx q[0];
rz(2.9045203) q[0];
rz(1.392662) q[1];
sx q[1];
rz(-1.9475513) q[1];
sx q[1];
rz(-2.1377835) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.061351731) q[0];
sx q[0];
rz(-2.6727242) q[0];
sx q[0];
rz(2.7659225) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6200203) q[2];
sx q[2];
rz(-0.78465377) q[2];
sx q[2];
rz(1.705738) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.3442093) q[1];
sx q[1];
rz(-0.64192574) q[1];
sx q[1];
rz(2.5688237) q[1];
rz(-pi) q[2];
rz(-2.5828894) q[3];
sx q[3];
rz(-1.9877583) q[3];
sx q[3];
rz(2.1897763) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6134593) q[2];
sx q[2];
rz(-1.5789092) q[2];
sx q[2];
rz(-0.51631874) q[2];
rz(0.83827072) q[3];
sx q[3];
rz(-1.6603989) q[3];
sx q[3];
rz(2.9576438) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8539921) q[0];
sx q[0];
rz(-0.28165278) q[0];
sx q[0];
rz(-0.91424346) q[0];
rz(2.5573348) q[1];
sx q[1];
rz(-1.7353568) q[1];
sx q[1];
rz(-2.2343238) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0153025) q[0];
sx q[0];
rz(-1.5526104) q[0];
sx q[0];
rz(-0.92211266) q[0];
rz(-pi) q[1];
rz(-0.24326218) q[2];
sx q[2];
rz(-0.62957669) q[2];
sx q[2];
rz(-2.7223029) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-3.0338194) q[1];
sx q[1];
rz(-0.93074742) q[1];
sx q[1];
rz(-0.99730753) q[1];
rz(-0.97213718) q[3];
sx q[3];
rz(-1.3479509) q[3];
sx q[3];
rz(3.0226662) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.82627901) q[2];
sx q[2];
rz(-1.7118914) q[2];
sx q[2];
rz(-0.86177525) q[2];
rz(-1.2365384) q[3];
sx q[3];
rz(-3.0573513) q[3];
sx q[3];
rz(-1.2110565) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5009907) q[0];
sx q[0];
rz(-0.7826829) q[0];
sx q[0];
rz(-1.030141) q[0];
rz(2.8252699) q[1];
sx q[1];
rz(-1.7681237) q[1];
sx q[1];
rz(0.28824678) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9194946) q[0];
sx q[0];
rz(-1.7159425) q[0];
sx q[0];
rz(-2.7371489) q[0];
rz(-pi) q[1];
rz(-0.64916237) q[2];
sx q[2];
rz(-1.1466951) q[2];
sx q[2];
rz(0.28605385) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.0287778) q[1];
sx q[1];
rz(-0.76422404) q[1];
sx q[1];
rz(-0.9420331) q[1];
rz(1.3449617) q[3];
sx q[3];
rz(-2.2299181) q[3];
sx q[3];
rz(0.82443217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2890702) q[2];
sx q[2];
rz(-1.4129637) q[2];
sx q[2];
rz(0.22656013) q[2];
rz(-1.4551) q[3];
sx q[3];
rz(-0.89613599) q[3];
sx q[3];
rz(-2.4494825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9914472) q[0];
sx q[0];
rz(-1.8658072) q[0];
sx q[0];
rz(2.5164497) q[0];
rz(1.217968) q[1];
sx q[1];
rz(-2.4324799) q[1];
sx q[1];
rz(-1.9295173) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8344515) q[0];
sx q[0];
rz(-1.092415) q[0];
sx q[0];
rz(-2.4605453) q[0];
x q[1];
rz(-1.3122968) q[2];
sx q[2];
rz(-2.2257559) q[2];
sx q[2];
rz(1.0884681) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.068262488) q[1];
sx q[1];
rz(-2.0113328) q[1];
sx q[1];
rz(3.0777626) q[1];
x q[2];
rz(2.3732568) q[3];
sx q[3];
rz(-2.7339122) q[3];
sx q[3];
rz(-0.74484315) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.7071699) q[2];
sx q[2];
rz(-1.7662798) q[2];
sx q[2];
rz(2.0909069) q[2];
rz(-2.5896942) q[3];
sx q[3];
rz(-0.69010186) q[3];
sx q[3];
rz(1.2845854) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62591775) q[0];
sx q[0];
rz(-2.3139625) q[0];
sx q[0];
rz(2.3232842) q[0];
rz(0.7863518) q[1];
sx q[1];
rz(-2.4813589) q[1];
sx q[1];
rz(-2.9207041) q[1];
rz(3.1037504) q[2];
sx q[2];
rz(-1.2076245) q[2];
sx q[2];
rz(-2.5943499) q[2];
rz(0.10765392) q[3];
sx q[3];
rz(-0.68895491) q[3];
sx q[3];
rz(0.72910492) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
