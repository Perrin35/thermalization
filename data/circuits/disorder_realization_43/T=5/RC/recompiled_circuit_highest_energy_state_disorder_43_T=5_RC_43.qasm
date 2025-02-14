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
rz(2.9146258) q[0];
sx q[0];
rz(-0.6114971) q[0];
sx q[0];
rz(0.55846941) q[0];
rz(0.59977579) q[1];
sx q[1];
rz(1.37473) q[1];
sx q[1];
rz(9.3461499) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0723018) q[0];
sx q[0];
rz(-0.8197533) q[0];
sx q[0];
rz(-0.79901521) q[0];
rz(-pi) q[1];
rz(-2.4628381) q[2];
sx q[2];
rz(-1.6683443) q[2];
sx q[2];
rz(0.3435979) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.47657644) q[1];
sx q[1];
rz(-1.652583) q[1];
sx q[1];
rz(0.25340124) q[1];
rz(1.2492248) q[3];
sx q[3];
rz(-1.235629) q[3];
sx q[3];
rz(2.9706819) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.2416396) q[2];
sx q[2];
rz(-3.0934379) q[2];
sx q[2];
rz(1.0345577) q[2];
rz(3.0311846) q[3];
sx q[3];
rz(-1.4805099) q[3];
sx q[3];
rz(2.835527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0050874) q[0];
sx q[0];
rz(-1.6965447) q[0];
sx q[0];
rz(1.9553631) q[0];
rz(-1.2677445) q[1];
sx q[1];
rz(-2.6823951) q[1];
sx q[1];
rz(1.3892106) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88436172) q[0];
sx q[0];
rz(-0.62998929) q[0];
sx q[0];
rz(1.0347741) q[0];
rz(-pi) q[1];
rz(-1.4591501) q[2];
sx q[2];
rz(-2.0497339) q[2];
sx q[2];
rz(-0.079733036) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.8858365) q[1];
sx q[1];
rz(-2.6955018) q[1];
sx q[1];
rz(0.85020937) q[1];
rz(-pi) q[2];
rz(3.0458353) q[3];
sx q[3];
rz(-2.1928722) q[3];
sx q[3];
rz(-2.2945537) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.17239751) q[2];
sx q[2];
rz(-1.1149422) q[2];
sx q[2];
rz(-1.7051075) q[2];
rz(-1.2596333) q[3];
sx q[3];
rz(-0.45537046) q[3];
sx q[3];
rz(3.1148541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
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
rz(2.7773975) q[0];
sx q[0];
rz(-1.3854249) q[0];
sx q[0];
rz(-1.9245603) q[0];
rz(-2.9265535) q[1];
sx q[1];
rz(-1.8937078) q[1];
sx q[1];
rz(-1.7020114) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0501409) q[0];
sx q[0];
rz(-0.9147075) q[0];
sx q[0];
rz(1.7277645) q[0];
rz(-1.7363988) q[2];
sx q[2];
rz(-0.51711997) q[2];
sx q[2];
rz(-1.0696326) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.176217) q[1];
sx q[1];
rz(-1.4439519) q[1];
sx q[1];
rz(-1.8017669) q[1];
x q[2];
rz(1.2974627) q[3];
sx q[3];
rz(-1.36051) q[3];
sx q[3];
rz(-3.1137636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.8000468) q[2];
sx q[2];
rz(-1.505082) q[2];
sx q[2];
rz(-0.81614196) q[2];
rz(-3.1331565) q[3];
sx q[3];
rz(-2.4237027) q[3];
sx q[3];
rz(1.5661904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4207882) q[0];
sx q[0];
rz(-1.1710465) q[0];
sx q[0];
rz(2.8724443) q[0];
rz(1.7587657) q[1];
sx q[1];
rz(-0.56126422) q[1];
sx q[1];
rz(-0.47163481) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4260226) q[0];
sx q[0];
rz(-1.7105118) q[0];
sx q[0];
rz(-2.2989397) q[0];
rz(0.083375562) q[2];
sx q[2];
rz(-2.7365587) q[2];
sx q[2];
rz(-2.0138559) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.1822575) q[1];
sx q[1];
rz(-1.8885837) q[1];
sx q[1];
rz(-2.5885568) q[1];
rz(-2.7876623) q[3];
sx q[3];
rz(-2.2056863) q[3];
sx q[3];
rz(-0.77087444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.65583324) q[2];
sx q[2];
rz(-1.7344319) q[2];
sx q[2];
rz(1.887623) q[2];
rz(0.29821011) q[3];
sx q[3];
rz(-1.8658274) q[3];
sx q[3];
rz(1.6037174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0373847) q[0];
sx q[0];
rz(-2.2315114) q[0];
sx q[0];
rz(-0.75697672) q[0];
rz(2.0825999) q[1];
sx q[1];
rz(-1.6964361) q[1];
sx q[1];
rz(-2.7991378) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.379102) q[0];
sx q[0];
rz(-0.44371334) q[0];
sx q[0];
rz(-2.0279473) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.191869) q[2];
sx q[2];
rz(-2.4034326) q[2];
sx q[2];
rz(-1.9110695) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.96794242) q[1];
sx q[1];
rz(-1.0999803) q[1];
sx q[1];
rz(1.1556666) q[1];
rz(-0.41564055) q[3];
sx q[3];
rz(-1.2977674) q[3];
sx q[3];
rz(1.331783) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.69514889) q[2];
sx q[2];
rz(-1.6764287) q[2];
sx q[2];
rz(2.5214419) q[2];
rz(-2.5946963) q[3];
sx q[3];
rz(-0.20709012) q[3];
sx q[3];
rz(-1.0774405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81724375) q[0];
sx q[0];
rz(-0.20340782) q[0];
sx q[0];
rz(-1.2212344) q[0];
rz(1.1034032) q[1];
sx q[1];
rz(-1.1627448) q[1];
sx q[1];
rz(2.3308636) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65064036) q[0];
sx q[0];
rz(-1.3155259) q[0];
sx q[0];
rz(-1.7351236) q[0];
rz(-2.2037494) q[2];
sx q[2];
rz(-1.019291) q[2];
sx q[2];
rz(-2.1126975) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.3103957) q[1];
sx q[1];
rz(-1.7094905) q[1];
sx q[1];
rz(3.0866361) q[1];
x q[2];
rz(1.9493413) q[3];
sx q[3];
rz(-1.8486128) q[3];
sx q[3];
rz(-1.5072106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.80060426) q[2];
sx q[2];
rz(-1.1977414) q[2];
sx q[2];
rz(1.1654589) q[2];
rz(0.10945877) q[3];
sx q[3];
rz(-2.1200924) q[3];
sx q[3];
rz(-0.060733184) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3133746) q[0];
sx q[0];
rz(-3.1309541) q[0];
sx q[0];
rz(-2.4995372) q[0];
rz(-2.8984046) q[1];
sx q[1];
rz(-0.41047341) q[1];
sx q[1];
rz(2.5004255) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0579906) q[0];
sx q[0];
rz(-0.80785368) q[0];
sx q[0];
rz(-2.3919748) q[0];
rz(2.7566822) q[2];
sx q[2];
rz(-0.78854568) q[2];
sx q[2];
rz(-0.53764082) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.6536246) q[1];
sx q[1];
rz(-1.2418126) q[1];
sx q[1];
rz(2.3423561) q[1];
rz(-pi) q[2];
rz(-2.083046) q[3];
sx q[3];
rz(-0.44971684) q[3];
sx q[3];
rz(-1.095677) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.0082561) q[2];
sx q[2];
rz(-1.9828321) q[2];
sx q[2];
rz(0.94375098) q[2];
rz(-2.4933695) q[3];
sx q[3];
rz(-0.74130487) q[3];
sx q[3];
rz(2.4751723) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42565313) q[0];
sx q[0];
rz(-2.3305927) q[0];
sx q[0];
rz(2.7124523) q[0];
rz(2.7375713) q[1];
sx q[1];
rz(-2.7262913) q[1];
sx q[1];
rz(-2.2228352) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0149834) q[0];
sx q[0];
rz(-1.0246687) q[0];
sx q[0];
rz(1.855143) q[0];
rz(2.0112923) q[2];
sx q[2];
rz(-0.28708946) q[2];
sx q[2];
rz(2.8237791) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.97087103) q[1];
sx q[1];
rz(-1.3525463) q[1];
sx q[1];
rz(1.8884601) q[1];
x q[2];
rz(-1.346896) q[3];
sx q[3];
rz(-2.8200275) q[3];
sx q[3];
rz(-2.0617776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.018844133) q[2];
sx q[2];
rz(-0.65905535) q[2];
sx q[2];
rz(-1.7913691) q[2];
rz(-2.8090779) q[3];
sx q[3];
rz(-0.88034383) q[3];
sx q[3];
rz(1.3688709) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3933082) q[0];
sx q[0];
rz(-2.4916593) q[0];
sx q[0];
rz(2.6848324) q[0];
rz(-1.4211753) q[1];
sx q[1];
rz(-1.8616118) q[1];
sx q[1];
rz(0.054904003) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0739912) q[0];
sx q[0];
rz(-1.5419863) q[0];
sx q[0];
rz(0.83782283) q[0];
x q[1];
rz(-0.41150062) q[2];
sx q[2];
rz(-1.9691281) q[2];
sx q[2];
rz(2.8702132) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.0753638) q[1];
sx q[1];
rz(-1.3540596) q[1];
sx q[1];
rz(-2.0677057) q[1];
rz(-2.6087481) q[3];
sx q[3];
rz(-2.5379764) q[3];
sx q[3];
rz(-0.6288022) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.2533337) q[2];
sx q[2];
rz(-0.71037018) q[2];
sx q[2];
rz(2.9448275) q[2];
rz(2.0678068) q[3];
sx q[3];
rz(-2.00627) q[3];
sx q[3];
rz(-0.73608583) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1573023) q[0];
sx q[0];
rz(-0.8534011) q[0];
sx q[0];
rz(2.2421457) q[0];
rz(2.8304214) q[1];
sx q[1];
rz(-1.2093465) q[1];
sx q[1];
rz(-1.4003632) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3547826) q[0];
sx q[0];
rz(-1.3616832) q[0];
sx q[0];
rz(2.1378921) q[0];
rz(-pi) q[1];
rz(2.0092404) q[2];
sx q[2];
rz(-0.52784) q[2];
sx q[2];
rz(2.6486625) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.33161641) q[1];
sx q[1];
rz(-0.73996937) q[1];
sx q[1];
rz(-0.59438594) q[1];
x q[2];
rz(-2.4254341) q[3];
sx q[3];
rz(-2.0422438) q[3];
sx q[3];
rz(0.52385073) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.060999) q[2];
sx q[2];
rz(-0.62005764) q[2];
sx q[2];
rz(1.9798123) q[2];
rz(2.0628085) q[3];
sx q[3];
rz(-0.99083841) q[3];
sx q[3];
rz(-2.259528) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5921191) q[0];
sx q[0];
rz(-1.5220806) q[0];
sx q[0];
rz(1.4694389) q[0];
rz(0.1334162) q[1];
sx q[1];
rz(-1.7620371) q[1];
sx q[1];
rz(-0.98099991) q[1];
rz(-1.5173923) q[2];
sx q[2];
rz(-1.5718062) q[2];
sx q[2];
rz(0.013769416) q[2];
rz(-1.6422938) q[3];
sx q[3];
rz(-1.0701978) q[3];
sx q[3];
rz(1.3842907) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
