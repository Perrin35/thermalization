OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.47473946) q[0];
sx q[0];
rz(-0.82959509) q[0];
sx q[0];
rz(0.15396804) q[0];
rz(0.83377588) q[1];
sx q[1];
rz(4.1339388) q[1];
sx q[1];
rz(9.0864656) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51988039) q[0];
sx q[0];
rz(-1.3268688) q[0];
sx q[0];
rz(1.2432616) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1233653) q[2];
sx q[2];
rz(-0.72927232) q[2];
sx q[2];
rz(2.4016618) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.073804341) q[1];
sx q[1];
rz(-1.6238302) q[1];
sx q[1];
rz(2.8552613) q[1];
x q[2];
rz(-2.0831574) q[3];
sx q[3];
rz(-1.5844987) q[3];
sx q[3];
rz(1.1289568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.14264318) q[2];
sx q[2];
rz(-2.8012186) q[2];
sx q[2];
rz(1.1738698) q[2];
rz(3.0657892) q[3];
sx q[3];
rz(-1.1444164) q[3];
sx q[3];
rz(3.048786) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8006111) q[0];
sx q[0];
rz(-2.0759463) q[0];
sx q[0];
rz(-3.0766292) q[0];
rz(2.5669572) q[1];
sx q[1];
rz(-2.7119633) q[1];
sx q[1];
rz(1.2423135) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3949462) q[0];
sx q[0];
rz(-1.442369) q[0];
sx q[0];
rz(-2.8917679) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9827036) q[2];
sx q[2];
rz(-0.94413589) q[2];
sx q[2];
rz(1.6275258) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.2989267) q[1];
sx q[1];
rz(-2.0999523) q[1];
sx q[1];
rz(1.3859205) q[1];
rz(-pi) q[2];
rz(-0.70988016) q[3];
sx q[3];
rz(-1.1871561) q[3];
sx q[3];
rz(-0.82304728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.80766455) q[2];
sx q[2];
rz(-2.0662722) q[2];
sx q[2];
rz(-2.5578965) q[2];
rz(0.57404533) q[3];
sx q[3];
rz(-1.1255001) q[3];
sx q[3];
rz(3.0103502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
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
rz(2.7211001) q[0];
sx q[0];
rz(-2.2476966) q[0];
sx q[0];
rz(0.72845355) q[0];
rz(-1.4942253) q[1];
sx q[1];
rz(-0.39847001) q[1];
sx q[1];
rz(-1.0167936) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3166312) q[0];
sx q[0];
rz(-1.5350071) q[0];
sx q[0];
rz(-3.0599942) q[0];
rz(-pi) q[1];
x q[1];
rz(0.74480199) q[2];
sx q[2];
rz(-0.6664657) q[2];
sx q[2];
rz(1.1643861) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.180417) q[1];
sx q[1];
rz(-2.3824586) q[1];
sx q[1];
rz(1.7708066) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2891644) q[3];
sx q[3];
rz(-1.6179807) q[3];
sx q[3];
rz(-0.55164528) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8016522) q[2];
sx q[2];
rz(-1.6093971) q[2];
sx q[2];
rz(2.0920848) q[2];
rz(-2.5028051) q[3];
sx q[3];
rz(-0.63126606) q[3];
sx q[3];
rz(-1.1857741) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3291572) q[0];
sx q[0];
rz(-1.2673459) q[0];
sx q[0];
rz(-1.4720434) q[0];
rz(-2.4064348) q[1];
sx q[1];
rz(-0.77886326) q[1];
sx q[1];
rz(-0.24681117) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1082552) q[0];
sx q[0];
rz(-2.4318683) q[0];
sx q[0];
rz(1.8002585) q[0];
x q[1];
rz(-0.47841448) q[2];
sx q[2];
rz(-1.7583499) q[2];
sx q[2];
rz(1.071655) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.81891638) q[1];
sx q[1];
rz(-2.7671742) q[1];
sx q[1];
rz(-2.9973292) q[1];
rz(1.1770583) q[3];
sx q[3];
rz(-0.65440946) q[3];
sx q[3];
rz(-2.8670093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.6639158) q[2];
sx q[2];
rz(-2.0229979) q[2];
sx q[2];
rz(-1.5412615) q[2];
rz(2.4345496) q[3];
sx q[3];
rz(-1.0975081) q[3];
sx q[3];
rz(0.88821205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8108869) q[0];
sx q[0];
rz(-0.72074497) q[0];
sx q[0];
rz(1.3274308) q[0];
rz(-1.5785626) q[1];
sx q[1];
rz(-2.6674318) q[1];
sx q[1];
rz(2.8932103) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36149704) q[0];
sx q[0];
rz(-2.1108315) q[0];
sx q[0];
rz(1.5310775) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7288754) q[2];
sx q[2];
rz(-1.386596) q[2];
sx q[2];
rz(-0.36518156) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.3779113) q[1];
sx q[1];
rz(-1.3994201) q[1];
sx q[1];
rz(0.59314368) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8304843) q[3];
sx q[3];
rz(-0.6725544) q[3];
sx q[3];
rz(-2.6274519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.7053232) q[2];
sx q[2];
rz(-2.1537809) q[2];
sx q[2];
rz(-0.79745897) q[2];
rz(-2.752839) q[3];
sx q[3];
rz(-0.60384408) q[3];
sx q[3];
rz(2.6388772) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1335063) q[0];
sx q[0];
rz(-3.0651423) q[0];
sx q[0];
rz(1.3457993) q[0];
rz(1.0812409) q[1];
sx q[1];
rz(-1.9045647) q[1];
sx q[1];
rz(3.016901) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7170982) q[0];
sx q[0];
rz(-1.4846804) q[0];
sx q[0];
rz(-0.17688246) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7331946) q[2];
sx q[2];
rz(-0.64986594) q[2];
sx q[2];
rz(-1.0330531) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.2248968) q[1];
sx q[1];
rz(-1.6273013) q[1];
sx q[1];
rz(-2.4042261) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3521306) q[3];
sx q[3];
rz(-1.339774) q[3];
sx q[3];
rz(-2.4276589) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.6283915) q[2];
sx q[2];
rz(-1.1867563) q[2];
sx q[2];
rz(-2.52264) q[2];
rz(1.0533054) q[3];
sx q[3];
rz(-2.9635933) q[3];
sx q[3];
rz(0.9296023) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58182794) q[0];
sx q[0];
rz(-1.3904089) q[0];
sx q[0];
rz(-2.0986309) q[0];
rz(2.6783121) q[1];
sx q[1];
rz(-2.0279341) q[1];
sx q[1];
rz(-1.0707062) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2830414) q[0];
sx q[0];
rz(-1.0867449) q[0];
sx q[0];
rz(-1.5714684) q[0];
rz(-pi) q[1];
rz(-1.4378909) q[2];
sx q[2];
rz(-1.226236) q[2];
sx q[2];
rz(0.29495707) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.067498265) q[1];
sx q[1];
rz(-2.3405582) q[1];
sx q[1];
rz(1.7098411) q[1];
rz(0.058810874) q[3];
sx q[3];
rz(-0.33274129) q[3];
sx q[3];
rz(-2.9174093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.7730007) q[2];
sx q[2];
rz(-0.7395145) q[2];
sx q[2];
rz(-2.8179742) q[2];
rz(-0.98179022) q[3];
sx q[3];
rz(-0.86172813) q[3];
sx q[3];
rz(-2.0543082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48802808) q[0];
sx q[0];
rz(-1.4166778) q[0];
sx q[0];
rz(1.9352242) q[0];
rz(-1.9288829) q[1];
sx q[1];
rz(-2.2884463) q[1];
sx q[1];
rz(2.1941197) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3377209) q[0];
sx q[0];
rz(-0.54715711) q[0];
sx q[0];
rz(-1.7659811) q[0];
rz(0.73080365) q[2];
sx q[2];
rz(-2.9060504) q[2];
sx q[2];
rz(1.8159602) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.265043) q[1];
sx q[1];
rz(-1.3339086) q[1];
sx q[1];
rz(-0.9898647) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.34668563) q[3];
sx q[3];
rz(-0.83804916) q[3];
sx q[3];
rz(-2.8447062) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5802713) q[2];
sx q[2];
rz(-1.7611971) q[2];
sx q[2];
rz(0.46978152) q[2];
rz(1.3011159) q[3];
sx q[3];
rz(-1.7062635) q[3];
sx q[3];
rz(2.8619213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
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
rz(0.25205055) q[0];
sx q[0];
rz(-2.7520576) q[0];
sx q[0];
rz(-1.8126194) q[0];
rz(0.7912311) q[1];
sx q[1];
rz(-2.8104517) q[1];
sx q[1];
rz(0.20283094) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3664353) q[0];
sx q[0];
rz(-0.046910722) q[0];
sx q[0];
rz(-0.58880083) q[0];
rz(-pi) q[1];
rz(1.780026) q[2];
sx q[2];
rz(-1.7893357) q[2];
sx q[2];
rz(-2.2184559) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.054246) q[1];
sx q[1];
rz(-3.0294703) q[1];
sx q[1];
rz(-1.1152769) q[1];
rz(-pi) q[2];
rz(-1.4167452) q[3];
sx q[3];
rz(-1.1136347) q[3];
sx q[3];
rz(-2.9726213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7245076) q[2];
sx q[2];
rz(-2.8736726) q[2];
sx q[2];
rz(0.99651304) q[2];
rz(2.7881682) q[3];
sx q[3];
rz(-0.74520183) q[3];
sx q[3];
rz(-2.3341808) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0614232) q[0];
sx q[0];
rz(-0.81289476) q[0];
sx q[0];
rz(0.18173519) q[0];
rz(0.043047992) q[1];
sx q[1];
rz(-2.4964066) q[1];
sx q[1];
rz(-2.8607686) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8107287) q[0];
sx q[0];
rz(-2.4952336) q[0];
sx q[0];
rz(-2.3848563) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7634723) q[2];
sx q[2];
rz(-1.780605) q[2];
sx q[2];
rz(-3.0378621) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.9210789) q[1];
sx q[1];
rz(-1.6069357) q[1];
sx q[1];
rz(1.5131348) q[1];
x q[2];
rz(2.9452768) q[3];
sx q[3];
rz(-1.2442949) q[3];
sx q[3];
rz(-1.452009) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3165555) q[2];
sx q[2];
rz(-1.8871769) q[2];
sx q[2];
rz(-2.5184856) q[2];
rz(1.0021707) q[3];
sx q[3];
rz(-1.3391756) q[3];
sx q[3];
rz(-0.56308693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4939209) q[0];
sx q[0];
rz(-1.5681842) q[0];
sx q[0];
rz(1.6012123) q[0];
rz(2.2676246) q[1];
sx q[1];
rz(-2.0762434) q[1];
sx q[1];
rz(0.11003065) q[1];
rz(0.80550823) q[2];
sx q[2];
rz(-2.0279573) q[2];
sx q[2];
rz(1.3062994) q[2];
rz(0.86250967) q[3];
sx q[3];
rz(-2.6211092) q[3];
sx q[3];
rz(0.72343788) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
