OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.83710837) q[0];
sx q[0];
rz(4.8298782) q[0];
sx q[0];
rz(9.7363135) q[0];
rz(2.7040634) q[1];
sx q[1];
rz(4.9649927) q[1];
sx q[1];
rz(8.8658219) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11249609) q[0];
sx q[0];
rz(-1.7100428) q[0];
sx q[0];
rz(-1.990156) q[0];
rz(3.121435) q[2];
sx q[2];
rz(-2.5841789) q[2];
sx q[2];
rz(-1.7928979) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.47145876) q[1];
sx q[1];
rz(-2.11073) q[1];
sx q[1];
rz(0.88175168) q[1];
x q[2];
rz(-1.8564838) q[3];
sx q[3];
rz(-1.1477071) q[3];
sx q[3];
rz(0.74564122) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3922334) q[2];
sx q[2];
rz(-1.8584676) q[2];
sx q[2];
rz(-2.5048845) q[2];
rz(-0.84896815) q[3];
sx q[3];
rz(-2.520112) q[3];
sx q[3];
rz(0.23392114) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7648776) q[0];
sx q[0];
rz(-0.24704084) q[0];
sx q[0];
rz(-2.9887181) q[0];
rz(-0.75694594) q[1];
sx q[1];
rz(-1.5544954) q[1];
sx q[1];
rz(2.1551932) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61030412) q[0];
sx q[0];
rz(-3.0695519) q[0];
sx q[0];
rz(2.5487367) q[0];
rz(-pi) q[1];
rz(2.5949391) q[2];
sx q[2];
rz(-0.84178998) q[2];
sx q[2];
rz(0.58568776) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.0609329) q[1];
sx q[1];
rz(-0.45127171) q[1];
sx q[1];
rz(0.56834759) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8996703) q[3];
sx q[3];
rz(-1.4423014) q[3];
sx q[3];
rz(-1.0458898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.5793005) q[2];
sx q[2];
rz(-1.9220756) q[2];
sx q[2];
rz(0.78312773) q[2];
rz(-3.1230208) q[3];
sx q[3];
rz(-1.637807) q[3];
sx q[3];
rz(-0.40772453) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(-1.0531533) q[0];
sx q[0];
rz(-2.8494371) q[0];
sx q[0];
rz(0.96167481) q[0];
rz(2.7812474) q[1];
sx q[1];
rz(-1.1018437) q[1];
sx q[1];
rz(0.12869421) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0592244) q[0];
sx q[0];
rz(-1.6186065) q[0];
sx q[0];
rz(-1.6533018) q[0];
rz(0.80758904) q[2];
sx q[2];
rz(-1.979504) q[2];
sx q[2];
rz(-1.7847716) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.2663181) q[1];
sx q[1];
rz(-2.6058063) q[1];
sx q[1];
rz(2.5712625) q[1];
rz(-pi) q[2];
x q[2];
rz(0.13462984) q[3];
sx q[3];
rz(-0.80727808) q[3];
sx q[3];
rz(1.5546297) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0814357) q[2];
sx q[2];
rz(-1.9443941) q[2];
sx q[2];
rz(1.241768) q[2];
rz(2.5545819) q[3];
sx q[3];
rz(-0.95654064) q[3];
sx q[3];
rz(-0.97755066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
rz(-0.63242763) q[0];
sx q[0];
rz(-2.2594663) q[0];
sx q[0];
rz(1.084491) q[0];
rz(-1.4831316) q[1];
sx q[1];
rz(-0.56743923) q[1];
sx q[1];
rz(0.09253563) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4400892) q[0];
sx q[0];
rz(-1.8532231) q[0];
sx q[0];
rz(0.28042067) q[0];
rz(-0.96524694) q[2];
sx q[2];
rz(-0.56376981) q[2];
sx q[2];
rz(0.077686003) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.9738237) q[1];
sx q[1];
rz(-1.2281706) q[1];
sx q[1];
rz(-1.509036) q[1];
rz(-2.9805698) q[3];
sx q[3];
rz(-2.2194214) q[3];
sx q[3];
rz(-0.30740689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.3080421) q[2];
sx q[2];
rz(-1.4414859) q[2];
sx q[2];
rz(0.33205024) q[2];
rz(1.0559233) q[3];
sx q[3];
rz(-0.27763405) q[3];
sx q[3];
rz(-0.61029303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32245359) q[0];
sx q[0];
rz(-1.8503014) q[0];
sx q[0];
rz(2.9300368) q[0];
rz(-1.8353204) q[1];
sx q[1];
rz(-1.2439367) q[1];
sx q[1];
rz(2.4938915) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54007285) q[0];
sx q[0];
rz(-1.7300907) q[0];
sx q[0];
rz(-1.457731) q[0];
x q[1];
rz(-1.4164045) q[2];
sx q[2];
rz(-1.2263745) q[2];
sx q[2];
rz(1.4167348) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9808637) q[1];
sx q[1];
rz(-0.54667066) q[1];
sx q[1];
rz(-1.1854118) q[1];
rz(-pi) q[2];
rz(3.1224498) q[3];
sx q[3];
rz(-1.5557516) q[3];
sx q[3];
rz(-2.6263833) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.5806879) q[2];
sx q[2];
rz(-0.40955341) q[2];
sx q[2];
rz(-0.69331759) q[2];
rz(-0.66926113) q[3];
sx q[3];
rz(-1.4368493) q[3];
sx q[3];
rz(-0.0049237331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82419056) q[0];
sx q[0];
rz(-0.22739246) q[0];
sx q[0];
rz(1.90907) q[0];
rz(-2.0690074) q[1];
sx q[1];
rz(-1.0718081) q[1];
sx q[1];
rz(0.17428621) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6264544) q[0];
sx q[0];
rz(-2.2720552) q[0];
sx q[0];
rz(-2.7324972) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2767378) q[2];
sx q[2];
rz(-1.0360498) q[2];
sx q[2];
rz(0.45644444) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.78813997) q[1];
sx q[1];
rz(-0.8359682) q[1];
sx q[1];
rz(1.047903) q[1];
rz(-pi) q[2];
rz(-3.0252881) q[3];
sx q[3];
rz(-1.8810086) q[3];
sx q[3];
rz(-2.5935964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.8217414) q[2];
sx q[2];
rz(-0.80703002) q[2];
sx q[2];
rz(-2.9439587) q[2];
rz(-2.8526784) q[3];
sx q[3];
rz(-0.87696004) q[3];
sx q[3];
rz(1.7355841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0777247) q[0];
sx q[0];
rz(-0.48982319) q[0];
sx q[0];
rz(2.9329964) q[0];
rz(2.1754307) q[1];
sx q[1];
rz(-1.9816793) q[1];
sx q[1];
rz(-1.6360412) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9150328) q[0];
sx q[0];
rz(-2.0953) q[0];
sx q[0];
rz(-2.4302308) q[0];
x q[1];
rz(3.0622919) q[2];
sx q[2];
rz(-2.7100025) q[2];
sx q[2];
rz(2.0183795) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.70430763) q[1];
sx q[1];
rz(-1.1204801) q[1];
sx q[1];
rz(-0.98547658) q[1];
x q[2];
rz(-0.286245) q[3];
sx q[3];
rz(-1.4713333) q[3];
sx q[3];
rz(0.092982987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.2839526) q[2];
sx q[2];
rz(-2.3985034) q[2];
sx q[2];
rz(0.097578438) q[2];
rz(1.7476667) q[3];
sx q[3];
rz(-1.7912309) q[3];
sx q[3];
rz(-2.424749) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7512648) q[0];
sx q[0];
rz(-1.8126235) q[0];
sx q[0];
rz(0.51399291) q[0];
rz(-0.12318525) q[1];
sx q[1];
rz(-0.24736483) q[1];
sx q[1];
rz(-0.93200144) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8228553) q[0];
sx q[0];
rz(-1.1294951) q[0];
sx q[0];
rz(-2.9916828) q[0];
rz(-1.2134238) q[2];
sx q[2];
rz(-0.66713152) q[2];
sx q[2];
rz(0.92598976) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.2145558) q[1];
sx q[1];
rz(-2.3408457) q[1];
sx q[1];
rz(0.75896778) q[1];
rz(-2.9429432) q[3];
sx q[3];
rz(-2.2774787) q[3];
sx q[3];
rz(-2.228565) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.0294068) q[2];
sx q[2];
rz(-1.1108578) q[2];
sx q[2];
rz(2.712148) q[2];
rz(1.2094234) q[3];
sx q[3];
rz(-2.7691787) q[3];
sx q[3];
rz(0.69303524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3747303) q[0];
sx q[0];
rz(-1.466789) q[0];
sx q[0];
rz(-1.3388348) q[0];
rz(0.70612899) q[1];
sx q[1];
rz(-1.8920205) q[1];
sx q[1];
rz(-2.1910117) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75854077) q[0];
sx q[0];
rz(-1.0785111) q[0];
sx q[0];
rz(-0.58347115) q[0];
x q[1];
rz(-1.9865932) q[2];
sx q[2];
rz(-1.5019413) q[2];
sx q[2];
rz(-1.5037675) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.7710167) q[1];
sx q[1];
rz(-2.5713213) q[1];
sx q[1];
rz(-1.2194521) q[1];
x q[2];
rz(-0.25456984) q[3];
sx q[3];
rz(-1.9707142) q[3];
sx q[3];
rz(-0.039777048) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.4799698) q[2];
sx q[2];
rz(-1.491549) q[2];
sx q[2];
rz(-0.48428145) q[2];
rz(0.92710036) q[3];
sx q[3];
rz(-1.8323332) q[3];
sx q[3];
rz(2.5203729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
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
rz(-2.9973688) q[0];
sx q[0];
rz(-0.08865083) q[0];
sx q[0];
rz(0.22928672) q[0];
rz(2.7067822) q[1];
sx q[1];
rz(-1.2229342) q[1];
sx q[1];
rz(0.71892175) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35480354) q[0];
sx q[0];
rz(-0.94593404) q[0];
sx q[0];
rz(-3.025269) q[0];
rz(-pi) q[1];
x q[1];
rz(0.11833338) q[2];
sx q[2];
rz(-2.1552857) q[2];
sx q[2];
rz(0.10533939) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.016991888) q[1];
sx q[1];
rz(-1.4233839) q[1];
sx q[1];
rz(0.37313811) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1925681) q[3];
sx q[3];
rz(-1.9927295) q[3];
sx q[3];
rz(2.4728647) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.7583313) q[2];
sx q[2];
rz(-1.9394082) q[2];
sx q[2];
rz(-2.3948005) q[2];
rz(2.2693999) q[3];
sx q[3];
rz(-2.3132497) q[3];
sx q[3];
rz(-2.1993568) q[3];
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
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9185716) q[0];
sx q[0];
rz(-1.3970319) q[0];
sx q[0];
rz(-1.3402517) q[0];
rz(-0.37721286) q[1];
sx q[1];
rz(-1.4706392) q[1];
sx q[1];
rz(-2.6249862) q[1];
rz(0.51858356) q[2];
sx q[2];
rz(-1.8443783) q[2];
sx q[2];
rz(1.5658866) q[2];
rz(-2.4242998) q[3];
sx q[3];
rz(-1.4581231) q[3];
sx q[3];
rz(0.067618528) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
