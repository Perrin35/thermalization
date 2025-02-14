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
rz(1.3156112) q[0];
sx q[0];
rz(5.4551107) q[0];
sx q[0];
rz(8.5760737) q[0];
rz(1.4915713) q[1];
sx q[1];
rz(-0.68869156) q[1];
sx q[1];
rz(-2.7149849) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85159447) q[0];
sx q[0];
rz(-1.6713872) q[0];
sx q[0];
rz(-1.6466584) q[0];
rz(0.83580221) q[2];
sx q[2];
rz(-2.9907132) q[2];
sx q[2];
rz(0.41656938) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.9200807) q[1];
sx q[1];
rz(-2.6903841) q[1];
sx q[1];
rz(-0.16772224) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3412104) q[3];
sx q[3];
rz(-1.1100148) q[3];
sx q[3];
rz(-0.15298259) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.2171057) q[2];
sx q[2];
rz(-1.3873528) q[2];
sx q[2];
rz(3.1021049) q[2];
rz(-1.0314137) q[3];
sx q[3];
rz(-2.1125427) q[3];
sx q[3];
rz(-1.903681) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6295488) q[0];
sx q[0];
rz(-1.5271674) q[0];
sx q[0];
rz(1.7426096) q[0];
rz(1.5139187) q[1];
sx q[1];
rz(-0.84295034) q[1];
sx q[1];
rz(-1.7248076) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91367224) q[0];
sx q[0];
rz(-2.3487072) q[0];
sx q[0];
rz(2.8280054) q[0];
x q[1];
rz(-2.3624647) q[2];
sx q[2];
rz(-1.3848403) q[2];
sx q[2];
rz(2.8715239) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.4128215) q[1];
sx q[1];
rz(-1.1755474) q[1];
sx q[1];
rz(-1.4229619) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6187158) q[3];
sx q[3];
rz(-1.9780434) q[3];
sx q[3];
rz(-1.4496692) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.9713355) q[2];
sx q[2];
rz(-2.038326) q[2];
sx q[2];
rz(-2.4309168) q[2];
rz(3.1273048) q[3];
sx q[3];
rz(-2.894214) q[3];
sx q[3];
rz(-1.3361196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9209442) q[0];
sx q[0];
rz(-0.9592239) q[0];
sx q[0];
rz(2.7253819) q[0];
rz(1.2843708) q[1];
sx q[1];
rz(-1.6651848) q[1];
sx q[1];
rz(0.31825569) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98781768) q[0];
sx q[0];
rz(-1.1889699) q[0];
sx q[0];
rz(2.8433958) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2936965) q[2];
sx q[2];
rz(-1.0506786) q[2];
sx q[2];
rz(-2.0682356) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.15662312) q[1];
sx q[1];
rz(-0.92069101) q[1];
sx q[1];
rz(-1.7414879) q[1];
x q[2];
rz(2.8937723) q[3];
sx q[3];
rz(-1.8704924) q[3];
sx q[3];
rz(-1.9639942) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.8754114) q[2];
sx q[2];
rz(-2.3616932) q[2];
sx q[2];
rz(1.2904588) q[2];
rz(-0.053704638) q[3];
sx q[3];
rz(-1.3196245) q[3];
sx q[3];
rz(1.7558338) q[3];
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
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5829492) q[0];
sx q[0];
rz(-0.68840331) q[0];
sx q[0];
rz(1.0097367) q[0];
rz(-0.91521493) q[1];
sx q[1];
rz(-1.5292294) q[1];
sx q[1];
rz(-2.7161652) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64176004) q[0];
sx q[0];
rz(-1.7374037) q[0];
sx q[0];
rz(-0.84979041) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.073510344) q[2];
sx q[2];
rz(-0.82907721) q[2];
sx q[2];
rz(0.69141594) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.3784619) q[1];
sx q[1];
rz(-1.9768856) q[1];
sx q[1];
rz(2.2927845) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5365547) q[3];
sx q[3];
rz(-0.16290191) q[3];
sx q[3];
rz(2.3690641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.90011826) q[2];
sx q[2];
rz(-1.5510617) q[2];
sx q[2];
rz(-0.11108622) q[2];
rz(-2.9491718) q[3];
sx q[3];
rz(-0.40168732) q[3];
sx q[3];
rz(2.4470611) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7531994) q[0];
sx q[0];
rz(-2.41112) q[0];
sx q[0];
rz(-2.2764192) q[0];
rz(-0.49258891) q[1];
sx q[1];
rz(-2.045423) q[1];
sx q[1];
rz(1.385484) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2619721) q[0];
sx q[0];
rz(-0.79085717) q[0];
sx q[0];
rz(1.7636931) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0086914) q[2];
sx q[2];
rz(-2.3302493) q[2];
sx q[2];
rz(-1.0824301) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.4069479) q[1];
sx q[1];
rz(-2.4666767) q[1];
sx q[1];
rz(2.4823443) q[1];
rz(-pi) q[2];
x q[2];
rz(2.319544) q[3];
sx q[3];
rz(-2.1072142) q[3];
sx q[3];
rz(0.82335237) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.2546786) q[2];
sx q[2];
rz(-2.6723537) q[2];
sx q[2];
rz(2.4349507) q[2];
rz(0.32736579) q[3];
sx q[3];
rz(-1.2923765) q[3];
sx q[3];
rz(2.4842026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3051598) q[0];
sx q[0];
rz(-1.7921472) q[0];
sx q[0];
rz(-2.7850372) q[0];
rz(0.56218475) q[1];
sx q[1];
rz(-2.0088582) q[1];
sx q[1];
rz(-0.98639375) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1323469) q[0];
sx q[0];
rz(-2.2725687) q[0];
sx q[0];
rz(-0.67108564) q[0];
x q[1];
rz(-0.34511995) q[2];
sx q[2];
rz(-1.2212996) q[2];
sx q[2];
rz(1.0687168) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.79365208) q[1];
sx q[1];
rz(-0.75885443) q[1];
sx q[1];
rz(0.037686613) q[1];
rz(-pi) q[2];
rz(-1.3273193) q[3];
sx q[3];
rz(-1.7469353) q[3];
sx q[3];
rz(2.1139262) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.3761313) q[2];
sx q[2];
rz(-0.71530801) q[2];
sx q[2];
rz(0.7005271) q[2];
rz(0.96945196) q[3];
sx q[3];
rz(-2.1078096) q[3];
sx q[3];
rz(0.9303003) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1932061) q[0];
sx q[0];
rz(-2.8128615) q[0];
sx q[0];
rz(-0.1513924) q[0];
rz(1.5104177) q[1];
sx q[1];
rz(-2.0163586) q[1];
sx q[1];
rz(1.6815394) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3356381) q[0];
sx q[0];
rz(-2.3344669) q[0];
sx q[0];
rz(-0.7025848) q[0];
rz(-pi) q[1];
rz(-1.0437772) q[2];
sx q[2];
rz(-1.8879381) q[2];
sx q[2];
rz(-0.70995599) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.5269176) q[1];
sx q[1];
rz(-1.9567031) q[1];
sx q[1];
rz(0.98395394) q[1];
x q[2];
rz(2.2112956) q[3];
sx q[3];
rz(-0.2411763) q[3];
sx q[3];
rz(2.0335787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.2363362) q[2];
sx q[2];
rz(-2.4589804) q[2];
sx q[2];
rz(2.7899817) q[2];
rz(0.44175276) q[3];
sx q[3];
rz(-2.2124002) q[3];
sx q[3];
rz(2.9898804) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.041895954) q[0];
sx q[0];
rz(-1.2956887) q[0];
sx q[0];
rz(-1.9805441) q[0];
rz(0.64758045) q[1];
sx q[1];
rz(-1.7627962) q[1];
sx q[1];
rz(-0.82239282) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0121545) q[0];
sx q[0];
rz(-2.4732051) q[0];
sx q[0];
rz(2.2447615) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6641409) q[2];
sx q[2];
rz(-2.3465354) q[2];
sx q[2];
rz(-2.8964122) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.5288189) q[1];
sx q[1];
rz(-0.83136035) q[1];
sx q[1];
rz(0.39503132) q[1];
x q[2];
rz(-2.0045207) q[3];
sx q[3];
rz(-1.2081283) q[3];
sx q[3];
rz(0.39946242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.9245727) q[2];
sx q[2];
rz(-1.1447039) q[2];
sx q[2];
rz(1.3059957) q[2];
rz(2.4235587) q[3];
sx q[3];
rz(-0.82457232) q[3];
sx q[3];
rz(-0.048862783) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43593916) q[0];
sx q[0];
rz(-1.0834563) q[0];
sx q[0];
rz(0.017024592) q[0];
rz(1.1964993) q[1];
sx q[1];
rz(-2.7462609) q[1];
sx q[1];
rz(-2.8065525) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3075313) q[0];
sx q[0];
rz(-2.1387707) q[0];
sx q[0];
rz(-3.1295883) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.64249383) q[2];
sx q[2];
rz(-1.9064925) q[2];
sx q[2];
rz(0.030980008) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.9043873) q[1];
sx q[1];
rz(-2.4438624) q[1];
sx q[1];
rz(-1.2751139) q[1];
x q[2];
rz(2.0552956) q[3];
sx q[3];
rz(-2.006975) q[3];
sx q[3];
rz(2.0791813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7177141) q[2];
sx q[2];
rz(-2.3953891) q[2];
sx q[2];
rz(-2.8625873) q[2];
rz(-1.3689857) q[3];
sx q[3];
rz(-1.5149346) q[3];
sx q[3];
rz(-2.2122808) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(1.6245215) q[0];
sx q[0];
rz(-0.14685024) q[0];
sx q[0];
rz(0.76706925) q[0];
rz(2.7686139) q[1];
sx q[1];
rz(-1.6423128) q[1];
sx q[1];
rz(0.17072089) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51119186) q[0];
sx q[0];
rz(-2.3053279) q[0];
sx q[0];
rz(-1.6096576) q[0];
rz(-1.7319948) q[2];
sx q[2];
rz(-1.9809763) q[2];
sx q[2];
rz(0.31482492) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.77826277) q[1];
sx q[1];
rz(-0.40503854) q[1];
sx q[1];
rz(-2.2447849) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.3091273) q[3];
sx q[3];
rz(-0.48764519) q[3];
sx q[3];
rz(-1.6548827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.79147044) q[2];
sx q[2];
rz(-0.56362027) q[2];
sx q[2];
rz(-3.138809) q[2];
rz(2.2400098) q[3];
sx q[3];
rz(-1.0222579) q[3];
sx q[3];
rz(-0.80488747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0881385) q[0];
sx q[0];
rz(-2.8202941) q[0];
sx q[0];
rz(1.1711076) q[0];
rz(-2.3114655) q[1];
sx q[1];
rz(-1.8142038) q[1];
sx q[1];
rz(1.289191) q[1];
rz(2.2466462) q[2];
sx q[2];
rz(-2.4104626) q[2];
sx q[2];
rz(1.8117803) q[2];
rz(0.90972539) q[3];
sx q[3];
rz(-0.88539888) q[3];
sx q[3];
rz(-1.51487) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
