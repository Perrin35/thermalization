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
rz(-0.19718117) q[0];
sx q[0];
rz(4.2137094) q[0];
sx q[0];
rz(11.022047) q[0];
rz(0.3577258) q[1];
sx q[1];
rz(-1.3074713) q[1];
sx q[1];
rz(0.37582418) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7505219) q[0];
sx q[0];
rz(-0.76996917) q[0];
sx q[0];
rz(-3.0379437) q[0];
rz(-pi) q[1];
rz(0.4514427) q[2];
sx q[2];
rz(-0.57381781) q[2];
sx q[2];
rz(0.30893477) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.5511939) q[1];
sx q[1];
rz(-0.85764581) q[1];
sx q[1];
rz(-2.3354095) q[1];
x q[2];
rz(2.8055111) q[3];
sx q[3];
rz(-2.1770526) q[3];
sx q[3];
rz(1.9466024) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.69349217) q[2];
sx q[2];
rz(-0.38303146) q[2];
sx q[2];
rz(-2.7310272) q[2];
rz(-2.6170714) q[3];
sx q[3];
rz(-2.9295242) q[3];
sx q[3];
rz(-1.1493523) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9583424) q[0];
sx q[0];
rz(-0.14160937) q[0];
sx q[0];
rz(-2.435922) q[0];
rz(-1.2880098) q[1];
sx q[1];
rz(-2.5649004) q[1];
sx q[1];
rz(-1.1340595) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81728888) q[0];
sx q[0];
rz(-1.9729483) q[0];
sx q[0];
rz(-2.3020846) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.34313029) q[2];
sx q[2];
rz(-1.3995106) q[2];
sx q[2];
rz(0.73395455) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.6780939) q[1];
sx q[1];
rz(-1.1927407) q[1];
sx q[1];
rz(2.3854756) q[1];
x q[2];
rz(1.39246) q[3];
sx q[3];
rz(-1.0127676) q[3];
sx q[3];
rz(-1.4205281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.80295339) q[2];
sx q[2];
rz(-2.6889375) q[2];
sx q[2];
rz(0.53042859) q[2];
rz(2.8545634) q[3];
sx q[3];
rz(-0.48873264) q[3];
sx q[3];
rz(-0.73634017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1257783) q[0];
sx q[0];
rz(-1.1017355) q[0];
sx q[0];
rz(1.0544448) q[0];
rz(-1.6619445) q[1];
sx q[1];
rz(-2.2190861) q[1];
sx q[1];
rz(2.0747917) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5402601) q[0];
sx q[0];
rz(-1.6076438) q[0];
sx q[0];
rz(-1.3829492) q[0];
rz(1.6629024) q[2];
sx q[2];
rz(-2.2450441) q[2];
sx q[2];
rz(2.9363652) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.2640564) q[1];
sx q[1];
rz(-1.1095558) q[1];
sx q[1];
rz(0.25564927) q[1];
rz(-pi) q[2];
rz(0.85632944) q[3];
sx q[3];
rz(-1.3420055) q[3];
sx q[3];
rz(-1.3135873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.58831424) q[2];
sx q[2];
rz(-1.30043) q[2];
sx q[2];
rz(-1.1021357) q[2];
rz(-2.6831902) q[3];
sx q[3];
rz(-1.5285243) q[3];
sx q[3];
rz(2.5100191) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.691953) q[0];
sx q[0];
rz(-2.3821558) q[0];
sx q[0];
rz(-0.48680437) q[0];
rz(-1.5961157) q[1];
sx q[1];
rz(-2.7524452) q[1];
sx q[1];
rz(2.5130689) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5768054) q[0];
sx q[0];
rz(-0.41713762) q[0];
sx q[0];
rz(1.23853) q[0];
rz(1.3447464) q[2];
sx q[2];
rz(-0.58494431) q[2];
sx q[2];
rz(1.9247885) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.13884096) q[1];
sx q[1];
rz(-1.3831487) q[1];
sx q[1];
rz(1.1739338) q[1];
rz(1.0221338) q[3];
sx q[3];
rz(-1.5096153) q[3];
sx q[3];
rz(-3.1312628) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.57807606) q[2];
sx q[2];
rz(-0.11512828) q[2];
sx q[2];
rz(2.3797177) q[2];
rz(0.21841194) q[3];
sx q[3];
rz(-0.29774791) q[3];
sx q[3];
rz(-1.0485605) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31576425) q[0];
sx q[0];
rz(-2.789412) q[0];
sx q[0];
rz(2.5370362) q[0];
rz(1.293659) q[1];
sx q[1];
rz(-2.3217521) q[1];
sx q[1];
rz(0.34436071) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71925795) q[0];
sx q[0];
rz(-0.08503332) q[0];
sx q[0];
rz(1.812395) q[0];
rz(0.62744572) q[2];
sx q[2];
rz(-1.2973669) q[2];
sx q[2];
rz(0.16579096) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1923848) q[1];
sx q[1];
rz(-2.4864962) q[1];
sx q[1];
rz(-0.75836436) q[1];
rz(-0.46724288) q[3];
sx q[3];
rz(-2.1911216) q[3];
sx q[3];
rz(2.1324699) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.0209501) q[2];
sx q[2];
rz(-0.63434333) q[2];
sx q[2];
rz(-0.20993385) q[2];
rz(-1.1076934) q[3];
sx q[3];
rz(-1.6898797) q[3];
sx q[3];
rz(0.2618739) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8629859) q[0];
sx q[0];
rz(-0.59026533) q[0];
sx q[0];
rz(3.0721831) q[0];
rz(0.033992652) q[1];
sx q[1];
rz(-2.4885204) q[1];
sx q[1];
rz(2.6896844) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44563189) q[0];
sx q[0];
rz(-1.9602141) q[0];
sx q[0];
rz(-1.3971863) q[0];
rz(-pi) q[1];
rz(2.4546409) q[2];
sx q[2];
rz(-1.0712306) q[2];
sx q[2];
rz(-2.9925008) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.44575366) q[1];
sx q[1];
rz(-1.5924982) q[1];
sx q[1];
rz(-0.0071956102) q[1];
rz(-pi) q[2];
x q[2];
rz(2.389598) q[3];
sx q[3];
rz(-2.7772336) q[3];
sx q[3];
rz(2.2235699) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.41428337) q[2];
sx q[2];
rz(-0.80237979) q[2];
sx q[2];
rz(1.9963973) q[2];
rz(2.1833873) q[3];
sx q[3];
rz(-1.2088135) q[3];
sx q[3];
rz(0.27829471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77022839) q[0];
sx q[0];
rz(-2.6537277) q[0];
sx q[0];
rz(0.86139739) q[0];
rz(0.96249145) q[1];
sx q[1];
rz(-2.1884514) q[1];
sx q[1];
rz(2.9846587) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3478137) q[0];
sx q[0];
rz(-2.2484178) q[0];
sx q[0];
rz(2.6676446) q[0];
rz(-1.292615) q[2];
sx q[2];
rz(-2.7713994) q[2];
sx q[2];
rz(-2.7687759) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.2409901) q[1];
sx q[1];
rz(-2.1442736) q[1];
sx q[1];
rz(-0.9830703) q[1];
rz(-pi) q[2];
rz(1.3488999) q[3];
sx q[3];
rz(-0.21824868) q[3];
sx q[3];
rz(1.143592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-3.018205) q[2];
sx q[2];
rz(-2.3473479) q[2];
sx q[2];
rz(2.0920853) q[2];
rz(0.47979313) q[3];
sx q[3];
rz(-0.19194651) q[3];
sx q[3];
rz(0.16201365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
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
rz(2.4341693) q[0];
sx q[0];
rz(-0.24288414) q[0];
sx q[0];
rz(0.49651399) q[0];
rz(-1.4855344) q[1];
sx q[1];
rz(-0.86692202) q[1];
sx q[1];
rz(0.022580126) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17538769) q[0];
sx q[0];
rz(-1.3577355) q[0];
sx q[0];
rz(-2.9609462) q[0];
rz(-pi) q[1];
rz(3.0737799) q[2];
sx q[2];
rz(-1.4275506) q[2];
sx q[2];
rz(0.80729173) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.7851614) q[1];
sx q[1];
rz(-1.3884434) q[1];
sx q[1];
rz(2.3556832) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1053505) q[3];
sx q[3];
rz(-1.4009985) q[3];
sx q[3];
rz(-1.2359197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.61985832) q[2];
sx q[2];
rz(-2.0010184) q[2];
sx q[2];
rz(0.015259585) q[2];
rz(0.36733019) q[3];
sx q[3];
rz(-2.6945249) q[3];
sx q[3];
rz(0.36146155) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58886445) q[0];
sx q[0];
rz(-0.22607729) q[0];
sx q[0];
rz(-0.077762522) q[0];
rz(3.0231158) q[1];
sx q[1];
rz(-2.8006554) q[1];
sx q[1];
rz(2.4854787) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8828147) q[0];
sx q[0];
rz(-1.4442632) q[0];
sx q[0];
rz(1.8765429) q[0];
x q[1];
rz(-0.17123789) q[2];
sx q[2];
rz(-1.3067553) q[2];
sx q[2];
rz(2.03351) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.95893491) q[1];
sx q[1];
rz(-1.9709341) q[1];
sx q[1];
rz(1.047664) q[1];
rz(2.3967152) q[3];
sx q[3];
rz(-0.69760579) q[3];
sx q[3];
rz(3.0234203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.6180827) q[2];
sx q[2];
rz(-1.9205576) q[2];
sx q[2];
rz(-0.54100424) q[2];
rz(-2.6909761) q[3];
sx q[3];
rz(-1.485598) q[3];
sx q[3];
rz(-0.97644067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4377874) q[0];
sx q[0];
rz(-1.7923651) q[0];
sx q[0];
rz(-3.034814) q[0];
rz(-1.5497426) q[1];
sx q[1];
rz(-2.6390862) q[1];
sx q[1];
rz(-1.3776616) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6984682) q[0];
sx q[0];
rz(-1.6148991) q[0];
sx q[0];
rz(-1.1209784) q[0];
rz(-pi) q[1];
rz(-1.7936034) q[2];
sx q[2];
rz(-1.2651332) q[2];
sx q[2];
rz(-0.7204186) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.71860524) q[1];
sx q[1];
rz(-1.8180625) q[1];
sx q[1];
rz(-2.894677) q[1];
x q[2];
rz(2.1414457) q[3];
sx q[3];
rz(-1.1091145) q[3];
sx q[3];
rz(0.0014071597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0926853) q[2];
sx q[2];
rz(-1.7203628) q[2];
sx q[2];
rz(0.29472026) q[2];
rz(-0.60485351) q[3];
sx q[3];
rz(-1.9113144) q[3];
sx q[3];
rz(-3.0321583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-2.7034364) q[0];
sx q[0];
rz(-1.3180757) q[0];
sx q[0];
rz(-0.99676589) q[0];
rz(-3.1238212) q[1];
sx q[1];
rz(-1.8780864) q[1];
sx q[1];
rz(-1.3395739) q[1];
rz(-0.80452917) q[2];
sx q[2];
rz(-1.1297516) q[2];
sx q[2];
rz(-0.17937112) q[2];
rz(-0.8682438) q[3];
sx q[3];
rz(-1.0414185) q[3];
sx q[3];
rz(0.99440893) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
