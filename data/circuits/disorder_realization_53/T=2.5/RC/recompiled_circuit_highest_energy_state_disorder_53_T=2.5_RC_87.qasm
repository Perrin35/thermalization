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
rz(1.9617957) q[0];
sx q[0];
rz(-2.0834162) q[0];
sx q[0];
rz(-2.7297468) q[0];
rz(-2.5208199) q[1];
sx q[1];
rz(-2.4359735) q[1];
sx q[1];
rz(0.128428) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1415048) q[0];
sx q[0];
rz(-0.94188085) q[0];
sx q[0];
rz(-1.5733243) q[0];
rz(-pi) q[1];
rz(1.0780222) q[2];
sx q[2];
rz(-1.2076305) q[2];
sx q[2];
rz(0.21445477) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.7196178) q[1];
sx q[1];
rz(-2.1901844) q[1];
sx q[1];
rz(-1.3307816) q[1];
x q[2];
rz(-1.7953231) q[3];
sx q[3];
rz(-2.7347964) q[3];
sx q[3];
rz(-1.612839) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.8992369) q[2];
sx q[2];
rz(-0.54348102) q[2];
sx q[2];
rz(0.71329722) q[2];
rz(-1.9244309) q[3];
sx q[3];
rz(-1.7141914) q[3];
sx q[3];
rz(-3.0349019) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84739581) q[0];
sx q[0];
rz(-1.0260181) q[0];
sx q[0];
rz(3.078939) q[0];
rz(-0.90473908) q[1];
sx q[1];
rz(-0.86974564) q[1];
sx q[1];
rz(3.0135801) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0580095) q[0];
sx q[0];
rz(-1.9045275) q[0];
sx q[0];
rz(0.34837848) q[0];
x q[1];
rz(-2.2163311) q[2];
sx q[2];
rz(-1.0224258) q[2];
sx q[2];
rz(-0.59522168) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.74740313) q[1];
sx q[1];
rz(-0.48626562) q[1];
sx q[1];
rz(-3.1267605) q[1];
rz(-pi) q[2];
x q[2];
rz(0.15897922) q[3];
sx q[3];
rz(-2.3975054) q[3];
sx q[3];
rz(-1.9682304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.5855828) q[2];
sx q[2];
rz(-1.2051008) q[2];
sx q[2];
rz(-2.0466059) q[2];
rz(2.5321142) q[3];
sx q[3];
rz(-2.1938117) q[3];
sx q[3];
rz(1.8402717) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0037917) q[0];
sx q[0];
rz(-0.68313685) q[0];
sx q[0];
rz(2.0004499) q[0];
rz(-1.2992651) q[1];
sx q[1];
rz(-1.3971993) q[1];
sx q[1];
rz(3.0019147) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8477301) q[0];
sx q[0];
rz(-2.8314857) q[0];
sx q[0];
rz(-2.1523813) q[0];
x q[1];
rz(-1.79931) q[2];
sx q[2];
rz(-1.7916792) q[2];
sx q[2];
rz(2.8877838) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.7157212) q[1];
sx q[1];
rz(-0.46498659) q[1];
sx q[1];
rz(2.6816508) q[1];
x q[2];
rz(-2.1909149) q[3];
sx q[3];
rz(-1.0178627) q[3];
sx q[3];
rz(-1.7861136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7695158) q[2];
sx q[2];
rz(-2.3774827) q[2];
sx q[2];
rz(0.3176983) q[2];
rz(2.5959173) q[3];
sx q[3];
rz(-0.92394865) q[3];
sx q[3];
rz(-1.8359449) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9724156) q[0];
sx q[0];
rz(-1.7298052) q[0];
sx q[0];
rz(0.2485982) q[0];
rz(1.2480805) q[1];
sx q[1];
rz(-2.0781519) q[1];
sx q[1];
rz(-0.6691106) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6921223) q[0];
sx q[0];
rz(-0.50296268) q[0];
sx q[0];
rz(0.19900222) q[0];
x q[1];
rz(-1.2641773) q[2];
sx q[2];
rz(-2.3117522) q[2];
sx q[2];
rz(1.4640984) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1966038) q[1];
sx q[1];
rz(-1.7763864) q[1];
sx q[1];
rz(0.14954549) q[1];
x q[2];
rz(2.6982299) q[3];
sx q[3];
rz(-2.0896455) q[3];
sx q[3];
rz(-2.8302397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.3847547) q[2];
sx q[2];
rz(-2.8110795) q[2];
sx q[2];
rz(3.030704) q[2];
rz(-0.97892654) q[3];
sx q[3];
rz(-1.7612235) q[3];
sx q[3];
rz(-3.0389649) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0961618) q[0];
sx q[0];
rz(-1.0935421) q[0];
sx q[0];
rz(-2.9683215) q[0];
rz(0.551956) q[1];
sx q[1];
rz(-2.5872784) q[1];
sx q[1];
rz(-0.90131235) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2103268) q[0];
sx q[0];
rz(-1.0823609) q[0];
sx q[0];
rz(0.68426056) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4129139) q[2];
sx q[2];
rz(-2.4167293) q[2];
sx q[2];
rz(-0.55654364) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.3696987) q[1];
sx q[1];
rz(-1.7189184) q[1];
sx q[1];
rz(0.019031899) q[1];
rz(-pi) q[2];
rz(-2.0088193) q[3];
sx q[3];
rz(-1.1896648) q[3];
sx q[3];
rz(2.1025859) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.99894714) q[2];
sx q[2];
rz(-1.877942) q[2];
sx q[2];
rz(0.20550814) q[2];
rz(-0.62028003) q[3];
sx q[3];
rz(-0.55448237) q[3];
sx q[3];
rz(1.8581355) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.810629) q[0];
sx q[0];
rz(-2.0331419) q[0];
sx q[0];
rz(-0.95322815) q[0];
rz(-0.78298059) q[1];
sx q[1];
rz(-2.6278261) q[1];
sx q[1];
rz(-0.89797529) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6655639) q[0];
sx q[0];
rz(-1.3485639) q[0];
sx q[0];
rz(2.1118947) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.4250323) q[2];
sx q[2];
rz(-1.9117182) q[2];
sx q[2];
rz(1.2630315) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6728357) q[1];
sx q[1];
rz(-1.7544946) q[1];
sx q[1];
rz(0.10803799) q[1];
x q[2];
rz(-2.9457568) q[3];
sx q[3];
rz(-1.2863428) q[3];
sx q[3];
rz(-1.0468512) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.7440685) q[2];
sx q[2];
rz(-0.70647883) q[2];
sx q[2];
rz(-1.2495329) q[2];
rz(1.9254098) q[3];
sx q[3];
rz(-1.774923) q[3];
sx q[3];
rz(-1.3528311) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3406496) q[0];
sx q[0];
rz(-0.018434374) q[0];
sx q[0];
rz(1.8311485) q[0];
rz(0.75366968) q[1];
sx q[1];
rz(-2.8050656) q[1];
sx q[1];
rz(-0.15187844) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.150249) q[0];
sx q[0];
rz(-1.4785355) q[0];
sx q[0];
rz(-1.4982066) q[0];
x q[1];
rz(-0.23252587) q[2];
sx q[2];
rz(-2.8029685) q[2];
sx q[2];
rz(1.9411638) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.79287) q[1];
sx q[1];
rz(-1.486549) q[1];
sx q[1];
rz(2.0920175) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8009284) q[3];
sx q[3];
rz(-2.4471183) q[3];
sx q[3];
rz(-2.9509773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.14022961) q[2];
sx q[2];
rz(-1.9915853) q[2];
sx q[2];
rz(1.4353732) q[2];
rz(0.6048778) q[3];
sx q[3];
rz(-1.6717654) q[3];
sx q[3];
rz(0.83355561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3816933) q[0];
sx q[0];
rz(-0.23867358) q[0];
sx q[0];
rz(-2.6651486) q[0];
rz(0.53642383) q[1];
sx q[1];
rz(-1.3342349) q[1];
sx q[1];
rz(0.36472067) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3140253) q[0];
sx q[0];
rz(-2.7430581) q[0];
sx q[0];
rz(2.3588595) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.1171072) q[2];
sx q[2];
rz(-2.6016781) q[2];
sx q[2];
rz(-1.2123666) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.1583511) q[1];
sx q[1];
rz(-0.73136202) q[1];
sx q[1];
rz(1.8476827) q[1];
rz(-pi) q[2];
x q[2];
rz(2.164806) q[3];
sx q[3];
rz(-1.5100118) q[3];
sx q[3];
rz(-1.9171627) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.81858188) q[2];
sx q[2];
rz(-1.3280758) q[2];
sx q[2];
rz(2.2200457) q[2];
rz(-1.3449097) q[3];
sx q[3];
rz(-1.4321046) q[3];
sx q[3];
rz(0.014605453) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0607249) q[0];
sx q[0];
rz(-0.23268555) q[0];
sx q[0];
rz(-3.0423855) q[0];
rz(-1.5946782) q[1];
sx q[1];
rz(-2.4989936) q[1];
sx q[1];
rz(-3.0329472) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0091054) q[0];
sx q[0];
rz(-2.0901839) q[0];
sx q[0];
rz(-0.035196134) q[0];
x q[1];
rz(2.2021239) q[2];
sx q[2];
rz(-1.924404) q[2];
sx q[2];
rz(-0.58009321) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.33230072) q[1];
sx q[1];
rz(-0.58590472) q[1];
sx q[1];
rz(-2.078474) q[1];
x q[2];
rz(-0.21317706) q[3];
sx q[3];
rz(-1.2090599) q[3];
sx q[3];
rz(-0.26644275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.8607311) q[2];
sx q[2];
rz(-0.18605575) q[2];
sx q[2];
rz(0.59530386) q[2];
rz(2.2286277) q[3];
sx q[3];
rz(-0.76415092) q[3];
sx q[3];
rz(-1.3179774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2309793) q[0];
sx q[0];
rz(-1.148104) q[0];
sx q[0];
rz(-0.30368152) q[0];
rz(0.23823711) q[1];
sx q[1];
rz(-0.99681011) q[1];
sx q[1];
rz(0.69746596) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8236341) q[0];
sx q[0];
rz(-2.1058361) q[0];
sx q[0];
rz(-1.3617424) q[0];
rz(-pi) q[1];
x q[1];
rz(0.35665193) q[2];
sx q[2];
rz(-1.4809297) q[2];
sx q[2];
rz(-1.9626475) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.94983236) q[1];
sx q[1];
rz(-1.0737541) q[1];
sx q[1];
rz(2.0326896) q[1];
x q[2];
rz(-0.26106278) q[3];
sx q[3];
rz(-1.1202284) q[3];
sx q[3];
rz(1.0999668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.4136228) q[2];
sx q[2];
rz(-2.0406849) q[2];
sx q[2];
rz(0.093416365) q[2];
rz(1.5824687) q[3];
sx q[3];
rz(-3.0714572) q[3];
sx q[3];
rz(-0.087433405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14417905) q[0];
sx q[0];
rz(-1.1638673) q[0];
sx q[0];
rz(-1.5735004) q[0];
rz(-2.6724124) q[1];
sx q[1];
rz(-2.4041685) q[1];
sx q[1];
rz(2.266177) q[1];
rz(-2.6826674) q[2];
sx q[2];
rz(-1.0815207) q[2];
sx q[2];
rz(-1.8455119) q[2];
rz(-0.31147891) q[3];
sx q[3];
rz(-2.3337915) q[3];
sx q[3];
rz(1.3561377) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
