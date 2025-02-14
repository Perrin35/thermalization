OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.8283451) q[0];
sx q[0];
rz(-0.77594835) q[0];
sx q[0];
rz(2.0394072) q[0];
rz(1.6232396) q[1];
sx q[1];
rz(2.0145388) q[1];
sx q[1];
rz(9.1993499) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.091087393) q[0];
sx q[0];
rz(-0.69982547) q[0];
sx q[0];
rz(-0.096629337) q[0];
rz(-pi) q[1];
rz(1.4233474) q[2];
sx q[2];
rz(-2.1814697) q[2];
sx q[2];
rz(1.5491279) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.3662339) q[1];
sx q[1];
rz(-0.45985301) q[1];
sx q[1];
rz(-2.2103106) q[1];
rz(-0.4938287) q[3];
sx q[3];
rz(-1.8803616) q[3];
sx q[3];
rz(2.1403596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.48308358) q[2];
sx q[2];
rz(-2.4675214) q[2];
sx q[2];
rz(-3.1209514) q[2];
rz(-0.189273) q[3];
sx q[3];
rz(-0.34957379) q[3];
sx q[3];
rz(-1.4741723) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78057688) q[0];
sx q[0];
rz(-0.78106946) q[0];
sx q[0];
rz(-2.1625157) q[0];
rz(-3.0304404) q[1];
sx q[1];
rz(-0.73768288) q[1];
sx q[1];
rz(2.0806064) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7951668) q[0];
sx q[0];
rz(-2.0316664) q[0];
sx q[0];
rz(2.5823442) q[0];
rz(-0.0060180863) q[2];
sx q[2];
rz(-0.13745795) q[2];
sx q[2];
rz(-2.6743367) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.5630305) q[1];
sx q[1];
rz(-2.6075224) q[1];
sx q[1];
rz(-0.41407163) q[1];
rz(-pi) q[2];
rz(-0.33943601) q[3];
sx q[3];
rz(-2.7314679) q[3];
sx q[3];
rz(-2.8462178) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.4336808) q[2];
sx q[2];
rz(-1.8209063) q[2];
sx q[2];
rz(2.743538) q[2];
rz(1.7178644) q[3];
sx q[3];
rz(-0.28702304) q[3];
sx q[3];
rz(-2.6557693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-2.5240391) q[0];
sx q[0];
rz(-2.6300639) q[0];
sx q[0];
rz(1.0523354) q[0];
rz(2.6984093) q[1];
sx q[1];
rz(-0.96142238) q[1];
sx q[1];
rz(0.65394941) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1285514) q[0];
sx q[0];
rz(-1.4440026) q[0];
sx q[0];
rz(-2.3953505) q[0];
rz(0.072418173) q[2];
sx q[2];
rz(-2.1927367) q[2];
sx q[2];
rz(-2.3427013) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.941085) q[1];
sx q[1];
rz(-1.9882747) q[1];
sx q[1];
rz(0.79065506) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7486568) q[3];
sx q[3];
rz(-1.3452969) q[3];
sx q[3];
rz(-2.5009843) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.7707278) q[2];
sx q[2];
rz(-0.76954049) q[2];
sx q[2];
rz(-0.23180836) q[2];
rz(1.2109463) q[3];
sx q[3];
rz(-1.1823267) q[3];
sx q[3];
rz(-0.32538357) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4732707) q[0];
sx q[0];
rz(-0.40437651) q[0];
sx q[0];
rz(0.37687287) q[0];
rz(-2.6301774) q[1];
sx q[1];
rz(-1.8835953) q[1];
sx q[1];
rz(2.2976141) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6451615) q[0];
sx q[0];
rz(-1.58824) q[0];
sx q[0];
rz(1.6806127) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8966214) q[2];
sx q[2];
rz(-0.60403901) q[2];
sx q[2];
rz(2.2789795) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.6214024) q[1];
sx q[1];
rz(-2.108413) q[1];
sx q[1];
rz(1.0034087) q[1];
rz(-pi) q[2];
rz(-3.1042905) q[3];
sx q[3];
rz(-2.6196194) q[3];
sx q[3];
rz(-0.69729174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.1692928) q[2];
sx q[2];
rz(-1.0147213) q[2];
sx q[2];
rz(-0.69804066) q[2];
rz(-1.1997148) q[3];
sx q[3];
rz(-0.90568393) q[3];
sx q[3];
rz(0.5996632) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19071628) q[0];
sx q[0];
rz(-0.27442014) q[0];
sx q[0];
rz(-3.0158667) q[0];
rz(1.0267286) q[1];
sx q[1];
rz(-1.3871565) q[1];
sx q[1];
rz(0.52621192) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65186319) q[0];
sx q[0];
rz(-1.4966156) q[0];
sx q[0];
rz(-0.47898888) q[0];
rz(-0.87398087) q[2];
sx q[2];
rz(-1.4451728) q[2];
sx q[2];
rz(-1.2619293) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.3354644) q[1];
sx q[1];
rz(-1.4683691) q[1];
sx q[1];
rz(-2.7598937) q[1];
rz(2.3461269) q[3];
sx q[3];
rz(-1.0364) q[3];
sx q[3];
rz(-0.55783844) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.0937664) q[2];
sx q[2];
rz(-0.83790773) q[2];
sx q[2];
rz(0.31025904) q[2];
rz(0.034916498) q[3];
sx q[3];
rz(-1.3860476) q[3];
sx q[3];
rz(-2.3518899) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-2.0652311) q[0];
sx q[0];
rz(-1.6595474) q[0];
sx q[0];
rz(-0.37102997) q[0];
rz(-0.064844355) q[1];
sx q[1];
rz(-2.3193181) q[1];
sx q[1];
rz(-1.8745905) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1112087) q[0];
sx q[0];
rz(-0.74978854) q[0];
sx q[0];
rz(2.0862338) q[0];
rz(-2.1999938) q[2];
sx q[2];
rz(-1.2824739) q[2];
sx q[2];
rz(2.0137816) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.53643287) q[1];
sx q[1];
rz(-0.39966941) q[1];
sx q[1];
rz(-0.088614806) q[1];
rz(-pi) q[2];
rz(-0.53474769) q[3];
sx q[3];
rz(-0.52899299) q[3];
sx q[3];
rz(0.76120602) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(3.0176257) q[2];
sx q[2];
rz(-1.664868) q[2];
sx q[2];
rz(1.6032479) q[2];
rz(2.7247143) q[3];
sx q[3];
rz(-2.6346801) q[3];
sx q[3];
rz(0.1703593) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-0.040319547) q[0];
sx q[0];
rz(-0.23389255) q[0];
sx q[0];
rz(-2.7129569) q[0];
rz(-1.1242695) q[1];
sx q[1];
rz(-1.5391896) q[1];
sx q[1];
rz(0.7695778) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3807629) q[0];
sx q[0];
rz(-1.5738945) q[0];
sx q[0];
rz(0.0018381434) q[0];
rz(-0.52346241) q[2];
sx q[2];
rz(-1.133073) q[2];
sx q[2];
rz(1.4283534) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.4691041) q[1];
sx q[1];
rz(-2.1616363) q[1];
sx q[1];
rz(-1.2149151) q[1];
x q[2];
rz(-1.0060723) q[3];
sx q[3];
rz(-2.5556563) q[3];
sx q[3];
rz(2.1028818) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.6577242) q[2];
sx q[2];
rz(-1.8572073) q[2];
sx q[2];
rz(-0.095495187) q[2];
rz(0.5419845) q[3];
sx q[3];
rz(-2.675455) q[3];
sx q[3];
rz(3.0123762) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99590456) q[0];
sx q[0];
rz(-2.2112084) q[0];
sx q[0];
rz(-3.0297739) q[0];
rz(1.3189141) q[1];
sx q[1];
rz(-1.2448064) q[1];
sx q[1];
rz(-1.8359312) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5163286) q[0];
sx q[0];
rz(-0.56076196) q[0];
sx q[0];
rz(-1.3365585) q[0];
x q[1];
rz(2.6839031) q[2];
sx q[2];
rz(-2.8297242) q[2];
sx q[2];
rz(1.0491187) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9010734) q[1];
sx q[1];
rz(-1.8522634) q[1];
sx q[1];
rz(2.752384) q[1];
x q[2];
rz(-1.9910802) q[3];
sx q[3];
rz(-2.0788361) q[3];
sx q[3];
rz(-0.86364323) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.8127415) q[2];
sx q[2];
rz(-1.988827) q[2];
sx q[2];
rz(2.0218772) q[2];
rz(-0.26851922) q[3];
sx q[3];
rz(-1.7374141) q[3];
sx q[3];
rz(-1.9560811) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28065228) q[0];
sx q[0];
rz(-0.76186162) q[0];
sx q[0];
rz(0.60894388) q[0];
rz(-0.87743419) q[1];
sx q[1];
rz(-1.0017706) q[1];
sx q[1];
rz(0.61224365) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5615879) q[0];
sx q[0];
rz(-1.1682434) q[0];
sx q[0];
rz(-2.3775565) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2127146) q[2];
sx q[2];
rz(-2.1109258) q[2];
sx q[2];
rz(0.071252099) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.22837199) q[1];
sx q[1];
rz(-0.84443426) q[1];
sx q[1];
rz(0.63131632) q[1];
x q[2];
rz(3.0123467) q[3];
sx q[3];
rz(-0.5448444) q[3];
sx q[3];
rz(-2.2428577) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0683384) q[2];
sx q[2];
rz(-1.5404258) q[2];
sx q[2];
rz(0.18745984) q[2];
rz(2.0295664) q[3];
sx q[3];
rz(-0.91809648) q[3];
sx q[3];
rz(-2.658127) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7336695) q[0];
sx q[0];
rz(-0.81403553) q[0];
sx q[0];
rz(0.2440051) q[0];
rz(-1.4834652) q[1];
sx q[1];
rz(-1.6226059) q[1];
sx q[1];
rz(2.3689178) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2940154) q[0];
sx q[0];
rz(-0.71888598) q[0];
sx q[0];
rz(-0.80475828) q[0];
rz(-1.7079321) q[2];
sx q[2];
rz(-1.8923645) q[2];
sx q[2];
rz(2.0069881) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.0341943) q[1];
sx q[1];
rz(-1.6813283) q[1];
sx q[1];
rz(-2.8128036) q[1];
rz(-pi) q[2];
rz(2.3546702) q[3];
sx q[3];
rz(-0.95810181) q[3];
sx q[3];
rz(-0.69132346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8030168) q[2];
sx q[2];
rz(-2.1803941) q[2];
sx q[2];
rz(-2.4147066) q[2];
rz(2.0342942) q[3];
sx q[3];
rz(-2.6328583) q[3];
sx q[3];
rz(1.0072964) q[3];
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
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0016639391) q[0];
sx q[0];
rz(-1.7353084) q[0];
sx q[0];
rz(-1.1575862) q[0];
rz(2.6955556) q[1];
sx q[1];
rz(-1.0855433) q[1];
sx q[1];
rz(-0.6378508) q[1];
rz(-0.28631306) q[2];
sx q[2];
rz(-2.7234017) q[2];
sx q[2];
rz(-1.0665733) q[2];
rz(0.31728716) q[3];
sx q[3];
rz(-2.696456) q[3];
sx q[3];
rz(1.0831931) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
