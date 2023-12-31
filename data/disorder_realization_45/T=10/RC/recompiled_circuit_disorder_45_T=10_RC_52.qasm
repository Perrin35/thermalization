OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.3061476) q[0];
sx q[0];
rz(-2.4581576) q[0];
sx q[0];
rz(-0.47877065) q[0];
rz(0.03102826) q[1];
sx q[1];
rz(5.1217084) q[1];
sx q[1];
rz(6.9245467) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7862579) q[0];
sx q[0];
rz(-1.9541652) q[0];
sx q[0];
rz(-1.3134365) q[0];
x q[1];
rz(0.37748572) q[2];
sx q[2];
rz(-1.971772) q[2];
sx q[2];
rz(-3.0149143) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.8909059) q[1];
sx q[1];
rz(-2.2911934) q[1];
sx q[1];
rz(-0.62178639) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.62946837) q[3];
sx q[3];
rz(-1.9068204) q[3];
sx q[3];
rz(1.1050841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.59149867) q[2];
sx q[2];
rz(-1.9248362) q[2];
sx q[2];
rz(-2.6634898) q[2];
rz(1.452662) q[3];
sx q[3];
rz(-2.0958459) q[3];
sx q[3];
rz(-12/(7*pi)) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1801382) q[0];
sx q[0];
rz(-0.77471662) q[0];
sx q[0];
rz(2.1226728) q[0];
rz(1.4787176) q[1];
sx q[1];
rz(-2.5264085) q[1];
sx q[1];
rz(-0.63308024) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7019254) q[0];
sx q[0];
rz(-0.81322008) q[0];
sx q[0];
rz(-1.4600091) q[0];
x q[1];
rz(-0.622153) q[2];
sx q[2];
rz(-2.2894147) q[2];
sx q[2];
rz(-2.6785786) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.9620348) q[1];
sx q[1];
rz(-1.1711367) q[1];
sx q[1];
rz(0.37079294) q[1];
x q[2];
rz(0.31140621) q[3];
sx q[3];
rz(-2.2283471) q[3];
sx q[3];
rz(2.1962375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.35976609) q[2];
sx q[2];
rz(-0.40010139) q[2];
sx q[2];
rz(-3.0241372) q[2];
rz(-2.84058) q[3];
sx q[3];
rz(-1.8071226) q[3];
sx q[3];
rz(-1.7787748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85787073) q[0];
sx q[0];
rz(-0.71325934) q[0];
sx q[0];
rz(3.0531847) q[0];
rz(1.1075426) q[1];
sx q[1];
rz(-0.91845599) q[1];
sx q[1];
rz(-2.450768) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3120964) q[0];
sx q[0];
rz(-1.4179686) q[0];
sx q[0];
rz(-1.6468847) q[0];
rz(-pi) q[1];
rz(1.3356528) q[2];
sx q[2];
rz(-2.3217839) q[2];
sx q[2];
rz(0.90607925) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.5882727) q[1];
sx q[1];
rz(-1.5447445) q[1];
sx q[1];
rz(1.7528898) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.74294528) q[3];
sx q[3];
rz(-2.6446614) q[3];
sx q[3];
rz(-0.039615354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.92418015) q[2];
sx q[2];
rz(-1.2574544) q[2];
sx q[2];
rz(3.0351191) q[2];
rz(1.4364093) q[3];
sx q[3];
rz(-0.36542106) q[3];
sx q[3];
rz(-1.6586554) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0687662) q[0];
sx q[0];
rz(-2.9077353) q[0];
sx q[0];
rz(0.52247125) q[0];
rz(-0.32896313) q[1];
sx q[1];
rz(-1.4986228) q[1];
sx q[1];
rz(-2.5879588) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.391511) q[0];
sx q[0];
rz(-0.9010074) q[0];
sx q[0];
rz(-2.1704587) q[0];
x q[1];
rz(-0.81981084) q[2];
sx q[2];
rz(-0.95633436) q[2];
sx q[2];
rz(-0.43624207) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.0370889) q[1];
sx q[1];
rz(-1.7040841) q[1];
sx q[1];
rz(-2.3549805) q[1];
rz(-pi) q[2];
x q[2];
rz(3.1133075) q[3];
sx q[3];
rz(-2.3445498) q[3];
sx q[3];
rz(-0.68238168) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.5840977) q[2];
sx q[2];
rz(-2.2258874) q[2];
sx q[2];
rz(-2.6468357) q[2];
rz(2.2385712) q[3];
sx q[3];
rz(-1.5938063) q[3];
sx q[3];
rz(1.8516699) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6506127) q[0];
sx q[0];
rz(-2.3968093) q[0];
sx q[0];
rz(-1.0850798) q[0];
rz(-2.1814573) q[1];
sx q[1];
rz(-1.1791869) q[1];
sx q[1];
rz(2.9575612) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9404011) q[0];
sx q[0];
rz(-1.860025) q[0];
sx q[0];
rz(1.4616696) q[0];
rz(0.63921914) q[2];
sx q[2];
rz(-2.6972065) q[2];
sx q[2];
rz(1.5612615) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6726471) q[1];
sx q[1];
rz(-0.99860672) q[1];
sx q[1];
rz(-2.1395626) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9156978) q[3];
sx q[3];
rz(-2.4028006) q[3];
sx q[3];
rz(1.6177288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.90157834) q[2];
sx q[2];
rz(-0.34873909) q[2];
sx q[2];
rz(-1.1408172) q[2];
rz(0.42282894) q[3];
sx q[3];
rz(-1.4774277) q[3];
sx q[3];
rz(-1.1269349) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35247701) q[0];
sx q[0];
rz(-1.1081835) q[0];
sx q[0];
rz(0.88678962) q[0];
rz(0.24041644) q[1];
sx q[1];
rz(-2.1227032) q[1];
sx q[1];
rz(-0.14850798) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7667023) q[0];
sx q[0];
rz(-0.74807157) q[0];
sx q[0];
rz(-2.9038249) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6806904) q[2];
sx q[2];
rz(-1.7981148) q[2];
sx q[2];
rz(2.0109039) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.099523274) q[1];
sx q[1];
rz(-1.4617141) q[1];
sx q[1];
rz(1.3101577) q[1];
x q[2];
rz(-2.4466483) q[3];
sx q[3];
rz(-2.5616025) q[3];
sx q[3];
rz(0.65141962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.285816) q[2];
sx q[2];
rz(-0.4824051) q[2];
sx q[2];
rz(-0.4883858) q[2];
rz(0.48940247) q[3];
sx q[3];
rz(-2.2198052) q[3];
sx q[3];
rz(-1.2667806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4483036) q[0];
sx q[0];
rz(-1.332809) q[0];
sx q[0];
rz(-0.81800246) q[0];
rz(-1.8891634) q[1];
sx q[1];
rz(-2.7493582) q[1];
sx q[1];
rz(2.0163527) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.274652) q[0];
sx q[0];
rz(-2.7526703) q[0];
sx q[0];
rz(-1.8453983) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9863425) q[2];
sx q[2];
rz(-1.8660188) q[2];
sx q[2];
rz(-1.0709907) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.4459671) q[1];
sx q[1];
rz(-1.8672767) q[1];
sx q[1];
rz(0.28709025) q[1];
x q[2];
rz(-2.4626477) q[3];
sx q[3];
rz(-0.79015398) q[3];
sx q[3];
rz(-1.3037579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0730878) q[2];
sx q[2];
rz(-0.98098522) q[2];
sx q[2];
rz(2.4970064) q[2];
rz(-1.4792431) q[3];
sx q[3];
rz(-0.95696604) q[3];
sx q[3];
rz(1.4054327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97061625) q[0];
sx q[0];
rz(-3.0727486) q[0];
sx q[0];
rz(-1.6059426) q[0];
rz(1.2212785) q[1];
sx q[1];
rz(-1.5131283) q[1];
sx q[1];
rz(2.1309526) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0518274) q[0];
sx q[0];
rz(-2.3065789) q[0];
sx q[0];
rz(2.9249973) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0534361) q[2];
sx q[2];
rz(-0.55609716) q[2];
sx q[2];
rz(0.14511395) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.1809363) q[1];
sx q[1];
rz(-1.3410543) q[1];
sx q[1];
rz(2.9147663) q[1];
rz(-pi) q[2];
rz(2.9270494) q[3];
sx q[3];
rz(-2.7326267) q[3];
sx q[3];
rz(1.4734801) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.5236139) q[2];
sx q[2];
rz(-0.31948677) q[2];
sx q[2];
rz(0.69407216) q[2];
rz(2.5726035) q[3];
sx q[3];
rz(-1.6945972) q[3];
sx q[3];
rz(-2.2038961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(-0.086031832) q[0];
sx q[0];
rz(-1.1431575) q[0];
sx q[0];
rz(-2.6468497) q[0];
rz(-2.5231979) q[1];
sx q[1];
rz(-1.4952375) q[1];
sx q[1];
rz(-3.0659952) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6583017) q[0];
sx q[0];
rz(-1.0169944) q[0];
sx q[0];
rz(2.7730586) q[0];
rz(-pi) q[1];
rz(1.7463023) q[2];
sx q[2];
rz(-1.5117466) q[2];
sx q[2];
rz(-0.67948558) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.96502111) q[1];
sx q[1];
rz(-2.2551943) q[1];
sx q[1];
rz(-1.7040764) q[1];
x q[2];
rz(-1.7528312) q[3];
sx q[3];
rz(-0.958003) q[3];
sx q[3];
rz(-0.36422563) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.8081234) q[2];
sx q[2];
rz(-1.418891) q[2];
sx q[2];
rz(-1.8010275) q[2];
rz(2.8373485) q[3];
sx q[3];
rz(-0.86383581) q[3];
sx q[3];
rz(-0.58469599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8463523) q[0];
sx q[0];
rz(-1.3051935) q[0];
sx q[0];
rz(2.9647968) q[0];
rz(1.8999752) q[1];
sx q[1];
rz(-2.5071564) q[1];
sx q[1];
rz(-0.44100824) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.037017578) q[0];
sx q[0];
rz(-1.3694166) q[0];
sx q[0];
rz(1.7778394) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4961692) q[2];
sx q[2];
rz(-1.4174889) q[2];
sx q[2];
rz(-0.15664936) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.1383789) q[1];
sx q[1];
rz(-2.5795476) q[1];
sx q[1];
rz(-2.2467062) q[1];
rz(-pi) q[2];
rz(0.40481683) q[3];
sx q[3];
rz(-0.66473085) q[3];
sx q[3];
rz(0.80617314) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.5796154) q[2];
sx q[2];
rz(-0.56869555) q[2];
sx q[2];
rz(-2.8397172) q[2];
rz(2.2484696) q[3];
sx q[3];
rz(-1.8538657) q[3];
sx q[3];
rz(-0.20475234) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4176035) q[0];
sx q[0];
rz(-1.2932734) q[0];
sx q[0];
rz(-1.4751157) q[0];
rz(0.026731116) q[1];
sx q[1];
rz(-1.4550799) q[1];
sx q[1];
rz(1.4310238) q[1];
rz(2.0886291) q[2];
sx q[2];
rz(-0.79635194) q[2];
sx q[2];
rz(1.2780381) q[2];
rz(1.0711014) q[3];
sx q[3];
rz(-2.069996) q[3];
sx q[3];
rz(1.3531006) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
