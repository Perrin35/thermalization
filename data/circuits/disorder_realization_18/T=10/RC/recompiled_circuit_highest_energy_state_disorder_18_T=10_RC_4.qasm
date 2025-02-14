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
rz(1.7914766) q[0];
sx q[0];
rz(-0.67576367) q[0];
sx q[0];
rz(-0.0078553353) q[0];
rz(0.76771843) q[1];
sx q[1];
rz(-1.5239198) q[1];
sx q[1];
rz(2.89892) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15930106) q[0];
sx q[0];
rz(-2.2758599) q[0];
sx q[0];
rz(0.16601913) q[0];
x q[1];
rz(-2.484876) q[2];
sx q[2];
rz(-2.5829331) q[2];
sx q[2];
rz(0.61121537) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3625177) q[1];
sx q[1];
rz(-1.9240161) q[1];
sx q[1];
rz(-1.5455957) q[1];
x q[2];
rz(-2.8318066) q[3];
sx q[3];
rz(-0.68947809) q[3];
sx q[3];
rz(0.30287095) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.0365888) q[2];
sx q[2];
rz(-1.236311) q[2];
sx q[2];
rz(0.71199065) q[2];
rz(-0.30501929) q[3];
sx q[3];
rz(-2.9192393) q[3];
sx q[3];
rz(-3.0960848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63118339) q[0];
sx q[0];
rz(-1.2774066) q[0];
sx q[0];
rz(1.3674059) q[0];
rz(3.101688) q[1];
sx q[1];
rz(-0.80723643) q[1];
sx q[1];
rz(-2.5423999) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0586226) q[0];
sx q[0];
rz(-2.6520067) q[0];
sx q[0];
rz(-2.8955196) q[0];
rz(-pi) q[1];
rz(0.56053253) q[2];
sx q[2];
rz(-0.90122242) q[2];
sx q[2];
rz(-0.59434429) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2476462) q[1];
sx q[1];
rz(-0.94280137) q[1];
sx q[1];
rz(1.1537854) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.68568315) q[3];
sx q[3];
rz(-1.7910379) q[3];
sx q[3];
rz(-1.9626647) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0590608) q[2];
sx q[2];
rz(-2.579687) q[2];
sx q[2];
rz(1.3085636) q[2];
rz(-3.1176873) q[3];
sx q[3];
rz(-1.6126361) q[3];
sx q[3];
rz(-1.5628975) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6858653) q[0];
sx q[0];
rz(-2.4754334) q[0];
sx q[0];
rz(2.300793) q[0];
rz(-2.1127545) q[1];
sx q[1];
rz(-0.40887555) q[1];
sx q[1];
rz(1.4604481) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.235229) q[0];
sx q[0];
rz(-1.1113941) q[0];
sx q[0];
rz(1.4109341) q[0];
rz(-pi) q[1];
rz(-1.9424136) q[2];
sx q[2];
rz(-2.016168) q[2];
sx q[2];
rz(-2.1567313) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.0277623) q[1];
sx q[1];
rz(-0.9200458) q[1];
sx q[1];
rz(-1.9939994) q[1];
rz(3.06006) q[3];
sx q[3];
rz(-2.53269) q[3];
sx q[3];
rz(2.1341153) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.2568405) q[2];
sx q[2];
rz(-1.931087) q[2];
sx q[2];
rz(2.9849226) q[2];
rz(-1.5001851) q[3];
sx q[3];
rz(-1.2431966) q[3];
sx q[3];
rz(-0.18178864) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0665862) q[0];
sx q[0];
rz(-1.8445419) q[0];
sx q[0];
rz(-0.080667607) q[0];
rz(2.0196041) q[1];
sx q[1];
rz(-0.64067084) q[1];
sx q[1];
rz(1.5871619) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4951404) q[0];
sx q[0];
rz(-2.3123154) q[0];
sx q[0];
rz(-0.089228169) q[0];
x q[1];
rz(1.6115613) q[2];
sx q[2];
rz(-1.6869378) q[2];
sx q[2];
rz(-2.9774506) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.5579368) q[1];
sx q[1];
rz(-2.2018345) q[1];
sx q[1];
rz(-0.93871745) q[1];
x q[2];
rz(-2.6770079) q[3];
sx q[3];
rz(-0.66028336) q[3];
sx q[3];
rz(-0.44251501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.88177219) q[2];
sx q[2];
rz(-1.0910923) q[2];
sx q[2];
rz(1.0238949) q[2];
rz(-2.3648868) q[3];
sx q[3];
rz(-1.9909765) q[3];
sx q[3];
rz(-1.0030494) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7316932) q[0];
sx q[0];
rz(-0.85407805) q[0];
sx q[0];
rz(-0.62752974) q[0];
rz(2.4420786) q[1];
sx q[1];
rz(-2.1414089) q[1];
sx q[1];
rz(2.2342009) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2058682) q[0];
sx q[0];
rz(-0.75399929) q[0];
sx q[0];
rz(-3.0831017) q[0];
x q[1];
rz(0.19030119) q[2];
sx q[2];
rz(-1.4138599) q[2];
sx q[2];
rz(-0.67853329) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.6945362) q[1];
sx q[1];
rz(-1.2061822) q[1];
sx q[1];
rz(-1.0178119) q[1];
rz(-1.6697407) q[3];
sx q[3];
rz(-2.6032748) q[3];
sx q[3];
rz(0.42699277) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0044535) q[2];
sx q[2];
rz(-1.8567825) q[2];
sx q[2];
rz(2.2966906) q[2];
rz(-0.6984624) q[3];
sx q[3];
rz(-0.74028492) q[3];
sx q[3];
rz(-2.3499878) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38248211) q[0];
sx q[0];
rz(-1.4610721) q[0];
sx q[0];
rz(1.8735029) q[0];
rz(1.5269439) q[1];
sx q[1];
rz(-2.0497649) q[1];
sx q[1];
rz(2.7472034) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6170664) q[0];
sx q[0];
rz(-1.5252648) q[0];
sx q[0];
rz(0.92638735) q[0];
rz(-pi) q[1];
x q[1];
rz(0.73494786) q[2];
sx q[2];
rz(-1.6326666) q[2];
sx q[2];
rz(2.8501373) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.1413381) q[1];
sx q[1];
rz(-0.91857498) q[1];
sx q[1];
rz(-1.5144996) q[1];
rz(-pi) q[2];
rz(-1.6485571) q[3];
sx q[3];
rz(-1.9304515) q[3];
sx q[3];
rz(0.80870562) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.7093198) q[2];
sx q[2];
rz(-3.0797112) q[2];
sx q[2];
rz(-0.56149948) q[2];
rz(-1.1310486) q[3];
sx q[3];
rz(-0.80315042) q[3];
sx q[3];
rz(-2.6530182) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5973709) q[0];
sx q[0];
rz(-2.9591296) q[0];
sx q[0];
rz(0.23319787) q[0];
rz(-0.11416642) q[1];
sx q[1];
rz(-0.79008594) q[1];
sx q[1];
rz(-2.9679969) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9812935) q[0];
sx q[0];
rz(-2.3620785) q[0];
sx q[0];
rz(-2.1008089) q[0];
rz(-pi) q[1];
x q[1];
rz(0.87145135) q[2];
sx q[2];
rz(-2.5837499) q[2];
sx q[2];
rz(-2.7047472) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.6919903) q[1];
sx q[1];
rz(-0.27562818) q[1];
sx q[1];
rz(-2.2066433) q[1];
x q[2];
rz(-1.5585445) q[3];
sx q[3];
rz(-1.3821332) q[3];
sx q[3];
rz(1.0808627) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.1058098) q[2];
sx q[2];
rz(-2.1828987) q[2];
sx q[2];
rz(2.2868273) q[2];
rz(-3.0586976) q[3];
sx q[3];
rz(-1.9708743) q[3];
sx q[3];
rz(0.054072592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4891124) q[0];
sx q[0];
rz(-1.880045) q[0];
sx q[0];
rz(1.1789119) q[0];
rz(1.0014125) q[1];
sx q[1];
rz(-1.8779571) q[1];
sx q[1];
rz(0.15403919) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6545367) q[0];
sx q[0];
rz(-0.77870071) q[0];
sx q[0];
rz(2.6683183) q[0];
rz(-0.70239046) q[2];
sx q[2];
rz(-0.99961126) q[2];
sx q[2];
rz(-1.0016425) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.6012965) q[1];
sx q[1];
rz(-2.5754693) q[1];
sx q[1];
rz(-0.98928063) q[1];
rz(-pi) q[2];
rz(-1.3074806) q[3];
sx q[3];
rz(-1.2905408) q[3];
sx q[3];
rz(1.6681886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.4607294) q[2];
sx q[2];
rz(-1.0793842) q[2];
sx q[2];
rz(0.66780773) q[2];
rz(1.4051416) q[3];
sx q[3];
rz(-1.7738155) q[3];
sx q[3];
rz(-0.63647979) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3619096) q[0];
sx q[0];
rz(-1.4754262) q[0];
sx q[0];
rz(2.5478126) q[0];
rz(1.8608015) q[1];
sx q[1];
rz(-2.180876) q[1];
sx q[1];
rz(-1.8437754) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0827328) q[0];
sx q[0];
rz(-1.4987336) q[0];
sx q[0];
rz(2.5565992) q[0];
rz(-2.5137481) q[2];
sx q[2];
rz(-1.825001) q[2];
sx q[2];
rz(-2.0275627) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.62420708) q[1];
sx q[1];
rz(-2.7605857) q[1];
sx q[1];
rz(-0.34492774) q[1];
rz(1.7421977) q[3];
sx q[3];
rz(-2.3706782) q[3];
sx q[3];
rz(2.7293929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.39591509) q[2];
sx q[2];
rz(-3.119645) q[2];
sx q[2];
rz(1.5094666) q[2];
rz(3.0294026) q[3];
sx q[3];
rz(-1.0271007) q[3];
sx q[3];
rz(1.7621015) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6701732) q[0];
sx q[0];
rz(-1.0059953) q[0];
sx q[0];
rz(-0.26563409) q[0];
rz(-0.98948014) q[1];
sx q[1];
rz(-1.8786636) q[1];
sx q[1];
rz(2.8094453) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0672537) q[0];
sx q[0];
rz(-1.566972) q[0];
sx q[0];
rz(-1.9042364) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0004326) q[2];
sx q[2];
rz(-0.62295914) q[2];
sx q[2];
rz(-1.0819544) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0451584) q[1];
sx q[1];
rz(-1.2306884) q[1];
sx q[1];
rz(3.0902658) q[1];
x q[2];
rz(-0.048236851) q[3];
sx q[3];
rz(-1.3833481) q[3];
sx q[3];
rz(1.6316044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.0997448) q[2];
sx q[2];
rz(-2.0691278) q[2];
sx q[2];
rz(2.9912046) q[2];
rz(-0.8485052) q[3];
sx q[3];
rz(-0.34981194) q[3];
sx q[3];
rz(-1.295804) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0583508) q[0];
sx q[0];
rz(-1.8235089) q[0];
sx q[0];
rz(1.9705809) q[0];
rz(-1.3311483) q[1];
sx q[1];
rz(-2.34453) q[1];
sx q[1];
rz(-1.1042368) q[1];
rz(1.642557) q[2];
sx q[2];
rz(-2.4007779) q[2];
sx q[2];
rz(3.130198) q[2];
rz(-2.4348197) q[3];
sx q[3];
rz(-1.5512244) q[3];
sx q[3];
rz(2.649818) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
