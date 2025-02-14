OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.3761223) q[0];
sx q[0];
rz(-1.0785311) q[0];
sx q[0];
rz(-1.3753139) q[0];
rz(3.5612192) q[1];
sx q[1];
rz(1.761275) q[1];
sx q[1];
rz(7.0786369) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8574816) q[0];
sx q[0];
rz(-1.6499551) q[0];
sx q[0];
rz(-0.37374951) q[0];
rz(-pi) q[1];
x q[1];
rz(0.4001053) q[2];
sx q[2];
rz(-1.4608619) q[2];
sx q[2];
rz(0.38450228) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.7401984) q[1];
sx q[1];
rz(-0.90448442) q[1];
sx q[1];
rz(1.8224276) q[1];
rz(3.0555435) q[3];
sx q[3];
rz(-1.3121735) q[3];
sx q[3];
rz(-1.8515967) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.85806435) q[2];
sx q[2];
rz(-0.29164803) q[2];
sx q[2];
rz(2.7533599) q[2];
rz(0.20902108) q[3];
sx q[3];
rz(-1.4744604) q[3];
sx q[3];
rz(-0.89116636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6727869) q[0];
sx q[0];
rz(-0.46872941) q[0];
sx q[0];
rz(0.76954532) q[0];
rz(-2.1107213) q[1];
sx q[1];
rz(-2.2831235) q[1];
sx q[1];
rz(1.4268202) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.010410943) q[0];
sx q[0];
rz(-2.7658434) q[0];
sx q[0];
rz(-0.6300169) q[0];
rz(1.6278817) q[2];
sx q[2];
rz(-2.2453614) q[2];
sx q[2];
rz(-3.0115173) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0116358) q[1];
sx q[1];
rz(-2.4677688) q[1];
sx q[1];
rz(-0.41566276) q[1];
rz(-pi) q[2];
rz(1.6321502) q[3];
sx q[3];
rz(-2.4262145) q[3];
sx q[3];
rz(2.3917907) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.1398805) q[2];
sx q[2];
rz(-0.44846815) q[2];
sx q[2];
rz(-1.2343538) q[2];
rz(-2.6185696) q[3];
sx q[3];
rz(-1.1250027) q[3];
sx q[3];
rz(-2.4856868) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6958385) q[0];
sx q[0];
rz(-3.1390751) q[0];
sx q[0];
rz(-0.54760951) q[0];
rz(-2.2721263) q[1];
sx q[1];
rz(-2.165386) q[1];
sx q[1];
rz(-1.4617823) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8903046) q[0];
sx q[0];
rz(-1.1532532) q[0];
sx q[0];
rz(1.1956716) q[0];
x q[1];
rz(2.8732576) q[2];
sx q[2];
rz(-0.78040022) q[2];
sx q[2];
rz(0.66106272) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.74341575) q[1];
sx q[1];
rz(-1.2733439) q[1];
sx q[1];
rz(-2.0259845) q[1];
rz(-pi) q[2];
rz(-2.2814155) q[3];
sx q[3];
rz(-1.0499991) q[3];
sx q[3];
rz(1.2788683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.72472858) q[2];
sx q[2];
rz(-0.29837307) q[2];
sx q[2];
rz(-0.63428632) q[2];
rz(1.0861081) q[3];
sx q[3];
rz(-2.6606798) q[3];
sx q[3];
rz(2.0079131) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5241765) q[0];
sx q[0];
rz(-2.9845181) q[0];
sx q[0];
rz(1.5217391) q[0];
rz(-2.1583083) q[1];
sx q[1];
rz(-2.0359928) q[1];
sx q[1];
rz(0.081667893) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8056792) q[0];
sx q[0];
rz(-1.4594074) q[0];
sx q[0];
rz(-2.1301444) q[0];
rz(-pi) q[1];
rz(3.1278243) q[2];
sx q[2];
rz(-1.2377316) q[2];
sx q[2];
rz(-0.027160732) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.5809243) q[1];
sx q[1];
rz(-2.4663975) q[1];
sx q[1];
rz(-2.2972754) q[1];
rz(-pi) q[2];
rz(3.087864) q[3];
sx q[3];
rz(-1.0815291) q[3];
sx q[3];
rz(0.51988822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.6433158) q[2];
sx q[2];
rz(-1.6643107) q[2];
sx q[2];
rz(-0.91900438) q[2];
rz(2.8194341) q[3];
sx q[3];
rz(-1.6291658) q[3];
sx q[3];
rz(2.7951516) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.084694013) q[0];
sx q[0];
rz(-2.0094805) q[0];
sx q[0];
rz(1.3380916) q[0];
rz(2.7390506) q[1];
sx q[1];
rz(-1.2502547) q[1];
sx q[1];
rz(1.58443) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8020222) q[0];
sx q[0];
rz(-2.5568107) q[0];
sx q[0];
rz(2.0680379) q[0];
x q[1];
rz(0.7942368) q[2];
sx q[2];
rz(-1.4068687) q[2];
sx q[2];
rz(-1.9203223) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.8338523) q[1];
sx q[1];
rz(-2.5870812) q[1];
sx q[1];
rz(-2.8608198) q[1];
rz(2.8731737) q[3];
sx q[3];
rz(-2.0257334) q[3];
sx q[3];
rz(2.0132298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.0255787) q[2];
sx q[2];
rz(-1.8639001) q[2];
sx q[2];
rz(-0.73520994) q[2];
rz(1.0544581) q[3];
sx q[3];
rz(-1.0337044) q[3];
sx q[3];
rz(-1.2247156) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9758107) q[0];
sx q[0];
rz(-2.7175792) q[0];
sx q[0];
rz(1.2531511) q[0];
rz(-1.5660628) q[1];
sx q[1];
rz(-1.1352481) q[1];
sx q[1];
rz(-2.6083686) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1178109) q[0];
sx q[0];
rz(-0.83988777) q[0];
sx q[0];
rz(2.8230142) q[0];
x q[1];
rz(1.7248956) q[2];
sx q[2];
rz(-2.0770501) q[2];
sx q[2];
rz(-2.7768898) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6937852) q[1];
sx q[1];
rz(-2.8485423) q[1];
sx q[1];
rz(-0.70290053) q[1];
rz(2.1142152) q[3];
sx q[3];
rz(-2.0680331) q[3];
sx q[3];
rz(-1.9774205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6553361) q[2];
sx q[2];
rz(-0.37898263) q[2];
sx q[2];
rz(2.3113225) q[2];
rz(-1.3116948) q[3];
sx q[3];
rz(-1.0887086) q[3];
sx q[3];
rz(2.8970498) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.526392) q[0];
sx q[0];
rz(-1.6561693) q[0];
sx q[0];
rz(-0.8335337) q[0];
rz(-2.1444881) q[1];
sx q[1];
rz(-1.3778967) q[1];
sx q[1];
rz(1.5901784) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1828022) q[0];
sx q[0];
rz(-1.1353462) q[0];
sx q[0];
rz(1.6531526) q[0];
x q[1];
rz(-1.0710414) q[2];
sx q[2];
rz(-2.4883371) q[2];
sx q[2];
rz(0.61966255) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.95495086) q[1];
sx q[1];
rz(-1.690505) q[1];
sx q[1];
rz(-0.11255944) q[1];
x q[2];
rz(1.3209913) q[3];
sx q[3];
rz(-2.0916998) q[3];
sx q[3];
rz(1.5977719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.3414063) q[2];
sx q[2];
rz(-1.032858) q[2];
sx q[2];
rz(-1.3979073) q[2];
rz(0.37311113) q[3];
sx q[3];
rz(-2.2171376) q[3];
sx q[3];
rz(1.0188518) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7031192) q[0];
sx q[0];
rz(-1.6699474) q[0];
sx q[0];
rz(-0.31914172) q[0];
rz(-1.9769662) q[1];
sx q[1];
rz(-1.9705557) q[1];
sx q[1];
rz(-3.107531) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0379763) q[0];
sx q[0];
rz(-1.2671794) q[0];
sx q[0];
rz(2.3453317) q[0];
x q[1];
rz(2.3715842) q[2];
sx q[2];
rz(-1.3427375) q[2];
sx q[2];
rz(-0.32021013) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.086255161) q[1];
sx q[1];
rz(-3.1172522) q[1];
sx q[1];
rz(1.1175977) q[1];
rz(-pi) q[2];
rz(-0.53197022) q[3];
sx q[3];
rz(-1.2571223) q[3];
sx q[3];
rz(1.6187291) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(3.0762735) q[2];
sx q[2];
rz(-1.7460145) q[2];
sx q[2];
rz(-0.070153959) q[2];
rz(-1.8326727) q[3];
sx q[3];
rz(-0.85246837) q[3];
sx q[3];
rz(-0.13516983) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.081130505) q[0];
sx q[0];
rz(-1.9996996) q[0];
sx q[0];
rz(-2.6666226) q[0];
rz(1.181107) q[1];
sx q[1];
rz(-2.2468552) q[1];
sx q[1];
rz(-0.917135) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1346996) q[0];
sx q[0];
rz(-1.0118766) q[0];
sx q[0];
rz(-1.2968282) q[0];
rz(-pi) q[1];
rz(-0.88240282) q[2];
sx q[2];
rz(-1.3254291) q[2];
sx q[2];
rz(0.74859017) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.9571984) q[1];
sx q[1];
rz(-2.9177114) q[1];
sx q[1];
rz(1.9489524) q[1];
x q[2];
rz(-2.6090129) q[3];
sx q[3];
rz(-0.6709698) q[3];
sx q[3];
rz(2.4155145) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.26943031) q[2];
sx q[2];
rz(-1.3698801) q[2];
sx q[2];
rz(-1.5838464) q[2];
rz(-0.15548429) q[3];
sx q[3];
rz(-1.0289861) q[3];
sx q[3];
rz(-2.1492929) q[3];
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
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6601722) q[0];
sx q[0];
rz(-1.6654797) q[0];
sx q[0];
rz(2.7760264) q[0];
rz(1.2008601) q[1];
sx q[1];
rz(-1.6114085) q[1];
sx q[1];
rz(1.0467122) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6088951) q[0];
sx q[0];
rz(-0.77975285) q[0];
sx q[0];
rz(-2.3391367) q[0];
rz(2.7002242) q[2];
sx q[2];
rz(-2.0725996) q[2];
sx q[2];
rz(2.8428915) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.2859058) q[1];
sx q[1];
rz(-1.0382723) q[1];
sx q[1];
rz(-2.3175815) q[1];
x q[2];
rz(-2.1165068) q[3];
sx q[3];
rz(-1.0574329) q[3];
sx q[3];
rz(1.8051763) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.6863579) q[2];
sx q[2];
rz(-1.1150259) q[2];
sx q[2];
rz(-2.4380747) q[2];
rz(-2.5840058) q[3];
sx q[3];
rz(-2.360354) q[3];
sx q[3];
rz(2.668146) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(-1.2932417) q[0];
sx q[0];
rz(-1.0307172) q[0];
sx q[0];
rz(0.95681554) q[0];
rz(-2.1059857) q[1];
sx q[1];
rz(-0.57301141) q[1];
sx q[1];
rz(1.116629) q[1];
rz(1.0445486) q[2];
sx q[2];
rz(-2.1696354) q[2];
sx q[2];
rz(-1.3641691) q[2];
rz(-2.9530408) q[3];
sx q[3];
rz(-0.22123541) q[3];
sx q[3];
rz(0.60047022) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
