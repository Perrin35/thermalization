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
rz(1.6990868) q[0];
sx q[0];
rz(-1.8449755) q[0];
sx q[0];
rz(-1.698864) q[0];
rz(-0.35248414) q[1];
sx q[1];
rz(-0.52630693) q[1];
sx q[1];
rz(-1.3679282) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5718846) q[0];
sx q[0];
rz(-1.8284599) q[0];
sx q[0];
rz(0.37057693) q[0];
rz(-pi) q[1];
x q[1];
rz(0.35570154) q[2];
sx q[2];
rz(-1.5354325) q[2];
sx q[2];
rz(-0.48330467) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.4579205) q[1];
sx q[1];
rz(-1.0506609) q[1];
sx q[1];
rz(-1.7207654) q[1];
rz(-pi) q[2];
x q[2];
rz(2.164293) q[3];
sx q[3];
rz(-2.233089) q[3];
sx q[3];
rz(0.73847258) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.0679396) q[2];
sx q[2];
rz(-1.0509793) q[2];
sx q[2];
rz(1.2037207) q[2];
rz(-2.1667571) q[3];
sx q[3];
rz(-2.1741512) q[3];
sx q[3];
rz(-1.2459285) q[3];
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
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75342733) q[0];
sx q[0];
rz(-2.0797256) q[0];
sx q[0];
rz(-2.2573096) q[0];
rz(0.70941225) q[1];
sx q[1];
rz(-1.8953036) q[1];
sx q[1];
rz(2.2394004) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1195148) q[0];
sx q[0];
rz(-1.365322) q[0];
sx q[0];
rz(-2.6809262) q[0];
rz(-pi) q[1];
rz(0.10295094) q[2];
sx q[2];
rz(-2.1753417) q[2];
sx q[2];
rz(-1.9905195) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.1556405) q[1];
sx q[1];
rz(-1.2725755) q[1];
sx q[1];
rz(1.8267426) q[1];
x q[2];
rz(3.0608474) q[3];
sx q[3];
rz(-1.7087548) q[3];
sx q[3];
rz(-1.6899275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.0238637) q[2];
sx q[2];
rz(-0.29568299) q[2];
sx q[2];
rz(1.6531061) q[2];
rz(1.214341) q[3];
sx q[3];
rz(-1.6818654) q[3];
sx q[3];
rz(-1.2795718) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0073256) q[0];
sx q[0];
rz(-1.7420344) q[0];
sx q[0];
rz(0.16635995) q[0];
rz(0.69794377) q[1];
sx q[1];
rz(-2.3614387) q[1];
sx q[1];
rz(-1.4770329) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.03334) q[0];
sx q[0];
rz(-1.1287924) q[0];
sx q[0];
rz(0.0033133034) q[0];
rz(-1.2899621) q[2];
sx q[2];
rz(-1.8321494) q[2];
sx q[2];
rz(-2.7203943) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.9117457) q[1];
sx q[1];
rz(-1.032879) q[1];
sx q[1];
rz(-0.53373806) q[1];
x q[2];
rz(-0.03607492) q[3];
sx q[3];
rz(-2.6429308) q[3];
sx q[3];
rz(2.7447678) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.4878238) q[2];
sx q[2];
rz(-1.1710125) q[2];
sx q[2];
rz(0.6146532) q[2];
rz(-1.4608308) q[3];
sx q[3];
rz(-1.5644282) q[3];
sx q[3];
rz(-3.0998668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5322658) q[0];
sx q[0];
rz(-1.8121239) q[0];
sx q[0];
rz(1.5186658) q[0];
rz(-0.80865771) q[1];
sx q[1];
rz(-1.2963908) q[1];
sx q[1];
rz(0.048695806) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0322389) q[0];
sx q[0];
rz(-1.236602) q[0];
sx q[0];
rz(2.4346057) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7596731) q[2];
sx q[2];
rz(-1.6798229) q[2];
sx q[2];
rz(-0.38709059) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.21203707) q[1];
sx q[1];
rz(-2.0494645) q[1];
sx q[1];
rz(2.5655866) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0297431) q[3];
sx q[3];
rz(-2.4623639) q[3];
sx q[3];
rz(1.3857057) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.101717) q[2];
sx q[2];
rz(-0.81284916) q[2];
sx q[2];
rz(1.2561049) q[2];
rz(3.0806372) q[3];
sx q[3];
rz(-0.90164369) q[3];
sx q[3];
rz(0.92756334) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
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
rz(-2.1561279) q[0];
sx q[0];
rz(-1.1277072) q[0];
sx q[0];
rz(-0.10966478) q[0];
rz(2.1127286) q[1];
sx q[1];
rz(-1.8548465) q[1];
sx q[1];
rz(-1.5922155) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4745195) q[0];
sx q[0];
rz(-1.1231849) q[0];
sx q[0];
rz(1.8846798) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1774753) q[2];
sx q[2];
rz(-2.5543382) q[2];
sx q[2];
rz(2.6850785) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.2031496) q[1];
sx q[1];
rz(-2.0324273) q[1];
sx q[1];
rz(1.0814971) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5331334) q[3];
sx q[3];
rz(-1.5364913) q[3];
sx q[3];
rz(-1.487285) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.5900383) q[2];
sx q[2];
rz(-2.8503214) q[2];
sx q[2];
rz(2.6596098) q[2];
rz(-0.41989741) q[3];
sx q[3];
rz(-1.0651383) q[3];
sx q[3];
rz(2.4317252) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9919306) q[0];
sx q[0];
rz(-2.2580632) q[0];
sx q[0];
rz(-2.4290207) q[0];
rz(2.271671) q[1];
sx q[1];
rz(-2.7674119) q[1];
sx q[1];
rz(0.023699997) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11940322) q[0];
sx q[0];
rz(-1.7339804) q[0];
sx q[0];
rz(-0.52363445) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.65448241) q[2];
sx q[2];
rz(-0.81771446) q[2];
sx q[2];
rz(0.47374642) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.5710793) q[1];
sx q[1];
rz(-0.16579443) q[1];
sx q[1];
rz(1.5618663) q[1];
x q[2];
rz(0.48051254) q[3];
sx q[3];
rz(-1.7578205) q[3];
sx q[3];
rz(-2.9399237) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.3370543) q[2];
sx q[2];
rz(-2.4961175) q[2];
sx q[2];
rz(-2.0013334) q[2];
rz(-1.1537457) q[3];
sx q[3];
rz(-1.9269582) q[3];
sx q[3];
rz(1.3212737) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.425728) q[0];
sx q[0];
rz(-1.3874929) q[0];
sx q[0];
rz(2.7737889) q[0];
rz(-1.6416719) q[1];
sx q[1];
rz(-2.2190614) q[1];
sx q[1];
rz(-1.0949162) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85830583) q[0];
sx q[0];
rz(-2.0525161) q[0];
sx q[0];
rz(2.879309) q[0];
x q[1];
rz(-1.7000654) q[2];
sx q[2];
rz(-2.1174927) q[2];
sx q[2];
rz(2.0055637) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.8290295) q[1];
sx q[1];
rz(-2.9051648) q[1];
sx q[1];
rz(-0.076024012) q[1];
x q[2];
rz(-1.1672375) q[3];
sx q[3];
rz(-1.5903683) q[3];
sx q[3];
rz(0.38364601) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0992004) q[2];
sx q[2];
rz(-1.7902057) q[2];
sx q[2];
rz(-0.055709906) q[2];
rz(-2.4646711) q[3];
sx q[3];
rz(-0.61546314) q[3];
sx q[3];
rz(-2.2630579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1894576) q[0];
sx q[0];
rz(-0.5624693) q[0];
sx q[0];
rz(-0.96543717) q[0];
rz(1.2646487) q[1];
sx q[1];
rz(-0.74600428) q[1];
sx q[1];
rz(0.3269349) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33990369) q[0];
sx q[0];
rz(-1.9205695) q[0];
sx q[0];
rz(-3.0236493) q[0];
x q[1];
rz(1.5599653) q[2];
sx q[2];
rz(-1.1248661) q[2];
sx q[2];
rz(1.2496834) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.2911243) q[1];
sx q[1];
rz(-1.7479436) q[1];
sx q[1];
rz(1.7825104) q[1];
rz(-pi) q[2];
rz(-3.0013004) q[3];
sx q[3];
rz(-0.45278087) q[3];
sx q[3];
rz(-2.0425532) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.7586729) q[2];
sx q[2];
rz(-0.70432538) q[2];
sx q[2];
rz(-0.79603377) q[2];
rz(3.054079) q[3];
sx q[3];
rz(-0.60860601) q[3];
sx q[3];
rz(-0.93761888) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7079119) q[0];
sx q[0];
rz(-1.3734564) q[0];
sx q[0];
rz(2.8117836) q[0];
rz(3.0682796) q[1];
sx q[1];
rz(-2.4344756) q[1];
sx q[1];
rz(2.0475533) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0596704) q[0];
sx q[0];
rz(-1.1003033) q[0];
sx q[0];
rz(1.5289896) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4118926) q[2];
sx q[2];
rz(-0.4022809) q[2];
sx q[2];
rz(2.2251468) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.7787012) q[1];
sx q[1];
rz(-2.5769632) q[1];
sx q[1];
rz(2.8933011) q[1];
x q[2];
rz(-1.6570143) q[3];
sx q[3];
rz(-0.43282193) q[3];
sx q[3];
rz(-2.509887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0426992) q[2];
sx q[2];
rz(-2.3697479) q[2];
sx q[2];
rz(-2.0043376) q[2];
rz(-0.62816652) q[3];
sx q[3];
rz(-2.7558432) q[3];
sx q[3];
rz(2.1342733) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6657418) q[0];
sx q[0];
rz(-0.51097149) q[0];
sx q[0];
rz(-1.092859) q[0];
rz(1.2650371) q[1];
sx q[1];
rz(-1.3053514) q[1];
sx q[1];
rz(1.0337894) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2898258) q[0];
sx q[0];
rz(-0.74301592) q[0];
sx q[0];
rz(0.8351589) q[0];
rz(-pi) q[1];
rz(0.25569465) q[2];
sx q[2];
rz(-0.77633023) q[2];
sx q[2];
rz(1.9071554) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.46230832) q[1];
sx q[1];
rz(-0.95272118) q[1];
sx q[1];
rz(-1.0584153) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3218811) q[3];
sx q[3];
rz(-0.14302172) q[3];
sx q[3];
rz(1.7211357) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.25896245) q[2];
sx q[2];
rz(-2.4179103) q[2];
sx q[2];
rz(2.9222729) q[2];
rz(0.80260459) q[3];
sx q[3];
rz(-1.8290627) q[3];
sx q[3];
rz(-2.820211) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6312859) q[0];
sx q[0];
rz(-1.8230556) q[0];
sx q[0];
rz(-2.0429116) q[0];
rz(-2.9639099) q[1];
sx q[1];
rz(-2.1251353) q[1];
sx q[1];
rz(2.9716117) q[1];
rz(-2.2732757) q[2];
sx q[2];
rz(-1.6196031) q[2];
sx q[2];
rz(-2.6283787) q[2];
rz(-3.0232676) q[3];
sx q[3];
rz(-1.3814303) q[3];
sx q[3];
rz(-0.43683972) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
