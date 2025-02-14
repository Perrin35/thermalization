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
rz(0.19652551) q[0];
sx q[0];
rz(-0.82634574) q[0];
sx q[0];
rz(0.91808051) q[0];
rz(3.0300568) q[1];
sx q[1];
rz(-1.7938951) q[1];
sx q[1];
rz(-1.5674051) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2409256) q[0];
sx q[0];
rz(-1.5673715) q[0];
sx q[0];
rz(-1.5802556) q[0];
rz(-0.55337535) q[2];
sx q[2];
rz(-2.5042483) q[2];
sx q[2];
rz(3.1269249) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.35167978) q[1];
sx q[1];
rz(-1.4654298) q[1];
sx q[1];
rz(-1.2672474) q[1];
rz(-1.449578) q[3];
sx q[3];
rz(-1.08076) q[3];
sx q[3];
rz(-0.13257172) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.427318) q[2];
sx q[2];
rz(-1.6552507) q[2];
sx q[2];
rz(1.2662668) q[2];
rz(1.0424987) q[3];
sx q[3];
rz(-1.6008585) q[3];
sx q[3];
rz(1.7460167) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.105044) q[0];
sx q[0];
rz(-2.5140913) q[0];
sx q[0];
rz(-0.096916048) q[0];
rz(0.80822432) q[1];
sx q[1];
rz(-2.7327635) q[1];
sx q[1];
rz(1.2867297) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5135749) q[0];
sx q[0];
rz(-2.6903408) q[0];
sx q[0];
rz(1.4892764) q[0];
rz(-2.419728) q[2];
sx q[2];
rz(-1.1099735) q[2];
sx q[2];
rz(-2.6913655) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.3073439) q[1];
sx q[1];
rz(-1.9527702) q[1];
sx q[1];
rz(0.81484199) q[1];
rz(-pi) q[2];
rz(0.80974354) q[3];
sx q[3];
rz(-0.47172037) q[3];
sx q[3];
rz(0.37744409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.6404932) q[2];
sx q[2];
rz(-2.2576136) q[2];
sx q[2];
rz(-0.38197771) q[2];
rz(-0.52465087) q[3];
sx q[3];
rz(-1.7347521) q[3];
sx q[3];
rz(2.9978571) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3749738) q[0];
sx q[0];
rz(-1.5856278) q[0];
sx q[0];
rz(2.4554456) q[0];
rz(1.067591) q[1];
sx q[1];
rz(-2.144404) q[1];
sx q[1];
rz(0.080726191) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2093182) q[0];
sx q[0];
rz(-2.9146412) q[0];
sx q[0];
rz(1.7448241) q[0];
x q[1];
rz(-0.65004316) q[2];
sx q[2];
rz(-1.0299242) q[2];
sx q[2];
rz(-0.75640565) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.9517652) q[1];
sx q[1];
rz(-1.8009342) q[1];
sx q[1];
rz(-1.5469502) q[1];
rz(-pi) q[2];
x q[2];
rz(0.78781294) q[3];
sx q[3];
rz(-1.2255049) q[3];
sx q[3];
rz(-2.6097176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.2637691) q[2];
sx q[2];
rz(-0.70085415) q[2];
sx q[2];
rz(-2.7027255) q[2];
rz(1.2980488) q[3];
sx q[3];
rz(-2.7545007) q[3];
sx q[3];
rz(2.7264061) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91442672) q[0];
sx q[0];
rz(-2.5974847) q[0];
sx q[0];
rz(-0.67657226) q[0];
rz(-3.0034972) q[1];
sx q[1];
rz(-1.0905677) q[1];
sx q[1];
rz(2.9409883) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.706079) q[0];
sx q[0];
rz(-1.9301751) q[0];
sx q[0];
rz(0.91889221) q[0];
rz(-2.8851231) q[2];
sx q[2];
rz(-0.64531981) q[2];
sx q[2];
rz(0.27094524) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.323206) q[1];
sx q[1];
rz(-2.1765059) q[1];
sx q[1];
rz(1.731864) q[1];
rz(0.5612026) q[3];
sx q[3];
rz(-1.4491014) q[3];
sx q[3];
rz(2.7336371) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.5033919) q[2];
sx q[2];
rz(-1.8663422) q[2];
sx q[2];
rz(-1.4208043) q[2];
rz(-0.11451379) q[3];
sx q[3];
rz(-1.6679461) q[3];
sx q[3];
rz(2.1421471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-0.81922174) q[0];
sx q[0];
rz(-2.350816) q[0];
sx q[0];
rz(0.97212273) q[0];
rz(-2.8179893) q[1];
sx q[1];
rz(-1.4515896) q[1];
sx q[1];
rz(-1.1824898) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0111368) q[0];
sx q[0];
rz(-0.19936518) q[0];
sx q[0];
rz(1.7749271) q[0];
x q[1];
rz(-0.57020541) q[2];
sx q[2];
rz(-1.5758762) q[2];
sx q[2];
rz(0.069381086) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2731635) q[1];
sx q[1];
rz(-1.4123962) q[1];
sx q[1];
rz(-3.0636154) q[1];
rz(-pi) q[2];
rz(0.72850169) q[3];
sx q[3];
rz(-0.71069169) q[3];
sx q[3];
rz(-1.817734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.16367308) q[2];
sx q[2];
rz(-1.7719496) q[2];
sx q[2];
rz(0.53691205) q[2];
rz(-0.72530693) q[3];
sx q[3];
rz(-3.091843) q[3];
sx q[3];
rz(0.79881001) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-0.86730114) q[0];
sx q[0];
rz(-2.0766356) q[0];
sx q[0];
rz(-1.4428447) q[0];
rz(-0.2105712) q[1];
sx q[1];
rz(-1.6981533) q[1];
sx q[1];
rz(2.494536) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66523401) q[0];
sx q[0];
rz(-1.826735) q[0];
sx q[0];
rz(0.38689918) q[0];
x q[1];
rz(-1.4616632) q[2];
sx q[2];
rz(-2.1246548) q[2];
sx q[2];
rz(2.6719029) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(3.0004859) q[1];
sx q[1];
rz(-2.5865411) q[1];
sx q[1];
rz(-0.6935814) q[1];
rz(-pi) q[2];
x q[2];
rz(1.296455) q[3];
sx q[3];
rz(-1.2452494) q[3];
sx q[3];
rz(2.5485768) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.0183384) q[2];
sx q[2];
rz(-1.0945357) q[2];
sx q[2];
rz(1.1798165) q[2];
rz(-1.8566462) q[3];
sx q[3];
rz(-0.40950567) q[3];
sx q[3];
rz(0.063260945) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1061851) q[0];
sx q[0];
rz(-1.6949061) q[0];
sx q[0];
rz(-2.2440198) q[0];
rz(-1.3342185) q[1];
sx q[1];
rz(-1.370627) q[1];
sx q[1];
rz(0.99475494) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.008381) q[0];
sx q[0];
rz(-1.8803839) q[0];
sx q[0];
rz(-0.82010834) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.45539664) q[2];
sx q[2];
rz(-1.0173305) q[2];
sx q[2];
rz(1.139251) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.23235944) q[1];
sx q[1];
rz(-2.5933752) q[1];
sx q[1];
rz(-2.1562804) q[1];
rz(3.1174421) q[3];
sx q[3];
rz(-1.2901879) q[3];
sx q[3];
rz(2.5606972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(3.0576386) q[2];
sx q[2];
rz(-1.8123764) q[2];
sx q[2];
rz(-0.29339054) q[2];
rz(3.0883664) q[3];
sx q[3];
rz(-2.2727727) q[3];
sx q[3];
rz(1.4330385) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5695802) q[0];
sx q[0];
rz(-1.7520289) q[0];
sx q[0];
rz(-2.8386175) q[0];
rz(1.4888034) q[1];
sx q[1];
rz(-1.3158512) q[1];
sx q[1];
rz(-0.66910076) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6953872) q[0];
sx q[0];
rz(-1.2406557) q[0];
sx q[0];
rz(2.6690525) q[0];
x q[1];
rz(1.2190518) q[2];
sx q[2];
rz(-1.3430077) q[2];
sx q[2];
rz(-2.5643189) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.41324612) q[1];
sx q[1];
rz(-1.9921682) q[1];
sx q[1];
rz(2.1284118) q[1];
rz(-pi) q[2];
rz(1.8101997) q[3];
sx q[3];
rz(-1.1576011) q[3];
sx q[3];
rz(-1.5904782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.1423433) q[2];
sx q[2];
rz(-1.8084904) q[2];
sx q[2];
rz(0.86714253) q[2];
rz(-2.0182746) q[3];
sx q[3];
rz(-2.2893548) q[3];
sx q[3];
rz(1.1423906) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2808696) q[0];
sx q[0];
rz(-2.6779802) q[0];
sx q[0];
rz(0.11775693) q[0];
rz(2.2640758) q[1];
sx q[1];
rz(-1.0209457) q[1];
sx q[1];
rz(-0.95474517) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3316393) q[0];
sx q[0];
rz(-2.3970986) q[0];
sx q[0];
rz(-2.5622423) q[0];
rz(-pi) q[1];
rz(0.81127848) q[2];
sx q[2];
rz(-0.64006058) q[2];
sx q[2];
rz(2.0612353) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.92312557) q[1];
sx q[1];
rz(-2.4187741) q[1];
sx q[1];
rz(1.2413526) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.84857984) q[3];
sx q[3];
rz(-0.89451087) q[3];
sx q[3];
rz(-3.0885895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.2473601) q[2];
sx q[2];
rz(-2.8017513) q[2];
sx q[2];
rz(-1.4840688) q[2];
rz(-1.6190716) q[3];
sx q[3];
rz(-1.7584453) q[3];
sx q[3];
rz(-1.9671666) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5840983) q[0];
sx q[0];
rz(-1.7099986) q[0];
sx q[0];
rz(0.75310055) q[0];
rz(-1.151471) q[1];
sx q[1];
rz(-2.617372) q[1];
sx q[1];
rz(1.5392736) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.879047) q[0];
sx q[0];
rz(-2.6521054) q[0];
sx q[0];
rz(-1.1734084) q[0];
rz(-pi) q[1];
rz(-0.87290092) q[2];
sx q[2];
rz(-1.1660327) q[2];
sx q[2];
rz(1.0412316) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.11048296) q[1];
sx q[1];
rz(-1.3785161) q[1];
sx q[1];
rz(-2.7483205) q[1];
x q[2];
rz(-0.33892858) q[3];
sx q[3];
rz(-2.458161) q[3];
sx q[3];
rz(0.11250699) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.3573542) q[2];
sx q[2];
rz(-2.9604762) q[2];
sx q[2];
rz(2.6757346) q[2];
rz(-2.0458131) q[3];
sx q[3];
rz(-2.1491137) q[3];
sx q[3];
rz(0.54212681) q[3];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.096238484) q[0];
sx q[0];
rz(-1.5662554) q[0];
sx q[0];
rz(1.5731496) q[0];
rz(2.0737598) q[1];
sx q[1];
rz(-2.7012431) q[1];
sx q[1];
rz(2.3213097) q[1];
rz(0.26366641) q[2];
sx q[2];
rz(-1.7963526) q[2];
sx q[2];
rz(0.11033478) q[2];
rz(-2.3798306) q[3];
sx q[3];
rz(-2.4488425) q[3];
sx q[3];
rz(-0.52656534) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
