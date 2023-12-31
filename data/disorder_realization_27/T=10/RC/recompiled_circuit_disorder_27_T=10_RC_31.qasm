OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(3.1383837) q[0];
sx q[0];
rz(-2.9870343) q[0];
sx q[0];
rz(2.4490693) q[0];
rz(-1.2094296) q[1];
sx q[1];
rz(-1.8930607) q[1];
sx q[1];
rz(-1.7564397) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9557578) q[0];
sx q[0];
rz(-1.251936) q[0];
sx q[0];
rz(-2.3715109) q[0];
rz(-pi) q[1];
rz(2.5897883) q[2];
sx q[2];
rz(-1.8066415) q[2];
sx q[2];
rz(-2.9896196) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.9618122) q[1];
sx q[1];
rz(-1.4084067) q[1];
sx q[1];
rz(-0.23602545) q[1];
rz(-pi) q[2];
rz(2.5296506) q[3];
sx q[3];
rz(-2.3903333) q[3];
sx q[3];
rz(-0.9179759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.2549071) q[2];
sx q[2];
rz(-0.79780769) q[2];
sx q[2];
rz(0.20516667) q[2];
rz(-0.77130476) q[3];
sx q[3];
rz(-0.78273928) q[3];
sx q[3];
rz(1.1024968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7339864) q[0];
sx q[0];
rz(-2.3953231) q[0];
sx q[0];
rz(0.45390391) q[0];
rz(-2.1167963) q[1];
sx q[1];
rz(-0.4075993) q[1];
sx q[1];
rz(1.9143547) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42392143) q[0];
sx q[0];
rz(-0.14980355) q[0];
sx q[0];
rz(-2.1013837) q[0];
rz(-2.0298376) q[2];
sx q[2];
rz(-1.0978205) q[2];
sx q[2];
rz(0.74726653) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.2482359) q[1];
sx q[1];
rz(-1.3125988) q[1];
sx q[1];
rz(-2.8398819) q[1];
x q[2];
rz(1.4144054) q[3];
sx q[3];
rz(-1.7692493) q[3];
sx q[3];
rz(0.95201492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.041302117) q[2];
sx q[2];
rz(-1.1854478) q[2];
sx q[2];
rz(-2.5741637) q[2];
rz(0.36519095) q[3];
sx q[3];
rz(-1.7285715) q[3];
sx q[3];
rz(-2.173483) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.658618) q[0];
sx q[0];
rz(-2.5768319) q[0];
sx q[0];
rz(0.89865249) q[0];
rz(0.99575106) q[1];
sx q[1];
rz(-1.5581222) q[1];
sx q[1];
rz(-2.8083037) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62045287) q[0];
sx q[0];
rz(-2.1846909) q[0];
sx q[0];
rz(-0.99434538) q[0];
rz(-pi) q[1];
rz(1.7477481) q[2];
sx q[2];
rz(-1.4743917) q[2];
sx q[2];
rz(2.5071438) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.2644314) q[1];
sx q[1];
rz(-2.1794771) q[1];
sx q[1];
rz(2.9560637) q[1];
rz(-pi) q[2];
rz(2.7379052) q[3];
sx q[3];
rz(-1.0567769) q[3];
sx q[3];
rz(-2.6727563) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.4553392) q[2];
sx q[2];
rz(-1.3505961) q[2];
sx q[2];
rz(1.1509482) q[2];
rz(-2.3006556) q[3];
sx q[3];
rz(-1.1154113) q[3];
sx q[3];
rz(-1.1545198) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6999917) q[0];
sx q[0];
rz(-1.5804407) q[0];
sx q[0];
rz(-0.68471318) q[0];
rz(-1.0355863) q[1];
sx q[1];
rz(-0.5077478) q[1];
sx q[1];
rz(1.9365786) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3048153) q[0];
sx q[0];
rz(-2.4495821) q[0];
sx q[0];
rz(-1.5185604) q[0];
rz(-pi) q[1];
rz(-0.92653805) q[2];
sx q[2];
rz(-1.7346003) q[2];
sx q[2];
rz(1.0545727) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.64441427) q[1];
sx q[1];
rz(-2.2481611) q[1];
sx q[1];
rz(2.7888984) q[1];
x q[2];
rz(0.86592309) q[3];
sx q[3];
rz(-1.5161627) q[3];
sx q[3];
rz(1.6161402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.8923607) q[2];
sx q[2];
rz(-1.7148596) q[2];
sx q[2];
rz(2.7704346) q[2];
rz(-1.4012198) q[3];
sx q[3];
rz(-2.4818082) q[3];
sx q[3];
rz(2.0223117) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.086833) q[0];
sx q[0];
rz(-2.355447) q[0];
sx q[0];
rz(-0.13312419) q[0];
rz(2.1482824) q[1];
sx q[1];
rz(-1.3860093) q[1];
sx q[1];
rz(-2.5865119) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.426794) q[0];
sx q[0];
rz(-1.5835276) q[0];
sx q[0];
rz(1.2554332) q[0];
rz(1.985717) q[2];
sx q[2];
rz(-1.6332111) q[2];
sx q[2];
rz(-1.2635363) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.3302147) q[1];
sx q[1];
rz(-1.1256309) q[1];
sx q[1];
rz(3.1203169) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5186148) q[3];
sx q[3];
rz(-1.6515886) q[3];
sx q[3];
rz(-2.3063456) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.83539) q[2];
sx q[2];
rz(-1.0058879) q[2];
sx q[2];
rz(0.13892697) q[2];
rz(-0.94240087) q[3];
sx q[3];
rz(-1.5706294) q[3];
sx q[3];
rz(-2.5901103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5979364) q[0];
sx q[0];
rz(-2.5718226) q[0];
sx q[0];
rz(2.561835) q[0];
rz(-3.014091) q[1];
sx q[1];
rz(-1.9516877) q[1];
sx q[1];
rz(1.6019843) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0149536) q[0];
sx q[0];
rz(-2.3987028) q[0];
sx q[0];
rz(-1.2680306) q[0];
rz(-pi) q[1];
rz(2.8069205) q[2];
sx q[2];
rz(-0.13813189) q[2];
sx q[2];
rz(-1.7931256) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.3547937) q[1];
sx q[1];
rz(-1.8537632) q[1];
sx q[1];
rz(1.1250886) q[1];
rz(-pi) q[2];
x q[2];
rz(0.29069889) q[3];
sx q[3];
rz(-0.906956) q[3];
sx q[3];
rz(-1.78474) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.55398983) q[2];
sx q[2];
rz(-0.25045276) q[2];
sx q[2];
rz(-0.26947752) q[2];
rz(-2.907471) q[3];
sx q[3];
rz(-0.52474371) q[3];
sx q[3];
rz(-0.060119303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44678974) q[0];
sx q[0];
rz(-0.61820522) q[0];
sx q[0];
rz(-2.457298) q[0];
rz(0.11958312) q[1];
sx q[1];
rz(-1.2922829) q[1];
sx q[1];
rz(-2.6228242) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51661922) q[0];
sx q[0];
rz(-2.3949351) q[0];
sx q[0];
rz(1.7685415) q[0];
rz(0.22612818) q[2];
sx q[2];
rz(-2.3580708) q[2];
sx q[2];
rz(-2.4353611) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.3541975) q[1];
sx q[1];
rz(-2.0320315) q[1];
sx q[1];
rz(0.72871491) q[1];
rz(-1.6325083) q[3];
sx q[3];
rz(-0.40137526) q[3];
sx q[3];
rz(0.20460953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.8873022) q[2];
sx q[2];
rz(-1.3672978) q[2];
sx q[2];
rz(-2.5781412) q[2];
rz(-3.0900132) q[3];
sx q[3];
rz(-1.1392461) q[3];
sx q[3];
rz(-2.1896867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1812487) q[0];
sx q[0];
rz(-0.5287756) q[0];
sx q[0];
rz(-1.7425591) q[0];
rz(0.78701204) q[1];
sx q[1];
rz(-2.0128638) q[1];
sx q[1];
rz(2.3972437) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23119152) q[0];
sx q[0];
rz(-1.4230799) q[0];
sx q[0];
rz(1.5157248) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0853945) q[2];
sx q[2];
rz(-1.4375763) q[2];
sx q[2];
rz(2.104987) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-3.0155322) q[1];
sx q[1];
rz(-2.0704401) q[1];
sx q[1];
rz(1.8263032) q[1];
rz(-pi) q[2];
rz(-2.6870071) q[3];
sx q[3];
rz(-1.5496407) q[3];
sx q[3];
rz(1.3190312) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.7156334) q[2];
sx q[2];
rz(-1.3954433) q[2];
sx q[2];
rz(2.5320833) q[2];
rz(-0.65731796) q[3];
sx q[3];
rz(-2.4980563) q[3];
sx q[3];
rz(2.8801584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9534) q[0];
sx q[0];
rz(-3.0472026) q[0];
sx q[0];
rz(-1.6375861) q[0];
rz(1.9001182) q[1];
sx q[1];
rz(-1.1499317) q[1];
sx q[1];
rz(-0.77493587) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.06668815) q[0];
sx q[0];
rz(-0.79830352) q[0];
sx q[0];
rz(2.128128) q[0];
rz(-pi) q[1];
rz(-2.0140892) q[2];
sx q[2];
rz(-1.4294251) q[2];
sx q[2];
rz(0.023250016) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.51674023) q[1];
sx q[1];
rz(-1.7588741) q[1];
sx q[1];
rz(2.6810758) q[1];
rz(2.9990066) q[3];
sx q[3];
rz(-1.3538059) q[3];
sx q[3];
rz(-2.2246974) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.9541786) q[2];
sx q[2];
rz(-2.9118907) q[2];
sx q[2];
rz(2.9837218) q[2];
rz(-1.9291417) q[3];
sx q[3];
rz(-2.082943) q[3];
sx q[3];
rz(-1.2780179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.071844) q[0];
sx q[0];
rz(-2.1691515) q[0];
sx q[0];
rz(2.9272595) q[0];
rz(-2.4841323) q[1];
sx q[1];
rz(-0.22413707) q[1];
sx q[1];
rz(-1.0459895) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8911274) q[0];
sx q[0];
rz(-1.9292826) q[0];
sx q[0];
rz(1.6842708) q[0];
rz(-pi) q[1];
rz(-1.5418391) q[2];
sx q[2];
rz(-0.98000295) q[2];
sx q[2];
rz(-0.40320413) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6007081) q[1];
sx q[1];
rz(-0.67544671) q[1];
sx q[1];
rz(2.1727968) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0797327) q[3];
sx q[3];
rz(-2.7354771) q[3];
sx q[3];
rz(-1.5861685) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.41632286) q[2];
sx q[2];
rz(-3.0815093) q[2];
sx q[2];
rz(-2.0521169) q[2];
rz(1.5754835) q[3];
sx q[3];
rz(-1.243306) q[3];
sx q[3];
rz(2.4889448) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8626704) q[0];
sx q[0];
rz(-0.60386064) q[0];
sx q[0];
rz(0.84558564) q[0];
rz(1.6090341) q[1];
sx q[1];
rz(-1.4668203) q[1];
sx q[1];
rz(-1.1046881) q[1];
rz(-2.7325148) q[2];
sx q[2];
rz(-1.2862051) q[2];
sx q[2];
rz(-0.4551879) q[2];
rz(1.6339176) q[3];
sx q[3];
rz(-2.1199385) q[3];
sx q[3];
rz(-1.4169823) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
