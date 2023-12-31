OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.5542334) q[0];
sx q[0];
rz(4.1242546) q[0];
sx q[0];
rz(10.186515) q[0];
rz(-2.3770483) q[1];
sx q[1];
rz(-1.0772871) q[1];
sx q[1];
rz(2.397937) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3418158) q[0];
sx q[0];
rz(-1.8108978) q[0];
sx q[0];
rz(0.067702985) q[0];
rz(-pi) q[1];
x q[1];
rz(0.15723575) q[2];
sx q[2];
rz(-0.6586282) q[2];
sx q[2];
rz(-0.14070357) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.8929157) q[1];
sx q[1];
rz(-0.49202737) q[1];
sx q[1];
rz(-0.9534652) q[1];
rz(-0.071804382) q[3];
sx q[3];
rz(-1.9339438) q[3];
sx q[3];
rz(-0.24955173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.7068229) q[2];
sx q[2];
rz(-2.7396024) q[2];
sx q[2];
rz(0.10786954) q[2];
rz(2.9989631) q[3];
sx q[3];
rz(-1.4167891) q[3];
sx q[3];
rz(-2.4690348) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.07664872) q[0];
sx q[0];
rz(-2.3363484) q[0];
sx q[0];
rz(0.23072492) q[0];
rz(1.8143066) q[1];
sx q[1];
rz(-2.4704411) q[1];
sx q[1];
rz(-0.040963106) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7777268) q[0];
sx q[0];
rz(-1.7330609) q[0];
sx q[0];
rz(-1.3217539) q[0];
rz(1.8789005) q[2];
sx q[2];
rz(-1.3034046) q[2];
sx q[2];
rz(-1.2958796) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7663824) q[1];
sx q[1];
rz(-2.3751405) q[1];
sx q[1];
rz(-2.4541928) q[1];
rz(-2.326194) q[3];
sx q[3];
rz(-2.6147463) q[3];
sx q[3];
rz(2.6032053) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.9791947) q[2];
sx q[2];
rz(-1.4731864) q[2];
sx q[2];
rz(2.5734148) q[2];
rz(-0.5125106) q[3];
sx q[3];
rz(-0.5957225) q[3];
sx q[3];
rz(-2.5797243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42999643) q[0];
sx q[0];
rz(-1.6226409) q[0];
sx q[0];
rz(2.6913753) q[0];
rz(1.2954767) q[1];
sx q[1];
rz(-1.1643012) q[1];
sx q[1];
rz(-0.67726642) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60904658) q[0];
sx q[0];
rz(-2.9653774) q[0];
sx q[0];
rz(-0.30058582) q[0];
rz(-pi) q[1];
x q[1];
rz(0.0063653221) q[2];
sx q[2];
rz(-1.5718939) q[2];
sx q[2];
rz(1.7249677) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.8202159) q[1];
sx q[1];
rz(-1.7019338) q[1];
sx q[1];
rz(-1.7507491) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4694091) q[3];
sx q[3];
rz(-2.5568092) q[3];
sx q[3];
rz(-0.72090805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.21851097) q[2];
sx q[2];
rz(-1.1542902) q[2];
sx q[2];
rz(-3.0947321) q[2];
rz(2.3299407) q[3];
sx q[3];
rz(-0.67957687) q[3];
sx q[3];
rz(2.9714382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1470404) q[0];
sx q[0];
rz(-3.0299598) q[0];
sx q[0];
rz(-3.0901093) q[0];
rz(-2.6507846) q[1];
sx q[1];
rz(-2.1665116) q[1];
sx q[1];
rz(-1.9690008) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20576142) q[0];
sx q[0];
rz(-1.7044221) q[0];
sx q[0];
rz(-0.66098102) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5355859) q[2];
sx q[2];
rz(-1.5835985) q[2];
sx q[2];
rz(2.9738604) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.40201515) q[1];
sx q[1];
rz(-0.47164729) q[1];
sx q[1];
rz(2.2775843) q[1];
x q[2];
rz(1.138932) q[3];
sx q[3];
rz(-1.5032839) q[3];
sx q[3];
rz(0.95089144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.2525758) q[2];
sx q[2];
rz(-0.93095195) q[2];
sx q[2];
rz(0.48689294) q[2];
rz(-1.9481109) q[3];
sx q[3];
rz(-0.32961696) q[3];
sx q[3];
rz(0.025432767) q[3];
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
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.837773) q[0];
sx q[0];
rz(-2.1210414) q[0];
sx q[0];
rz(-1.3254962) q[0];
rz(2.5505113) q[1];
sx q[1];
rz(-0.68060827) q[1];
sx q[1];
rz(-2.4868734) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7692208) q[0];
sx q[0];
rz(-1.1290316) q[0];
sx q[0];
rz(1.9584993) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.40277092) q[2];
sx q[2];
rz(-2.6908814) q[2];
sx q[2];
rz(-0.41926256) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.66203413) q[1];
sx q[1];
rz(-0.91842945) q[1];
sx q[1];
rz(1.0865092) q[1];
x q[2];
rz(-1.4773024) q[3];
sx q[3];
rz(-1.9545385) q[3];
sx q[3];
rz(-1.5941217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.6902265) q[2];
sx q[2];
rz(-1.2008685) q[2];
sx q[2];
rz(0.5711242) q[2];
rz(-0.58978224) q[3];
sx q[3];
rz(-0.46001205) q[3];
sx q[3];
rz(3.0736198) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(2.1968483) q[0];
sx q[0];
rz(-1.9510883) q[0];
sx q[0];
rz(2.8425472) q[0];
rz(1.8213182) q[1];
sx q[1];
rz(-0.25779217) q[1];
sx q[1];
rz(1.6437795) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4371571) q[0];
sx q[0];
rz(-2.2309982) q[0];
sx q[0];
rz(-0.066083834) q[0];
rz(-0.97700714) q[2];
sx q[2];
rz(-1.0736246) q[2];
sx q[2];
rz(1.1360816) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.99243473) q[1];
sx q[1];
rz(-2.5870393) q[1];
sx q[1];
rz(-1.0013594) q[1];
x q[2];
rz(1.8201581) q[3];
sx q[3];
rz(-1.2994088) q[3];
sx q[3];
rz(-2.4621778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.6254639) q[2];
sx q[2];
rz(-1.7136145) q[2];
sx q[2];
rz(-3.1141172) q[2];
rz(0.52250683) q[3];
sx q[3];
rz(-0.80111879) q[3];
sx q[3];
rz(-0.81645042) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7291173) q[0];
sx q[0];
rz(-2.0232047) q[0];
sx q[0];
rz(3.0274042) q[0];
rz(2.1633637) q[1];
sx q[1];
rz(-2.5949635) q[1];
sx q[1];
rz(-0.79089975) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7366911) q[0];
sx q[0];
rz(-1.6653403) q[0];
sx q[0];
rz(-2.1928284) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2355455) q[2];
sx q[2];
rz(-2.2860048) q[2];
sx q[2];
rz(-0.4244948) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.8154527) q[1];
sx q[1];
rz(-0.37323144) q[1];
sx q[1];
rz(-1.3532072) q[1];
x q[2];
rz(-2.7780276) q[3];
sx q[3];
rz(-2.225038) q[3];
sx q[3];
rz(-0.082106575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.2699282) q[2];
sx q[2];
rz(-1.1130788) q[2];
sx q[2];
rz(0.86501914) q[2];
rz(-0.67251742) q[3];
sx q[3];
rz(-2.0130242) q[3];
sx q[3];
rz(-0.095656693) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4928116) q[0];
sx q[0];
rz(-2.6219941) q[0];
sx q[0];
rz(2.6742324) q[0];
rz(2.6043747) q[1];
sx q[1];
rz(-0.98058128) q[1];
sx q[1];
rz(-0.25407243) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8763037) q[0];
sx q[0];
rz(-2.9348001) q[0];
sx q[0];
rz(-2.8121023) q[0];
rz(-pi) q[1];
rz(-1.5723096) q[2];
sx q[2];
rz(-1.5760033) q[2];
sx q[2];
rz(1.4323915) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.2368187) q[1];
sx q[1];
rz(-0.33729759) q[1];
sx q[1];
rz(-2.7493613) q[1];
rz(-pi) q[2];
rz(0.11717511) q[3];
sx q[3];
rz(-0.60248884) q[3];
sx q[3];
rz(0.52091) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.6415928) q[2];
sx q[2];
rz(-2.3708512) q[2];
sx q[2];
rz(0.33561486) q[2];
rz(-2.8149758) q[3];
sx q[3];
rz(-0.89544046) q[3];
sx q[3];
rz(1.7232822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0582054) q[0];
sx q[0];
rz(-2.931262) q[0];
sx q[0];
rz(-2.0876419) q[0];
rz(2.9846233) q[1];
sx q[1];
rz(-1.7228246) q[1];
sx q[1];
rz(-0.98186791) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4440585) q[0];
sx q[0];
rz(-0.42675787) q[0];
sx q[0];
rz(-2.8696637) q[0];
rz(3.0332546) q[2];
sx q[2];
rz(-2.1370558) q[2];
sx q[2];
rz(2.1926751) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.30584221) q[1];
sx q[1];
rz(-1.7562477) q[1];
sx q[1];
rz(-3.0830543) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3217661) q[3];
sx q[3];
rz(-0.8753652) q[3];
sx q[3];
rz(1.673656) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2173569) q[2];
sx q[2];
rz(-1.057426) q[2];
sx q[2];
rz(0.99739933) q[2];
rz(-0.50619566) q[3];
sx q[3];
rz(-0.96118569) q[3];
sx q[3];
rz(-3.1072646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(0.85703325) q[0];
sx q[0];
rz(-2.4301346) q[0];
sx q[0];
rz(-2.6249264) q[0];
rz(-2.754028) q[1];
sx q[1];
rz(-1.060408) q[1];
sx q[1];
rz(0.26836747) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3049406) q[0];
sx q[0];
rz(-0.9965082) q[0];
sx q[0];
rz(-2.4313297) q[0];
x q[1];
rz(2.1568314) q[2];
sx q[2];
rz(-1.1981989) q[2];
sx q[2];
rz(-2.7013456) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.25071535) q[1];
sx q[1];
rz(-2.4836575) q[1];
sx q[1];
rz(0.4472181) q[1];
rz(1.9645421) q[3];
sx q[3];
rz(-2.5265794) q[3];
sx q[3];
rz(-2.7494489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.9364075) q[2];
sx q[2];
rz(-0.78339094) q[2];
sx q[2];
rz(-2.5777585) q[2];
rz(1.1901723) q[3];
sx q[3];
rz(-0.95919132) q[3];
sx q[3];
rz(-2.4998375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.664809) q[0];
sx q[0];
rz(-1.2510779) q[0];
sx q[0];
rz(-1.0673987) q[0];
rz(1.8019567) q[1];
sx q[1];
rz(-1.7032774) q[1];
sx q[1];
rz(1.3443321) q[1];
rz(2.6396991) q[2];
sx q[2];
rz(-0.5770275) q[2];
sx q[2];
rz(1.8736476) q[2];
rz(-0.56223829) q[3];
sx q[3];
rz(-0.36459618) q[3];
sx q[3];
rz(0.96095745) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
