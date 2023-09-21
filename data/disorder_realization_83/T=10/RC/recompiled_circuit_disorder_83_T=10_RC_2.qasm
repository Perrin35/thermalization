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
rz(-2.1589307) q[0];
sx q[0];
rz(0.76173705) q[0];
rz(-2.3770483) q[1];
sx q[1];
rz(-1.0772871) q[1];
sx q[1];
rz(-0.74365562) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0640472) q[0];
sx q[0];
rz(-2.8923058) q[0];
sx q[0];
rz(-1.8403948) q[0];
rz(1.6913937) q[2];
sx q[2];
rz(-2.2199124) q[2];
sx q[2];
rz(2.803034) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.570959) q[1];
sx q[1];
rz(-1.966241) q[1];
sx q[1];
rz(-0.30084893) q[1];
rz(-pi) q[2];
x q[2];
rz(0.071804382) q[3];
sx q[3];
rz(-1.9339438) q[3];
sx q[3];
rz(0.24955173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.4347697) q[2];
sx q[2];
rz(-2.7396024) q[2];
sx q[2];
rz(-3.0337231) q[2];
rz(0.14262959) q[3];
sx q[3];
rz(-1.4167891) q[3];
sx q[3];
rz(2.4690348) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
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
rz(-2.9108677) q[0];
rz(1.8143066) q[1];
sx q[1];
rz(-0.67115152) q[1];
sx q[1];
rz(0.040963106) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3687392) q[0];
sx q[0];
rz(-2.8452747) q[0];
sx q[0];
rz(-2.1570737) q[0];
rz(0.27994056) q[2];
sx q[2];
rz(-1.2739812) q[2];
sx q[2];
rz(2.9505626) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.47536182) q[1];
sx q[1];
rz(-1.0050217) q[1];
sx q[1];
rz(2.1192141) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9713692) q[3];
sx q[3];
rz(-1.2188606) q[3];
sx q[3];
rz(-0.34917253) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.162398) q[2];
sx q[2];
rz(-1.6684063) q[2];
sx q[2];
rz(2.5734148) q[2];
rz(0.5125106) q[3];
sx q[3];
rz(-0.5957225) q[3];
sx q[3];
rz(2.5797243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7115962) q[0];
sx q[0];
rz(-1.6226409) q[0];
sx q[0];
rz(2.6913753) q[0];
rz(-1.846116) q[1];
sx q[1];
rz(-1.9772915) q[1];
sx q[1];
rz(0.67726642) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8836425) q[0];
sx q[0];
rz(-1.6227239) q[0];
sx q[0];
rz(0.16846637) q[0];
x q[1];
rz(-1.5696987) q[2];
sx q[2];
rz(-1.564431) q[2];
sx q[2];
rz(0.15417834) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.22563572) q[1];
sx q[1];
rz(-1.7491873) q[1];
sx q[1];
rz(-0.13326463) q[1];
rz(-pi) q[2];
rz(-2.4694091) q[3];
sx q[3];
rz(-0.58478343) q[3];
sx q[3];
rz(0.72090805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.21851097) q[2];
sx q[2];
rz(-1.9873025) q[2];
sx q[2];
rz(3.0947321) q[2];
rz(0.81165195) q[3];
sx q[3];
rz(-2.4620158) q[3];
sx q[3];
rz(-0.17015447) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1470404) q[0];
sx q[0];
rz(-0.11163286) q[0];
sx q[0];
rz(3.0901093) q[0];
rz(-0.4908081) q[1];
sx q[1];
rz(-0.97508109) q[1];
sx q[1];
rz(-1.9690008) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9358312) q[0];
sx q[0];
rz(-1.7044221) q[0];
sx q[0];
rz(-2.4806116) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.1191191) q[2];
sx q[2];
rz(-0.60612504) q[2];
sx q[2];
rz(-1.7200574) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.1662894) q[1];
sx q[1];
rz(-1.2180093) q[1];
sx q[1];
rz(-0.3198448) q[1];
rz(3.067279) q[3];
sx q[3];
rz(-2.001611) q[3];
sx q[3];
rz(2.5527692) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.88901687) q[2];
sx q[2];
rz(-0.93095195) q[2];
sx q[2];
rz(2.6546997) q[2];
rz(-1.9481109) q[3];
sx q[3];
rz(-2.8119757) q[3];
sx q[3];
rz(-0.025432767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.837773) q[0];
sx q[0];
rz(-1.0205512) q[0];
sx q[0];
rz(1.3254962) q[0];
rz(2.5505113) q[1];
sx q[1];
rz(-0.68060827) q[1];
sx q[1];
rz(-2.4868734) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39010534) q[0];
sx q[0];
rz(-2.562398) q[0];
sx q[0];
rz(2.4672227) q[0];
rz(-2.7388217) q[2];
sx q[2];
rz(-0.45071128) q[2];
sx q[2];
rz(-0.41926256) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.4795585) q[1];
sx q[1];
rz(-2.2231632) q[1];
sx q[1];
rz(1.0865092) q[1];
x q[2];
rz(-2.914364) q[3];
sx q[3];
rz(-2.7471746) q[3];
sx q[3];
rz(1.3487181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.6902265) q[2];
sx q[2];
rz(-1.9407242) q[2];
sx q[2];
rz(0.5711242) q[2];
rz(0.58978224) q[3];
sx q[3];
rz(-2.6815806) q[3];
sx q[3];
rz(-0.067972876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.94474435) q[0];
sx q[0];
rz(-1.9510883) q[0];
sx q[0];
rz(-0.29904547) q[0];
rz(1.3202745) q[1];
sx q[1];
rz(-0.25779217) q[1];
sx q[1];
rz(-1.6437795) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2346674) q[0];
sx q[0];
rz(-1.6229796) q[0];
sx q[0];
rz(-2.232057) q[0];
rz(-0.97700714) q[2];
sx q[2];
rz(-1.0736246) q[2];
sx q[2];
rz(1.1360816) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.0648246) q[1];
sx q[1];
rz(-1.8586564) q[1];
sx q[1];
rz(-2.0516146) q[1];
rz(2.4160556) q[3];
sx q[3];
rz(-0.36645884) q[3];
sx q[3];
rz(-0.08034245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.51612878) q[2];
sx q[2];
rz(-1.7136145) q[2];
sx q[2];
rz(-0.027475474) q[2];
rz(-2.6190858) q[3];
sx q[3];
rz(-2.3404739) q[3];
sx q[3];
rz(0.81645042) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41247535) q[0];
sx q[0];
rz(-1.1183879) q[0];
sx q[0];
rz(0.11418848) q[0];
rz(2.1633637) q[1];
sx q[1];
rz(-0.54662919) q[1];
sx q[1];
rz(0.79089975) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40490155) q[0];
sx q[0];
rz(-1.6653403) q[0];
sx q[0];
rz(0.94876429) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3979264) q[2];
sx q[2];
rz(-1.8218092) q[2];
sx q[2];
rz(-0.92168346) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.32614) q[1];
sx q[1];
rz(-0.37323144) q[1];
sx q[1];
rz(1.7883854) q[1];
rz(2.2579455) q[3];
sx q[3];
rz(-1.8568608) q[3];
sx q[3];
rz(-1.7162232) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.2699282) q[2];
sx q[2];
rz(-1.1130788) q[2];
sx q[2];
rz(2.2765735) q[2];
rz(-2.4690752) q[3];
sx q[3];
rz(-2.0130242) q[3];
sx q[3];
rz(0.095656693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6487811) q[0];
sx q[0];
rz(-2.6219941) q[0];
sx q[0];
rz(-2.6742324) q[0];
rz(0.53721792) q[1];
sx q[1];
rz(-0.98058128) q[1];
sx q[1];
rz(0.25407243) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8763037) q[0];
sx q[0];
rz(-0.2067925) q[0];
sx q[0];
rz(-0.32949038) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8587564) q[2];
sx q[2];
rz(-3.1361702) q[2];
sx q[2];
rz(1.1495513) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.904774) q[1];
sx q[1];
rz(-2.8042951) q[1];
sx q[1];
rz(-0.39223139) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.11717511) q[3];
sx q[3];
rz(-0.60248884) q[3];
sx q[3];
rz(-0.52091) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.49999985) q[2];
sx q[2];
rz(-2.3708512) q[2];
sx q[2];
rz(-2.8059778) q[2];
rz(-0.32661682) q[3];
sx q[3];
rz(-0.89544046) q[3];
sx q[3];
rz(-1.7232822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0582054) q[0];
sx q[0];
rz(-2.931262) q[0];
sx q[0];
rz(2.0876419) q[0];
rz(0.15696934) q[1];
sx q[1];
rz(-1.7228246) q[1];
sx q[1];
rz(0.98186791) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1218402) q[0];
sx q[0];
rz(-1.4593908) q[0];
sx q[0];
rz(0.41282546) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4023196) q[2];
sx q[2];
rz(-2.5661764) q[2];
sx q[2];
rz(-1.9926496) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.2541483) q[1];
sx q[1];
rz(-1.5132628) q[1];
sx q[1];
rz(-1.7565586) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.28717678) q[3];
sx q[3];
rz(-0.73162006) q[3];
sx q[3];
rz(-1.8457796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.92423576) q[2];
sx q[2];
rz(-2.0841667) q[2];
sx q[2];
rz(0.99739933) q[2];
rz(2.635397) q[3];
sx q[3];
rz(-0.96118569) q[3];
sx q[3];
rz(-3.1072646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2845594) q[0];
sx q[0];
rz(-0.71145809) q[0];
sx q[0];
rz(-2.6249264) q[0];
rz(0.38756469) q[1];
sx q[1];
rz(-1.060408) q[1];
sx q[1];
rz(0.26836747) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4386913) q[0];
sx q[0];
rz(-0.99150204) q[0];
sx q[0];
rz(-2.2772574) q[0];
rz(-pi) q[1];
rz(0.95558138) q[2];
sx q[2];
rz(-0.68253839) q[2];
sx q[2];
rz(0.62894422) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.8908773) q[1];
sx q[1];
rz(-0.6579352) q[1];
sx q[1];
rz(2.6943745) q[1];
x q[2];
rz(-1.9645421) q[3];
sx q[3];
rz(-0.6150133) q[3];
sx q[3];
rz(-2.7494489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.2051852) q[2];
sx q[2];
rz(-0.78339094) q[2];
sx q[2];
rz(0.56383413) q[2];
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
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47678369) q[0];
sx q[0];
rz(-1.2510779) q[0];
sx q[0];
rz(-1.0673987) q[0];
rz(1.339636) q[1];
sx q[1];
rz(-1.4383153) q[1];
sx q[1];
rz(-1.7972606) q[1];
rz(-0.50189353) q[2];
sx q[2];
rz(-0.5770275) q[2];
sx q[2];
rz(1.8736476) q[2];
rz(-1.7715122) q[3];
sx q[3];
rz(-1.8772535) q[3];
sx q[3];
rz(1.5542961) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
