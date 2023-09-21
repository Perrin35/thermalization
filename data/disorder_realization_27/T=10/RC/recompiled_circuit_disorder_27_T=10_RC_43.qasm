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
rz(1.9321631) q[1];
sx q[1];
rz(-1.2485319) q[1];
sx q[1];
rz(1.7564397) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4614852) q[0];
sx q[0];
rz(-2.2930817) q[0];
sx q[0];
rz(-2.0018342) q[0];
rz(2.7117549) q[2];
sx q[2];
rz(-2.5463383) q[2];
sx q[2];
rz(2.0855479) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.17978046) q[1];
sx q[1];
rz(-1.7331859) q[1];
sx q[1];
rz(-0.23602545) q[1];
x q[2];
rz(-2.0632283) q[3];
sx q[3];
rz(-2.1636117) q[3];
sx q[3];
rz(-1.4584695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.8866855) q[2];
sx q[2];
rz(-0.79780769) q[2];
sx q[2];
rz(-2.936426) q[2];
rz(-2.3702879) q[3];
sx q[3];
rz(-0.78273928) q[3];
sx q[3];
rz(2.0390959) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40760621) q[0];
sx q[0];
rz(-0.74626958) q[0];
sx q[0];
rz(-2.6876887) q[0];
rz(-2.1167963) q[1];
sx q[1];
rz(-2.7339934) q[1];
sx q[1];
rz(-1.9143547) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95943806) q[0];
sx q[0];
rz(-1.4417138) q[0];
sx q[0];
rz(-0.076230787) q[0];
rz(2.0298376) q[2];
sx q[2];
rz(-2.0437721) q[2];
sx q[2];
rz(-2.3943261) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.36453516) q[1];
sx q[1];
rz(-0.39452663) q[1];
sx q[1];
rz(-0.72655861) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7271872) q[3];
sx q[3];
rz(-1.3723433) q[3];
sx q[3];
rz(-2.1895777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.041302117) q[2];
sx q[2];
rz(-1.9561448) q[2];
sx q[2];
rz(2.5741637) q[2];
rz(0.36519095) q[3];
sx q[3];
rz(-1.4130211) q[3];
sx q[3];
rz(2.173483) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48297468) q[0];
sx q[0];
rz(-2.5768319) q[0];
sx q[0];
rz(-2.2429402) q[0];
rz(-0.99575106) q[1];
sx q[1];
rz(-1.5581222) q[1];
sx q[1];
rz(-0.333289) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8329187) q[0];
sx q[0];
rz(-1.1090288) q[0];
sx q[0];
rz(2.4426016) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0684418) q[2];
sx q[2];
rz(-2.9403254) q[2];
sx q[2];
rz(1.7114491) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0888575) q[1];
sx q[1];
rz(-2.5087068) q[1];
sx q[1];
rz(-1.8295374) q[1];
rz(-pi) q[2];
rz(-2.1786147) q[3];
sx q[3];
rz(-0.64219785) q[3];
sx q[3];
rz(-1.1841139) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.4553392) q[2];
sx q[2];
rz(-1.7909966) q[2];
sx q[2];
rz(-1.1509482) q[2];
rz(-0.84093705) q[3];
sx q[3];
rz(-1.1154113) q[3];
sx q[3];
rz(-1.9870728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6999917) q[0];
sx q[0];
rz(-1.561152) q[0];
sx q[0];
rz(-0.68471318) q[0];
rz(1.0355863) q[1];
sx q[1];
rz(-0.5077478) q[1];
sx q[1];
rz(-1.9365786) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22574632) q[0];
sx q[0];
rz(-1.5374743) q[0];
sx q[0];
rz(-0.87945625) q[0];
x q[1];
rz(2.2150546) q[2];
sx q[2];
rz(-1.7346003) q[2];
sx q[2];
rz(1.0545727) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.1531096) q[1];
sx q[1];
rz(-1.298269) q[1];
sx q[1];
rz(0.86221282) q[1];
rz(-pi) q[2];
rz(-0.071675008) q[3];
sx q[3];
rz(-2.274401) q[3];
sx q[3];
rz(-3.049831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.24923199) q[2];
sx q[2];
rz(-1.7148596) q[2];
sx q[2];
rz(-2.7704346) q[2];
rz(-1.4012198) q[3];
sx q[3];
rz(-0.6597844) q[3];
sx q[3];
rz(-2.0223117) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.054759653) q[0];
sx q[0];
rz(-0.78614569) q[0];
sx q[0];
rz(0.13312419) q[0];
rz(2.1482824) q[1];
sx q[1];
rz(-1.3860093) q[1];
sx q[1];
rz(0.55508074) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0017437) q[0];
sx q[0];
rz(-1.2554597) q[0];
sx q[0];
rz(-3.128201) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1558756) q[2];
sx q[2];
rz(-1.5083815) q[2];
sx q[2];
rz(-1.2635363) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3913369) q[1];
sx q[1];
rz(-1.5515944) q[1];
sx q[1];
rz(1.125543) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.57226945) q[3];
sx q[3];
rz(-0.096147691) q[3];
sx q[3];
rz(-2.8807246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.30620265) q[2];
sx q[2];
rz(-1.0058879) q[2];
sx q[2];
rz(0.13892697) q[2];
rz(2.1991918) q[3];
sx q[3];
rz(-1.5706294) q[3];
sx q[3];
rz(0.55148235) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5979364) q[0];
sx q[0];
rz(-2.5718226) q[0];
sx q[0];
rz(0.57975769) q[0];
rz(3.014091) q[1];
sx q[1];
rz(-1.189905) q[1];
sx q[1];
rz(-1.5396083) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72551661) q[0];
sx q[0];
rz(-0.86891876) q[0];
sx q[0];
rz(-0.26728018) q[0];
rz(-pi) q[1];
rz(1.6164262) q[2];
sx q[2];
rz(-1.4403733) q[2];
sx q[2];
rz(1.6861196) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.74486596) q[1];
sx q[1];
rz(-2.6187881) q[1];
sx q[1];
rz(0.97739873) q[1];
rz(0.88605373) q[3];
sx q[3];
rz(-1.7985117) q[3];
sx q[3];
rz(0.39623228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.5876028) q[2];
sx q[2];
rz(-0.25045276) q[2];
sx q[2];
rz(2.8721151) q[2];
rz(-0.23412165) q[3];
sx q[3];
rz(-0.52474371) q[3];
sx q[3];
rz(0.060119303) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6948029) q[0];
sx q[0];
rz(-0.61820522) q[0];
sx q[0];
rz(-2.457298) q[0];
rz(-0.11958312) q[1];
sx q[1];
rz(-1.8493098) q[1];
sx q[1];
rz(0.51876846) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9414026) q[0];
sx q[0];
rz(-1.4369643) q[0];
sx q[0];
rz(-2.3076513) q[0];
x q[1];
rz(0.77063009) q[2];
sx q[2];
rz(-1.7297041) q[2];
sx q[2];
rz(-0.7030013) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.24592933) q[1];
sx q[1];
rz(-2.3024125) q[1];
sx q[1];
rz(-0.64114665) q[1];
x q[2];
rz(0.026168907) q[3];
sx q[3];
rz(-1.1702288) q[3];
sx q[3];
rz(-2.8699584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8873022) q[2];
sx q[2];
rz(-1.7742949) q[2];
sx q[2];
rz(-2.5781412) q[2];
rz(-3.0900132) q[3];
sx q[3];
rz(-1.1392461) q[3];
sx q[3];
rz(0.95190597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96034399) q[0];
sx q[0];
rz(-0.5287756) q[0];
sx q[0];
rz(1.3990336) q[0];
rz(-0.78701204) q[1];
sx q[1];
rz(-1.1287289) q[1];
sx q[1];
rz(2.3972437) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58957419) q[0];
sx q[0];
rz(-0.15757832) q[0];
sx q[0];
rz(-2.7872859) q[0];
rz(-pi) q[1];
rz(-2.9888399) q[2];
sx q[2];
rz(-1.061201) q[2];
sx q[2];
rz(0.60915142) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.12606049) q[1];
sx q[1];
rz(-2.0704401) q[1];
sx q[1];
rz(1.3152895) q[1];
x q[2];
rz(-0.45458557) q[3];
sx q[3];
rz(-1.5496407) q[3];
sx q[3];
rz(-1.3190312) q[3];
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
rz(-0.60950935) q[2];
rz(-2.4842747) q[3];
sx q[3];
rz(-0.64353639) q[3];
sx q[3];
rz(-0.26143423) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9534) q[0];
sx q[0];
rz(-0.094390079) q[0];
sx q[0];
rz(-1.6375861) q[0];
rz(-1.9001182) q[1];
sx q[1];
rz(-1.1499317) q[1];
sx q[1];
rz(0.77493587) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.06668815) q[0];
sx q[0];
rz(-0.79830352) q[0];
sx q[0];
rz(1.0134646) q[0];
x q[1];
rz(0.1562642) q[2];
sx q[2];
rz(-2.0093577) q[2];
sx q[2];
rz(1.4807448) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.96156582) q[1];
sx q[1];
rz(-1.119009) q[1];
sx q[1];
rz(-1.3614484) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9990066) q[3];
sx q[3];
rz(-1.3538059) q[3];
sx q[3];
rz(-2.2246974) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.187414) q[2];
sx q[2];
rz(-2.9118907) q[2];
sx q[2];
rz(-0.15787086) q[2];
rz(1.212451) q[3];
sx q[3];
rz(-2.082943) q[3];
sx q[3];
rz(-1.2780179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-2.0697486) q[0];
sx q[0];
rz(-2.1691515) q[0];
sx q[0];
rz(0.21433314) q[0];
rz(-2.4841323) q[1];
sx q[1];
rz(-0.22413707) q[1];
sx q[1];
rz(2.0956031) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5770618) q[0];
sx q[0];
rz(-2.7663167) q[0];
sx q[0];
rz(0.29348404) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5506053) q[2];
sx q[2];
rz(-1.5467484) q[2];
sx q[2];
rz(1.1514593) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.6007081) q[1];
sx q[1];
rz(-0.67544671) q[1];
sx q[1];
rz(-2.1727968) q[1];
x q[2];
rz(1.5973741) q[3];
sx q[3];
rz(-1.9760895) q[3];
sx q[3];
rz(1.4881031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7252698) q[2];
sx q[2];
rz(-3.0815093) q[2];
sx q[2];
rz(1.0894758) q[2];
rz(-1.5661092) q[3];
sx q[3];
rz(-1.243306) q[3];
sx q[3];
rz(-0.65264788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8626704) q[0];
sx q[0];
rz(-2.537732) q[0];
sx q[0];
rz(-2.296007) q[0];
rz(-1.5325585) q[1];
sx q[1];
rz(-1.4668203) q[1];
sx q[1];
rz(-1.1046881) q[1];
rz(0.63411843) q[2];
sx q[2];
rz(-0.49370439) q[2];
sx q[2];
rz(1.6903071) q[2];
rz(0.10272051) q[3];
sx q[3];
rz(-0.552388) q[3];
sx q[3];
rz(-1.296476) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];