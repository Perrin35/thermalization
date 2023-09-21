OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.3553319) q[0];
sx q[0];
rz(0.064602764) q[0];
sx q[0];
rz(6.3048007) q[0];
rz(2.1463483) q[1];
sx q[1];
rz(-1.8145476) q[1];
sx q[1];
rz(-1.8099161) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74237139) q[0];
sx q[0];
rz(-2.135015) q[0];
sx q[0];
rz(-2.7447634) q[0];
rz(-pi) q[1];
rz(-0.10279074) q[2];
sx q[2];
rz(-1.3445026) q[2];
sx q[2];
rz(0.13470995) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5035489) q[1];
sx q[1];
rz(-0.81144864) q[1];
sx q[1];
rz(-1.7502977) q[1];
rz(-pi) q[2];
x q[2];
rz(0.33820037) q[3];
sx q[3];
rz(-0.97465289) q[3];
sx q[3];
rz(-0.81830922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.8109479) q[2];
sx q[2];
rz(-1.6487048) q[2];
sx q[2];
rz(2.5374106) q[2];
rz(-2.1172681) q[3];
sx q[3];
rz(-1.9842792) q[3];
sx q[3];
rz(-2.0143051) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8169096) q[0];
sx q[0];
rz(-3.1284101) q[0];
sx q[0];
rz(1.0634134) q[0];
rz(0.88513199) q[1];
sx q[1];
rz(-1.5848031) q[1];
sx q[1];
rz(-0.0016454776) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9379399) q[0];
sx q[0];
rz(-1.5601888) q[0];
sx q[0];
rz(1.525735) q[0];
x q[1];
rz(-2.4992141) q[2];
sx q[2];
rz(-1.3628236) q[2];
sx q[2];
rz(0.013052879) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.21316321) q[1];
sx q[1];
rz(-2.3191116) q[1];
sx q[1];
rz(-0.33555062) q[1];
rz(0.1699154) q[3];
sx q[3];
rz(-1.8222457) q[3];
sx q[3];
rz(-0.91472317) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.985618) q[2];
sx q[2];
rz(-1.5094455) q[2];
sx q[2];
rz(2.4334811) q[2];
rz(1.0937141) q[3];
sx q[3];
rz(-1.3392859) q[3];
sx q[3];
rz(0.99350199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.920632) q[0];
sx q[0];
rz(-1.7376124) q[0];
sx q[0];
rz(0.8272585) q[0];
rz(-0.0050841252) q[1];
sx q[1];
rz(-1.9283483) q[1];
sx q[1];
rz(-1.089383) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0076865772) q[0];
sx q[0];
rz(-1.0431164) q[0];
sx q[0];
rz(-2.6469127) q[0];
rz(1.9956279) q[2];
sx q[2];
rz(-1.0152738) q[2];
sx q[2];
rz(-2.5512763) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.3315862) q[1];
sx q[1];
rz(-0.3948822) q[1];
sx q[1];
rz(1.7942795) q[1];
rz(-1.009922) q[3];
sx q[3];
rz(-2.2787333) q[3];
sx q[3];
rz(-2.5885955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.0777145) q[2];
sx q[2];
rz(-0.87746799) q[2];
sx q[2];
rz(2.2568978) q[2];
rz(1.1832773) q[3];
sx q[3];
rz(-1.8235455) q[3];
sx q[3];
rz(-1.6515091) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
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
rz(-1.5723715) q[0];
sx q[0];
rz(-2.5771993) q[0];
sx q[0];
rz(-0.72682056) q[0];
rz(0.80365333) q[1];
sx q[1];
rz(-1.0875965) q[1];
sx q[1];
rz(0.35983905) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7774178) q[0];
sx q[0];
rz(-2.4628277) q[0];
sx q[0];
rz(-1.8666408) q[0];
rz(-pi) q[1];
rz(-2.4904576) q[2];
sx q[2];
rz(-1.480181) q[2];
sx q[2];
rz(0.091094253) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.29405669) q[1];
sx q[1];
rz(-1.1514074) q[1];
sx q[1];
rz(1.9860752) q[1];
rz(-pi) q[2];
rz(-2.2877341) q[3];
sx q[3];
rz(-0.44970185) q[3];
sx q[3];
rz(-0.85292294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.5741253) q[2];
sx q[2];
rz(-1.987477) q[2];
sx q[2];
rz(2.3763669) q[2];
rz(2.3848173) q[3];
sx q[3];
rz(-0.59195834) q[3];
sx q[3];
rz(-2.9005907) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1446447) q[0];
sx q[0];
rz(-2.6648271) q[0];
sx q[0];
rz(-1.7255406) q[0];
rz(-0.36711806) q[1];
sx q[1];
rz(-1.3857625) q[1];
sx q[1];
rz(2.1062772) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.028776289) q[0];
sx q[0];
rz(-1.1338286) q[0];
sx q[0];
rz(-2.2321738) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5315227) q[2];
sx q[2];
rz(-1.0796667) q[2];
sx q[2];
rz(-0.66205762) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.6742179) q[1];
sx q[1];
rz(-1.4754703) q[1];
sx q[1];
rz(-1.6256888) q[1];
rz(-pi) q[2];
x q[2];
rz(0.64621332) q[3];
sx q[3];
rz(-0.34020243) q[3];
sx q[3];
rz(-2.7681805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.3570024) q[2];
sx q[2];
rz(-1.4214397) q[2];
sx q[2];
rz(-0.76888293) q[2];
rz(-2.8055577) q[3];
sx q[3];
rz(-2.3612645) q[3];
sx q[3];
rz(1.6736354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1795905) q[0];
sx q[0];
rz(-2.5080894) q[0];
sx q[0];
rz(1.996421) q[0];
rz(-2.0369453) q[1];
sx q[1];
rz(-1.2303753) q[1];
sx q[1];
rz(-2.9343658) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4743054) q[0];
sx q[0];
rz(-1.4315839) q[0];
sx q[0];
rz(-2.5179203) q[0];
rz(-pi) q[1];
rz(2.0358762) q[2];
sx q[2];
rz(-1.0208703) q[2];
sx q[2];
rz(-1.9351026) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8689649) q[1];
sx q[1];
rz(-2.7902863) q[1];
sx q[1];
rz(-0.79877324) q[1];
x q[2];
rz(1.3979785) q[3];
sx q[3];
rz(-0.76068766) q[3];
sx q[3];
rz(0.93483227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.3383639) q[2];
sx q[2];
rz(-0.88023606) q[2];
sx q[2];
rz(2.4198789) q[2];
rz(-1.2747814) q[3];
sx q[3];
rz(-1.5694247) q[3];
sx q[3];
rz(-0.26003626) q[3];
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
rz(2.8969144) q[0];
sx q[0];
rz(-1.5667863) q[0];
sx q[0];
rz(-0.73202837) q[0];
rz(-3.1320944) q[1];
sx q[1];
rz(-0.54832012) q[1];
sx q[1];
rz(0.2917372) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4371944) q[0];
sx q[0];
rz(-1.5196374) q[0];
sx q[0];
rz(-2.722446) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8740011) q[2];
sx q[2];
rz(-1.6723987) q[2];
sx q[2];
rz(-0.50819699) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.2597255) q[1];
sx q[1];
rz(-1.7181267) q[1];
sx q[1];
rz(-0.30585652) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.52280207) q[3];
sx q[3];
rz(-1.3005946) q[3];
sx q[3];
rz(-3.0772046) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.8774524) q[2];
sx q[2];
rz(-1.143127) q[2];
sx q[2];
rz(-2.3042802) q[2];
rz(-1.9705747) q[3];
sx q[3];
rz(-1.5191017) q[3];
sx q[3];
rz(-3.0795735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0983122) q[0];
sx q[0];
rz(-2.4536112) q[0];
sx q[0];
rz(-0.6638546) q[0];
rz(3.035416) q[1];
sx q[1];
rz(-2.5352434) q[1];
sx q[1];
rz(0.95867872) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3631358) q[0];
sx q[0];
rz(-1.0582557) q[0];
sx q[0];
rz(0.97495671) q[0];
rz(-pi) q[1];
rz(2.1088976) q[2];
sx q[2];
rz(-1.8128464) q[2];
sx q[2];
rz(0.48697105) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.465185) q[1];
sx q[1];
rz(-2.0945815) q[1];
sx q[1];
rz(-3.0182059) q[1];
rz(-pi) q[2];
x q[2];
rz(0.58416768) q[3];
sx q[3];
rz(-1.6730047) q[3];
sx q[3];
rz(-2.971873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.133193) q[2];
sx q[2];
rz(-1.9063176) q[2];
sx q[2];
rz(2.2224902) q[2];
rz(1.3778) q[3];
sx q[3];
rz(-1.2104687) q[3];
sx q[3];
rz(-1.1834043) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8063426) q[0];
sx q[0];
rz(-2.2882473) q[0];
sx q[0];
rz(-2.7046955) q[0];
rz(-0.70029744) q[1];
sx q[1];
rz(-1.4236139) q[1];
sx q[1];
rz(0.46554309) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3532928) q[0];
sx q[0];
rz(-2.4294937) q[0];
sx q[0];
rz(0.13029356) q[0];
rz(1.2836254) q[2];
sx q[2];
rz(-0.42771491) q[2];
sx q[2];
rz(-2.5123793) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.000396) q[1];
sx q[1];
rz(-1.9138412) q[1];
sx q[1];
rz(2.6617906) q[1];
rz(0.97165473) q[3];
sx q[3];
rz(-0.42704901) q[3];
sx q[3];
rz(-3.0640102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.7312701) q[2];
sx q[2];
rz(-2.6022544) q[2];
sx q[2];
rz(1.643606) q[2];
rz(0.20478976) q[3];
sx q[3];
rz(-2.1361165) q[3];
sx q[3];
rz(-2.8695316) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64514226) q[0];
sx q[0];
rz(-2.5800939) q[0];
sx q[0];
rz(2.0196594) q[0];
rz(2.3902068) q[1];
sx q[1];
rz(-1.6151927) q[1];
sx q[1];
rz(-0.5823935) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96888992) q[0];
sx q[0];
rz(-0.70969289) q[0];
sx q[0];
rz(1.885528) q[0];
rz(1.8495314) q[2];
sx q[2];
rz(-0.50694743) q[2];
sx q[2];
rz(-2.7915733) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.45422428) q[1];
sx q[1];
rz(-0.2174938) q[1];
sx q[1];
rz(0.99517676) q[1];
rz(-1.5428316) q[3];
sx q[3];
rz(-1.7879322) q[3];
sx q[3];
rz(2.12487) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.0104684) q[2];
sx q[2];
rz(-0.98103395) q[2];
sx q[2];
rz(-0.76254145) q[2];
rz(-1.4108346) q[3];
sx q[3];
rz(-0.91791955) q[3];
sx q[3];
rz(-1.5677174) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0653771) q[0];
sx q[0];
rz(-2.1552754) q[0];
sx q[0];
rz(-1.7393204) q[0];
rz(-1.8021884) q[1];
sx q[1];
rz(-1.4504455) q[1];
sx q[1];
rz(1.6557678) q[1];
rz(-0.81417685) q[2];
sx q[2];
rz(-1.7780751) q[2];
sx q[2];
rz(1.8572394) q[2];
rz(2.0740261) q[3];
sx q[3];
rz(-1.9964841) q[3];
sx q[3];
rz(-1.869429) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
