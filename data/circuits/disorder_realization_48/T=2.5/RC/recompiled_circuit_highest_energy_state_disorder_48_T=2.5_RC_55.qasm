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
rz(-0.35204044) q[0];
sx q[0];
rz(-0.8249324) q[0];
sx q[0];
rz(0.53476778) q[0];
rz(0.98400247) q[1];
sx q[1];
rz(-2.6360631) q[1];
sx q[1];
rz(1.2619789) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6184632) q[0];
sx q[0];
rz(-0.97181706) q[0];
sx q[0];
rz(1.8904786) q[0];
rz(0.87902576) q[2];
sx q[2];
rz(-2.3874385) q[2];
sx q[2];
rz(1.1471495) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.76319219) q[1];
sx q[1];
rz(-1.0240558) q[1];
sx q[1];
rz(1.631447) q[1];
rz(-2.6747236) q[3];
sx q[3];
rz(-2.7109466) q[3];
sx q[3];
rz(2.1809585) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.24903211) q[2];
sx q[2];
rz(-0.51237115) q[2];
sx q[2];
rz(1.6585635) q[2];
rz(2.856971) q[3];
sx q[3];
rz(-0.62323815) q[3];
sx q[3];
rz(2.0774138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2738709) q[0];
sx q[0];
rz(-0.72672788) q[0];
sx q[0];
rz(-3.100585) q[0];
rz(1.2076591) q[1];
sx q[1];
rz(-2.9137847) q[1];
sx q[1];
rz(1.6809195) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2940833) q[0];
sx q[0];
rz(-2.9929711) q[0];
sx q[0];
rz(0.20047184) q[0];
x q[1];
rz(-0.56378491) q[2];
sx q[2];
rz(-0.85104686) q[2];
sx q[2];
rz(-1.189718) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.25819689) q[1];
sx q[1];
rz(-0.25942311) q[1];
sx q[1];
rz(1.4197465) q[1];
rz(-1.8051265) q[3];
sx q[3];
rz(-0.70584345) q[3];
sx q[3];
rz(0.88283759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.66551208) q[2];
sx q[2];
rz(-1.6802695) q[2];
sx q[2];
rz(1.3903138) q[2];
rz(0.0811854) q[3];
sx q[3];
rz(-2.472671) q[3];
sx q[3];
rz(1.4359052) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3944893) q[0];
sx q[0];
rz(-3.1293226) q[0];
sx q[0];
rz(0.61007208) q[0];
rz(1.3129781) q[1];
sx q[1];
rz(-2.4459631) q[1];
sx q[1];
rz(-2.5909766) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0673772) q[0];
sx q[0];
rz(-2.2226637) q[0];
sx q[0];
rz(2.0280736) q[0];
x q[1];
rz(2.4464452) q[2];
sx q[2];
rz(-0.41054976) q[2];
sx q[2];
rz(-1.6352959) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.6414601) q[1];
sx q[1];
rz(-1.271903) q[1];
sx q[1];
rz(-2.9445093) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7254225) q[3];
sx q[3];
rz(-0.75115132) q[3];
sx q[3];
rz(-0.42645833) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.47784558) q[2];
sx q[2];
rz(-2.9691594) q[2];
sx q[2];
rz(-2.0752068) q[2];
rz(1.1586698) q[3];
sx q[3];
rz(-1.7585157) q[3];
sx q[3];
rz(3.0995479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31024194) q[0];
sx q[0];
rz(-0.61315918) q[0];
sx q[0];
rz(-0.9170652) q[0];
rz(-0.65545583) q[1];
sx q[1];
rz(-0.90924811) q[1];
sx q[1];
rz(-0.38958946) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8389908) q[0];
sx q[0];
rz(-2.9377958) q[0];
sx q[0];
rz(-1.9857282) q[0];
rz(-pi) q[1];
rz(-2.2558446) q[2];
sx q[2];
rz(-2.0221524) q[2];
sx q[2];
rz(-1.3944172) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.0181737) q[1];
sx q[1];
rz(-1.5999854) q[1];
sx q[1];
rz(1.1694714) q[1];
x q[2];
rz(0.84221496) q[3];
sx q[3];
rz(-2.4024253) q[3];
sx q[3];
rz(0.30876866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.3098844) q[2];
sx q[2];
rz(-0.94253057) q[2];
sx q[2];
rz(1.3523098) q[2];
rz(-0.74812198) q[3];
sx q[3];
rz(-1.3040521) q[3];
sx q[3];
rz(1.5266533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8378976) q[0];
sx q[0];
rz(-0.56981531) q[0];
sx q[0];
rz(1.8001528) q[0];
rz(-0.98126423) q[1];
sx q[1];
rz(-1.860268) q[1];
sx q[1];
rz(2.431638) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34559743) q[0];
sx q[0];
rz(-1.5899204) q[0];
sx q[0];
rz(-2.8157793) q[0];
rz(1.2125489) q[2];
sx q[2];
rz(-0.60645559) q[2];
sx q[2];
rz(2.0403634) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6483575) q[1];
sx q[1];
rz(-1.4005737) q[1];
sx q[1];
rz(-1.3551718) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.59779915) q[3];
sx q[3];
rz(-1.6493268) q[3];
sx q[3];
rz(-0.30818916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.183341) q[2];
sx q[2];
rz(-0.7603344) q[2];
sx q[2];
rz(-0.90671268) q[2];
rz(-0.18236154) q[3];
sx q[3];
rz(-1.7014818) q[3];
sx q[3];
rz(-0.85842925) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7957423) q[0];
sx q[0];
rz(-2.1868571) q[0];
sx q[0];
rz(-0.20172754) q[0];
rz(1.3564823) q[1];
sx q[1];
rz(-0.67142612) q[1];
sx q[1];
rz(-1.1511525) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9485474) q[0];
sx q[0];
rz(-0.47848114) q[0];
sx q[0];
rz(-1.1043666) q[0];
rz(-0.60672673) q[2];
sx q[2];
rz(-1.9579685) q[2];
sx q[2];
rz(0.78578709) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.56111873) q[1];
sx q[1];
rz(-0.72787933) q[1];
sx q[1];
rz(0.62904686) q[1];
rz(2.6309507) q[3];
sx q[3];
rz(-1.770442) q[3];
sx q[3];
rz(0.24195237) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.98728937) q[2];
sx q[2];
rz(-0.58649784) q[2];
sx q[2];
rz(-0.34745535) q[2];
rz(2.3197428) q[3];
sx q[3];
rz(-1.3653267) q[3];
sx q[3];
rz(1.52389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8007941) q[0];
sx q[0];
rz(-2.4947385) q[0];
sx q[0];
rz(1.1232173) q[0];
rz(0.42876354) q[1];
sx q[1];
rz(-2.2518497) q[1];
sx q[1];
rz(-2.9926328) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87951794) q[0];
sx q[0];
rz(-1.5929758) q[0];
sx q[0];
rz(-2.0878077) q[0];
rz(-pi) q[1];
rz(2.5507418) q[2];
sx q[2];
rz(-0.35195165) q[2];
sx q[2];
rz(-2.1972434) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.625714) q[1];
sx q[1];
rz(-0.82059723) q[1];
sx q[1];
rz(1.7734115) q[1];
x q[2];
rz(-0.034986939) q[3];
sx q[3];
rz(-1.1652428) q[3];
sx q[3];
rz(1.4407033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1296156) q[2];
sx q[2];
rz(-0.35795438) q[2];
sx q[2];
rz(-0.32148662) q[2];
rz(-1.1164411) q[3];
sx q[3];
rz(-0.8588841) q[3];
sx q[3];
rz(-1.727625) q[3];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1534934) q[0];
sx q[0];
rz(-1.7970002) q[0];
sx q[0];
rz(-0.022911428) q[0];
rz(-1.5494391) q[1];
sx q[1];
rz(-1.0433082) q[1];
sx q[1];
rz(1.9291838) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5488729) q[0];
sx q[0];
rz(-0.61036068) q[0];
sx q[0];
rz(-1.1038853) q[0];
rz(-pi) q[1];
rz(0.83974482) q[2];
sx q[2];
rz(-1.296412) q[2];
sx q[2];
rz(1.7626761) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.072619) q[1];
sx q[1];
rz(-2.6311462) q[1];
sx q[1];
rz(-0.65565311) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6818265) q[3];
sx q[3];
rz(-2.4250406) q[3];
sx q[3];
rz(2.6727303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.93078485) q[2];
sx q[2];
rz(-0.19889861) q[2];
sx q[2];
rz(2.1442294) q[2];
rz(1.2584244) q[3];
sx q[3];
rz(-1.1910028) q[3];
sx q[3];
rz(1.2303801) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
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
rz(-2.9835984) q[0];
sx q[0];
rz(-0.65584922) q[0];
sx q[0];
rz(0.47437814) q[0];
rz(-2.5899218) q[1];
sx q[1];
rz(-2.3730979) q[1];
sx q[1];
rz(-2.4840568) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5112572) q[0];
sx q[0];
rz(-1.5310107) q[0];
sx q[0];
rz(-1.5471094) q[0];
rz(-pi) q[1];
rz(-1.1649407) q[2];
sx q[2];
rz(-1.7271583) q[2];
sx q[2];
rz(-0.49977068) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.20831693) q[1];
sx q[1];
rz(-2.7728656) q[1];
sx q[1];
rz(1.5255724) q[1];
rz(-2.1980638) q[3];
sx q[3];
rz(-0.17788685) q[3];
sx q[3];
rz(-1.8454144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.9115596) q[2];
sx q[2];
rz(-0.40731373) q[2];
sx q[2];
rz(-1.8302906) q[2];
rz(0.71410549) q[3];
sx q[3];
rz(-1.2823391) q[3];
sx q[3];
rz(-0.56984058) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1353726) q[0];
sx q[0];
rz(-2.1923809) q[0];
sx q[0];
rz(0.60829341) q[0];
rz(-2.3237806) q[1];
sx q[1];
rz(-1.1759956) q[1];
sx q[1];
rz(-2.3086595) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.333182) q[0];
sx q[0];
rz(-1.6224553) q[0];
sx q[0];
rz(-0.42148659) q[0];
x q[1];
rz(1.8274177) q[2];
sx q[2];
rz(-1.875281) q[2];
sx q[2];
rz(-1.2198795) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3994308) q[1];
sx q[1];
rz(-0.99282904) q[1];
sx q[1];
rz(1.7294243) q[1];
x q[2];
rz(-1.6021483) q[3];
sx q[3];
rz(-0.97797457) q[3];
sx q[3];
rz(-1.5167936) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.9587162) q[2];
sx q[2];
rz(-2.7887838) q[2];
sx q[2];
rz(0.58615169) q[2];
rz(0.90306774) q[3];
sx q[3];
rz(-2.51913) q[3];
sx q[3];
rz(-3.0557475) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2387954) q[0];
sx q[0];
rz(-2.0113404) q[0];
sx q[0];
rz(2.197862) q[0];
rz(-2.0441652) q[1];
sx q[1];
rz(-0.25249093) q[1];
sx q[1];
rz(2.9846334) q[1];
rz(-3.0022754) q[2];
sx q[2];
rz(-1.9092154) q[2];
sx q[2];
rz(-0.42001324) q[2];
rz(-1.0020574) q[3];
sx q[3];
rz(-0.78227038) q[3];
sx q[3];
rz(-3.0113358) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
