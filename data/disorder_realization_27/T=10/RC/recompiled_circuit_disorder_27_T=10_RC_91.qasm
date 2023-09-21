OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-3.1448016) q[0];
sx q[0];
rz(0.15455833) q[0];
sx q[0];
rz(6.9757087) q[0];
rz(1.9321631) q[1];
sx q[1];
rz(-1.2485319) q[1];
sx q[1];
rz(1.7564397) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4614852) q[0];
sx q[0];
rz(-0.84851096) q[0];
sx q[0];
rz(1.1397584) q[0];
rz(-pi) q[1];
rz(-1.8458457) q[2];
sx q[2];
rz(-1.0359456) q[2];
sx q[2];
rz(-1.5616852) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4298809) q[1];
sx q[1];
rz(-1.3379339) q[1];
sx q[1];
rz(-1.7377322) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4888943) q[3];
sx q[3];
rz(-1.1678809) q[3];
sx q[3];
rz(-0.17890113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2549071) q[2];
sx q[2];
rz(-0.79780769) q[2];
sx q[2];
rz(-2.936426) q[2];
rz(0.77130476) q[3];
sx q[3];
rz(-0.78273928) q[3];
sx q[3];
rz(-1.1024968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40760621) q[0];
sx q[0];
rz(-2.3953231) q[0];
sx q[0];
rz(-2.6876887) q[0];
rz(1.0247963) q[1];
sx q[1];
rz(-0.4075993) q[1];
sx q[1];
rz(-1.227238) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42392143) q[0];
sx q[0];
rz(-2.9917891) q[0];
sx q[0];
rz(-1.040209) q[0];
rz(0.71364673) q[2];
sx q[2];
rz(-0.64672856) q[2];
sx q[2];
rz(-1.5681859) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.2482359) q[1];
sx q[1];
rz(-1.3125988) q[1];
sx q[1];
rz(-2.8398819) q[1];
rz(-pi) q[2];
rz(-1.4144054) q[3];
sx q[3];
rz(-1.7692493) q[3];
sx q[3];
rz(2.1895777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(3.1002905) q[2];
sx q[2];
rz(-1.1854478) q[2];
sx q[2];
rz(-2.5741637) q[2];
rz(2.7764017) q[3];
sx q[3];
rz(-1.4130211) q[3];
sx q[3];
rz(-2.173483) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.658618) q[0];
sx q[0];
rz(-0.56476074) q[0];
sx q[0];
rz(2.2429402) q[0];
rz(0.99575106) q[1];
sx q[1];
rz(-1.5581222) q[1];
sx q[1];
rz(-2.8083037) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22526564) q[0];
sx q[0];
rz(-0.81575459) q[0];
sx q[0];
rz(2.4832721) q[0];
x q[1];
rz(0.097924175) q[2];
sx q[2];
rz(-1.3946748) q[2];
sx q[2];
rz(-0.91913659) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.7283199) q[1];
sx q[1];
rz(-1.4188758) q[1];
sx q[1];
rz(-2.1876213) q[1];
rz(-pi) q[2];
rz(-0.40368747) q[3];
sx q[3];
rz(-1.0567769) q[3];
sx q[3];
rz(-2.6727563) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.68625346) q[2];
sx q[2];
rz(-1.7909966) q[2];
sx q[2];
rz(1.1509482) q[2];
rz(-2.3006556) q[3];
sx q[3];
rz(-2.0261814) q[3];
sx q[3];
rz(1.1545198) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
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
rz(-2.6338449) q[1];
sx q[1];
rz(-1.9365786) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8367774) q[0];
sx q[0];
rz(-2.4495821) q[0];
sx q[0];
rz(1.5185604) q[0];
rz(-0.20385216) q[2];
sx q[2];
rz(-2.2050369) q[2];
sx q[2];
rz(-0.39433345) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.64441427) q[1];
sx q[1];
rz(-2.2481611) q[1];
sx q[1];
rz(0.35269423) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6550001) q[3];
sx q[3];
rz(-2.4349672) q[3];
sx q[3];
rz(-3.1228309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.24923199) q[2];
sx q[2];
rz(-1.426733) q[2];
sx q[2];
rz(2.7704346) q[2];
rz(-1.7403729) q[3];
sx q[3];
rz(-2.4818082) q[3];
sx q[3];
rz(1.1192809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.054759653) q[0];
sx q[0];
rz(-0.78614569) q[0];
sx q[0];
rz(3.0084685) q[0];
rz(2.1482824) q[1];
sx q[1];
rz(-1.7555833) q[1];
sx q[1];
rz(2.5865119) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0017437) q[0];
sx q[0];
rz(-1.2554597) q[0];
sx q[0];
rz(-3.128201) q[0];
x q[1];
rz(-1.985717) q[2];
sx q[2];
rz(-1.6332111) q[2];
sx q[2];
rz(1.2635363) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.86075393) q[1];
sx q[1];
rz(-2.6959531) q[1];
sx q[1];
rz(-1.526236) q[1];
x q[2];
rz(-1.5186148) q[3];
sx q[3];
rz(-1.6515886) q[3];
sx q[3];
rz(2.3063456) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.30620265) q[2];
sx q[2];
rz(-2.1357048) q[2];
sx q[2];
rz(-0.13892697) q[2];
rz(-2.1991918) q[3];
sx q[3];
rz(-1.5709632) q[3];
sx q[3];
rz(0.55148235) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5979364) q[0];
sx q[0];
rz(-0.56977001) q[0];
sx q[0];
rz(2.561835) q[0];
rz(-0.12750164) q[1];
sx q[1];
rz(-1.9516877) q[1];
sx q[1];
rz(1.5396083) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72551661) q[0];
sx q[0];
rz(-2.2726739) q[0];
sx q[0];
rz(0.26728018) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8069205) q[2];
sx q[2];
rz(-3.0034608) q[2];
sx q[2];
rz(1.3484671) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.3967267) q[1];
sx q[1];
rz(-2.6187881) q[1];
sx q[1];
rz(-2.1641939) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9220011) q[3];
sx q[3];
rz(-0.71577365) q[3];
sx q[3];
rz(-2.2367665) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.55398983) q[2];
sx q[2];
rz(-2.8911399) q[2];
sx q[2];
rz(0.26947752) q[2];
rz(-2.907471) q[3];
sx q[3];
rz(-0.52474371) q[3];
sx q[3];
rz(-0.060119303) q[3];
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
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44678974) q[0];
sx q[0];
rz(-2.5233874) q[0];
sx q[0];
rz(0.68429464) q[0];
rz(-3.0220095) q[1];
sx q[1];
rz(-1.2922829) q[1];
sx q[1];
rz(0.51876846) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25012384) q[0];
sx q[0];
rz(-0.84202535) q[0];
sx q[0];
rz(0.17983371) q[0];
rz(0.22612818) q[2];
sx q[2];
rz(-2.3580708) q[2];
sx q[2];
rz(-2.4353611) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.59473945) q[1];
sx q[1];
rz(-0.93187983) q[1];
sx q[1];
rz(2.1584595) q[1];
rz(-1.9714868) q[3];
sx q[3];
rz(-1.5948933) q[3];
sx q[3];
rz(-1.8322242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.2542904) q[2];
sx q[2];
rz(-1.7742949) q[2];
sx q[2];
rz(0.56345144) q[2];
rz(-0.051579483) q[3];
sx q[3];
rz(-2.0023465) q[3];
sx q[3];
rz(-2.1896867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96034399) q[0];
sx q[0];
rz(-0.5287756) q[0];
sx q[0];
rz(1.7425591) q[0];
rz(-0.78701204) q[1];
sx q[1];
rz(-1.1287289) q[1];
sx q[1];
rz(-0.74434892) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5520185) q[0];
sx q[0];
rz(-0.15757832) q[0];
sx q[0];
rz(-0.35430674) q[0];
rz(2.9888399) q[2];
sx q[2];
rz(-2.0803917) q[2];
sx q[2];
rz(-2.5324412) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8213615) q[1];
sx q[1];
rz(-1.7944971) q[1];
sx q[1];
rz(-2.6279468) q[1];
rz(-0.45458557) q[3];
sx q[3];
rz(-1.5919519) q[3];
sx q[3];
rz(-1.8225614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.7156334) q[2];
sx q[2];
rz(-1.3954433) q[2];
sx q[2];
rz(-0.60950935) q[2];
rz(2.4842747) q[3];
sx q[3];
rz(-0.64353639) q[3];
sx q[3];
rz(-2.8801584) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1881926) q[0];
sx q[0];
rz(-0.094390079) q[0];
sx q[0];
rz(-1.6375861) q[0];
rz(1.9001182) q[1];
sx q[1];
rz(-1.1499317) q[1];
sx q[1];
rz(2.3666568) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0938213) q[0];
sx q[0];
rz(-1.1822961) q[0];
sx q[0];
rz(-0.85431487) q[0];
x q[1];
rz(1.2504134) q[2];
sx q[2];
rz(-0.46386007) q[2];
sx q[2];
rz(1.3055717) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.6248524) q[1];
sx q[1];
rz(-1.3827185) q[1];
sx q[1];
rz(2.6810758) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.14258607) q[3];
sx q[3];
rz(-1.3538059) q[3];
sx q[3];
rz(0.91689527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.9541786) q[2];
sx q[2];
rz(-2.9118907) q[2];
sx q[2];
rz(-2.9837218) q[2];
rz(1.212451) q[3];
sx q[3];
rz(-1.0586497) q[3];
sx q[3];
rz(-1.8635748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-1.071844) q[0];
sx q[0];
rz(-0.97244111) q[0];
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
rz(2.8911274) q[0];
sx q[0];
rz(-1.9292826) q[0];
sx q[0];
rz(-1.6842708) q[0];
rz(-pi) q[1];
rz(0.59098737) q[2];
sx q[2];
rz(-1.5948442) q[2];
sx q[2];
rz(1.1514593) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.5221094) q[1];
sx q[1];
rz(-1.932718) q[1];
sx q[1];
rz(-2.1543909) q[1];
x q[2];
rz(0.06185992) q[3];
sx q[3];
rz(-2.7354771) q[3];
sx q[3];
rz(1.5861685) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.41632286) q[2];
sx q[2];
rz(-3.0815093) q[2];
sx q[2];
rz(2.0521169) q[2];
rz(-1.5754835) q[3];
sx q[3];
rz(-1.8982866) q[3];
sx q[3];
rz(-0.65264788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2789223) q[0];
sx q[0];
rz(-2.537732) q[0];
sx q[0];
rz(-2.296007) q[0];
rz(1.5325585) q[1];
sx q[1];
rz(-1.6747723) q[1];
sx q[1];
rz(2.0369045) q[1];
rz(-2.5074742) q[2];
sx q[2];
rz(-0.49370439) q[2];
sx q[2];
rz(1.6903071) q[2];
rz(-1.6339176) q[3];
sx q[3];
rz(-1.0216542) q[3];
sx q[3];
rz(1.7246104) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];