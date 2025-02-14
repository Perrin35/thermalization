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
rz(0.32349411) q[0];
sx q[0];
rz(-0.20322023) q[0];
sx q[0];
rz(0.30015552) q[0];
rz(-0.35429859) q[1];
sx q[1];
rz(-1.9072671) q[1];
sx q[1];
rz(0.77114463) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82282372) q[0];
sx q[0];
rz(-1.7576801) q[0];
sx q[0];
rz(-2.9479638) q[0];
rz(2.6080898) q[2];
sx q[2];
rz(-1.9868104) q[2];
sx q[2];
rz(2.5737284) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.59471666) q[1];
sx q[1];
rz(-1.1961814) q[1];
sx q[1];
rz(1.131011) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0047831) q[3];
sx q[3];
rz(-1.5094637) q[3];
sx q[3];
rz(-0.086333806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.99042088) q[2];
sx q[2];
rz(-0.46151084) q[2];
sx q[2];
rz(-2.0578461) q[2];
rz(0.57461965) q[3];
sx q[3];
rz(-1.5465522) q[3];
sx q[3];
rz(0.77742022) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3144749) q[0];
sx q[0];
rz(-0.62905335) q[0];
sx q[0];
rz(2.7431059) q[0];
rz(-1.9972948) q[1];
sx q[1];
rz(-2.5356348) q[1];
sx q[1];
rz(0.95432895) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1935611) q[0];
sx q[0];
rz(-1.5753509) q[0];
sx q[0];
rz(-1.5702374) q[0];
rz(-pi) q[1];
rz(1.6923087) q[2];
sx q[2];
rz(-0.6612311) q[2];
sx q[2];
rz(2.0362771) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1283537) q[1];
sx q[1];
rz(-1.8132105) q[1];
sx q[1];
rz(2.1552312) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4407473) q[3];
sx q[3];
rz(-2.4748445) q[3];
sx q[3];
rz(1.4481883) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.9053456) q[2];
sx q[2];
rz(-2.6409642) q[2];
sx q[2];
rz(3.0312669) q[2];
rz(-2.3162383) q[3];
sx q[3];
rz(-2.3767411) q[3];
sx q[3];
rz(-2.4013405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1684619) q[0];
sx q[0];
rz(-0.2178807) q[0];
sx q[0];
rz(2.2595898) q[0];
rz(-2.1130256) q[1];
sx q[1];
rz(-2.5879637) q[1];
sx q[1];
rz(0.91241589) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4140625) q[0];
sx q[0];
rz(-1.4986218) q[0];
sx q[0];
rz(-1.1185924) q[0];
x q[1];
rz(-0.71850168) q[2];
sx q[2];
rz(-1.6616491) q[2];
sx q[2];
rz(-1.529734) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.6710224) q[1];
sx q[1];
rz(-2.467944) q[1];
sx q[1];
rz(-1.0224875) q[1];
rz(2.8863593) q[3];
sx q[3];
rz(-1.1369656) q[3];
sx q[3];
rz(-1.1800769) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.87938157) q[2];
sx q[2];
rz(-2.8935581) q[2];
sx q[2];
rz(-0.81825078) q[2];
rz(-0.36000559) q[3];
sx q[3];
rz(-1.8428948) q[3];
sx q[3];
rz(-0.037809614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70571947) q[0];
sx q[0];
rz(-2.191045) q[0];
sx q[0];
rz(1.1608423) q[0];
rz(-2.1707161) q[1];
sx q[1];
rz(-0.91578805) q[1];
sx q[1];
rz(-0.11347778) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.023983) q[0];
sx q[0];
rz(-1.9736088) q[0];
sx q[0];
rz(-2.4871123) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9275437) q[2];
sx q[2];
rz(-0.037790701) q[2];
sx q[2];
rz(2.1434857) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.8673528) q[1];
sx q[1];
rz(-0.93207658) q[1];
sx q[1];
rz(-0.8886877) q[1];
rz(-pi) q[2];
rz(0.63469074) q[3];
sx q[3];
rz(-1.029976) q[3];
sx q[3];
rz(1.5405003) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.11161441) q[2];
sx q[2];
rz(-1.3361715) q[2];
sx q[2];
rz(2.2335936) q[2];
rz(-0.32634431) q[3];
sx q[3];
rz(-0.55485266) q[3];
sx q[3];
rz(-2.7197796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(0.57080764) q[0];
sx q[0];
rz(-0.76199216) q[0];
sx q[0];
rz(-2.1245891) q[0];
rz(0.48655888) q[1];
sx q[1];
rz(-2.1162972) q[1];
sx q[1];
rz(-3.0466381) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61750283) q[0];
sx q[0];
rz(-1.4727778) q[0];
sx q[0];
rz(1.4630415) q[0];
rz(-pi) q[1];
rz(2.9385511) q[2];
sx q[2];
rz(-2.3066694) q[2];
sx q[2];
rz(-2.1020213) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.50852114) q[1];
sx q[1];
rz(-1.8446559) q[1];
sx q[1];
rz(2.7809596) q[1];
x q[2];
rz(2.8311555) q[3];
sx q[3];
rz(-1.2708029) q[3];
sx q[3];
rz(2.7315549) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.93593705) q[2];
sx q[2];
rz(-0.63465261) q[2];
sx q[2];
rz(2.671799) q[2];
rz(-0.69822407) q[3];
sx q[3];
rz(-1.2263115) q[3];
sx q[3];
rz(-0.32018143) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6093589) q[0];
sx q[0];
rz(-0.73901743) q[0];
sx q[0];
rz(-0.21482378) q[0];
rz(-0.23669446) q[1];
sx q[1];
rz(-2.5645945) q[1];
sx q[1];
rz(-2.7223041) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40060166) q[0];
sx q[0];
rz(-0.090734847) q[0];
sx q[0];
rz(1.0450793) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0299267) q[2];
sx q[2];
rz(-1.8395784) q[2];
sx q[2];
rz(2.2044685) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.57945573) q[1];
sx q[1];
rz(-0.39135763) q[1];
sx q[1];
rz(1.0117024) q[1];
x q[2];
rz(2.9131439) q[3];
sx q[3];
rz(-0.72586021) q[3];
sx q[3];
rz(1.5184234) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.99570167) q[2];
sx q[2];
rz(-3.006101) q[2];
sx q[2];
rz(0.091391407) q[2];
rz(0.35178301) q[3];
sx q[3];
rz(-2.3942949) q[3];
sx q[3];
rz(-0.465213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2202989) q[0];
sx q[0];
rz(-1.969368) q[0];
sx q[0];
rz(1.175513) q[0];
rz(-0.57918864) q[1];
sx q[1];
rz(-0.80668956) q[1];
sx q[1];
rz(-2.9520292) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2611147) q[0];
sx q[0];
rz(-1.7715329) q[0];
sx q[0];
rz(-1.6693382) q[0];
rz(-pi) q[1];
rz(0.86935027) q[2];
sx q[2];
rz(-1.6371718) q[2];
sx q[2];
rz(-0.36176031) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.71566641) q[1];
sx q[1];
rz(-2.3235011) q[1];
sx q[1];
rz(-3.1311684) q[1];
x q[2];
rz(1.0503672) q[3];
sx q[3];
rz(-1.6230021) q[3];
sx q[3];
rz(-2.5676651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.2438948) q[2];
sx q[2];
rz(-1.8420668) q[2];
sx q[2];
rz(2.6356836) q[2];
rz(-0.77887744) q[3];
sx q[3];
rz(-0.39611045) q[3];
sx q[3];
rz(1.2946543) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8706354) q[0];
sx q[0];
rz(-1.0519692) q[0];
sx q[0];
rz(-1.9860995) q[0];
rz(3.1404066) q[1];
sx q[1];
rz(-2.3540034) q[1];
sx q[1];
rz(-2.7364065) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7195994) q[0];
sx q[0];
rz(-1.1047537) q[0];
sx q[0];
rz(-2.5525722) q[0];
rz(-2.9878163) q[2];
sx q[2];
rz(-1.9755873) q[2];
sx q[2];
rz(0.60196934) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.5444774) q[1];
sx q[1];
rz(-1.4946723) q[1];
sx q[1];
rz(0.80634574) q[1];
rz(-1.5022329) q[3];
sx q[3];
rz(-1.4430832) q[3];
sx q[3];
rz(-2.8721916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.323641) q[2];
sx q[2];
rz(-1.1672856) q[2];
sx q[2];
rz(-2.6200068) q[2];
rz(0.096605435) q[3];
sx q[3];
rz(-2.1229027) q[3];
sx q[3];
rz(-0.49083138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38799649) q[0];
sx q[0];
rz(-0.8803519) q[0];
sx q[0];
rz(0.11465797) q[0];
rz(-1.8278587) q[1];
sx q[1];
rz(-1.6856245) q[1];
sx q[1];
rz(-3.0065261) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0226267) q[0];
sx q[0];
rz(-2.8014136) q[0];
sx q[0];
rz(-0.69891286) q[0];
rz(-pi) q[1];
rz(0.39188633) q[2];
sx q[2];
rz(-1.8183437) q[2];
sx q[2];
rz(2.4154169) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.8986021) q[1];
sx q[1];
rz(-2.1861316) q[1];
sx q[1];
rz(-0.19886417) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0681048) q[3];
sx q[3];
rz(-1.013143) q[3];
sx q[3];
rz(-0.36493767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.8257398) q[2];
sx q[2];
rz(-1.9261381) q[2];
sx q[2];
rz(-0.94557554) q[2];
rz(-2.8816176) q[3];
sx q[3];
rz(-1.0228478) q[3];
sx q[3];
rz(2.0135349) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.05397) q[0];
sx q[0];
rz(-1.0902417) q[0];
sx q[0];
rz(1.4029652) q[0];
rz(-2.2258017) q[1];
sx q[1];
rz(-0.61640888) q[1];
sx q[1];
rz(-0.91555196) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4175891) q[0];
sx q[0];
rz(-1.6678828) q[0];
sx q[0];
rz(-1.4396981) q[0];
rz(-pi) q[1];
rz(1.0491285) q[2];
sx q[2];
rz(-1.4031488) q[2];
sx q[2];
rz(-2.0190792) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1656395) q[1];
sx q[1];
rz(-1.136186) q[1];
sx q[1];
rz(-2.6967816) q[1];
rz(-pi) q[2];
rz(1.7677018) q[3];
sx q[3];
rz(-2.422296) q[3];
sx q[3];
rz(-0.48433892) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.12330684) q[2];
sx q[2];
rz(-2.4985963) q[2];
sx q[2];
rz(2.7389738) q[2];
rz(1.9649402) q[3];
sx q[3];
rz(-2.8549356) q[3];
sx q[3];
rz(0.76796842) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1422414) q[0];
sx q[0];
rz(-2.169123) q[0];
sx q[0];
rz(2.1219667) q[0];
rz(2.453852) q[1];
sx q[1];
rz(-2.3383457) q[1];
sx q[1];
rz(1.8170423) q[1];
rz(-0.8310995) q[2];
sx q[2];
rz(-0.88242006) q[2];
sx q[2];
rz(-0.77280828) q[2];
rz(1.3152033) q[3];
sx q[3];
rz(-0.9556613) q[3];
sx q[3];
rz(-1.2355604) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
