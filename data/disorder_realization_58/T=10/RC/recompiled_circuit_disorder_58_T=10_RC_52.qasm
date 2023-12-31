OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.20237246) q[0];
sx q[0];
rz(-2.7352754) q[0];
sx q[0];
rz(2.321474) q[0];
rz(-0.36110538) q[1];
sx q[1];
rz(-2.5087924) q[1];
sx q[1];
rz(-2.3109205) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9554064) q[0];
sx q[0];
rz(-2.6468228) q[0];
sx q[0];
rz(-2.4525989) q[0];
rz(-2.0482424) q[2];
sx q[2];
rz(-2.2508143) q[2];
sx q[2];
rz(0.29104656) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.7100916) q[1];
sx q[1];
rz(-0.7845062) q[1];
sx q[1];
rz(-1.1169408) q[1];
rz(-1.9740231) q[3];
sx q[3];
rz(-1.9864169) q[3];
sx q[3];
rz(-0.88413873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.43710199) q[2];
sx q[2];
rz(-2.7097242) q[2];
sx q[2];
rz(3.0214156) q[2];
rz(-1.1581356) q[3];
sx q[3];
rz(-1.7430256) q[3];
sx q[3];
rz(2.2226298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(3.0607818) q[0];
sx q[0];
rz(-1.3588384) q[0];
sx q[0];
rz(-0.91180116) q[0];
rz(0.78951019) q[1];
sx q[1];
rz(-2.1531838) q[1];
sx q[1];
rz(-2.8149014) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9625898) q[0];
sx q[0];
rz(-1.0674745) q[0];
sx q[0];
rz(-2.3110564) q[0];
x q[1];
rz(-2.0147444) q[2];
sx q[2];
rz(-1.4269281) q[2];
sx q[2];
rz(2.4614046) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.82636278) q[1];
sx q[1];
rz(-2.3781207) q[1];
sx q[1];
rz(-0.96674322) q[1];
rz(-pi) q[2];
x q[2];
rz(0.88300206) q[3];
sx q[3];
rz(-1.5044754) q[3];
sx q[3];
rz(1.6460713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.5102753) q[2];
sx q[2];
rz(-0.77284208) q[2];
sx q[2];
rz(1.8544244) q[2];
rz(3.0316947) q[3];
sx q[3];
rz(-1.4108312) q[3];
sx q[3];
rz(1.7597642) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54863769) q[0];
sx q[0];
rz(-2.4023963) q[0];
sx q[0];
rz(-0.32989311) q[0];
rz(2.864481) q[1];
sx q[1];
rz(-1.3169293) q[1];
sx q[1];
rz(-1.057391) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1304504) q[0];
sx q[0];
rz(-1.5913977) q[0];
sx q[0];
rz(-1.9111454) q[0];
rz(-pi) q[1];
rz(0.51809394) q[2];
sx q[2];
rz(-2.7535541) q[2];
sx q[2];
rz(-2.5990017) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.33144618) q[1];
sx q[1];
rz(-0.76008893) q[1];
sx q[1];
rz(2.0600832) q[1];
rz(-pi) q[2];
x q[2];
rz(0.40773817) q[3];
sx q[3];
rz(-0.45555112) q[3];
sx q[3];
rz(-1.7689442) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3588336) q[2];
sx q[2];
rz(-1.1153355) q[2];
sx q[2];
rz(-1.7791629) q[2];
rz(0.6247012) q[3];
sx q[3];
rz(-1.0995068) q[3];
sx q[3];
rz(0.81956285) q[3];
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
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7841004) q[0];
sx q[0];
rz(-2.6153013) q[0];
sx q[0];
rz(-2.6065361) q[0];
rz(2.0013981) q[1];
sx q[1];
rz(-1.3277206) q[1];
sx q[1];
rz(-0.16539703) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7968984) q[0];
sx q[0];
rz(-1.3860774) q[0];
sx q[0];
rz(-3.0152507) q[0];
rz(-pi) q[1];
rz(-2.303316) q[2];
sx q[2];
rz(-2.8998313) q[2];
sx q[2];
rz(-1.4571112) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.71388984) q[1];
sx q[1];
rz(-1.0298924) q[1];
sx q[1];
rz(1.9546024) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9966647) q[3];
sx q[3];
rz(-0.59755675) q[3];
sx q[3];
rz(1.4262517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.39607221) q[2];
sx q[2];
rz(-1.8053651) q[2];
sx q[2];
rz(2.9361434) q[2];
rz(-2.0139587) q[3];
sx q[3];
rz(-1.1536359) q[3];
sx q[3];
rz(-1.2566465) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66185343) q[0];
sx q[0];
rz(-0.91402188) q[0];
sx q[0];
rz(2.7868295) q[0];
rz(-1.154249) q[1];
sx q[1];
rz(-2.216979) q[1];
sx q[1];
rz(-2.9096471) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1922798) q[0];
sx q[0];
rz(-0.67040196) q[0];
sx q[0];
rz(2.288726) q[0];
rz(-pi) q[1];
rz(2.4460692) q[2];
sx q[2];
rz(-1.6791653) q[2];
sx q[2];
rz(-2.7930789) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.0143118) q[1];
sx q[1];
rz(-2.8297272) q[1];
sx q[1];
rz(-0.77906268) q[1];
x q[2];
rz(2.7730745) q[3];
sx q[3];
rz(-1.5889865) q[3];
sx q[3];
rz(-0.25572488) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0014687) q[2];
sx q[2];
rz(-1.8073558) q[2];
sx q[2];
rz(-1.9011964) q[2];
rz(-0.59605789) q[3];
sx q[3];
rz(-1.8362935) q[3];
sx q[3];
rz(1.8381455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0444788) q[0];
sx q[0];
rz(-0.070274027) q[0];
sx q[0];
rz(2.9388359) q[0];
rz(-2.1525106) q[1];
sx q[1];
rz(-1.6977856) q[1];
sx q[1];
rz(-2.1441377) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8231682) q[0];
sx q[0];
rz(-2.423375) q[0];
sx q[0];
rz(-1.3658001) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4100254) q[2];
sx q[2];
rz(-1.288207) q[2];
sx q[2];
rz(-2.3716795) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.84038489) q[1];
sx q[1];
rz(-1.9921229) q[1];
sx q[1];
rz(0.2627443) q[1];
rz(2.0632083) q[3];
sx q[3];
rz(-0.78740722) q[3];
sx q[3];
rz(-1.9099964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.22770195) q[2];
sx q[2];
rz(-1.1662377) q[2];
sx q[2];
rz(2.690199) q[2];
rz(2.732892) q[3];
sx q[3];
rz(-1.6025851) q[3];
sx q[3];
rz(-1.9394978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8354427) q[0];
sx q[0];
rz(-2.7375484) q[0];
sx q[0];
rz(-2.5174482) q[0];
rz(-1.6250601) q[1];
sx q[1];
rz(-2.8846278) q[1];
sx q[1];
rz(0.61378941) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81543024) q[0];
sx q[0];
rz(-0.76061941) q[0];
sx q[0];
rz(-1.2994231) q[0];
rz(-pi) q[1];
rz(-2.9250547) q[2];
sx q[2];
rz(-0.29007402) q[2];
sx q[2];
rz(-1.1494344) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.1544513) q[1];
sx q[1];
rz(-0.11949355) q[1];
sx q[1];
rz(-2.6277072) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8011439) q[3];
sx q[3];
rz(-0.47603546) q[3];
sx q[3];
rz(1.3766833) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.43341407) q[2];
sx q[2];
rz(-2.2257979) q[2];
sx q[2];
rz(-1.1748574) q[2];
rz(0.60837778) q[3];
sx q[3];
rz(-1.404168) q[3];
sx q[3];
rz(1.7181989) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2414918) q[0];
sx q[0];
rz(-1.2986203) q[0];
sx q[0];
rz(-0.4883782) q[0];
rz(1.5178559) q[1];
sx q[1];
rz(-1.7428215) q[1];
sx q[1];
rz(-0.98446313) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0122089) q[0];
sx q[0];
rz(-1.0249656) q[0];
sx q[0];
rz(-1.6145541) q[0];
rz(-pi) q[1];
rz(1.9756873) q[2];
sx q[2];
rz(-0.84120175) q[2];
sx q[2];
rz(1.0798432) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.8967646) q[1];
sx q[1];
rz(-0.36537376) q[1];
sx q[1];
rz(0.016896292) q[1];
rz(-pi) q[2];
rz(-2.5823309) q[3];
sx q[3];
rz(-0.54702938) q[3];
sx q[3];
rz(2.3032041) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4349334) q[2];
sx q[2];
rz(-1.443053) q[2];
sx q[2];
rz(2.6064176) q[2];
rz(-2.0914071) q[3];
sx q[3];
rz(-1.340056) q[3];
sx q[3];
rz(2.1323269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.500279) q[0];
sx q[0];
rz(-2.2135493) q[0];
sx q[0];
rz(-2.4556659) q[0];
rz(-0.39086875) q[1];
sx q[1];
rz(-0.94376826) q[1];
sx q[1];
rz(2.2156782) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6305011) q[0];
sx q[0];
rz(-1.7112268) q[0];
sx q[0];
rz(1.7846084) q[0];
x q[1];
rz(-2.3472896) q[2];
sx q[2];
rz(-2.1199193) q[2];
sx q[2];
rz(-1.4505475) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.8030501) q[1];
sx q[1];
rz(-2.0085137) q[1];
sx q[1];
rz(0.16191698) q[1];
x q[2];
rz(-1.4521452) q[3];
sx q[3];
rz(-1.7735529) q[3];
sx q[3];
rz(-2.6076536) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.19568504) q[2];
sx q[2];
rz(-2.6699799) q[2];
sx q[2];
rz(-2.6055028) q[2];
rz(-2.7219971) q[3];
sx q[3];
rz(-2.5511238) q[3];
sx q[3];
rz(1.6528116) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0855899) q[0];
sx q[0];
rz(-0.35878006) q[0];
sx q[0];
rz(2.7465903) q[0];
rz(1.5123873) q[1];
sx q[1];
rz(-1.869166) q[1];
sx q[1];
rz(2.1283456) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.917406) q[0];
sx q[0];
rz(-1.4316443) q[0];
sx q[0];
rz(-1.3361206) q[0];
x q[1];
rz(3.0229438) q[2];
sx q[2];
rz(-2.7855706) q[2];
sx q[2];
rz(-2.7729386) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0638748) q[1];
sx q[1];
rz(-0.43242726) q[1];
sx q[1];
rz(-2.9090911) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4471531) q[3];
sx q[3];
rz(-0.89530066) q[3];
sx q[3];
rz(1.465786) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.44832486) q[2];
sx q[2];
rz(-1.5193181) q[2];
sx q[2];
rz(-0.75941336) q[2];
rz(1.7761207) q[3];
sx q[3];
rz(-1.1080192) q[3];
sx q[3];
rz(2.571648) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41291819) q[0];
sx q[0];
rz(-1.8771794) q[0];
sx q[0];
rz(-1.9369453) q[0];
rz(-2.5683174) q[1];
sx q[1];
rz(-1.9075867) q[1];
sx q[1];
rz(1.8214068) q[1];
rz(2.0157094) q[2];
sx q[2];
rz(-1.5013668) q[2];
sx q[2];
rz(-2.8253386) q[2];
rz(3.0782386) q[3];
sx q[3];
rz(-2.2676716) q[3];
sx q[3];
rz(0.36361658) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
