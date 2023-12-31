OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.5946755) q[0];
sx q[0];
rz(-1.0008873) q[0];
sx q[0];
rz(-0.21240182) q[0];
rz(-2.4266333) q[1];
sx q[1];
rz(-0.7874878) q[1];
sx q[1];
rz(1.8600872) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7377388) q[0];
sx q[0];
rz(-1.3699342) q[0];
sx q[0];
rz(1.2234729) q[0];
rz(-pi) q[1];
rz(1.6002866) q[2];
sx q[2];
rz(-2.2631096) q[2];
sx q[2];
rz(-2.8147547) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.233853) q[1];
sx q[1];
rz(-2.6820161) q[1];
sx q[1];
rz(-1.185226) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1164066) q[3];
sx q[3];
rz(-0.14305275) q[3];
sx q[3];
rz(1.5798626) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.32221258) q[2];
sx q[2];
rz(-0.11495049) q[2];
sx q[2];
rz(-1.5167351) q[2];
rz(0.24762282) q[3];
sx q[3];
rz(-1.473749) q[3];
sx q[3];
rz(2.4860399) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4441929) q[0];
sx q[0];
rz(-1.4780761) q[0];
sx q[0];
rz(-0.15629388) q[0];
rz(-0.63931757) q[1];
sx q[1];
rz(-1.0126065) q[1];
sx q[1];
rz(1.6527536) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50870897) q[0];
sx q[0];
rz(-1.7299621) q[0];
sx q[0];
rz(2.5863618) q[0];
rz(-pi) q[1];
rz(0.70398753) q[2];
sx q[2];
rz(-2.1805602) q[2];
sx q[2];
rz(-3.0733382) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.7557775) q[1];
sx q[1];
rz(-2.1375244) q[1];
sx q[1];
rz(-2.4789026) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5403455) q[3];
sx q[3];
rz(-1.7073365) q[3];
sx q[3];
rz(2.2583435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.2911825) q[2];
sx q[2];
rz(-3.0427142) q[2];
sx q[2];
rz(-2.8448811) q[2];
rz(-0.3324278) q[3];
sx q[3];
rz(-1.6231096) q[3];
sx q[3];
rz(-1.9077574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5082224) q[0];
sx q[0];
rz(-1.8828266) q[0];
sx q[0];
rz(0.53031522) q[0];
rz(0.5439533) q[1];
sx q[1];
rz(-1.1214316) q[1];
sx q[1];
rz(1.189032) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2169164) q[0];
sx q[0];
rz(-2.8042256) q[0];
sx q[0];
rz(-0.49305537) q[0];
rz(-pi) q[1];
x q[1];
rz(0.15709917) q[2];
sx q[2];
rz(-0.25114775) q[2];
sx q[2];
rz(0.37065333) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6991068) q[1];
sx q[1];
rz(-0.92256303) q[1];
sx q[1];
rz(0.68251619) q[1];
rz(-pi) q[2];
rz(-2.6044106) q[3];
sx q[3];
rz(-2.3677285) q[3];
sx q[3];
rz(2.8230132) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.236078) q[2];
sx q[2];
rz(-1.4730075) q[2];
sx q[2];
rz(-1.8900324) q[2];
rz(2.6873612) q[3];
sx q[3];
rz(-2.7082704) q[3];
sx q[3];
rz(2.7864578) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-0.77804756) q[0];
sx q[0];
rz(-0.30968928) q[0];
sx q[0];
rz(1.3409412) q[0];
rz(-0.60607934) q[1];
sx q[1];
rz(-0.93901912) q[1];
sx q[1];
rz(-1.1846503) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.143648) q[0];
sx q[0];
rz(-2.8638683) q[0];
sx q[0];
rz(0.68302897) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0834288) q[2];
sx q[2];
rz(-1.9789654) q[2];
sx q[2];
rz(-2.4525016) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.9222316) q[1];
sx q[1];
rz(-2.1582099) q[1];
sx q[1];
rz(0.14031336) q[1];
rz(1.2042449) q[3];
sx q[3];
rz(-1.1712211) q[3];
sx q[3];
rz(1.9384345) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.9081395) q[2];
sx q[2];
rz(-0.88580695) q[2];
sx q[2];
rz(-0.63956368) q[2];
rz(0.8574287) q[3];
sx q[3];
rz(-1.144751) q[3];
sx q[3];
rz(2.6517984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63067591) q[0];
sx q[0];
rz(-2.5998721) q[0];
sx q[0];
rz(-2.6211991) q[0];
rz(0.10294542) q[1];
sx q[1];
rz(-1.0504477) q[1];
sx q[1];
rz(-2.4127814) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2136912) q[0];
sx q[0];
rz(-1.7753082) q[0];
sx q[0];
rz(-1.6708899) q[0];
rz(-pi) q[1];
rz(-0.16291933) q[2];
sx q[2];
rz(-0.85328057) q[2];
sx q[2];
rz(-1.7999072) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.5506958) q[1];
sx q[1];
rz(-1.3662852) q[1];
sx q[1];
rz(1.2414316) q[1];
rz(-2.1843188) q[3];
sx q[3];
rz(-1.1144708) q[3];
sx q[3];
rz(-2.7077655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.0985428) q[2];
sx q[2];
rz(-2.7816732) q[2];
sx q[2];
rz(-2.4318802) q[2];
rz(-2.5085311) q[3];
sx q[3];
rz(-1.9934306) q[3];
sx q[3];
rz(3.0098797) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1028041) q[0];
sx q[0];
rz(-0.10237256) q[0];
sx q[0];
rz(-0.88371712) q[0];
rz(-0.9206413) q[1];
sx q[1];
rz(-1.0620774) q[1];
sx q[1];
rz(-0.19764915) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5716498) q[0];
sx q[0];
rz(-2.2075704) q[0];
sx q[0];
rz(-2.3417579) q[0];
x q[1];
rz(-3.0350787) q[2];
sx q[2];
rz(-1.4918367) q[2];
sx q[2];
rz(-2.4626412) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.27053988) q[1];
sx q[1];
rz(-2.6846243) q[1];
sx q[1];
rz(0.096710042) q[1];
rz(-pi) q[2];
rz(0.71182735) q[3];
sx q[3];
rz(-0.89755745) q[3];
sx q[3];
rz(2.2664546) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.027823042) q[2];
sx q[2];
rz(-1.7980857) q[2];
sx q[2];
rz(-3.0659884) q[2];
rz(-1.6522225) q[3];
sx q[3];
rz(-0.39947709) q[3];
sx q[3];
rz(-2.3125125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7224715) q[0];
sx q[0];
rz(-2.2561181) q[0];
sx q[0];
rz(-3.0995195) q[0];
rz(1.178859) q[1];
sx q[1];
rz(-1.1584287) q[1];
sx q[1];
rz(1.1694318) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7021381) q[0];
sx q[0];
rz(-1.5041782) q[0];
sx q[0];
rz(1.735795) q[0];
rz(0.39016907) q[2];
sx q[2];
rz(-1.8346268) q[2];
sx q[2];
rz(-0.90261501) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.52160727) q[1];
sx q[1];
rz(-1.0343401) q[1];
sx q[1];
rz(-1.9636088) q[1];
rz(1.5285989) q[3];
sx q[3];
rz(-1.9584667) q[3];
sx q[3];
rz(1.2927511) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.823267) q[2];
sx q[2];
rz(-1.1151168) q[2];
sx q[2];
rz(0.014483359) q[2];
rz(1.5197586) q[3];
sx q[3];
rz(-1.2967916) q[3];
sx q[3];
rz(-0.55317944) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13977215) q[0];
sx q[0];
rz(-2.791239) q[0];
sx q[0];
rz(-0.56234223) q[0];
rz(1.7100122) q[1];
sx q[1];
rz(-1.7446012) q[1];
sx q[1];
rz(2.6838578) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6477752) q[0];
sx q[0];
rz(-1.4871948) q[0];
sx q[0];
rz(-1.093822) q[0];
rz(-pi) q[1];
rz(-0.87403239) q[2];
sx q[2];
rz(-2.0119785) q[2];
sx q[2];
rz(2.2639084) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.22842912) q[1];
sx q[1];
rz(-2.404681) q[1];
sx q[1];
rz(-1.5494898) q[1];
x q[2];
rz(2.0972341) q[3];
sx q[3];
rz(-1.8705771) q[3];
sx q[3];
rz(1.6696904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.83595014) q[2];
sx q[2];
rz(-0.93985158) q[2];
sx q[2];
rz(-2.8590554) q[2];
rz(2.1841168) q[3];
sx q[3];
rz(-1.7493533) q[3];
sx q[3];
rz(2.6628475) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95865059) q[0];
sx q[0];
rz(-2.4562953) q[0];
sx q[0];
rz(-1.7392993) q[0];
rz(0.69933403) q[1];
sx q[1];
rz(-1.8383086) q[1];
sx q[1];
rz(1.8797849) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14784797) q[0];
sx q[0];
rz(-1.1435259) q[0];
sx q[0];
rz(-1.2742313) q[0];
x q[1];
rz(0.65323921) q[2];
sx q[2];
rz(-1.7368894) q[2];
sx q[2];
rz(1.7476029) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.1101493) q[1];
sx q[1];
rz(-1.823339) q[1];
sx q[1];
rz(-0.7048216) q[1];
rz(-pi) q[2];
rz(-0.52957876) q[3];
sx q[3];
rz(-1.7486608) q[3];
sx q[3];
rz(-1.9999268) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.70301473) q[2];
sx q[2];
rz(-1.5038749) q[2];
sx q[2];
rz(-1.6837439) q[2];
rz(2.4750989) q[3];
sx q[3];
rz(-1.7137182) q[3];
sx q[3];
rz(2.3898747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89501971) q[0];
sx q[0];
rz(-2.331215) q[0];
sx q[0];
rz(-1.1599468) q[0];
rz(0.054140422) q[1];
sx q[1];
rz(-1.663798) q[1];
sx q[1];
rz(1.0704401) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4186231) q[0];
sx q[0];
rz(-1.7374991) q[0];
sx q[0];
rz(0.82794257) q[0];
rz(-pi) q[1];
x q[1];
rz(1.530933) q[2];
sx q[2];
rz(-2.3764052) q[2];
sx q[2];
rz(2.6305692) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.4813303) q[1];
sx q[1];
rz(-0.80033014) q[1];
sx q[1];
rz(-2.7923613) q[1];
rz(-0.6108547) q[3];
sx q[3];
rz(-2.8419552) q[3];
sx q[3];
rz(-0.50869298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.76876172) q[2];
sx q[2];
rz(-2.5925345) q[2];
sx q[2];
rz(1.5268415) q[2];
rz(0.57957831) q[3];
sx q[3];
rz(-1.3468346) q[3];
sx q[3];
rz(-0.65210623) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6282745) q[0];
sx q[0];
rz(-1.5341298) q[0];
sx q[0];
rz(-2.4232724) q[0];
rz(2.1075481) q[1];
sx q[1];
rz(-0.68987344) q[1];
sx q[1];
rz(-0.72763163) q[1];
rz(3.1038531) q[2];
sx q[2];
rz(-2.3024617) q[2];
sx q[2];
rz(3.0842177) q[2];
rz(-0.53229971) q[3];
sx q[3];
rz(-1.1333864) q[3];
sx q[3];
rz(-0.65896853) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
