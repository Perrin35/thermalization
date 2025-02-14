OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.76685846) q[0];
sx q[0];
rz(-0.91349608) q[0];
sx q[0];
rz(1.7537533) q[0];
rz(-2.9517382) q[1];
sx q[1];
rz(-1.318765) q[1];
sx q[1];
rz(2.8554822) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.090534276) q[0];
sx q[0];
rz(-1.5947475) q[0];
sx q[0];
rz(-0.46047237) q[0];
rz(-0.95979664) q[2];
sx q[2];
rz(-1.5033443) q[2];
sx q[2];
rz(-2.5416201) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.82569796) q[1];
sx q[1];
rz(-1.8301536) q[1];
sx q[1];
rz(0.56866912) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0864618) q[3];
sx q[3];
rz(-1.1305439) q[3];
sx q[3];
rz(-0.80894473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-3.0949969) q[2];
sx q[2];
rz(-1.4904212) q[2];
sx q[2];
rz(-2.4563834) q[2];
rz(1.4677216) q[3];
sx q[3];
rz(-2.0592212) q[3];
sx q[3];
rz(0.26836029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0005223) q[0];
sx q[0];
rz(-1.5474316) q[0];
sx q[0];
rz(-2.1970314) q[0];
rz(-0.073307723) q[1];
sx q[1];
rz(-0.83275515) q[1];
sx q[1];
rz(-1.8836969) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2680141) q[0];
sx q[0];
rz(-1.1167545) q[0];
sx q[0];
rz(0.94768967) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7750568) q[2];
sx q[2];
rz(-2.0957207) q[2];
sx q[2];
rz(-1.7634942) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.98698178) q[1];
sx q[1];
rz(-1.715259) q[1];
sx q[1];
rz(-0.39244907) q[1];
rz(-pi) q[2];
rz(3.1044534) q[3];
sx q[3];
rz(-1.5243153) q[3];
sx q[3];
rz(-0.4872191) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8153317) q[2];
sx q[2];
rz(-1.3911284) q[2];
sx q[2];
rz(-1.0247256) q[2];
rz(-2.4552086) q[3];
sx q[3];
rz(-0.73203433) q[3];
sx q[3];
rz(-2.2288403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3150385) q[0];
sx q[0];
rz(-1.2737561) q[0];
sx q[0];
rz(0.81845534) q[0];
rz(1.9326899) q[1];
sx q[1];
rz(-1.5070288) q[1];
sx q[1];
rz(-2.2817629) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6085109) q[0];
sx q[0];
rz(-1.2405335) q[0];
sx q[0];
rz(-2.1950366) q[0];
rz(2.6099714) q[2];
sx q[2];
rz(-2.3477656) q[2];
sx q[2];
rz(2.9693111) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.3834741) q[1];
sx q[1];
rz(-0.84953298) q[1];
sx q[1];
rz(0.63754941) q[1];
rz(3.0257175) q[3];
sx q[3];
rz(-1.795035) q[3];
sx q[3];
rz(-1.7622926) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.39765057) q[2];
sx q[2];
rz(-2.8068779) q[2];
sx q[2];
rz(-0.115455) q[2];
rz(-2.7502821) q[3];
sx q[3];
rz(-2.1152928) q[3];
sx q[3];
rz(-1.1439884) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7855969) q[0];
sx q[0];
rz(-1.1859256) q[0];
sx q[0];
rz(-2.9010229) q[0];
rz(-1.7754414) q[1];
sx q[1];
rz(-1.6484478) q[1];
sx q[1];
rz(1.2496525) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8927596) q[0];
sx q[0];
rz(-0.9380149) q[0];
sx q[0];
rz(-2.8248293) q[0];
x q[1];
rz(1.0860822) q[2];
sx q[2];
rz(-2.4668985) q[2];
sx q[2];
rz(-2.7759068) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.0084997) q[1];
sx q[1];
rz(-1.805882) q[1];
sx q[1];
rz(-2.4753768) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0765188) q[3];
sx q[3];
rz(-2.1902962) q[3];
sx q[3];
rz(2.2175585) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.25632855) q[2];
sx q[2];
rz(-0.54577959) q[2];
sx q[2];
rz(1.3679999) q[2];
rz(0.3041501) q[3];
sx q[3];
rz(-1.5658028) q[3];
sx q[3];
rz(-0.69653851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91365415) q[0];
sx q[0];
rz(-2.9253687) q[0];
sx q[0];
rz(-0.53713334) q[0];
rz(-0.74983239) q[1];
sx q[1];
rz(-2.0593819) q[1];
sx q[1];
rz(-0.73712635) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1196302) q[0];
sx q[0];
rz(-0.11521327) q[0];
sx q[0];
rz(-2.553279) q[0];
rz(1.7456878) q[2];
sx q[2];
rz(-1.6564257) q[2];
sx q[2];
rz(2.0599287) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.3143718) q[1];
sx q[1];
rz(-1.4431568) q[1];
sx q[1];
rz(-2.8807544) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6629822) q[3];
sx q[3];
rz(-1.8927729) q[3];
sx q[3];
rz(-2.3809402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.6970814) q[2];
sx q[2];
rz(-2.7883174) q[2];
sx q[2];
rz(-0.14661655) q[2];
rz(-1.692449) q[3];
sx q[3];
rz(-1.861898) q[3];
sx q[3];
rz(-0.52391887) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61831063) q[0];
sx q[0];
rz(-1.016541) q[0];
sx q[0];
rz(-1.7200394) q[0];
rz(1.8824185) q[1];
sx q[1];
rz(-2.4651395) q[1];
sx q[1];
rz(2.7377985) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97966563) q[0];
sx q[0];
rz(-2.296519) q[0];
sx q[0];
rz(2.7620074) q[0];
x q[1];
rz(-1.7558891) q[2];
sx q[2];
rz(-1.3178692) q[2];
sx q[2];
rz(1.3651939) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.4686069) q[1];
sx q[1];
rz(-0.98277175) q[1];
sx q[1];
rz(2.468416) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.29098068) q[3];
sx q[3];
rz(-0.70937362) q[3];
sx q[3];
rz(0.13252276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8518565) q[2];
sx q[2];
rz(-1.7278262) q[2];
sx q[2];
rz(2.9750321) q[2];
rz(-1.5038495) q[3];
sx q[3];
rz(-0.47347355) q[3];
sx q[3];
rz(3.04305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.903776) q[0];
sx q[0];
rz(-0.38919583) q[0];
sx q[0];
rz(0.87752262) q[0];
rz(-0.96218836) q[1];
sx q[1];
rz(-1.1233556) q[1];
sx q[1];
rz(-2.7872564) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0577045) q[0];
sx q[0];
rz(-1.5885) q[0];
sx q[0];
rz(2.903865) q[0];
rz(-1.7298752) q[2];
sx q[2];
rz(-2.9949778) q[2];
sx q[2];
rz(-0.40381144) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.52088941) q[1];
sx q[1];
rz(-2.1701309) q[1];
sx q[1];
rz(-2.4769737) q[1];
x q[2];
rz(0.50994344) q[3];
sx q[3];
rz(-1.5048319) q[3];
sx q[3];
rz(-0.13291453) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3366036) q[2];
sx q[2];
rz(-2.4038834) q[2];
sx q[2];
rz(-2.0965516) q[2];
rz(0.40237829) q[3];
sx q[3];
rz(-1.1674403) q[3];
sx q[3];
rz(-1.4761402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7104915) q[0];
sx q[0];
rz(-0.11678188) q[0];
sx q[0];
rz(0.85103273) q[0];
rz(-0.35375133) q[1];
sx q[1];
rz(-1.7314792) q[1];
sx q[1];
rz(-2.2869349) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2609278) q[0];
sx q[0];
rz(-1.3890084) q[0];
sx q[0];
rz(-0.28566912) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.607367) q[2];
sx q[2];
rz(-2.8534128) q[2];
sx q[2];
rz(2.2792918) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.15957279) q[1];
sx q[1];
rz(-2.5347485) q[1];
sx q[1];
rz(2.8446376) q[1];
rz(-1.6519671) q[3];
sx q[3];
rz(-1.1690835) q[3];
sx q[3];
rz(-2.4430526) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9510368) q[2];
sx q[2];
rz(-2.4943116) q[2];
sx q[2];
rz(0.62038842) q[2];
rz(0.22647151) q[3];
sx q[3];
rz(-1.3969996) q[3];
sx q[3];
rz(-2.5605719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-1.9613551) q[0];
sx q[0];
rz(-1.9397475) q[0];
sx q[0];
rz(3.1264547) q[0];
rz(2.119428) q[1];
sx q[1];
rz(-0.41003761) q[1];
sx q[1];
rz(1.2000363) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2101636) q[0];
sx q[0];
rz(-0.72648049) q[0];
sx q[0];
rz(-1.8379712) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6300417) q[2];
sx q[2];
rz(-0.97840913) q[2];
sx q[2];
rz(0.5991583) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.85644279) q[1];
sx q[1];
rz(-2.0822444) q[1];
sx q[1];
rz(0.034354547) q[1];
rz(-pi) q[2];
rz(-1.0608114) q[3];
sx q[3];
rz(-2.0037162) q[3];
sx q[3];
rz(1.0266538) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6156561) q[2];
sx q[2];
rz(-2.2286712) q[2];
sx q[2];
rz(-2.9021662) q[2];
rz(0.058852363) q[3];
sx q[3];
rz(-2.3299496) q[3];
sx q[3];
rz(-0.27590251) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1926159) q[0];
sx q[0];
rz(-1.2139576) q[0];
sx q[0];
rz(2.1666727) q[0];
rz(2.4661567) q[1];
sx q[1];
rz(-2.2223739) q[1];
sx q[1];
rz(-0.98178896) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47768053) q[0];
sx q[0];
rz(-2.4943611) q[0];
sx q[0];
rz(-1.1301229) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7300225) q[2];
sx q[2];
rz(-1.6333505) q[2];
sx q[2];
rz(2.6175509) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.2416545) q[1];
sx q[1];
rz(-1.3119427) q[1];
sx q[1];
rz(-2.3090825) q[1];
x q[2];
rz(1.124106) q[3];
sx q[3];
rz(-1.6047995) q[3];
sx q[3];
rz(1.407935) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.56741095) q[2];
sx q[2];
rz(-2.1705706) q[2];
sx q[2];
rz(-1.2249472) q[2];
rz(2.775906) q[3];
sx q[3];
rz(-2.1920125) q[3];
sx q[3];
rz(-2.001781) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0443307) q[0];
sx q[0];
rz(-1.8987569) q[0];
sx q[0];
rz(1.4415997) q[0];
rz(3.0615831) q[1];
sx q[1];
rz(-1.7623822) q[1];
sx q[1];
rz(-0.0026127876) q[1];
rz(-1.710113) q[2];
sx q[2];
rz(-2.0124544) q[2];
sx q[2];
rz(-2.6228946) q[2];
rz(-3.0946685) q[3];
sx q[3];
rz(-0.45294807) q[3];
sx q[3];
rz(3.1035085) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
