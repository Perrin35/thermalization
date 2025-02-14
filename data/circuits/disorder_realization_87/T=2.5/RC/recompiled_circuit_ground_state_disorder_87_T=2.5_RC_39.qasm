OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.7423695) q[0];
sx q[0];
rz(2.7288781) q[0];
sx q[0];
rz(5.4469845) q[0];
rz(0.77464453) q[1];
sx q[1];
rz(-0.062008468) q[1];
sx q[1];
rz(-1.0083415) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22577408) q[0];
sx q[0];
rz(-1.649722) q[0];
sx q[0];
rz(1.945334) q[0];
x q[1];
rz(-2.600619) q[2];
sx q[2];
rz(-0.43593513) q[2];
sx q[2];
rz(1.3672787) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-3.0550067) q[1];
sx q[1];
rz(-1.0957967) q[1];
sx q[1];
rz(1.7473298) q[1];
x q[2];
rz(-1.7640231) q[3];
sx q[3];
rz(-1.591103) q[3];
sx q[3];
rz(1.8373674) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.3081554) q[2];
sx q[2];
rz(-2.1980632) q[2];
sx q[2];
rz(0.65307871) q[2];
rz(-0.11450442) q[3];
sx q[3];
rz(-1.7479892) q[3];
sx q[3];
rz(2.0792927) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
rz(-2.7750875) q[0];
sx q[0];
rz(-2.1108284) q[0];
sx q[0];
rz(-1.7979701) q[0];
rz(1.1603181) q[1];
sx q[1];
rz(-0.41028816) q[1];
sx q[1];
rz(0.62058273) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68408332) q[0];
sx q[0];
rz(-0.39827049) q[0];
sx q[0];
rz(-1.0578367) q[0];
x q[1];
rz(-1.1158022) q[2];
sx q[2];
rz(-2.4206851) q[2];
sx q[2];
rz(-1.4269331) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.81765905) q[1];
sx q[1];
rz(-2.226737) q[1];
sx q[1];
rz(-2.0327912) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.28264265) q[3];
sx q[3];
rz(-0.97393721) q[3];
sx q[3];
rz(2.7925036) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.6985942) q[2];
sx q[2];
rz(-1.6220762) q[2];
sx q[2];
rz(0.2300187) q[2];
rz(-1.5551785) q[3];
sx q[3];
rz(-1.9929726) q[3];
sx q[3];
rz(2.8659099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8860633) q[0];
sx q[0];
rz(-0.38917381) q[0];
sx q[0];
rz(1.080876) q[0];
rz(0.64322645) q[1];
sx q[1];
rz(-0.79935646) q[1];
sx q[1];
rz(3.0562775) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9935308) q[0];
sx q[0];
rz(-0.9328273) q[0];
sx q[0];
rz(0.83420269) q[0];
rz(-pi) q[1];
rz(2.5695989) q[2];
sx q[2];
rz(-1.3912918) q[2];
sx q[2];
rz(-3.1310022) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.7300295) q[1];
sx q[1];
rz(-1.4394898) q[1];
sx q[1];
rz(2.9741785) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5253228) q[3];
sx q[3];
rz(-2.1915428) q[3];
sx q[3];
rz(0.57102247) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.2049415) q[2];
sx q[2];
rz(-3.1320269) q[2];
sx q[2];
rz(-0.19634518) q[2];
rz(0.28228545) q[3];
sx q[3];
rz(-1.117027) q[3];
sx q[3];
rz(1.0452247) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1333756) q[0];
sx q[0];
rz(-2.5947925) q[0];
sx q[0];
rz(2.7561482) q[0];
rz(1.4133981) q[1];
sx q[1];
rz(-1.9368659) q[1];
sx q[1];
rz(-0.24994303) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88943931) q[0];
sx q[0];
rz(-1.1595386) q[0];
sx q[0];
rz(0.90775872) q[0];
rz(-1.6310299) q[2];
sx q[2];
rz(-1.6585729) q[2];
sx q[2];
rz(-0.46030948) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.6648207) q[1];
sx q[1];
rz(-2.6198434) q[1];
sx q[1];
rz(2.797965) q[1];
rz(2.945845) q[3];
sx q[3];
rz(-2.3382814) q[3];
sx q[3];
rz(-2.3565528) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.46127737) q[2];
sx q[2];
rz(-1.2948493) q[2];
sx q[2];
rz(2.9052022) q[2];
rz(-1.8810898) q[3];
sx q[3];
rz(-0.25006306) q[3];
sx q[3];
rz(2.4642956) q[3];
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
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33609718) q[0];
sx q[0];
rz(-1.6039055) q[0];
sx q[0];
rz(-2.6564823) q[0];
rz(-1.1803892) q[1];
sx q[1];
rz(-1.5472658) q[1];
sx q[1];
rz(-0.83650437) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77959767) q[0];
sx q[0];
rz(-1.9266202) q[0];
sx q[0];
rz(1.9076288) q[0];
rz(2.2274687) q[2];
sx q[2];
rz(-2.6466922) q[2];
sx q[2];
rz(0.47593853) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.066557601) q[1];
sx q[1];
rz(-1.6187985) q[1];
sx q[1];
rz(-2.2407299) q[1];
rz(-pi) q[2];
rz(-1.5077293) q[3];
sx q[3];
rz(-1.9708037) q[3];
sx q[3];
rz(-2.0162752) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.89852077) q[2];
sx q[2];
rz(-2.409755) q[2];
sx q[2];
rz(-0.76954976) q[2];
rz(-0.92929333) q[3];
sx q[3];
rz(-1.1309036) q[3];
sx q[3];
rz(-0.23269674) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9418075) q[0];
sx q[0];
rz(-0.48518825) q[0];
sx q[0];
rz(0.76882452) q[0];
rz(1.3517514) q[1];
sx q[1];
rz(-0.47305802) q[1];
sx q[1];
rz(0.88968712) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79737299) q[0];
sx q[0];
rz(-1.780172) q[0];
sx q[0];
rz(-1.1816417) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.0013498505) q[2];
sx q[2];
rz(-1.6151516) q[2];
sx q[2];
rz(2.7173017) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.5701712) q[1];
sx q[1];
rz(-1.4246877) q[1];
sx q[1];
rz(0.04695462) q[1];
rz(-2.6906312) q[3];
sx q[3];
rz(-1.4672674) q[3];
sx q[3];
rz(-0.24117392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.6512904) q[2];
sx q[2];
rz(-1.7040665) q[2];
sx q[2];
rz(-1.5865405) q[2];
rz(-1.8709315) q[3];
sx q[3];
rz(-2.7226166) q[3];
sx q[3];
rz(-3.0095625) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1726058) q[0];
sx q[0];
rz(-2.0033328) q[0];
sx q[0];
rz(1.1440811) q[0];
rz(-2.3106958) q[1];
sx q[1];
rz(-1.7832719) q[1];
sx q[1];
rz(-2.3544618) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89881508) q[0];
sx q[0];
rz(-1.4393596) q[0];
sx q[0];
rz(-1.7064894) q[0];
rz(0.93573715) q[2];
sx q[2];
rz(-2.3496186) q[2];
sx q[2];
rz(2.9284649) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.28879189) q[1];
sx q[1];
rz(-2.0388985) q[1];
sx q[1];
rz(1.0352943) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1803341) q[3];
sx q[3];
rz(-1.045908) q[3];
sx q[3];
rz(-2.386664) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1181011) q[2];
sx q[2];
rz(-2.1716437) q[2];
sx q[2];
rz(-0.033163158) q[2];
rz(-1.5983332) q[3];
sx q[3];
rz(-2.4256746) q[3];
sx q[3];
rz(1.5020812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24932662) q[0];
sx q[0];
rz(-1.5557657) q[0];
sx q[0];
rz(-1.7484885) q[0];
rz(-0.45685592) q[1];
sx q[1];
rz(-1.1445878) q[1];
sx q[1];
rz(-2.0410062) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83525676) q[0];
sx q[0];
rz(-3.0655711) q[0];
sx q[0];
rz(3.118909) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4462561) q[2];
sx q[2];
rz(-2.0954663) q[2];
sx q[2];
rz(2.6161043) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7707899) q[1];
sx q[1];
rz(-2.2375771) q[1];
sx q[1];
rz(3.0984224) q[1];
rz(-pi) q[2];
x q[2];
rz(0.3914137) q[3];
sx q[3];
rz(-1.7000704) q[3];
sx q[3];
rz(-1.5864652) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.042772375) q[2];
sx q[2];
rz(-2.6440812) q[2];
sx q[2];
rz(1.5808606) q[2];
rz(-0.76357311) q[3];
sx q[3];
rz(-1.6288501) q[3];
sx q[3];
rz(1.0360576) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40912691) q[0];
sx q[0];
rz(-2.3865073) q[0];
sx q[0];
rz(2.8501046) q[0];
rz(1.2358707) q[1];
sx q[1];
rz(-1.9449077) q[1];
sx q[1];
rz(0.12289563) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9366347) q[0];
sx q[0];
rz(-2.7911148) q[0];
sx q[0];
rz(0.15021439) q[0];
rz(-pi) q[1];
rz(2.9192186) q[2];
sx q[2];
rz(-1.8463328) q[2];
sx q[2];
rz(-2.0982519) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.025561573) q[1];
sx q[1];
rz(-1.4605547) q[1];
sx q[1];
rz(-1.7787537) q[1];
x q[2];
rz(-3.1355643) q[3];
sx q[3];
rz(-1.7154126) q[3];
sx q[3];
rz(-0.58352375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.62873944) q[2];
sx q[2];
rz(-1.7491919) q[2];
sx q[2];
rz(1.0635618) q[2];
rz(-0.17744803) q[3];
sx q[3];
rz(-2.1833503) q[3];
sx q[3];
rz(1.0183081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7040831) q[0];
sx q[0];
rz(-2.6884485) q[0];
sx q[0];
rz(-1.7105239) q[0];
rz(0.60072947) q[1];
sx q[1];
rz(-0.44980106) q[1];
sx q[1];
rz(-2.5240555) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3356614) q[0];
sx q[0];
rz(-1.7863723) q[0];
sx q[0];
rz(-0.12713253) q[0];
rz(-pi) q[1];
rz(-2.6170116) q[2];
sx q[2];
rz(-0.89898212) q[2];
sx q[2];
rz(-1.2685405) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.85166288) q[1];
sx q[1];
rz(-1.549822) q[1];
sx q[1];
rz(-2.1037654) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.57281877) q[3];
sx q[3];
rz(-1.8316934) q[3];
sx q[3];
rz(3.0134921) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.90126976) q[2];
sx q[2];
rz(-1.0003041) q[2];
sx q[2];
rz(2.1642302) q[2];
rz(-1.0824925) q[3];
sx q[3];
rz(-0.87477028) q[3];
sx q[3];
rz(-3.0860743) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33901535) q[0];
sx q[0];
rz(-0.76887283) q[0];
sx q[0];
rz(-0.73706891) q[0];
rz(0.045724178) q[1];
sx q[1];
rz(-0.8174236) q[1];
sx q[1];
rz(1.5541706) q[1];
rz(2.634766) q[2];
sx q[2];
rz(-1.097737) q[2];
sx q[2];
rz(-3.1180827) q[2];
rz(-1.7405199) q[3];
sx q[3];
rz(-2.2206497) q[3];
sx q[3];
rz(2.7322265) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
