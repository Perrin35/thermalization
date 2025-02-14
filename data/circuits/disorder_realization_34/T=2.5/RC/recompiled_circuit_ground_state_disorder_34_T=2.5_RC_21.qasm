OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.3736149) q[0];
sx q[0];
rz(-0.64426214) q[0];
sx q[0];
rz(-0.86384073) q[0];
rz(1.5732425) q[1];
sx q[1];
rz(-1.9116126) q[1];
sx q[1];
rz(0.65043515) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38017958) q[0];
sx q[0];
rz(-2.5963547) q[0];
sx q[0];
rz(-2.7392967) q[0];
rz(2.8312248) q[2];
sx q[2];
rz(-1.5987607) q[2];
sx q[2];
rz(-1.8486277) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.1623154) q[1];
sx q[1];
rz(-2.3188496) q[1];
sx q[1];
rz(2.2869799) q[1];
x q[2];
rz(-0.69919805) q[3];
sx q[3];
rz(-0.16575925) q[3];
sx q[3];
rz(2.0729271) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.9170561) q[2];
sx q[2];
rz(-1.7660331) q[2];
sx q[2];
rz(-1.7898233) q[2];
rz(-2.3228862) q[3];
sx q[3];
rz(-1.4592146) q[3];
sx q[3];
rz(1.5514577) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3025892) q[0];
sx q[0];
rz(-0.67794696) q[0];
sx q[0];
rz(2.5658521) q[0];
rz(0.88306824) q[1];
sx q[1];
rz(-1.6022316) q[1];
sx q[1];
rz(-1.8227122) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7144189) q[0];
sx q[0];
rz(-1.321081) q[0];
sx q[0];
rz(0.7575794) q[0];
rz(-pi) q[1];
rz(2.3568866) q[2];
sx q[2];
rz(-1.5532576) q[2];
sx q[2];
rz(-2.6477046) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-3.0134351) q[1];
sx q[1];
rz(-0.87608209) q[1];
sx q[1];
rz(2.6753268) q[1];
x q[2];
rz(-2.8208597) q[3];
sx q[3];
rz(-1.7400636) q[3];
sx q[3];
rz(1.2812268) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.64112249) q[2];
sx q[2];
rz(-1.0176071) q[2];
sx q[2];
rz(2.9024331) q[2];
rz(-0.33254361) q[3];
sx q[3];
rz(-1.6859237) q[3];
sx q[3];
rz(1.0913764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34514937) q[0];
sx q[0];
rz(-2.394016) q[0];
sx q[0];
rz(2.8969966) q[0];
rz(0.39464513) q[1];
sx q[1];
rz(-2.5375745) q[1];
sx q[1];
rz(1.9371202) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24893471) q[0];
sx q[0];
rz(-1.5818074) q[0];
sx q[0];
rz(1.5599773) q[0];
rz(-pi) q[1];
rz(2.0599686) q[2];
sx q[2];
rz(-1.2252639) q[2];
sx q[2];
rz(-0.45128497) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.2358346) q[1];
sx q[1];
rz(-1.8595962) q[1];
sx q[1];
rz(-2.1101084) q[1];
rz(-pi) q[2];
rz(2.0709023) q[3];
sx q[3];
rz(-1.9311889) q[3];
sx q[3];
rz(-2.5640324) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.7977153) q[2];
sx q[2];
rz(-0.98261967) q[2];
sx q[2];
rz(-2.9003411) q[2];
rz(-1.7734211) q[3];
sx q[3];
rz(-1.3894812) q[3];
sx q[3];
rz(2.1696137) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5695246) q[0];
sx q[0];
rz(-1.854874) q[0];
sx q[0];
rz(-2.0401814) q[0];
rz(-1.4934348) q[1];
sx q[1];
rz(-0.28683528) q[1];
sx q[1];
rz(2.1410087) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34383044) q[0];
sx q[0];
rz(-1.6407818) q[0];
sx q[0];
rz(1.153341) q[0];
x q[1];
rz(2.1270833) q[2];
sx q[2];
rz(-1.5624701) q[2];
sx q[2];
rz(2.3187092) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.1305589) q[1];
sx q[1];
rz(-0.10422626) q[1];
sx q[1];
rz(0.24918814) q[1];
rz(2.9977418) q[3];
sx q[3];
rz(-1.1531623) q[3];
sx q[3];
rz(-0.9506433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.7857886) q[2];
sx q[2];
rz(-1.641909) q[2];
sx q[2];
rz(0.22264063) q[2];
rz(1.6629793) q[3];
sx q[3];
rz(-0.40545774) q[3];
sx q[3];
rz(1.2610029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-1.3123689) q[0];
sx q[0];
rz(-1.1382505) q[0];
sx q[0];
rz(-2.9421222) q[0];
rz(-1.624674) q[1];
sx q[1];
rz(-1.3360887) q[1];
sx q[1];
rz(2.6768501) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7620649) q[0];
sx q[0];
rz(-1.7025954) q[0];
sx q[0];
rz(1.2799311) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5114485) q[2];
sx q[2];
rz(-1.7918824) q[2];
sx q[2];
rz(-0.57002744) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.66942635) q[1];
sx q[1];
rz(-1.8490044) q[1];
sx q[1];
rz(2.8406155) q[1];
rz(-pi) q[2];
rz(3.1384077) q[3];
sx q[3];
rz(-0.69744686) q[3];
sx q[3];
rz(1.3435329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(3.0988079) q[2];
sx q[2];
rz(-1.3036737) q[2];
sx q[2];
rz(0.25351563) q[2];
rz(0.86067307) q[3];
sx q[3];
rz(-2.6200675) q[3];
sx q[3];
rz(-1.1389987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-3.0694224) q[0];
sx q[0];
rz(-1.5266029) q[0];
sx q[0];
rz(1.7813064) q[0];
rz(2.5379429) q[1];
sx q[1];
rz(-2.344548) q[1];
sx q[1];
rz(-0.48702494) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1384772) q[0];
sx q[0];
rz(-2.2238408) q[0];
sx q[0];
rz(1.6621291) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6658044) q[2];
sx q[2];
rz(-2.3072647) q[2];
sx q[2];
rz(-0.67259865) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0700775) q[1];
sx q[1];
rz(-1.9299475) q[1];
sx q[1];
rz(2.9643494) q[1];
rz(-0.2305687) q[3];
sx q[3];
rz(-1.9631223) q[3];
sx q[3];
rz(2.2986064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.6765678) q[2];
sx q[2];
rz(-1.0176696) q[2];
sx q[2];
rz(-2.6505995) q[2];
rz(-2.6077014) q[3];
sx q[3];
rz(-0.54755727) q[3];
sx q[3];
rz(3.0965613) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8462867) q[0];
sx q[0];
rz(-2.5745109) q[0];
sx q[0];
rz(1.2714161) q[0];
rz(-0.92453399) q[1];
sx q[1];
rz(-1.144616) q[1];
sx q[1];
rz(1.0251934) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36350016) q[0];
sx q[0];
rz(-2.1937167) q[0];
sx q[0];
rz(-1.8353375) q[0];
x q[1];
rz(1.2751638) q[2];
sx q[2];
rz(-2.4951093) q[2];
sx q[2];
rz(1.531284) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.89910903) q[1];
sx q[1];
rz(-2.5216853) q[1];
sx q[1];
rz(1.0044881) q[1];
rz(-2.4348172) q[3];
sx q[3];
rz(-1.7590344) q[3];
sx q[3];
rz(0.95043236) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.2123432) q[2];
sx q[2];
rz(-2.847147) q[2];
sx q[2];
rz(2.219131) q[2];
rz(-2.755002) q[3];
sx q[3];
rz(-1.8542733) q[3];
sx q[3];
rz(-1.6283901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59231049) q[0];
sx q[0];
rz(-3.1070502) q[0];
sx q[0];
rz(0.70796815) q[0];
rz(-2.6225846) q[1];
sx q[1];
rz(-0.52878562) q[1];
sx q[1];
rz(-2.4755898) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3579263) q[0];
sx q[0];
rz(-2.5979498) q[0];
sx q[0];
rz(-1.3915791) q[0];
rz(-pi) q[1];
rz(3.0254499) q[2];
sx q[2];
rz(-1.2525038) q[2];
sx q[2];
rz(1.4333985) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4909926) q[1];
sx q[1];
rz(-1.566399) q[1];
sx q[1];
rz(-2.472509) q[1];
rz(-pi) q[2];
rz(1.610393) q[3];
sx q[3];
rz(-0.61568135) q[3];
sx q[3];
rz(0.99056584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.4697326) q[2];
sx q[2];
rz(-1.1598776) q[2];
sx q[2];
rz(1.0603909) q[2];
rz(-1.4276069) q[3];
sx q[3];
rz(-1.5771882) q[3];
sx q[3];
rz(2.3744627) q[3];
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
rz(-1.0107182) q[0];
sx q[0];
rz(-2.3227782) q[0];
sx q[0];
rz(2.4406216) q[0];
rz(-2.2531807) q[1];
sx q[1];
rz(-2.7263434) q[1];
sx q[1];
rz(-1.7078687) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0232359) q[0];
sx q[0];
rz(-2.6161086) q[0];
sx q[0];
rz(-2.9246287) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5979684) q[2];
sx q[2];
rz(-1.9765696) q[2];
sx q[2];
rz(-2.9494065) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-3.0317185) q[1];
sx q[1];
rz(-1.5496792) q[1];
sx q[1];
rz(-1.5892522) q[1];
rz(-2.4495226) q[3];
sx q[3];
rz(-1.0502397) q[3];
sx q[3];
rz(2.6259921) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.1514757) q[2];
sx q[2];
rz(-1.2274123) q[2];
sx q[2];
rz(-1.6647476) q[2];
rz(0.19674033) q[3];
sx q[3];
rz(-1.1367831) q[3];
sx q[3];
rz(-0.93973947) q[3];
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
rz(1.5340586) q[0];
sx q[0];
rz(-1.5787831) q[0];
sx q[0];
rz(-1.1210572) q[0];
rz(2.6634482) q[1];
sx q[1];
rz(-2.4604764) q[1];
sx q[1];
rz(-0.71472439) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5276565) q[0];
sx q[0];
rz(-0.52866565) q[0];
sx q[0];
rz(0.40724004) q[0];
rz(1.0214731) q[2];
sx q[2];
rz(-1.2397814) q[2];
sx q[2];
rz(-2.2084055) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.0352382) q[1];
sx q[1];
rz(-1.0627305) q[1];
sx q[1];
rz(0.50915995) q[1];
rz(0.86831324) q[3];
sx q[3];
rz(-0.83783093) q[3];
sx q[3];
rz(2.365685) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.62632442) q[2];
sx q[2];
rz(-1.3268665) q[2];
sx q[2];
rz(1.932762) q[2];
rz(-1.442046) q[3];
sx q[3];
rz(-2.2550826) q[3];
sx q[3];
rz(0.065571688) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6616853) q[0];
sx q[0];
rz(-1.8330782) q[0];
sx q[0];
rz(-0.080396419) q[0];
rz(-2.5936364) q[1];
sx q[1];
rz(-2.4520271) q[1];
sx q[1];
rz(1.3507623) q[1];
rz(1.1153658) q[2];
sx q[2];
rz(-0.51389384) q[2];
sx q[2];
rz(0.92432164) q[2];
rz(-0.73765909) q[3];
sx q[3];
rz(-2.548703) q[3];
sx q[3];
rz(-0.92787837) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
