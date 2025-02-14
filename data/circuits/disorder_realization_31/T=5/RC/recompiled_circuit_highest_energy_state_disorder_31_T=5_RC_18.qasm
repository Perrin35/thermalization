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
rz(-1.8588869) q[0];
sx q[0];
rz(-1.016322) q[0];
sx q[0];
rz(-1.8860201) q[0];
rz(-1.3610871) q[1];
sx q[1];
rz(-0.95870107) q[1];
sx q[1];
rz(-2.5348742) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62846216) q[0];
sx q[0];
rz(-2.2431122) q[0];
sx q[0];
rz(0.99654128) q[0];
rz(-pi) q[1];
rz(3.1090281) q[2];
sx q[2];
rz(-1.3059907) q[2];
sx q[2];
rz(0.84003583) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.1944988) q[1];
sx q[1];
rz(-1.7003807) q[1];
sx q[1];
rz(-1.938185) q[1];
rz(-pi) q[2];
x q[2];
rz(0.16096551) q[3];
sx q[3];
rz(-1.6897298) q[3];
sx q[3];
rz(1.87784) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.9386193) q[2];
sx q[2];
rz(-1.1569269) q[2];
sx q[2];
rz(0.17091664) q[2];
rz(-1.1275229) q[3];
sx q[3];
rz(-0.5564965) q[3];
sx q[3];
rz(-3.0881622) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4776066) q[0];
sx q[0];
rz(-0.32653102) q[0];
sx q[0];
rz(-2.4531181) q[0];
rz(1.7636048) q[1];
sx q[1];
rz(-1.0208027) q[1];
sx q[1];
rz(-3.0000906) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3162075) q[0];
sx q[0];
rz(-1.0757462) q[0];
sx q[0];
rz(-2.0284257) q[0];
rz(-pi) q[1];
rz(0.42948855) q[2];
sx q[2];
rz(-1.4802684) q[2];
sx q[2];
rz(1.2538101) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.32530871) q[1];
sx q[1];
rz(-0.46093309) q[1];
sx q[1];
rz(0.67616762) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7670632) q[3];
sx q[3];
rz(-0.39159039) q[3];
sx q[3];
rz(2.151769) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.7576393) q[2];
sx q[2];
rz(-0.44718224) q[2];
sx q[2];
rz(-3.1128856) q[2];
rz(-2.7590397) q[3];
sx q[3];
rz(-1.5111204) q[3];
sx q[3];
rz(0.20146519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89556995) q[0];
sx q[0];
rz(-1.0423132) q[0];
sx q[0];
rz(-0.37164715) q[0];
rz(2.9314575) q[1];
sx q[1];
rz(-2.7919283) q[1];
sx q[1];
rz(1.9336611) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18603095) q[0];
sx q[0];
rz(-0.81350858) q[0];
sx q[0];
rz(0.71788089) q[0];
x q[1];
rz(1.5919331) q[2];
sx q[2];
rz(-1.653228) q[2];
sx q[2];
rz(1.7475413) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.4000798) q[1];
sx q[1];
rz(-2.1472405) q[1];
sx q[1];
rz(2.4680572) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.34545107) q[3];
sx q[3];
rz(-1.5323223) q[3];
sx q[3];
rz(2.1360122) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5189884) q[2];
sx q[2];
rz(-3.0486139) q[2];
sx q[2];
rz(-0.63817111) q[2];
rz(-2.7291164) q[3];
sx q[3];
rz(-1.4700438) q[3];
sx q[3];
rz(-2.0392058) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.07688044) q[0];
sx q[0];
rz(-0.91016722) q[0];
sx q[0];
rz(-1.3077211) q[0];
rz(-0.67963302) q[1];
sx q[1];
rz(-0.72703397) q[1];
sx q[1];
rz(-1.9576498) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9480431) q[0];
sx q[0];
rz(-1.8947161) q[0];
sx q[0];
rz(0.69635038) q[0];
rz(-pi) q[1];
rz(-1.6865191) q[2];
sx q[2];
rz(-1.8050131) q[2];
sx q[2];
rz(-0.17143347) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.7412926) q[1];
sx q[1];
rz(-1.6389936) q[1];
sx q[1];
rz(2.791021) q[1];
rz(-pi) q[2];
rz(-0.77149646) q[3];
sx q[3];
rz(-0.1813387) q[3];
sx q[3];
rz(-2.4490631) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(3.1188941) q[2];
sx q[2];
rz(-0.1476295) q[2];
sx q[2];
rz(-1.0559319) q[2];
rz(0.28327709) q[3];
sx q[3];
rz(-1.8746459) q[3];
sx q[3];
rz(-1.383708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.088575514) q[0];
sx q[0];
rz(-0.96240369) q[0];
sx q[0];
rz(-1.7330633) q[0];
rz(-2.9119496) q[1];
sx q[1];
rz(-0.85669986) q[1];
sx q[1];
rz(-1.9502669) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9339122) q[0];
sx q[0];
rz(-1.4290853) q[0];
sx q[0];
rz(-1.974154) q[0];
rz(-pi) q[1];
rz(3.0158305) q[2];
sx q[2];
rz(-1.874141) q[2];
sx q[2];
rz(1.4810918) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8247485) q[1];
sx q[1];
rz(-1.1543546) q[1];
sx q[1];
rz(-1.29711) q[1];
x q[2];
rz(0.91893371) q[3];
sx q[3];
rz(-1.100538) q[3];
sx q[3];
rz(1.7277499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4952937) q[2];
sx q[2];
rz(-2.7806492) q[2];
sx q[2];
rz(1.4781282) q[2];
rz(-0.16119257) q[3];
sx q[3];
rz(-1.6094004) q[3];
sx q[3];
rz(2.351779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5021055) q[0];
sx q[0];
rz(-0.4244856) q[0];
sx q[0];
rz(-0.91833997) q[0];
rz(0.25686747) q[1];
sx q[1];
rz(-0.99595064) q[1];
sx q[1];
rz(1.1161944) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74280069) q[0];
sx q[0];
rz(-1.0519807) q[0];
sx q[0];
rz(2.1543845) q[0];
rz(1.6565422) q[2];
sx q[2];
rz(-1.298279) q[2];
sx q[2];
rz(-2.9086824) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.5600561) q[1];
sx q[1];
rz(-1.7651084) q[1];
sx q[1];
rz(-0.24340731) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5258216) q[3];
sx q[3];
rz(-0.29524657) q[3];
sx q[3];
rz(2.8930882) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9548698) q[2];
sx q[2];
rz(-0.8064417) q[2];
sx q[2];
rz(-0.40697971) q[2];
rz(0.61378971) q[3];
sx q[3];
rz(-1.611462) q[3];
sx q[3];
rz(2.6500402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.75519049) q[0];
sx q[0];
rz(-0.84119216) q[0];
sx q[0];
rz(0.92591539) q[0];
rz(-0.47261247) q[1];
sx q[1];
rz(-2.0874529) q[1];
sx q[1];
rz(-0.47017631) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42379728) q[0];
sx q[0];
rz(-1.4873624) q[0];
sx q[0];
rz(-2.3071737) q[0];
rz(-2.0704001) q[2];
sx q[2];
rz(-0.70935574) q[2];
sx q[2];
rz(-1.5200523) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.8452085) q[1];
sx q[1];
rz(-1.607158) q[1];
sx q[1];
rz(-2.3567289) q[1];
x q[2];
rz(2.3005465) q[3];
sx q[3];
rz(-1.5435092) q[3];
sx q[3];
rz(-0.59852615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.65972313) q[2];
sx q[2];
rz(-2.0707776) q[2];
sx q[2];
rz(-2.3107963) q[2];
rz(-3.0610436) q[3];
sx q[3];
rz(-2.060067) q[3];
sx q[3];
rz(-0.95675937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.047121) q[0];
sx q[0];
rz(-1.6166649) q[0];
sx q[0];
rz(-2.9834874) q[0];
rz(-2.415601) q[1];
sx q[1];
rz(-0.43015614) q[1];
sx q[1];
rz(0.37011883) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86240101) q[0];
sx q[0];
rz(-1.6625615) q[0];
sx q[0];
rz(2.807995) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0453395) q[2];
sx q[2];
rz(-2.2703746) q[2];
sx q[2];
rz(-2.2945905) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7230795) q[1];
sx q[1];
rz(-1.1795768) q[1];
sx q[1];
rz(-1.1488951) q[1];
rz(-2.1922969) q[3];
sx q[3];
rz(-2.5007476) q[3];
sx q[3];
rz(-0.31539886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.98081723) q[2];
sx q[2];
rz(-1.2951916) q[2];
sx q[2];
rz(0.13713169) q[2];
rz(1.49617) q[3];
sx q[3];
rz(-2.4479595) q[3];
sx q[3];
rz(2.6579198) q[3];
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
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5982323) q[0];
sx q[0];
rz(-1.0022663) q[0];
sx q[0];
rz(-0.12635669) q[0];
rz(2.5670746) q[1];
sx q[1];
rz(-2.2210329) q[1];
sx q[1];
rz(-0.7472907) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3355646) q[0];
sx q[0];
rz(-1.3845516) q[0];
sx q[0];
rz(-1.1769017) q[0];
rz(-pi) q[1];
rz(-2.4392082) q[2];
sx q[2];
rz(-1.1051854) q[2];
sx q[2];
rz(-1.7543751) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9979228) q[1];
sx q[1];
rz(-1.2005245) q[1];
sx q[1];
rz(1.0940432) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7036311) q[3];
sx q[3];
rz(-1.3107302) q[3];
sx q[3];
rz(-2.4822134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.0390465) q[2];
sx q[2];
rz(-1.1209844) q[2];
sx q[2];
rz(-0.54753629) q[2];
rz(1.840379) q[3];
sx q[3];
rz(-1.1373212) q[3];
sx q[3];
rz(3.014452) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8319156) q[0];
sx q[0];
rz(-1.3818106) q[0];
sx q[0];
rz(0.58050138) q[0];
rz(0.25513395) q[1];
sx q[1];
rz(-0.83108035) q[1];
sx q[1];
rz(1.1855804) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0609293) q[0];
sx q[0];
rz(-0.97968819) q[0];
sx q[0];
rz(-2.3800058) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9337148) q[2];
sx q[2];
rz(-1.5693773) q[2];
sx q[2];
rz(-0.05677536) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.23981389) q[1];
sx q[1];
rz(-1.1019152) q[1];
sx q[1];
rz(-1.4711051) q[1];
rz(-pi) q[2];
rz(-2.2907652) q[3];
sx q[3];
rz(-1.730207) q[3];
sx q[3];
rz(-0.042378332) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.7287951) q[2];
sx q[2];
rz(-0.86270657) q[2];
sx q[2];
rz(2.1346788) q[2];
rz(-1.5236731) q[3];
sx q[3];
rz(-0.98098743) q[3];
sx q[3];
rz(1.1716918) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11378743) q[0];
sx q[0];
rz(-1.2600949) q[0];
sx q[0];
rz(-2.8366198) q[0];
rz(-1.8995359) q[1];
sx q[1];
rz(-1.2868953) q[1];
sx q[1];
rz(-2.3931265) q[1];
rz(1.9340877) q[2];
sx q[2];
rz(-2.6836094) q[2];
sx q[2];
rz(-1.2140254) q[2];
rz(2.2119568) q[3];
sx q[3];
rz(-0.84992483) q[3];
sx q[3];
rz(0.58569943) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
