OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.6239983) q[0];
sx q[0];
rz(2.6156293) q[0];
sx q[0];
rz(9.2064657) q[0];
rz(-1.6649618) q[1];
sx q[1];
rz(2.6638439) q[1];
sx q[1];
rz(12.01241) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4628578) q[0];
sx q[0];
rz(-2.4241872) q[0];
sx q[0];
rz(-0.82490246) q[0];
rz(-pi) q[1];
rz(3.0402115) q[2];
sx q[2];
rz(-0.53115618) q[2];
sx q[2];
rz(2.6918333) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.7377555) q[1];
sx q[1];
rz(-2.3970876) q[1];
sx q[1];
rz(-1.8831403) q[1];
rz(-pi) q[2];
rz(-1.1470818) q[3];
sx q[3];
rz(-2.0689575) q[3];
sx q[3];
rz(-0.26396449) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.7654984) q[2];
sx q[2];
rz(-2.8480397) q[2];
sx q[2];
rz(1.2676839) q[2];
rz(-2.3119161) q[3];
sx q[3];
rz(-1.6206348) q[3];
sx q[3];
rz(2.7177496) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.053452881) q[0];
sx q[0];
rz(-1.3351048) q[0];
sx q[0];
rz(1.394519) q[0];
rz(-0.8949737) q[1];
sx q[1];
rz(-2.112969) q[1];
sx q[1];
rz(1.5590394) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.549986) q[0];
sx q[0];
rz(-0.9261407) q[0];
sx q[0];
rz(-2.2695973) q[0];
rz(-pi) q[1];
rz(2.2649647) q[2];
sx q[2];
rz(-2.6149984) q[2];
sx q[2];
rz(-1.2132298) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.5176508) q[1];
sx q[1];
rz(-1.883916) q[1];
sx q[1];
rz(0.70065686) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1760494) q[3];
sx q[3];
rz(-2.6355834) q[3];
sx q[3];
rz(3.0288896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.47722694) q[2];
sx q[2];
rz(-1.5218647) q[2];
sx q[2];
rz(3.1108372) q[2];
rz(2.6886046) q[3];
sx q[3];
rz(-0.22946295) q[3];
sx q[3];
rz(0.17791137) q[3];
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
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7396616) q[0];
sx q[0];
rz(-3.0406096) q[0];
sx q[0];
rz(-2.3133551) q[0];
rz(-3.0939057) q[1];
sx q[1];
rz(-0.86367718) q[1];
sx q[1];
rz(-1.2275068) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0704884) q[0];
sx q[0];
rz(-1.0336735) q[0];
sx q[0];
rz(2.0893196) q[0];
x q[1];
rz(0.21593185) q[2];
sx q[2];
rz(-1.1732444) q[2];
sx q[2];
rz(2.2550607) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.8125898) q[1];
sx q[1];
rz(-2.0101476) q[1];
sx q[1];
rz(2.8529608) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5121721) q[3];
sx q[3];
rz(-0.72949648) q[3];
sx q[3];
rz(0.54704715) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(3.0487655) q[2];
sx q[2];
rz(-0.91492492) q[2];
sx q[2];
rz(1.1009334) q[2];
rz(-3.1277505) q[3];
sx q[3];
rz(-1.3661386) q[3];
sx q[3];
rz(0.83465105) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5562627) q[0];
sx q[0];
rz(-1.7381373) q[0];
sx q[0];
rz(-0.334326) q[0];
rz(0.73257929) q[1];
sx q[1];
rz(-0.94894797) q[1];
sx q[1];
rz(-0.11071959) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3363655) q[0];
sx q[0];
rz(-0.76167652) q[0];
sx q[0];
rz(3.1195927) q[0];
x q[1];
rz(1.6640856) q[2];
sx q[2];
rz(-1.7356725) q[2];
sx q[2];
rz(-3.0704569) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.2313614) q[1];
sx q[1];
rz(-2.1276603) q[1];
sx q[1];
rz(-0.75480144) q[1];
x q[2];
rz(-2.3663051) q[3];
sx q[3];
rz(-1.9248157) q[3];
sx q[3];
rz(-0.86590761) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.4577786) q[2];
sx q[2];
rz(-0.28202287) q[2];
sx q[2];
rz(1.7337743) q[2];
rz(-1.9862566) q[3];
sx q[3];
rz(-1.9852394) q[3];
sx q[3];
rz(0.19481625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4487576) q[0];
sx q[0];
rz(-2.223707) q[0];
sx q[0];
rz(-2.1441929) q[0];
rz(1.6150486) q[1];
sx q[1];
rz(-2.5020182) q[1];
sx q[1];
rz(-1.4195199) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90516312) q[0];
sx q[0];
rz(-1.7337203) q[0];
sx q[0];
rz(1.7927367) q[0];
x q[1];
rz(2.3676374) q[2];
sx q[2];
rz(-1.4593235) q[2];
sx q[2];
rz(-0.74711266) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.43806009) q[1];
sx q[1];
rz(-0.72444455) q[1];
sx q[1];
rz(-2.0634335) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6035242) q[3];
sx q[3];
rz(-1.447289) q[3];
sx q[3];
rz(1.1203958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.851696) q[2];
sx q[2];
rz(-0.09446129) q[2];
sx q[2];
rz(0.32290253) q[2];
rz(2.0276535) q[3];
sx q[3];
rz(-2.0506004) q[3];
sx q[3];
rz(-0.41306257) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36393976) q[0];
sx q[0];
rz(-2.8033065) q[0];
sx q[0];
rz(-2.3349578) q[0];
rz(-0.58397645) q[1];
sx q[1];
rz(-1.1176502) q[1];
sx q[1];
rz(-1.1899828) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89878824) q[0];
sx q[0];
rz(-1.4477647) q[0];
sx q[0];
rz(1.9210299) q[0];
rz(-pi) q[1];
rz(2.7940668) q[2];
sx q[2];
rz(-0.56171562) q[2];
sx q[2];
rz(2.7507741) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5197503) q[1];
sx q[1];
rz(-1.519366) q[1];
sx q[1];
rz(-1.4019045) q[1];
x q[2];
rz(-1.7403931) q[3];
sx q[3];
rz(-0.84753321) q[3];
sx q[3];
rz(0.71454988) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.40818647) q[2];
sx q[2];
rz(-1.5032282) q[2];
sx q[2];
rz(-2.9564986) q[2];
rz(1.5271651) q[3];
sx q[3];
rz(-1.4027169) q[3];
sx q[3];
rz(0.68814284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9516893) q[0];
sx q[0];
rz(-0.95174319) q[0];
sx q[0];
rz(0.21251799) q[0];
rz(-2.3566133) q[1];
sx q[1];
rz(-1.6467983) q[1];
sx q[1];
rz(-2.3627538) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23379743) q[0];
sx q[0];
rz(-1.0186695) q[0];
sx q[0];
rz(2.2566081) q[0];
x q[1];
rz(0.58765192) q[2];
sx q[2];
rz(-1.8853123) q[2];
sx q[2];
rz(1.3101729) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.6952371) q[1];
sx q[1];
rz(-1.3160222) q[1];
sx q[1];
rz(0.22344113) q[1];
rz(-pi) q[2];
rz(-2.1206585) q[3];
sx q[3];
rz(-2.3519197) q[3];
sx q[3];
rz(2.0146973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.6276041) q[2];
sx q[2];
rz(-0.49391654) q[2];
sx q[2];
rz(0.40840515) q[2];
rz(-2.4397395) q[3];
sx q[3];
rz(-0.93188325) q[3];
sx q[3];
rz(-1.2592038) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7931165) q[0];
sx q[0];
rz(-1.2154673) q[0];
sx q[0];
rz(3.075573) q[0];
rz(-1.5090212) q[1];
sx q[1];
rz(-2.0965818) q[1];
sx q[1];
rz(-2.1836233) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5292266) q[0];
sx q[0];
rz(-0.54915308) q[0];
sx q[0];
rz(-0.56580881) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3249469) q[2];
sx q[2];
rz(-1.1743675) q[2];
sx q[2];
rz(-2.3289837) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.9878475) q[1];
sx q[1];
rz(-1.8829131) q[1];
sx q[1];
rz(-2.4448443) q[1];
rz(-pi) q[2];
x q[2];
rz(0.10029467) q[3];
sx q[3];
rz(-1.7236606) q[3];
sx q[3];
rz(-2.4748442) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.3079188) q[2];
sx q[2];
rz(-1.6180399) q[2];
sx q[2];
rz(-1.7306805) q[2];
rz(0.69502568) q[3];
sx q[3];
rz(-1.6618988) q[3];
sx q[3];
rz(-3.1237349) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0787635) q[0];
sx q[0];
rz(-1.48209) q[0];
sx q[0];
rz(-2.1516946) q[0];
rz(0.46317378) q[1];
sx q[1];
rz(-1.4521234) q[1];
sx q[1];
rz(-1.90082) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1027066) q[0];
sx q[0];
rz(-0.30065824) q[0];
sx q[0];
rz(1.7666398) q[0];
x q[1];
rz(-0.2026338) q[2];
sx q[2];
rz(-1.1362193) q[2];
sx q[2];
rz(3.0960577) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.74490101) q[1];
sx q[1];
rz(-2.7518175) q[1];
sx q[1];
rz(0.74662965) q[1];
rz(-pi) q[2];
x q[2];
rz(0.95043358) q[3];
sx q[3];
rz(-1.1945121) q[3];
sx q[3];
rz(-2.2971414) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.8076294) q[2];
sx q[2];
rz(-1.7619851) q[2];
sx q[2];
rz(0.71869746) q[2];
rz(-1.9888318) q[3];
sx q[3];
rz(-1.7199687) q[3];
sx q[3];
rz(-2.1271472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54866791) q[0];
sx q[0];
rz(-2.7935226) q[0];
sx q[0];
rz(0.80192178) q[0];
rz(-2.0536664) q[1];
sx q[1];
rz(-1.3366924) q[1];
sx q[1];
rz(-0.71802872) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82336003) q[0];
sx q[0];
rz(-1.2498858) q[0];
sx q[0];
rz(2.4145187) q[0];
rz(-0.74943351) q[2];
sx q[2];
rz(-2.4735056) q[2];
sx q[2];
rz(-1.6783448) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.57783571) q[1];
sx q[1];
rz(-0.33591649) q[1];
sx q[1];
rz(1.729924) q[1];
rz(2.2949785) q[3];
sx q[3];
rz(-0.9808971) q[3];
sx q[3];
rz(2.5835832) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.7903018) q[2];
sx q[2];
rz(-2.9238034) q[2];
sx q[2];
rz(-2.6289319) q[2];
rz(0.74448186) q[3];
sx q[3];
rz(-2.1148041) q[3];
sx q[3];
rz(-2.3775533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22513334) q[0];
sx q[0];
rz(-1.2703348) q[0];
sx q[0];
rz(-1.2402007) q[0];
rz(2.4254639) q[1];
sx q[1];
rz(-0.54812535) q[1];
sx q[1];
rz(-2.5352238) q[1];
rz(1.345558) q[2];
sx q[2];
rz(-1.6788531) q[2];
sx q[2];
rz(-2.0137871) q[2];
rz(1.2848787) q[3];
sx q[3];
rz(-0.34964041) q[3];
sx q[3];
rz(1.6007363) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
