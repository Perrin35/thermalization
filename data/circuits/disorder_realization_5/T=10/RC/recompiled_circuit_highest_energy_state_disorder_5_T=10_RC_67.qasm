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
rz(-1.9755826) q[0];
sx q[0];
rz(-2.2628885) q[0];
sx q[0];
rz(2.9455844) q[0];
rz(-1.4930383) q[1];
sx q[1];
rz(-0.9382481) q[1];
sx q[1];
rz(0.87747639) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24841079) q[0];
sx q[0];
rz(-2.0201004) q[0];
sx q[0];
rz(-0.23907678) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2219561) q[2];
sx q[2];
rz(-1.9720805) q[2];
sx q[2];
rz(-0.87043412) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.2171195) q[1];
sx q[1];
rz(-1.8397325) q[1];
sx q[1];
rz(-2.0321789) q[1];
x q[2];
rz(-1.0935173) q[3];
sx q[3];
rz(-1.358489) q[3];
sx q[3];
rz(-0.72960723) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8087372) q[2];
sx q[2];
rz(-0.70781195) q[2];
sx q[2];
rz(0.037192496) q[2];
rz(-2.228179) q[3];
sx q[3];
rz(-0.98228684) q[3];
sx q[3];
rz(1.5090401) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(2.7971802) q[0];
sx q[0];
rz(-2.0672815) q[0];
sx q[0];
rz(-0.2221701) q[0];
rz(0.19368681) q[1];
sx q[1];
rz(-1.2682468) q[1];
sx q[1];
rz(2.8817315) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79767449) q[0];
sx q[0];
rz(-1.1278544) q[0];
sx q[0];
rz(-2.8721149) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1104576) q[2];
sx q[2];
rz(-0.20763131) q[2];
sx q[2];
rz(-2.5467444) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.7227962) q[1];
sx q[1];
rz(-1.7823311) q[1];
sx q[1];
rz(1.7593033) q[1];
x q[2];
rz(2.6333359) q[3];
sx q[3];
rz(-0.82672366) q[3];
sx q[3];
rz(-1.9342157) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3746609) q[2];
sx q[2];
rz(-0.58284512) q[2];
sx q[2];
rz(0.15677162) q[2];
rz(-0.80592704) q[3];
sx q[3];
rz(-2.0511878) q[3];
sx q[3];
rz(0.52483112) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90461516) q[0];
sx q[0];
rz(-0.1122864) q[0];
sx q[0];
rz(1.9554546) q[0];
rz(0.063749464) q[1];
sx q[1];
rz(-1.6155764) q[1];
sx q[1];
rz(-2.729111) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97226364) q[0];
sx q[0];
rz(-2.0893851) q[0];
sx q[0];
rz(-0.56904582) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2507452) q[2];
sx q[2];
rz(-1.6932704) q[2];
sx q[2];
rz(-2.9735931) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7802935) q[1];
sx q[1];
rz(-0.34443203) q[1];
sx q[1];
rz(-2.5028178) q[1];
rz(-2.2687421) q[3];
sx q[3];
rz(-0.79668364) q[3];
sx q[3];
rz(-2.7484675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.56846642) q[2];
sx q[2];
rz(-2.7012479) q[2];
sx q[2];
rz(-1.1661412) q[2];
rz(2.3447573) q[3];
sx q[3];
rz(-1.7723869) q[3];
sx q[3];
rz(-0.78127512) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9154938) q[0];
sx q[0];
rz(-1.6330566) q[0];
sx q[0];
rz(-1.728212) q[0];
rz(2.6426897) q[1];
sx q[1];
rz(-1.3830802) q[1];
sx q[1];
rz(-1.8477207) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2436566) q[0];
sx q[0];
rz(-1.5439646) q[0];
sx q[0];
rz(1.4856337) q[0];
rz(1.8982688) q[2];
sx q[2];
rz(-1.6427543) q[2];
sx q[2];
rz(1.5775934) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.7728987) q[1];
sx q[1];
rz(-3.0795018) q[1];
sx q[1];
rz(-0.05120488) q[1];
rz(3.0741229) q[3];
sx q[3];
rz(-1.5088135) q[3];
sx q[3];
rz(0.27674473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.52900195) q[2];
sx q[2];
rz(-2.4606073) q[2];
sx q[2];
rz(2.9717818) q[2];
rz(0.42810193) q[3];
sx q[3];
rz(-1.7120818) q[3];
sx q[3];
rz(-2.9812109) q[3];
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
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47996461) q[0];
sx q[0];
rz(-2.7492838) q[0];
sx q[0];
rz(-0.83537927) q[0];
rz(-0.63940489) q[1];
sx q[1];
rz(-2.4765922) q[1];
sx q[1];
rz(2.0506052) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3314455) q[0];
sx q[0];
rz(-1.5584599) q[0];
sx q[0];
rz(-0.69633616) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6409482) q[2];
sx q[2];
rz(-2.8029059) q[2];
sx q[2];
rz(-3.088986) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.6907566) q[1];
sx q[1];
rz(-0.95162205) q[1];
sx q[1];
rz(0.49003933) q[1];
rz(1.0398846) q[3];
sx q[3];
rz(-0.58947488) q[3];
sx q[3];
rz(-2.7949751) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.368448) q[2];
sx q[2];
rz(-1.1870563) q[2];
sx q[2];
rz(-3.1375569) q[2];
rz(-1.076738) q[3];
sx q[3];
rz(-0.69759798) q[3];
sx q[3];
rz(2.9546886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3183597) q[0];
sx q[0];
rz(-1.7802745) q[0];
sx q[0];
rz(2.5079492) q[0];
rz(2.7404495) q[1];
sx q[1];
rz(-2.432343) q[1];
sx q[1];
rz(2.4493682) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5790235) q[0];
sx q[0];
rz(-0.85183598) q[0];
sx q[0];
rz(-1.0114874) q[0];
x q[1];
rz(0.31816407) q[2];
sx q[2];
rz(-0.79472322) q[2];
sx q[2];
rz(0.77407167) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-3.0226188) q[1];
sx q[1];
rz(-1.8198208) q[1];
sx q[1];
rz(-0.98414661) q[1];
rz(-1.7217466) q[3];
sx q[3];
rz(-1.7541837) q[3];
sx q[3];
rz(-2.2269785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.9358518) q[2];
sx q[2];
rz(-2.3571099) q[2];
sx q[2];
rz(1.9007696) q[2];
rz(-1.9567018) q[3];
sx q[3];
rz(-1.840206) q[3];
sx q[3];
rz(1.5968116) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1320837) q[0];
sx q[0];
rz(-1.8385831) q[0];
sx q[0];
rz(1.3707772) q[0];
rz(2.7235183) q[1];
sx q[1];
rz(-1.6590174) q[1];
sx q[1];
rz(-0.48666993) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70145973) q[0];
sx q[0];
rz(-2.1281689) q[0];
sx q[0];
rz(-0.11451829) q[0];
rz(-pi) q[1];
x q[1];
rz(0.37481793) q[2];
sx q[2];
rz(-2.0897802) q[2];
sx q[2];
rz(-1.5216375) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.0508065) q[1];
sx q[1];
rz(-1.7722844) q[1];
sx q[1];
rz(2.9606716) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.671516) q[3];
sx q[3];
rz(-1.5288908) q[3];
sx q[3];
rz(2.3500729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.70249867) q[2];
sx q[2];
rz(-2.0926026) q[2];
sx q[2];
rz(0.19719633) q[2];
rz(0.50944734) q[3];
sx q[3];
rz(-2.8809437) q[3];
sx q[3];
rz(2.313224) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2800804) q[0];
sx q[0];
rz(-2.0674288) q[0];
sx q[0];
rz(-1.7751088) q[0];
rz(0.81661433) q[1];
sx q[1];
rz(-2.0520703) q[1];
sx q[1];
rz(-2.8588967) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1408843) q[0];
sx q[0];
rz(-1.750803) q[0];
sx q[0];
rz(-2.1732844) q[0];
rz(-pi) q[1];
rz(2.8545998) q[2];
sx q[2];
rz(-2.8176135) q[2];
sx q[2];
rz(0.92008495) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.41425175) q[1];
sx q[1];
rz(-1.7569572) q[1];
sx q[1];
rz(-2.258979) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6201309) q[3];
sx q[3];
rz(-1.2288501) q[3];
sx q[3];
rz(-1.6638883) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.1064328) q[2];
sx q[2];
rz(-1.0059493) q[2];
sx q[2];
rz(2.1584568) q[2];
rz(1.2784917) q[3];
sx q[3];
rz(-2.216414) q[3];
sx q[3];
rz(0.92888752) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(0.66016692) q[0];
sx q[0];
rz(-1.3496512) q[0];
sx q[0];
rz(-2.8283258) q[0];
rz(2.616864) q[1];
sx q[1];
rz(-2.8786761) q[1];
sx q[1];
rz(1.8892586) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6853297) q[0];
sx q[0];
rz(-2.0976817) q[0];
sx q[0];
rz(1.4312126) q[0];
rz(2.1519442) q[2];
sx q[2];
rz(-2.2454442) q[2];
sx q[2];
rz(-2.0591813) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5572667) q[1];
sx q[1];
rz(-1.0385822) q[1];
sx q[1];
rz(2.5126712) q[1];
rz(-pi) q[2];
rz(1.657293) q[3];
sx q[3];
rz(-0.75631006) q[3];
sx q[3];
rz(-2.5027517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.19920443) q[2];
sx q[2];
rz(-2.51666) q[2];
sx q[2];
rz(1.7601298) q[2];
rz(-0.30113164) q[3];
sx q[3];
rz(-0.98265177) q[3];
sx q[3];
rz(-2.7605831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.23406601) q[0];
sx q[0];
rz(-1.0028239) q[0];
sx q[0];
rz(1.7940849) q[0];
rz(2.2566336) q[1];
sx q[1];
rz(-2.5159409) q[1];
sx q[1];
rz(1.398483) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8475474) q[0];
sx q[0];
rz(-1.1963524) q[0];
sx q[0];
rz(-3.1301296) q[0];
rz(-pi) q[1];
rz(-2.9580206) q[2];
sx q[2];
rz(-2.2677448) q[2];
sx q[2];
rz(1.1810034) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.3264638) q[1];
sx q[1];
rz(-0.088610709) q[1];
sx q[1];
rz(2.7068702) q[1];
x q[2];
rz(1.0085671) q[3];
sx q[3];
rz(-1.3763714) q[3];
sx q[3];
rz(-1.0510707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.053190319) q[2];
sx q[2];
rz(-1.8209499) q[2];
sx q[2];
rz(2.1957446) q[2];
rz(2.451402) q[3];
sx q[3];
rz(-1.918957) q[3];
sx q[3];
rz(1.3010196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.239554) q[0];
sx q[0];
rz(-1.6033462) q[0];
sx q[0];
rz(-0.70815804) q[0];
rz(3.0802849) q[1];
sx q[1];
rz(-2.5262482) q[1];
sx q[1];
rz(1.6940438) q[1];
rz(0.270025) q[2];
sx q[2];
rz(-0.30120987) q[2];
sx q[2];
rz(0.78200151) q[2];
rz(-0.072515247) q[3];
sx q[3];
rz(-1.391165) q[3];
sx q[3];
rz(2.388849) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
