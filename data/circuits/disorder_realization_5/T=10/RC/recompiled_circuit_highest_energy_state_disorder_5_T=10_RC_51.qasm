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
rz(1.6485543) q[1];
sx q[1];
rz(4.0798408) q[1];
sx q[1];
rz(8.5473016) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2169133) q[0];
sx q[0];
rz(-1.7857505) q[0];
sx q[0];
rz(-2.0314905) q[0];
x q[1];
rz(1.9196366) q[2];
sx q[2];
rz(-1.9720805) q[2];
sx q[2];
rz(2.2711585) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.92447313) q[1];
sx q[1];
rz(-1.8397325) q[1];
sx q[1];
rz(-1.1094138) q[1];
rz(1.0935173) q[3];
sx q[3];
rz(-1.7831037) q[3];
sx q[3];
rz(2.4119854) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.33285546) q[2];
sx q[2];
rz(-2.4337807) q[2];
sx q[2];
rz(0.037192496) q[2];
rz(-2.228179) q[3];
sx q[3];
rz(-2.1593058) q[3];
sx q[3];
rz(-1.5090401) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7971802) q[0];
sx q[0];
rz(-2.0672815) q[0];
sx q[0];
rz(0.2221701) q[0];
rz(-0.19368681) q[1];
sx q[1];
rz(-1.2682468) q[1];
sx q[1];
rz(0.2598612) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79767449) q[0];
sx q[0];
rz(-2.0137383) q[0];
sx q[0];
rz(2.8721149) q[0];
x q[1];
rz(-1.3919984) q[2];
sx q[2];
rz(-1.4646718) q[2];
sx q[2];
rz(0.44580844) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.4603468) q[1];
sx q[1];
rz(-0.28239861) q[1];
sx q[1];
rz(-2.4241134) q[1];
rz(-pi) q[2];
rz(2.3823795) q[3];
sx q[3];
rz(-1.2046283) q[3];
sx q[3];
rz(3.1389583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.3746609) q[2];
sx q[2];
rz(-0.58284512) q[2];
sx q[2];
rz(0.15677162) q[2];
rz(-0.80592704) q[3];
sx q[3];
rz(-1.0904049) q[3];
sx q[3];
rz(2.6167615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
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
rz(0.90461516) q[0];
sx q[0];
rz(-0.1122864) q[0];
sx q[0];
rz(1.1861381) q[0];
rz(0.063749464) q[1];
sx q[1];
rz(-1.6155764) q[1];
sx q[1];
rz(-2.729111) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8500687) q[0];
sx q[0];
rz(-1.0837892) q[0];
sx q[0];
rz(2.1662232) q[0];
rz(2.2507452) q[2];
sx q[2];
rz(-1.6932704) q[2];
sx q[2];
rz(2.9735931) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.7802935) q[1];
sx q[1];
rz(-2.7971606) q[1];
sx q[1];
rz(2.5028178) q[1];
rz(2.5600912) q[3];
sx q[3];
rz(-2.1505754) q[3];
sx q[3];
rz(1.2691154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.5731262) q[2];
sx q[2];
rz(-2.7012479) q[2];
sx q[2];
rz(1.1661412) q[2];
rz(2.3447573) q[3];
sx q[3];
rz(-1.3692057) q[3];
sx q[3];
rz(-2.3603175) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9154938) q[0];
sx q[0];
rz(-1.508536) q[0];
sx q[0];
rz(-1.4133806) q[0];
rz(-0.49890292) q[1];
sx q[1];
rz(-1.7585124) q[1];
sx q[1];
rz(-1.2938719) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2436566) q[0];
sx q[0];
rz(-1.597628) q[0];
sx q[0];
rz(-1.4856337) q[0];
x q[1];
rz(1.7912553) q[2];
sx q[2];
rz(-0.33500698) q[2];
sx q[2];
rz(-2.9398244) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.41999757) q[1];
sx q[1];
rz(-1.5087869) q[1];
sx q[1];
rz(1.5739784) q[1];
x q[2];
rz(2.3975375) q[3];
sx q[3];
rz(-0.091587154) q[3];
sx q[3];
rz(1.1055783) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.6125907) q[2];
sx q[2];
rz(-2.4606073) q[2];
sx q[2];
rz(-0.16981086) q[2];
rz(-2.7134907) q[3];
sx q[3];
rz(-1.4295108) q[3];
sx q[3];
rz(2.9812109) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.661628) q[0];
sx q[0];
rz(-2.7492838) q[0];
sx q[0];
rz(0.83537927) q[0];
rz(-0.63940489) q[1];
sx q[1];
rz(-2.4765922) q[1];
sx q[1];
rz(-1.0909874) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3314455) q[0];
sx q[0];
rz(-1.5584599) q[0];
sx q[0];
rz(-2.4452565) q[0];
x q[1];
rz(-0.29971896) q[2];
sx q[2];
rz(-1.4106361) q[2];
sx q[2];
rz(-1.0417787) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.94746339) q[1];
sx q[1];
rz(-0.76912472) q[1];
sx q[1];
rz(-2.1544653) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0939319) q[3];
sx q[3];
rz(-1.2854648) q[3];
sx q[3];
rz(2.3714575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.368448) q[2];
sx q[2];
rz(-1.1870563) q[2];
sx q[2];
rz(-0.0040357987) q[2];
rz(-1.076738) q[3];
sx q[3];
rz(-0.69759798) q[3];
sx q[3];
rz(2.9546886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3183597) q[0];
sx q[0];
rz(-1.3613181) q[0];
sx q[0];
rz(2.5079492) q[0];
rz(2.7404495) q[1];
sx q[1];
rz(-0.70924962) q[1];
sx q[1];
rz(-2.4493682) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5790235) q[0];
sx q[0];
rz(-0.85183598) q[0];
sx q[0];
rz(-2.1301053) q[0];
x q[1];
rz(-1.2622617) q[2];
sx q[2];
rz(-0.82595982) q[2];
sx q[2];
rz(2.8070297) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.2894404) q[1];
sx q[1];
rz(-1.0045144) q[1];
sx q[1];
rz(0.29636611) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7217466) q[3];
sx q[3];
rz(-1.387409) q[3];
sx q[3];
rz(0.91461411) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.20574084) q[2];
sx q[2];
rz(-2.3571099) q[2];
sx q[2];
rz(-1.240823) q[2];
rz(-1.9567018) q[3];
sx q[3];
rz(-1.3013867) q[3];
sx q[3];
rz(-1.5968116) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1320837) q[0];
sx q[0];
rz(-1.8385831) q[0];
sx q[0];
rz(-1.3707772) q[0];
rz(-0.41807434) q[1];
sx q[1];
rz(-1.6590174) q[1];
sx q[1];
rz(-0.48666993) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70145973) q[0];
sx q[0];
rz(-2.1281689) q[0];
sx q[0];
rz(-3.0270744) q[0];
rz(-1.0202706) q[2];
sx q[2];
rz(-1.8943059) q[2];
sx q[2];
rz(-0.24187096) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.8312953) q[1];
sx q[1];
rz(-2.8716209) q[1];
sx q[1];
rz(-0.84862535) q[1];
rz(-pi) q[2];
rz(-1.6177931) q[3];
sx q[3];
rz(-1.1011657) q[3];
sx q[3];
rz(0.75799537) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.70249867) q[2];
sx q[2];
rz(-1.0489901) q[2];
sx q[2];
rz(2.9443963) q[2];
rz(0.50944734) q[3];
sx q[3];
rz(-0.26064894) q[3];
sx q[3];
rz(0.82836866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2800804) q[0];
sx q[0];
rz(-1.0741638) q[0];
sx q[0];
rz(1.3664838) q[0];
rz(-0.81661433) q[1];
sx q[1];
rz(-1.0895224) q[1];
sx q[1];
rz(-2.8588967) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44756457) q[0];
sx q[0];
rz(-2.1622133) q[0];
sx q[0];
rz(-2.924218) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8545998) q[2];
sx q[2];
rz(-2.8176135) q[2];
sx q[2];
rz(-2.2215077) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.3779439) q[1];
sx q[1];
rz(-2.4326583) q[1];
sx q[1];
rz(1.8590742) q[1];
rz(-pi) q[2];
rz(1.9603086) q[3];
sx q[3];
rz(-2.0593025) q[3];
sx q[3];
rz(-0.28340411) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(3.1064328) q[2];
sx q[2];
rz(-2.1356434) q[2];
sx q[2];
rz(2.1584568) q[2];
rz(-1.8631009) q[3];
sx q[3];
rz(-0.92517868) q[3];
sx q[3];
rz(-0.92888752) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66016692) q[0];
sx q[0];
rz(-1.7919414) q[0];
sx q[0];
rz(2.8283258) q[0];
rz(-2.616864) q[1];
sx q[1];
rz(-2.8786761) q[1];
sx q[1];
rz(1.252334) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9565283) q[0];
sx q[0];
rz(-1.6913497) q[0];
sx q[0];
rz(-2.610449) q[0];
rz(-pi) q[1];
x q[1];
rz(0.98964844) q[2];
sx q[2];
rz(-2.2454442) q[2];
sx q[2];
rz(-1.0824114) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.5572667) q[1];
sx q[1];
rz(-1.0385822) q[1];
sx q[1];
rz(0.6289215) q[1];
rz(-pi) q[2];
rz(-2.325237) q[3];
sx q[3];
rz(-1.6301148) q[3];
sx q[3];
rz(-0.99494464) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.9423882) q[2];
sx q[2];
rz(-0.62493268) q[2];
sx q[2];
rz(1.7601298) q[2];
rz(2.840461) q[3];
sx q[3];
rz(-0.98265177) q[3];
sx q[3];
rz(0.38100955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9075266) q[0];
sx q[0];
rz(-1.0028239) q[0];
sx q[0];
rz(1.3475077) q[0];
rz(2.2566336) q[1];
sx q[1];
rz(-2.5159409) q[1];
sx q[1];
rz(-1.7431097) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32537726) q[0];
sx q[0];
rz(-2.7669816) q[0];
sx q[0];
rz(1.5416359) q[0];
rz(-pi) q[1];
x q[1];
rz(0.86549564) q[2];
sx q[2];
rz(-1.7112321) q[2];
sx q[2];
rz(2.8704134) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.2513548) q[1];
sx q[1];
rz(-1.4904463) q[1];
sx q[1];
rz(-1.608196) q[1];
rz(0.22866727) q[3];
sx q[3];
rz(-2.1211984) q[3];
sx q[3];
rz(2.500734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.053190319) q[2];
sx q[2];
rz(-1.3206427) q[2];
sx q[2];
rz(0.94584805) q[2];
rz(-2.451402) q[3];
sx q[3];
rz(-1.2226356) q[3];
sx q[3];
rz(-1.8405731) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-2.239554) q[0];
sx q[0];
rz(-1.6033462) q[0];
sx q[0];
rz(-0.70815804) q[0];
rz(-3.0802849) q[1];
sx q[1];
rz(-0.61534449) q[1];
sx q[1];
rz(-1.4475488) q[1];
rz(-0.29091111) q[2];
sx q[2];
rz(-1.6500191) q[2];
sx q[2];
rz(2.0943841) q[2];
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
