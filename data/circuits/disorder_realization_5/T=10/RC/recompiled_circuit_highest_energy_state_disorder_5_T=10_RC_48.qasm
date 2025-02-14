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
rz(1.16601) q[0];
sx q[0];
rz(2.2628885) q[0];
sx q[0];
rz(9.2287697) q[0];
rz(-1.4930383) q[1];
sx q[1];
rz(-0.9382481) q[1];
sx q[1];
rz(0.87747639) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9246793) q[0];
sx q[0];
rz(-1.3558421) q[0];
sx q[0];
rz(-2.0314905) q[0];
rz(2.7174905) q[2];
sx q[2];
rz(-1.2507032) q[2];
sx q[2];
rz(2.3001075) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.6266105) q[1];
sx q[1];
rz(-1.1272073) q[1];
sx q[1];
rz(-0.29859467) q[1];
rz(-pi) q[2];
rz(-1.1320513) q[3];
sx q[3];
rz(-2.6225704) q[3];
sx q[3];
rz(-1.2281017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.8087372) q[2];
sx q[2];
rz(-2.4337807) q[2];
sx q[2];
rz(-3.1044002) q[2];
rz(2.228179) q[3];
sx q[3];
rz(-2.1593058) q[3];
sx q[3];
rz(-1.6325525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7971802) q[0];
sx q[0];
rz(-2.0672815) q[0];
sx q[0];
rz(0.2221701) q[0];
rz(0.19368681) q[1];
sx q[1];
rz(-1.8733459) q[1];
sx q[1];
rz(0.2598612) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3439182) q[0];
sx q[0];
rz(-1.1278544) q[0];
sx q[0];
rz(-2.8721149) q[0];
x q[1];
rz(1.031135) q[2];
sx q[2];
rz(-2.9339613) q[2];
sx q[2];
rz(-0.59484824) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.7227962) q[1];
sx q[1];
rz(-1.3592615) q[1];
sx q[1];
rz(1.3822894) q[1];
rz(-pi) q[2];
rz(-2.0570898) q[3];
sx q[3];
rz(-2.2688365) q[3];
sx q[3];
rz(-1.2459038) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.3746609) q[2];
sx q[2];
rz(-2.5587475) q[2];
sx q[2];
rz(-0.15677162) q[2];
rz(-0.80592704) q[3];
sx q[3];
rz(-1.0904049) q[3];
sx q[3];
rz(-0.52483112) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2369775) q[0];
sx q[0];
rz(-0.1122864) q[0];
sx q[0];
rz(1.1861381) q[0];
rz(3.0778432) q[1];
sx q[1];
rz(-1.5260162) q[1];
sx q[1];
rz(-2.729111) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.169329) q[0];
sx q[0];
rz(-2.0893851) q[0];
sx q[0];
rz(0.56904582) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3774728) q[2];
sx q[2];
rz(-2.4524322) q[2];
sx q[2];
rz(-1.2528407) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8349065) q[1];
sx q[1];
rz(-1.2962771) q[1];
sx q[1];
rz(-1.7815018) q[1];
rz(0.90610151) q[3];
sx q[3];
rz(-2.0482488) q[3];
sx q[3];
rz(-2.4942644) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.5731262) q[2];
sx q[2];
rz(-0.44034475) q[2];
sx q[2];
rz(-1.1661412) q[2];
rz(0.79683534) q[3];
sx q[3];
rz(-1.7723869) q[3];
sx q[3];
rz(0.78127512) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2260988) q[0];
sx q[0];
rz(-1.508536) q[0];
sx q[0];
rz(-1.728212) q[0];
rz(-2.6426897) q[1];
sx q[1];
rz(-1.7585124) q[1];
sx q[1];
rz(-1.8477207) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5099611) q[0];
sx q[0];
rz(-0.089279739) q[0];
sx q[0];
rz(1.2651612) q[0];
rz(-1.2433238) q[2];
sx q[2];
rz(-1.6427543) q[2];
sx q[2];
rz(1.5775934) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.1509959) q[1];
sx q[1];
rz(-1.5739723) q[1];
sx q[1];
rz(3.0795829) q[1];
rz(-pi) q[2];
rz(-1.5086725) q[3];
sx q[3];
rz(-1.5034564) q[3];
sx q[3];
rz(-1.289866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.52900195) q[2];
sx q[2];
rz(-0.68098536) q[2];
sx q[2];
rz(0.16981086) q[2];
rz(0.42810193) q[3];
sx q[3];
rz(-1.7120818) q[3];
sx q[3];
rz(-2.9812109) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47996461) q[0];
sx q[0];
rz(-0.39230883) q[0];
sx q[0];
rz(-2.3062134) q[0];
rz(2.5021878) q[1];
sx q[1];
rz(-0.6650005) q[1];
sx q[1];
rz(1.0909874) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8101472) q[0];
sx q[0];
rz(-1.5584599) q[0];
sx q[0];
rz(-0.69633616) q[0];
x q[1];
rz(1.7382938) q[2];
sx q[2];
rz(-1.27503) q[2];
sx q[2];
rz(0.57826051) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6907566) q[1];
sx q[1];
rz(-0.95162205) q[1];
sx q[1];
rz(-0.49003933) q[1];
rz(-pi) q[2];
rz(-1.0476607) q[3];
sx q[3];
rz(-1.8561279) q[3];
sx q[3];
rz(-0.77013515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.368448) q[2];
sx q[2];
rz(-1.9545363) q[2];
sx q[2];
rz(3.1375569) q[2];
rz(-1.076738) q[3];
sx q[3];
rz(-0.69759798) q[3];
sx q[3];
rz(-0.18690404) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8232329) q[0];
sx q[0];
rz(-1.3613181) q[0];
sx q[0];
rz(-2.5079492) q[0];
rz(2.7404495) q[1];
sx q[1];
rz(-2.432343) q[1];
sx q[1];
rz(-0.69222442) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8190097) q[0];
sx q[0];
rz(-0.87912175) q[0];
sx q[0];
rz(-2.5965967) q[0];
rz(-2.8234286) q[2];
sx q[2];
rz(-0.79472322) q[2];
sx q[2];
rz(0.77407167) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.8068829) q[1];
sx q[1];
rz(-0.63155424) q[1];
sx q[1];
rz(-1.1401661) q[1];
x q[2];
rz(-2.9561437) q[3];
sx q[3];
rz(-1.4223961) q[3];
sx q[3];
rz(-2.4576791) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.9358518) q[2];
sx q[2];
rz(-0.78448272) q[2];
sx q[2];
rz(1.240823) q[2];
rz(-1.1848909) q[3];
sx q[3];
rz(-1.3013867) q[3];
sx q[3];
rz(-1.544781) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-2.1320837) q[0];
sx q[0];
rz(-1.8385831) q[0];
sx q[0];
rz(-1.3707772) q[0];
rz(2.7235183) q[1];
sx q[1];
rz(-1.6590174) q[1];
sx q[1];
rz(-0.48666993) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3330227) q[0];
sx q[0];
rz(-1.6679224) q[0];
sx q[0];
rz(1.0104695) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.37481793) q[2];
sx q[2];
rz(-2.0897802) q[2];
sx q[2];
rz(-1.6199552) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.0907862) q[1];
sx q[1];
rz(-1.3693083) q[1];
sx q[1];
rz(-0.18092107) q[1];
x q[2];
rz(-0.47007665) q[3];
sx q[3];
rz(-1.6127018) q[3];
sx q[3];
rz(2.3500729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.70249867) q[2];
sx q[2];
rz(-1.0489901) q[2];
sx q[2];
rz(-2.9443963) q[2];
rz(0.50944734) q[3];
sx q[3];
rz(-0.26064894) q[3];
sx q[3];
rz(-2.313224) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2800804) q[0];
sx q[0];
rz(-1.0741638) q[0];
sx q[0];
rz(-1.7751088) q[0];
rz(-2.3249783) q[1];
sx q[1];
rz(-1.0895224) q[1];
sx q[1];
rz(-0.28269592) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44756457) q[0];
sx q[0];
rz(-2.1622133) q[0];
sx q[0];
rz(-0.21737463) q[0];
rz(-pi) q[1];
rz(2.8300072) q[2];
sx q[2];
rz(-1.4805613) q[2];
sx q[2];
rz(-2.2180706) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.3779439) q[1];
sx q[1];
rz(-2.4326583) q[1];
sx q[1];
rz(-1.8590742) q[1];
rz(-pi) q[2];
rz(-0.52146179) q[3];
sx q[3];
rz(-1.2288501) q[3];
sx q[3];
rz(1.6638883) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.1064328) q[2];
sx q[2];
rz(-1.0059493) q[2];
sx q[2];
rz(-2.1584568) q[2];
rz(1.8631009) q[3];
sx q[3];
rz(-2.216414) q[3];
sx q[3];
rz(-0.92888752) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66016692) q[0];
sx q[0];
rz(-1.3496512) q[0];
sx q[0];
rz(-2.8283258) q[0];
rz(-0.52472862) q[1];
sx q[1];
rz(-0.26291651) q[1];
sx q[1];
rz(1.252334) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9577873) q[0];
sx q[0];
rz(-2.5982214) q[0];
sx q[0];
rz(2.9068391) q[0];
x q[1];
rz(2.5400794) q[2];
sx q[2];
rz(-0.85962112) q[2];
sx q[2];
rz(-1.892923) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.5460427) q[1];
sx q[1];
rz(-0.79985207) q[1];
sx q[1];
rz(0.78592664) q[1];
rz(-pi) q[2];
x q[2];
rz(1.657293) q[3];
sx q[3];
rz(-0.75631006) q[3];
sx q[3];
rz(0.63884097) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.9423882) q[2];
sx q[2];
rz(-2.51666) q[2];
sx q[2];
rz(1.7601298) q[2];
rz(2.840461) q[3];
sx q[3];
rz(-2.1589409) q[3];
sx q[3];
rz(2.7605831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23406601) q[0];
sx q[0];
rz(-2.1387687) q[0];
sx q[0];
rz(1.7940849) q[0];
rz(0.8849591) q[1];
sx q[1];
rz(-0.62565175) q[1];
sx q[1];
rz(-1.7431097) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8475474) q[0];
sx q[0];
rz(-1.9452403) q[0];
sx q[0];
rz(3.1301296) q[0];
rz(1.3560881) q[2];
sx q[2];
rz(-2.4248059) q[2];
sx q[2];
rz(1.4625664) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.32244476) q[1];
sx q[1];
rz(-1.5335173) q[1];
sx q[1];
rz(-0.080406043) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.22866727) q[3];
sx q[3];
rz(-2.1211984) q[3];
sx q[3];
rz(0.64085863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(3.0884023) q[2];
sx q[2];
rz(-1.3206427) q[2];
sx q[2];
rz(0.94584805) q[2];
rz(-0.69019067) q[3];
sx q[3];
rz(-1.918957) q[3];
sx q[3];
rz(1.3010196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.239554) q[0];
sx q[0];
rz(-1.5382465) q[0];
sx q[0];
rz(2.4334346) q[0];
rz(0.061307727) q[1];
sx q[1];
rz(-0.61534449) q[1];
sx q[1];
rz(-1.4475488) q[1];
rz(-2.8715677) q[2];
sx q[2];
rz(-0.30120987) q[2];
sx q[2];
rz(0.78200151) q[2];
rz(3.0690774) q[3];
sx q[3];
rz(-1.391165) q[3];
sx q[3];
rz(2.388849) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
