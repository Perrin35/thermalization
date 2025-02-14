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
rz(-0.87870413) q[0];
sx q[0];
rz(0.19600828) q[0];
rz(-1.4930383) q[1];
sx q[1];
rz(-0.9382481) q[1];
sx q[1];
rz(-2.2641163) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24841079) q[0];
sx q[0];
rz(-2.0201004) q[0];
sx q[0];
rz(-2.9025159) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2219561) q[2];
sx q[2];
rz(-1.1695122) q[2];
sx q[2];
rz(-0.87043412) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2171195) q[1];
sx q[1];
rz(-1.8397325) q[1];
sx q[1];
rz(-2.0321789) q[1];
x q[2];
rz(-0.23807293) q[3];
sx q[3];
rz(-1.1050944) q[3];
sx q[3];
rz(0.7326441) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.8087372) q[2];
sx q[2];
rz(-2.4337807) q[2];
sx q[2];
rz(3.1044002) q[2];
rz(0.91341364) q[3];
sx q[3];
rz(-2.1593058) q[3];
sx q[3];
rz(1.6325525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7971802) q[0];
sx q[0];
rz(-1.0743112) q[0];
sx q[0];
rz(-2.9194226) q[0];
rz(0.19368681) q[1];
sx q[1];
rz(-1.2682468) q[1];
sx q[1];
rz(2.8817315) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3701137) q[0];
sx q[0];
rz(-2.6278017) q[0];
sx q[0];
rz(-2.0822078) q[0];
x q[1];
rz(1.7495943) q[2];
sx q[2];
rz(-1.6769209) q[2];
sx q[2];
rz(-0.44580844) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.4603468) q[1];
sx q[1];
rz(-0.28239861) q[1];
sx q[1];
rz(0.71747924) q[1];
x q[2];
rz(-1.0845029) q[3];
sx q[3];
rz(-0.87275617) q[3];
sx q[3];
rz(1.8956888) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.3746609) q[2];
sx q[2];
rz(-2.5587475) q[2];
sx q[2];
rz(-0.15677162) q[2];
rz(2.3356656) q[3];
sx q[3];
rz(-2.0511878) q[3];
sx q[3];
rz(-2.6167615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2369775) q[0];
sx q[0];
rz(-3.0293063) q[0];
sx q[0];
rz(-1.1861381) q[0];
rz(-0.063749464) q[1];
sx q[1];
rz(-1.5260162) q[1];
sx q[1];
rz(-2.729111) q[1];
sx q[2];
rz(-pi/2) q[2];
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
rz(-0.15699129) q[2];
sx q[2];
rz(-0.89688939) q[2];
sx q[2];
rz(-1.5012596) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.3220249) q[1];
sx q[1];
rz(-1.7734999) q[1];
sx q[1];
rz(2.8611819) q[1];
rz(-pi) q[2];
rz(2.2354911) q[3];
sx q[3];
rz(-1.0933439) q[3];
sx q[3];
rz(-2.4942644) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.5731262) q[2];
sx q[2];
rz(-2.7012479) q[2];
sx q[2];
rz(1.1661412) q[2];
rz(-2.3447573) q[3];
sx q[3];
rz(-1.3692057) q[3];
sx q[3];
rz(2.3603175) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9154938) q[0];
sx q[0];
rz(-1.508536) q[0];
sx q[0];
rz(1.4133806) q[0];
rz(-0.49890292) q[1];
sx q[1];
rz(-1.7585124) q[1];
sx q[1];
rz(-1.2938719) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5099611) q[0];
sx q[0];
rz(-0.089279739) q[0];
sx q[0];
rz(-1.2651612) q[0];
rz(-pi) q[1];
rz(-1.3503374) q[2];
sx q[2];
rz(-2.8065857) q[2];
sx q[2];
rz(2.9398244) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.41999757) q[1];
sx q[1];
rz(-1.5087869) q[1];
sx q[1];
rz(1.5739784) q[1];
x q[2];
rz(-3.0741229) q[3];
sx q[3];
rz(-1.6327792) q[3];
sx q[3];
rz(-2.8648479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6125907) q[2];
sx q[2];
rz(-0.68098536) q[2];
sx q[2];
rz(2.9717818) q[2];
rz(-2.7134907) q[3];
sx q[3];
rz(-1.7120818) q[3];
sx q[3];
rz(-2.9812109) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.661628) q[0];
sx q[0];
rz(-0.39230883) q[0];
sx q[0];
rz(-0.83537927) q[0];
rz(-0.63940489) q[1];
sx q[1];
rz(-2.4765922) q[1];
sx q[1];
rz(-1.0909874) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22459659) q[0];
sx q[0];
rz(-0.69642717) q[0];
sx q[0];
rz(3.1223608) q[0];
x q[1];
rz(2.6409482) q[2];
sx q[2];
rz(-2.8029059) q[2];
sx q[2];
rz(3.088986) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.18026946) q[1];
sx q[1];
rz(-1.1774447) q[1];
sx q[1];
rz(2.250227) q[1];
rz(-1.0476607) q[3];
sx q[3];
rz(-1.8561279) q[3];
sx q[3];
rz(2.3714575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.7731446) q[2];
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
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(1.8232329) q[0];
sx q[0];
rz(-1.7802745) q[0];
sx q[0];
rz(-2.5079492) q[0];
rz(-2.7404495) q[1];
sx q[1];
rz(-0.70924962) q[1];
sx q[1];
rz(2.4493682) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5244103) q[0];
sx q[0];
rz(-1.9815191) q[0];
sx q[0];
rz(-0.80142148) q[0];
rz(-pi) q[1];
rz(-2.8234286) q[2];
sx q[2];
rz(-2.3468694) q[2];
sx q[2];
rz(-0.77407167) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.8068829) q[1];
sx q[1];
rz(-2.5100384) q[1];
sx q[1];
rz(2.0014265) q[1];
rz(-1.419846) q[3];
sx q[3];
rz(-1.387409) q[3];
sx q[3];
rz(0.91461411) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.9358518) q[2];
sx q[2];
rz(-2.3571099) q[2];
sx q[2];
rz(1.240823) q[2];
rz(1.1848909) q[3];
sx q[3];
rz(-1.840206) q[3];
sx q[3];
rz(-1.544781) q[3];
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
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.009509) q[0];
sx q[0];
rz(-1.8385831) q[0];
sx q[0];
rz(1.3707772) q[0];
rz(0.41807434) q[1];
sx q[1];
rz(-1.4825753) q[1];
sx q[1];
rz(-0.48666993) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91557568) q[0];
sx q[0];
rz(-0.56779438) q[0];
sx q[0];
rz(1.3894807) q[0];
x q[1];
rz(-0.37481793) q[2];
sx q[2];
rz(-2.0897802) q[2];
sx q[2];
rz(-1.6199552) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.44341893) q[1];
sx q[1];
rz(-1.3935745) q[1];
sx q[1];
rz(-1.366057) q[1];
rz(-1.5237996) q[3];
sx q[3];
rz(-1.1011657) q[3];
sx q[3];
rz(2.3835973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.70249867) q[2];
sx q[2];
rz(-2.0926026) q[2];
sx q[2];
rz(2.9443963) q[2];
rz(-0.50944734) q[3];
sx q[3];
rz(-2.8809437) q[3];
sx q[3];
rz(0.82836866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
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
rz(-1.2800804) q[0];
sx q[0];
rz(-1.0741638) q[0];
sx q[0];
rz(1.7751088) q[0];
rz(-2.3249783) q[1];
sx q[1];
rz(-2.0520703) q[1];
sx q[1];
rz(-2.8588967) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1408843) q[0];
sx q[0];
rz(-1.3907897) q[0];
sx q[0];
rz(-2.1732844) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6655695) q[2];
sx q[2];
rz(-1.2605209) q[2];
sx q[2];
rz(2.5233334) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.7636488) q[1];
sx q[1];
rz(-0.70893439) q[1];
sx q[1];
rz(-1.2825185) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9603086) q[3];
sx q[3];
rz(-1.0822902) q[3];
sx q[3];
rz(-0.28340411) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.035159811) q[2];
sx q[2];
rz(-1.0059493) q[2];
sx q[2];
rz(-2.1584568) q[2];
rz(1.2784917) q[3];
sx q[3];
rz(-0.92517868) q[3];
sx q[3];
rz(-0.92888752) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66016692) q[0];
sx q[0];
rz(-1.7919414) q[0];
sx q[0];
rz(-2.8283258) q[0];
rz(-2.616864) q[1];
sx q[1];
rz(-0.26291651) q[1];
sx q[1];
rz(-1.252334) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6853297) q[0];
sx q[0];
rz(-1.043911) q[0];
sx q[0];
rz(-1.4312126) q[0];
x q[1];
rz(0.76340126) q[2];
sx q[2];
rz(-1.127725) q[2];
sx q[2];
rz(3.0425261) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5572667) q[1];
sx q[1];
rz(-1.0385822) q[1];
sx q[1];
rz(0.6289215) q[1];
rz(-pi) q[2];
rz(-3.0602686) q[3];
sx q[3];
rz(-0.81800753) q[3];
sx q[3];
rz(2.6214056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.9423882) q[2];
sx q[2];
rz(-0.62493268) q[2];
sx q[2];
rz(1.7601298) q[2];
rz(-2.840461) q[3];
sx q[3];
rz(-2.1589409) q[3];
sx q[3];
rz(0.38100955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9075266) q[0];
sx q[0];
rz(-2.1387687) q[0];
sx q[0];
rz(1.3475077) q[0];
rz(2.2566336) q[1];
sx q[1];
rz(-2.5159409) q[1];
sx q[1];
rz(1.398483) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32537726) q[0];
sx q[0];
rz(-0.37461108) q[0];
sx q[0];
rz(1.5999567) q[0];
x q[1];
rz(-0.18357205) q[2];
sx q[2];
rz(-2.2677448) q[2];
sx q[2];
rz(-1.1810034) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.2513548) q[1];
sx q[1];
rz(-1.4904463) q[1];
sx q[1];
rz(1.608196) q[1];
rz(-pi) q[2];
rz(-0.22866727) q[3];
sx q[3];
rz(-2.1211984) q[3];
sx q[3];
rz(-2.500734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.053190319) q[2];
sx q[2];
rz(-1.8209499) q[2];
sx q[2];
rz(2.1957446) q[2];
rz(-0.69019067) q[3];
sx q[3];
rz(-1.2226356) q[3];
sx q[3];
rz(-1.3010196) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.239554) q[0];
sx q[0];
rz(-1.5382465) q[0];
sx q[0];
rz(2.4334346) q[0];
rz(-0.061307727) q[1];
sx q[1];
rz(-2.5262482) q[1];
sx q[1];
rz(1.6940438) q[1];
rz(0.270025) q[2];
sx q[2];
rz(-0.30120987) q[2];
sx q[2];
rz(0.78200151) q[2];
rz(-1.1911627) q[3];
sx q[3];
rz(-0.1935696) q[3];
sx q[3];
rz(-0.36657666) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
