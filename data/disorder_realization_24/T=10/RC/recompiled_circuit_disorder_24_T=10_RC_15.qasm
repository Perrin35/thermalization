OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.55968094) q[0];
sx q[0];
rz(-1.0868602) q[0];
sx q[0];
rz(-1.342919) q[0];
rz(7.8328447) q[1];
sx q[1];
rz(3.0631493) q[1];
sx q[1];
rz(6.9258239) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4896048) q[0];
sx q[0];
rz(-1.5986686) q[0];
sx q[0];
rz(0.53919381) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.3354934) q[2];
sx q[2];
rz(-0.53288424) q[2];
sx q[2];
rz(0.83836183) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.027015162) q[1];
sx q[1];
rz(-1.6903965) q[1];
sx q[1];
rz(0.12160614) q[1];
rz(-pi) q[2];
rz(-0.23407614) q[3];
sx q[3];
rz(-1.6062859) q[3];
sx q[3];
rz(1.2897829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.47444433) q[2];
sx q[2];
rz(-2.4953304) q[2];
sx q[2];
rz(-1.2791963) q[2];
rz(2.4228418) q[3];
sx q[3];
rz(-1.5703392) q[3];
sx q[3];
rz(-2.8180502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
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
rz(-0.46368018) q[0];
sx q[0];
rz(-1.7872515) q[0];
sx q[0];
rz(0.98051488) q[0];
rz(-2.9837043) q[1];
sx q[1];
rz(-1.6442559) q[1];
sx q[1];
rz(0.78871361) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6533587) q[0];
sx q[0];
rz(-1.4569067) q[0];
sx q[0];
rz(-1.4466404) q[0];
rz(-pi) q[1];
rz(-3.0208203) q[2];
sx q[2];
rz(-0.25794068) q[2];
sx q[2];
rz(1.0735807) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9849557) q[1];
sx q[1];
rz(-0.64132323) q[1];
sx q[1];
rz(2.6167469) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.73123587) q[3];
sx q[3];
rz(-0.9160708) q[3];
sx q[3];
rz(2.5996641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.366189) q[2];
sx q[2];
rz(-1.0752233) q[2];
sx q[2];
rz(-2.9555087) q[2];
rz(-0.6535334) q[3];
sx q[3];
rz(-1.3862405) q[3];
sx q[3];
rz(0.11793605) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9300951) q[0];
sx q[0];
rz(-2.1989172) q[0];
sx q[0];
rz(-0.63013664) q[0];
rz(-0.12763003) q[1];
sx q[1];
rz(-2.4928513) q[1];
sx q[1];
rz(2.4198467) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5979413) q[0];
sx q[0];
rz(-1.3184034) q[0];
sx q[0];
rz(-2.4734205) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1100936) q[2];
sx q[2];
rz(-2.2113872) q[2];
sx q[2];
rz(-3.1415423) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5210515) q[1];
sx q[1];
rz(-2.6958145) q[1];
sx q[1];
rz(3.0522703) q[1];
x q[2];
rz(0.94354043) q[3];
sx q[3];
rz(-1.2949847) q[3];
sx q[3];
rz(-0.73345473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.3815986) q[2];
sx q[2];
rz(-0.95278946) q[2];
sx q[2];
rz(0.90908137) q[2];
rz(-0.51820731) q[3];
sx q[3];
rz(-1.0639233) q[3];
sx q[3];
rz(-0.28809965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
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
rz(-0.58406126) q[0];
sx q[0];
rz(-2.4052305) q[0];
sx q[0];
rz(3.0990565) q[0];
rz(-2.361239) q[1];
sx q[1];
rz(-0.50023729) q[1];
sx q[1];
rz(0.76400486) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1799058) q[0];
sx q[0];
rz(-0.67877239) q[0];
sx q[0];
rz(1.5508482) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9742083) q[2];
sx q[2];
rz(-2.4902654) q[2];
sx q[2];
rz(-1.960388) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.3180247) q[1];
sx q[1];
rz(-1.2122224) q[1];
sx q[1];
rz(1.093868) q[1];
x q[2];
rz(-0.66391151) q[3];
sx q[3];
rz(-1.1949364) q[3];
sx q[3];
rz(1.0192878) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.2970695) q[2];
sx q[2];
rz(-1.2244747) q[2];
sx q[2];
rz(-2.6085473) q[2];
rz(2.879203) q[3];
sx q[3];
rz(-0.636594) q[3];
sx q[3];
rz(0.46245241) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9168636) q[0];
sx q[0];
rz(-2.841195) q[0];
sx q[0];
rz(-1.8977144) q[0];
rz(0.22661701) q[1];
sx q[1];
rz(-0.77426568) q[1];
sx q[1];
rz(-0.40333834) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6809083) q[0];
sx q[0];
rz(-1.7366689) q[0];
sx q[0];
rz(2.9537863) q[0];
rz(-pi) q[1];
rz(2.7062347) q[2];
sx q[2];
rz(-0.73409664) q[2];
sx q[2];
rz(1.3370607) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.12521872) q[1];
sx q[1];
rz(-0.55129904) q[1];
sx q[1];
rz(0.2972879) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0395398) q[3];
sx q[3];
rz(-0.96635339) q[3];
sx q[3];
rz(1.8464551) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.32101813) q[2];
sx q[2];
rz(-1.0514739) q[2];
sx q[2];
rz(-0.78249758) q[2];
rz(-1.1123505) q[3];
sx q[3];
rz(-2.7045043) q[3];
sx q[3];
rz(-1.0085683) q[3];
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
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.430442) q[0];
sx q[0];
rz(-0.40200457) q[0];
sx q[0];
rz(-0.45561403) q[0];
rz(0.09952155) q[1];
sx q[1];
rz(-2.0030622) q[1];
sx q[1];
rz(3.1351556) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.422056) q[0];
sx q[0];
rz(-0.7932084) q[0];
sx q[0];
rz(-0.35736812) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3077277) q[2];
sx q[2];
rz(-0.38882133) q[2];
sx q[2];
rz(-0.013465492) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.16033123) q[1];
sx q[1];
rz(-2.6604974) q[1];
sx q[1];
rz(1.0465924) q[1];
x q[2];
rz(0.75140679) q[3];
sx q[3];
rz(-0.91114984) q[3];
sx q[3];
rz(-0.71099647) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.36879888) q[2];
sx q[2];
rz(-1.5051179) q[2];
sx q[2];
rz(2.2951365) q[2];
rz(2.1438697) q[3];
sx q[3];
rz(-2.3322451) q[3];
sx q[3];
rz(-1.458582) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0434175) q[0];
sx q[0];
rz(-0.88413969) q[0];
sx q[0];
rz(-3.085882) q[0];
rz(-2.3588691) q[1];
sx q[1];
rz(-1.1306154) q[1];
sx q[1];
rz(1.1605211) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6877277) q[0];
sx q[0];
rz(-1.223432) q[0];
sx q[0];
rz(1.0930644) q[0];
rz(-pi) q[1];
rz(-2.0918526) q[2];
sx q[2];
rz(-1.9869291) q[2];
sx q[2];
rz(3.1396438) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.50642636) q[1];
sx q[1];
rz(-2.4392358) q[1];
sx q[1];
rz(-0.96794767) q[1];
rz(-2.7619744) q[3];
sx q[3];
rz(-1.9306722) q[3];
sx q[3];
rz(0.66756638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.065585) q[2];
sx q[2];
rz(-2.2202754) q[2];
sx q[2];
rz(-2.356142) q[2];
rz(-2.3857332) q[3];
sx q[3];
rz(-1.8463585) q[3];
sx q[3];
rz(0.41539645) q[3];
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
rz(-1.8287559) q[0];
sx q[0];
rz(-3.0506595) q[0];
sx q[0];
rz(-1.1428517) q[0];
rz(1.3061334) q[1];
sx q[1];
rz(-1.9530692) q[1];
sx q[1];
rz(-0.41608861) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0932255) q[0];
sx q[0];
rz(-1.0972293) q[0];
sx q[0];
rz(-0.22305365) q[0];
rz(-pi) q[1];
rz(-2.0390688) q[2];
sx q[2];
rz(-1.4434012) q[2];
sx q[2];
rz(2.1786736) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.0033274) q[1];
sx q[1];
rz(-1.4107553) q[1];
sx q[1];
rz(2.7368967) q[1];
rz(-pi) q[2];
rz(3.0307426) q[3];
sx q[3];
rz(-2.0409611) q[3];
sx q[3];
rz(-1.449031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.1064421) q[2];
sx q[2];
rz(-1.4539377) q[2];
sx q[2];
rz(0.32315928) q[2];
rz(-0.20251003) q[3];
sx q[3];
rz(-1.2812307) q[3];
sx q[3];
rz(0.57730738) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.768196) q[0];
sx q[0];
rz(-0.62921262) q[0];
sx q[0];
rz(-1.9966104) q[0];
rz(-1.1960944) q[1];
sx q[1];
rz(-0.15592608) q[1];
sx q[1];
rz(-0.51913613) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7767169) q[0];
sx q[0];
rz(-1.4380699) q[0];
sx q[0];
rz(-1.1908635) q[0];
rz(-1.8066508) q[2];
sx q[2];
rz(-1.354885) q[2];
sx q[2];
rz(-2.4049135) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.5362894) q[1];
sx q[1];
rz(-1.6603866) q[1];
sx q[1];
rz(2.3593966) q[1];
rz(-pi) q[2];
x q[2];
rz(0.2378283) q[3];
sx q[3];
rz(-1.7416218) q[3];
sx q[3];
rz(-1.2549787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.8273932) q[2];
sx q[2];
rz(-1.2131571) q[2];
sx q[2];
rz(0.040977565) q[2];
rz(2.273902) q[3];
sx q[3];
rz(-0.65338457) q[3];
sx q[3];
rz(-0.51122558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39559078) q[0];
sx q[0];
rz(-2.262291) q[0];
sx q[0];
rz(-1.5270773) q[0];
rz(-1.7136259) q[1];
sx q[1];
rz(-2.7797996) q[1];
sx q[1];
rz(-1.4987) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.98793) q[0];
sx q[0];
rz(-1.4953574) q[0];
sx q[0];
rz(3.1213785) q[0];
rz(-pi) q[1];
rz(1.8685568) q[2];
sx q[2];
rz(-0.16212633) q[2];
sx q[2];
rz(-0.33725421) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.1415256) q[1];
sx q[1];
rz(-0.86874092) q[1];
sx q[1];
rz(2.6932004) q[1];
x q[2];
rz(2.5999971) q[3];
sx q[3];
rz(-2.9062727) q[3];
sx q[3];
rz(3.1129587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.3830118) q[2];
sx q[2];
rz(-2.0643533) q[2];
sx q[2];
rz(2.995058) q[2];
rz(-0.81418973) q[3];
sx q[3];
rz(-2.345293) q[3];
sx q[3];
rz(-2.7248785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9183337) q[0];
sx q[0];
rz(-1.2467361) q[0];
sx q[0];
rz(0.99714739) q[0];
rz(1.2659484) q[1];
sx q[1];
rz(-2.1389778) q[1];
sx q[1];
rz(-1.9139342) q[1];
rz(-1.3314432) q[2];
sx q[2];
rz(-1.6806921) q[2];
sx q[2];
rz(0.42721911) q[2];
rz(2.416715) q[3];
sx q[3];
rz(-1.0387883) q[3];
sx q[3];
rz(-3.0742857) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
