OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.37031072) q[0];
sx q[0];
rz(4.1376576) q[0];
sx q[0];
rz(7.1538038) q[0];
rz(2.119996) q[1];
sx q[1];
rz(-2.8586913) q[1];
sx q[1];
rz(0.14970782) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71404845) q[0];
sx q[0];
rz(-2.2872426) q[0];
sx q[0];
rz(2.2358405) q[0];
rz(2.1120464) q[2];
sx q[2];
rz(-2.0143348) q[2];
sx q[2];
rz(2.8461547) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.183061) q[1];
sx q[1];
rz(-1.6089988) q[1];
sx q[1];
rz(1.7512291) q[1];
x q[2];
rz(1.2655067) q[3];
sx q[3];
rz(-0.57170924) q[3];
sx q[3];
rz(-1.3523462) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.25847882) q[2];
sx q[2];
rz(-0.09720619) q[2];
sx q[2];
rz(2.4374938) q[2];
rz(0.95300931) q[3];
sx q[3];
rz(-0.97218958) q[3];
sx q[3];
rz(-1.7378418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54884058) q[0];
sx q[0];
rz(-1.5177746) q[0];
sx q[0];
rz(2.5090704) q[0];
rz(-2.6951492) q[1];
sx q[1];
rz(-1.4182785) q[1];
sx q[1];
rz(0.65223637) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3632293) q[0];
sx q[0];
rz(-3.1104381) q[0];
sx q[0];
rz(2.7345783) q[0];
x q[1];
rz(0.89276887) q[2];
sx q[2];
rz(-0.36913482) q[2];
sx q[2];
rz(2.6004651) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.0374239) q[1];
sx q[1];
rz(-1.2836604) q[1];
sx q[1];
rz(1.8038521) q[1];
rz(2.2964301) q[3];
sx q[3];
rz(-0.80544986) q[3];
sx q[3];
rz(-1.0172539) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.5474881) q[2];
sx q[2];
rz(-1.0141806) q[2];
sx q[2];
rz(-1.1616421) q[2];
rz(1.9836327) q[3];
sx q[3];
rz(-2.0722814) q[3];
sx q[3];
rz(-1.6903711) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0911672) q[0];
sx q[0];
rz(-2.7624891) q[0];
sx q[0];
rz(-0.85025775) q[0];
rz(0.49750528) q[1];
sx q[1];
rz(-1.9583227) q[1];
sx q[1];
rz(1.7920378) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.772086) q[0];
sx q[0];
rz(-1.3319351) q[0];
sx q[0];
rz(2.4580965) q[0];
x q[1];
rz(-1.6367958) q[2];
sx q[2];
rz(-1.608035) q[2];
sx q[2];
rz(-2.6413692) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.459356) q[1];
sx q[1];
rz(-2.5767234) q[1];
sx q[1];
rz(0.45046803) q[1];
rz(-pi) q[2];
rz(-0.0098185929) q[3];
sx q[3];
rz(-1.9968642) q[3];
sx q[3];
rz(-2.1992418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.1094018) q[2];
sx q[2];
rz(-1.2357864) q[2];
sx q[2];
rz(0.90399495) q[2];
rz(-0.30113014) q[3];
sx q[3];
rz(-1.3826933) q[3];
sx q[3];
rz(-1.3999456) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6501453) q[0];
sx q[0];
rz(-0.97147816) q[0];
sx q[0];
rz(-1.7310463) q[0];
rz(0.63181216) q[1];
sx q[1];
rz(-1.8099064) q[1];
sx q[1];
rz(-3.1052123) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.13658) q[0];
sx q[0];
rz(-0.67626017) q[0];
sx q[0];
rz(-1.9146634) q[0];
x q[1];
rz(1.8965917) q[2];
sx q[2];
rz(-2.1278283) q[2];
sx q[2];
rz(-1.7011124) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.034678) q[1];
sx q[1];
rz(-2.2639096) q[1];
sx q[1];
rz(-2.9702859) q[1];
rz(-pi) q[2];
rz(-2.7500238) q[3];
sx q[3];
rz(-2.8885926) q[3];
sx q[3];
rz(-2.5391425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.92695421) q[2];
sx q[2];
rz(-1.4321233) q[2];
sx q[2];
rz(-1.1882163) q[2];
rz(-2.4711117) q[3];
sx q[3];
rz(-1.2256349) q[3];
sx q[3];
rz(2.5454583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4398414) q[0];
sx q[0];
rz(-1.8816467) q[0];
sx q[0];
rz(0.16648509) q[0];
rz(-2.3855551) q[1];
sx q[1];
rz(-2.2557204) q[1];
sx q[1];
rz(0.23434815) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4317961) q[0];
sx q[0];
rz(-2.095982) q[0];
sx q[0];
rz(-1.964142) q[0];
rz(-pi) q[1];
rz(-3.0460065) q[2];
sx q[2];
rz(-2.0628953) q[2];
sx q[2];
rz(1.0620067) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.4489331) q[1];
sx q[1];
rz(-1.4217136) q[1];
sx q[1];
rz(-0.74933185) q[1];
x q[2];
rz(-3.0626416) q[3];
sx q[3];
rz(-1.8936994) q[3];
sx q[3];
rz(-2.3372646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.4328737) q[2];
sx q[2];
rz(-1.9659698) q[2];
sx q[2];
rz(2.9210572) q[2];
rz(-0.43705127) q[3];
sx q[3];
rz(-2.1190937) q[3];
sx q[3];
rz(0.76550686) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7261312) q[0];
sx q[0];
rz(-0.3188062) q[0];
sx q[0];
rz(2.3244526) q[0];
rz(-0.56610402) q[1];
sx q[1];
rz(-1.348446) q[1];
sx q[1];
rz(1.1436499) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2957942) q[0];
sx q[0];
rz(-2.7685071) q[0];
sx q[0];
rz(0.11008115) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8311062) q[2];
sx q[2];
rz(-2.227265) q[2];
sx q[2];
rz(1.7973763) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.027187849) q[1];
sx q[1];
rz(-1.4804375) q[1];
sx q[1];
rz(2.0143709) q[1];
x q[2];
rz(-1.2915217) q[3];
sx q[3];
rz(-1.0696628) q[3];
sx q[3];
rz(0.53370332) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.75366655) q[2];
sx q[2];
rz(-1.5619229) q[2];
sx q[2];
rz(-0.4894408) q[2];
rz(-0.22805452) q[3];
sx q[3];
rz(-1.8830048) q[3];
sx q[3];
rz(2.7155546) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(1.56931) q[0];
sx q[0];
rz(-0.64240488) q[0];
sx q[0];
rz(1.8547159) q[0];
rz(2.4781748) q[1];
sx q[1];
rz(-1.5692915) q[1];
sx q[1];
rz(1.2333262) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39499261) q[0];
sx q[0];
rz(-1.9548423) q[0];
sx q[0];
rz(0.011944255) q[0];
rz(2.2912824) q[2];
sx q[2];
rz(-1.0734953) q[2];
sx q[2];
rz(-0.26956272) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.26112939) q[1];
sx q[1];
rz(-1.5848586) q[1];
sx q[1];
rz(-0.1111828) q[1];
rz(-pi) q[2];
rz(-0.4864278) q[3];
sx q[3];
rz(-1.2573811) q[3];
sx q[3];
rz(-3.0048971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-3.104091) q[2];
sx q[2];
rz(-0.23510322) q[2];
sx q[2];
rz(1.0160149) q[2];
rz(3.0715023) q[3];
sx q[3];
rz(-1.2066119) q[3];
sx q[3];
rz(2.0751374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
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
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12524097) q[0];
sx q[0];
rz(-0.68269435) q[0];
sx q[0];
rz(1.4455147) q[0];
rz(-0.21487543) q[1];
sx q[1];
rz(-2.3863249) q[1];
sx q[1];
rz(1.8833556) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94851516) q[0];
sx q[0];
rz(-0.45398871) q[0];
sx q[0];
rz(1.8817188) q[0];
rz(-pi) q[1];
rz(2.2082981) q[2];
sx q[2];
rz(-2.1142695) q[2];
sx q[2];
rz(0.18804929) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.82842365) q[1];
sx q[1];
rz(-2.1143882) q[1];
sx q[1];
rz(1.6924752) q[1];
rz(-2.2860252) q[3];
sx q[3];
rz(-1.1547525) q[3];
sx q[3];
rz(-0.95203979) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.68226472) q[2];
sx q[2];
rz(-0.1846281) q[2];
sx q[2];
rz(-2.5788467) q[2];
rz(2.941926) q[3];
sx q[3];
rz(-0.66185799) q[3];
sx q[3];
rz(1.1220042) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
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
rz(0.11583081) q[0];
sx q[0];
rz(-1.04302) q[0];
sx q[0];
rz(-0.2510221) q[0];
rz(0.42731467) q[1];
sx q[1];
rz(-1.2298093) q[1];
sx q[1];
rz(-3.1138611) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4902892) q[0];
sx q[0];
rz(-0.9285183) q[0];
sx q[0];
rz(0.62255967) q[0];
rz(-0.26720033) q[2];
sx q[2];
rz(-0.88951096) q[2];
sx q[2];
rz(0.53182488) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.7293538) q[1];
sx q[1];
rz(-0.60935271) q[1];
sx q[1];
rz(-2.9669697) q[1];
x q[2];
rz(3.1157007) q[3];
sx q[3];
rz(-1.504244) q[3];
sx q[3];
rz(0.47749146) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.015908265) q[2];
sx q[2];
rz(-2.1618844) q[2];
sx q[2];
rz(-1.998418) q[2];
rz(0.14287359) q[3];
sx q[3];
rz(-1.6294799) q[3];
sx q[3];
rz(-2.2843602) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7509572) q[0];
sx q[0];
rz(-1.9366783) q[0];
sx q[0];
rz(2.6877158) q[0];
rz(-0.67165309) q[1];
sx q[1];
rz(-1.6891054) q[1];
sx q[1];
rz(-0.25751105) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89850241) q[0];
sx q[0];
rz(-1.9542964) q[0];
sx q[0];
rz(2.6570508) q[0];
x q[1];
rz(-2.5970039) q[2];
sx q[2];
rz(-1.1368183) q[2];
sx q[2];
rz(1.4101654) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.9092907) q[1];
sx q[1];
rz(-1.5850987) q[1];
sx q[1];
rz(-2.4453352) q[1];
rz(-0.24210838) q[3];
sx q[3];
rz(-2.0383516) q[3];
sx q[3];
rz(-0.20413354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8132849) q[2];
sx q[2];
rz(-1.4537145) q[2];
sx q[2];
rz(-0.60662398) q[2];
rz(2.666752) q[3];
sx q[3];
rz(-2.1947221) q[3];
sx q[3];
rz(-2.301208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3474779) q[0];
sx q[0];
rz(-1.5119727) q[0];
sx q[0];
rz(2.150362) q[0];
rz(0.22656245) q[1];
sx q[1];
rz(-1.4145874) q[1];
sx q[1];
rz(0.58691595) q[1];
rz(-2.6735641) q[2];
sx q[2];
rz(-1.2266908) q[2];
sx q[2];
rz(-2.5897349) q[2];
rz(-1.0874891) q[3];
sx q[3];
rz(-1.3369505) q[3];
sx q[3];
rz(-1.7575775) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
