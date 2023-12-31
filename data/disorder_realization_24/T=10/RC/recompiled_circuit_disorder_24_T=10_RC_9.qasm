OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.5819117) q[0];
sx q[0];
rz(-2.0547325) q[0];
sx q[0];
rz(1.342919) q[0];
rz(7.8328447) q[1];
sx q[1];
rz(3.0631493) q[1];
sx q[1];
rz(6.9258239) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6519878) q[0];
sx q[0];
rz(-1.5986686) q[0];
sx q[0];
rz(-0.53919381) q[0];
x q[1];
rz(-1.7625916) q[2];
sx q[2];
rz(-1.0704874) q[2];
sx q[2];
rz(1.9185916) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.5832311) q[1];
sx q[1];
rz(-1.4500631) q[1];
sx q[1];
rz(-1.6912778) q[1];
x q[2];
rz(-1.60728) q[3];
sx q[3];
rz(-1.3368703) q[3];
sx q[3];
rz(2.8521188) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.47444433) q[2];
sx q[2];
rz(-0.64626226) q[2];
sx q[2];
rz(1.2791963) q[2];
rz(2.4228418) q[3];
sx q[3];
rz(-1.5712534) q[3];
sx q[3];
rz(-0.32354245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(0.46368018) q[0];
sx q[0];
rz(-1.7872515) q[0];
sx q[0];
rz(-0.98051488) q[0];
rz(0.15788831) q[1];
sx q[1];
rz(-1.4973367) q[1];
sx q[1];
rz(2.352879) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.096743874) q[0];
sx q[0];
rz(-1.4474488) q[0];
sx q[0];
rz(-0.11476536) q[0];
x q[1];
rz(3.0208203) q[2];
sx q[2];
rz(-2.883652) q[2];
sx q[2];
rz(-2.068012) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.6106818) q[1];
sx q[1];
rz(-1.0265988) q[1];
sx q[1];
rz(-1.2128085) q[1];
rz(-pi) q[2];
rz(2.2867145) q[3];
sx q[3];
rz(-2.2023871) q[3];
sx q[3];
rz(-0.43254334) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.366189) q[2];
sx q[2];
rz(-2.0663694) q[2];
sx q[2];
rz(2.9555087) q[2];
rz(-2.4880593) q[3];
sx q[3];
rz(-1.7553522) q[3];
sx q[3];
rz(-3.0236566) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9300951) q[0];
sx q[0];
rz(-0.94267541) q[0];
sx q[0];
rz(-0.63013664) q[0];
rz(3.0139626) q[1];
sx q[1];
rz(-2.4928513) q[1];
sx q[1];
rz(2.4198467) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9741309) q[0];
sx q[0];
rz(-0.92739096) q[0];
sx q[0];
rz(-1.8882303) q[0];
rz(-pi) q[1];
rz(-0.69408728) q[2];
sx q[2];
rz(-1.9352479) q[2];
sx q[2];
rz(1.8592161) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.5210515) q[1];
sx q[1];
rz(-2.6958145) q[1];
sx q[1];
rz(-3.0522703) q[1];
x q[2];
rz(-2.1980522) q[3];
sx q[3];
rz(-1.846608) q[3];
sx q[3];
rz(-2.4081379) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.3815986) q[2];
sx q[2];
rz(-2.1888032) q[2];
sx q[2];
rz(-2.2325113) q[2];
rz(0.51820731) q[3];
sx q[3];
rz(-1.0639233) q[3];
sx q[3];
rz(0.28809965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5575314) q[0];
sx q[0];
rz(-0.73636213) q[0];
sx q[0];
rz(0.042536143) q[0];
rz(0.78035367) q[1];
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
rz(-2.9873136) q[0];
sx q[0];
rz(-0.89218441) q[0];
sx q[0];
rz(-0.016088386) q[0];
x q[1];
rz(1.4444703) q[2];
sx q[2];
rz(-2.2115123) q[2];
sx q[2];
rz(1.390552) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.15061513) q[1];
sx q[1];
rz(-0.58826485) q[1];
sx q[1];
rz(0.88612835) q[1];
rz(-2.0352827) q[3];
sx q[3];
rz(-0.96040695) q[3];
sx q[3];
rz(0.27184091) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.8445231) q[2];
sx q[2];
rz(-1.917118) q[2];
sx q[2];
rz(-0.53304535) q[2];
rz(-2.879203) q[3];
sx q[3];
rz(-2.5049987) q[3];
sx q[3];
rz(-2.6791402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2247291) q[0];
sx q[0];
rz(-2.841195) q[0];
sx q[0];
rz(-1.2438783) q[0];
rz(0.22661701) q[1];
sx q[1];
rz(-0.77426568) q[1];
sx q[1];
rz(-0.40333834) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.078743155) q[0];
sx q[0];
rz(-1.3855977) q[0];
sx q[0];
rz(-1.4020105) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.68571217) q[2];
sx q[2];
rz(-1.2843686) q[2];
sx q[2];
rz(0.098766947) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.9215645) q[1];
sx q[1];
rz(-2.0953396) q[1];
sx q[1];
rz(-1.7490053) q[1];
x q[2];
rz(0.63286085) q[3];
sx q[3];
rz(-0.78213464) q[3];
sx q[3];
rz(-0.49304214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8205745) q[2];
sx q[2];
rz(-2.0901188) q[2];
sx q[2];
rz(-2.3590951) q[2];
rz(-1.1123505) q[3];
sx q[3];
rz(-2.7045043) q[3];
sx q[3];
rz(2.1330244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7111506) q[0];
sx q[0];
rz(-0.40200457) q[0];
sx q[0];
rz(0.45561403) q[0];
rz(3.0420711) q[1];
sx q[1];
rz(-2.0030622) q[1];
sx q[1];
rz(-3.1351556) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71953668) q[0];
sx q[0];
rz(-0.7932084) q[0];
sx q[0];
rz(2.7842245) q[0];
rz(-pi) q[1];
rz(1.3077277) q[2];
sx q[2];
rz(-0.38882133) q[2];
sx q[2];
rz(-3.1281272) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.2575063) q[1];
sx q[1];
rz(-1.8045366) q[1];
sx q[1];
rz(-1.1463548) q[1];
rz(-pi) q[2];
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
rz(-2.7727938) q[2];
sx q[2];
rz(-1.6364748) q[2];
sx q[2];
rz(-0.84645611) q[2];
rz(2.1438697) q[3];
sx q[3];
rz(-0.80934757) q[3];
sx q[3];
rz(1.458582) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.098175123) q[0];
sx q[0];
rz(-2.257453) q[0];
sx q[0];
rz(-3.085882) q[0];
rz(2.3588691) q[1];
sx q[1];
rz(-1.1306154) q[1];
sx q[1];
rz(-1.1605211) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6768778) q[0];
sx q[0];
rz(-2.5589295) q[0];
sx q[0];
rz(-2.23784) q[0];
rz(-pi) q[1];
rz(-0.84476446) q[2];
sx q[2];
rz(-2.4869707) q[2];
sx q[2];
rz(0.95915937) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.5481422) q[1];
sx q[1];
rz(-1.9458117) q[1];
sx q[1];
rz(-0.96199357) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.79302391) q[3];
sx q[3];
rz(-2.6245955) q[3];
sx q[3];
rz(-1.5152064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.065585) q[2];
sx q[2];
rz(-0.92131725) q[2];
sx q[2];
rz(0.78545061) q[2];
rz(-0.75585946) q[3];
sx q[3];
rz(-1.8463585) q[3];
sx q[3];
rz(2.7261962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8287559) q[0];
sx q[0];
rz(-0.090933196) q[0];
sx q[0];
rz(-1.1428517) q[0];
rz(1.8354592) q[1];
sx q[1];
rz(-1.1885234) q[1];
sx q[1];
rz(2.725504) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0932255) q[0];
sx q[0];
rz(-2.0443633) q[0];
sx q[0];
rz(0.22305365) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.14256723) q[2];
sx q[2];
rz(-1.1066184) q[2];
sx q[2];
rz(2.4695421) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.0033274) q[1];
sx q[1];
rz(-1.7308373) q[1];
sx q[1];
rz(0.40469594) q[1];
rz(-pi) q[2];
x q[2];
rz(0.11085005) q[3];
sx q[3];
rz(-2.0409611) q[3];
sx q[3];
rz(1.449031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.1064421) q[2];
sx q[2];
rz(-1.687655) q[2];
sx q[2];
rz(-2.8184334) q[2];
rz(-2.9390826) q[3];
sx q[3];
rz(-1.860362) q[3];
sx q[3];
rz(-2.5642853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.768196) q[0];
sx q[0];
rz(-0.62921262) q[0];
sx q[0];
rz(1.9966104) q[0];
rz(-1.9454983) q[1];
sx q[1];
rz(-0.15592608) q[1];
sx q[1];
rz(-2.6224565) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7767169) q[0];
sx q[0];
rz(-1.7035228) q[0];
sx q[0];
rz(-1.9507292) q[0];
rz(-pi) q[1];
rz(0.22186188) q[2];
sx q[2];
rz(-1.3405181) q[2];
sx q[2];
rz(2.3589173) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.12431006) q[1];
sx q[1];
rz(-0.78622422) q[1];
sx q[1];
rz(3.0148274) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9037644) q[3];
sx q[3];
rz(-1.7416218) q[3];
sx q[3];
rz(-1.886614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.3141994) q[2];
sx q[2];
rz(-1.9284356) q[2];
sx q[2];
rz(-3.1006151) q[2];
rz(-2.273902) q[3];
sx q[3];
rz(-0.65338457) q[3];
sx q[3];
rz(-2.6303671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39559078) q[0];
sx q[0];
rz(-0.87930167) q[0];
sx q[0];
rz(1.5270773) q[0];
rz(-1.7136259) q[1];
sx q[1];
rz(-0.36179301) q[1];
sx q[1];
rz(1.4987) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7258543) q[0];
sx q[0];
rz(-3.0634974) q[0];
sx q[0];
rz(-1.8321091) q[0];
rz(-pi) q[1];
rz(-1.8685568) q[2];
sx q[2];
rz(-0.16212633) q[2];
sx q[2];
rz(0.33725421) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.1415256) q[1];
sx q[1];
rz(-0.86874092) q[1];
sx q[1];
rz(-0.4483923) q[1];
x q[2];
rz(0.54159553) q[3];
sx q[3];
rz(-2.9062727) q[3];
sx q[3];
rz(-3.1129587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.3830118) q[2];
sx q[2];
rz(-1.0772394) q[2];
sx q[2];
rz(0.14653462) q[2];
rz(-2.3274029) q[3];
sx q[3];
rz(-2.345293) q[3];
sx q[3];
rz(-0.41671419) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
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
rz(1.3314432) q[2];
sx q[2];
rz(-1.4609006) q[2];
sx q[2];
rz(-2.7143735) q[2];
rz(-0.90445789) q[3];
sx q[3];
rz(-0.96257985) q[3];
sx q[3];
rz(-1.0812159) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
