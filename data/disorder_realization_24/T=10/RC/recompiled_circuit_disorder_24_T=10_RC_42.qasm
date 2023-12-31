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
rz(-1.5919332) q[1];
sx q[1];
rz(-3.0631493) q[1];
sx q[1];
rz(0.6426386) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6519878) q[0];
sx q[0];
rz(-1.5986686) q[0];
sx q[0];
rz(-0.53919381) q[0];
x q[1];
rz(1.7625916) q[2];
sx q[2];
rz(-2.0711053) q[2];
sx q[2];
rz(1.9185916) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.5583615) q[1];
sx q[1];
rz(-1.4500631) q[1];
sx q[1];
rz(-1.6912778) q[1];
x q[2];
rz(-1.5343127) q[3];
sx q[3];
rz(-1.3368703) q[3];
sx q[3];
rz(-2.8521188) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.6671483) q[2];
sx q[2];
rz(-0.64626226) q[2];
sx q[2];
rz(1.8623964) q[2];
rz(0.71875087) q[3];
sx q[3];
rz(-1.5703392) q[3];
sx q[3];
rz(-0.32354245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46368018) q[0];
sx q[0];
rz(-1.3543411) q[0];
sx q[0];
rz(2.1610778) q[0];
rz(-2.9837043) q[1];
sx q[1];
rz(-1.6442559) q[1];
sx q[1];
rz(-2.352879) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4855027) q[0];
sx q[0];
rz(-2.9733109) q[0];
sx q[0];
rz(2.3165354) q[0];
rz(-pi) q[1];
rz(-0.12077232) q[2];
sx q[2];
rz(-2.883652) q[2];
sx q[2];
rz(-2.068012) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.5309108) q[1];
sx q[1];
rz(-1.0265988) q[1];
sx q[1];
rz(1.9287841) q[1];
rz(0.85487811) q[3];
sx q[3];
rz(-2.2023871) q[3];
sx q[3];
rz(-2.7090493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.7754037) q[2];
sx q[2];
rz(-2.0663694) q[2];
sx q[2];
rz(0.186084) q[2];
rz(2.4880593) q[3];
sx q[3];
rz(-1.3862405) q[3];
sx q[3];
rz(-3.0236566) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2114975) q[0];
sx q[0];
rz(-2.1989172) q[0];
sx q[0];
rz(-0.63013664) q[0];
rz(-0.12763003) q[1];
sx q[1];
rz(-0.6487414) q[1];
sx q[1];
rz(0.72174597) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8080224) q[0];
sx q[0];
rz(-2.434242) q[0];
sx q[0];
rz(-2.7471514) q[0];
x q[1];
rz(1.1100936) q[2];
sx q[2];
rz(-0.93020541) q[2];
sx q[2];
rz(5.0355807e-05) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.6205412) q[1];
sx q[1];
rz(-2.6958145) q[1];
sx q[1];
rz(0.089322395) q[1];
x q[2];
rz(2.0201163) q[3];
sx q[3];
rz(-0.67766261) q[3];
sx q[3];
rz(2.6636056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.7599941) q[2];
sx q[2];
rz(-0.95278946) q[2];
sx q[2];
rz(-0.90908137) q[2];
rz(-2.6233853) q[3];
sx q[3];
rz(-2.0776694) q[3];
sx q[3];
rz(2.853493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(2.5575314) q[0];
sx q[0];
rz(-2.4052305) q[0];
sx q[0];
rz(3.0990565) q[0];
rz(-0.78035367) q[1];
sx q[1];
rz(-2.6413554) q[1];
sx q[1];
rz(0.76400486) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9616868) q[0];
sx q[0];
rz(-2.4628203) q[0];
sx q[0];
rz(-1.5508482) q[0];
rz(-0.1673844) q[2];
sx q[2];
rz(-2.4902654) q[2];
sx q[2];
rz(-1.1812047) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.2149787) q[1];
sx q[1];
rz(-1.1264631) q[1];
sx q[1];
rz(-2.742393) q[1];
rz(-pi) q[2];
x q[2];
rz(0.66391151) q[3];
sx q[3];
rz(-1.1949364) q[3];
sx q[3];
rz(2.1223048) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.2970695) q[2];
sx q[2];
rz(-1.917118) q[2];
sx q[2];
rz(2.6085473) q[2];
rz(0.26238966) q[3];
sx q[3];
rz(-2.5049987) q[3];
sx q[3];
rz(0.46245241) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
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
rz(-1.9168636) q[0];
sx q[0];
rz(-2.841195) q[0];
sx q[0];
rz(1.8977144) q[0];
rz(-0.22661701) q[1];
sx q[1];
rz(-0.77426568) q[1];
sx q[1];
rz(-2.7382543) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0628495) q[0];
sx q[0];
rz(-1.7559949) q[0];
sx q[0];
rz(1.7395822) q[0];
rz(-pi) q[1];
rz(0.68571217) q[2];
sx q[2];
rz(-1.8572241) q[2];
sx q[2];
rz(-3.0428257) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(3.0163739) q[1];
sx q[1];
rz(-2.5902936) q[1];
sx q[1];
rz(0.2972879) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0395398) q[3];
sx q[3];
rz(-2.1752393) q[3];
sx q[3];
rz(-1.8464551) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.8205745) q[2];
sx q[2];
rz(-1.0514739) q[2];
sx q[2];
rz(0.78249758) q[2];
rz(-2.0292422) q[3];
sx q[3];
rz(-2.7045043) q[3];
sx q[3];
rz(1.0085683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7111506) q[0];
sx q[0];
rz(-2.7395881) q[0];
sx q[0];
rz(0.45561403) q[0];
rz(0.09952155) q[1];
sx q[1];
rz(-2.0030622) q[1];
sx q[1];
rz(3.1351556) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59505263) q[0];
sx q[0];
rz(-1.8227302) q[0];
sx q[0];
rz(2.3810054) q[0];
x q[1];
rz(1.8338649) q[2];
sx q[2];
rz(-2.7527713) q[2];
sx q[2];
rz(0.013465492) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.9812614) q[1];
sx q[1];
rz(-0.48109522) q[1];
sx q[1];
rz(-2.0950003) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3859343) q[3];
sx q[3];
rz(-1.0010127) q[3];
sx q[3];
rz(2.8017686) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.7727938) q[2];
sx q[2];
rz(-1.6364748) q[2];
sx q[2];
rz(0.84645611) q[2];
rz(0.99772292) q[3];
sx q[3];
rz(-0.80934757) q[3];
sx q[3];
rz(-1.458582) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
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
rz(-1.9810716) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6768778) q[0];
sx q[0];
rz(-2.5589295) q[0];
sx q[0];
rz(2.23784) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0497401) q[2];
sx q[2];
rz(-1.9869291) q[2];
sx q[2];
rz(-3.1396438) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.50642636) q[1];
sx q[1];
rz(-0.70235683) q[1];
sx q[1];
rz(0.96794767) q[1];
x q[2];
rz(-1.9556932) q[3];
sx q[3];
rz(-1.2166096) q[3];
sx q[3];
rz(-2.3779496) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.065585) q[2];
sx q[2];
rz(-0.92131725) q[2];
sx q[2];
rz(0.78545061) q[2];
rz(2.3857332) q[3];
sx q[3];
rz(-1.2952341) q[3];
sx q[3];
rz(-2.7261962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3128368) q[0];
sx q[0];
rz(-0.090933196) q[0];
sx q[0];
rz(-1.998741) q[0];
rz(1.8354592) q[1];
sx q[1];
rz(-1.1885234) q[1];
sx q[1];
rz(2.725504) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0932255) q[0];
sx q[0];
rz(-1.0972293) q[0];
sx q[0];
rz(0.22305365) q[0];
x q[1];
rz(2.0390688) q[2];
sx q[2];
rz(-1.4434012) q[2];
sx q[2];
rz(0.96291908) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.5059698) q[1];
sx q[1];
rz(-1.9700248) q[1];
sx q[1];
rz(1.3969621) q[1];
rz(-pi) q[2];
rz(-2.0434521) q[3];
sx q[3];
rz(-1.4720159) q[3];
sx q[3];
rz(0.071382513) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.0351506) q[2];
sx q[2];
rz(-1.4539377) q[2];
sx q[2];
rz(-0.32315928) q[2];
rz(2.9390826) q[3];
sx q[3];
rz(-1.860362) q[3];
sx q[3];
rz(2.5642853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37339661) q[0];
sx q[0];
rz(-2.51238) q[0];
sx q[0];
rz(-1.9966104) q[0];
rz(-1.9454983) q[1];
sx q[1];
rz(-2.9856666) q[1];
sx q[1];
rz(2.6224565) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11408344) q[0];
sx q[0];
rz(-2.7402096) q[0];
sx q[0];
rz(1.2252349) q[0];
rz(-pi) q[1];
rz(1.8066508) q[2];
sx q[2];
rz(-1.7867076) q[2];
sx q[2];
rz(0.73667919) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.0874333) q[1];
sx q[1];
rz(-2.3490153) q[1];
sx q[1];
rz(1.6967609) q[1];
rz(-pi) q[2];
rz(1.395123) q[3];
sx q[3];
rz(-1.3364949) q[3];
sx q[3];
rz(-0.3570041) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.3141994) q[2];
sx q[2];
rz(-1.9284356) q[2];
sx q[2];
rz(3.1006151) q[2];
rz(0.86769062) q[3];
sx q[3];
rz(-0.65338457) q[3];
sx q[3];
rz(0.51122558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39559078) q[0];
sx q[0];
rz(-2.262291) q[0];
sx q[0];
rz(-1.6145153) q[0];
rz(-1.4279667) q[1];
sx q[1];
rz(-2.7797996) q[1];
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
x q[1];
rz(3.0936436) q[2];
sx q[2];
rz(-1.7257294) q[2];
sx q[2];
rz(-2.5028554) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1415256) q[1];
sx q[1];
rz(-0.86874092) q[1];
sx q[1];
rz(0.4483923) q[1];
rz(-pi) q[2];
rz(0.54159553) q[3];
sx q[3];
rz(-0.23532) q[3];
sx q[3];
rz(3.1129587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.75858086) q[2];
sx q[2];
rz(-1.0772394) q[2];
sx q[2];
rz(-2.995058) q[2];
rz(2.3274029) q[3];
sx q[3];
rz(-2.345293) q[3];
sx q[3];
rz(-2.7248785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9183337) q[0];
sx q[0];
rz(-1.8948566) q[0];
sx q[0];
rz(-2.1444453) q[0];
rz(-1.2659484) q[1];
sx q[1];
rz(-1.0026149) q[1];
sx q[1];
rz(1.2276585) q[1];
rz(-2.0064034) q[2];
sx q[2];
rz(-0.26293593) q[2];
sx q[2];
rz(1.5756366) q[2];
rz(-0.72600611) q[3];
sx q[3];
rz(-0.86961679) q[3];
sx q[3];
rz(1.1179954) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
