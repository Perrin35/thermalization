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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6519878) q[0];
sx q[0];
rz(-1.5986686) q[0];
sx q[0];
rz(-0.53919381) q[0];
rz(0.3354934) q[2];
sx q[2];
rz(-0.53288424) q[2];
sx q[2];
rz(2.3032308) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.027015162) q[1];
sx q[1];
rz(-1.6903965) q[1];
sx q[1];
rz(-3.0199865) q[1];
rz(-pi) q[2];
x q[2];
rz(0.1518941) q[3];
sx q[3];
rz(-0.23670247) q[3];
sx q[3];
rz(3.0083002) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.6671483) q[2];
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
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6779125) q[0];
sx q[0];
rz(-1.3543411) q[0];
sx q[0];
rz(0.98051488) q[0];
rz(0.15788831) q[1];
sx q[1];
rz(-1.6442559) q[1];
sx q[1];
rz(0.78871361) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.096743874) q[0];
sx q[0];
rz(-1.4474488) q[0];
sx q[0];
rz(3.0268273) q[0];
x q[1];
rz(1.5390225) q[2];
sx q[2];
rz(-1.3147768) q[2];
sx q[2];
rz(-2.1928744) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.6106818) q[1];
sx q[1];
rz(-1.0265988) q[1];
sx q[1];
rz(1.2128085) q[1];
x q[2];
rz(-2.3716281) q[3];
sx q[3];
rz(-1.0125481) q[3];
sx q[3];
rz(1.6127197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.7754037) q[2];
sx q[2];
rz(-2.0663694) q[2];
sx q[2];
rz(0.186084) q[2];
rz(-2.4880593) q[3];
sx q[3];
rz(-1.3862405) q[3];
sx q[3];
rz(3.0236566) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9300951) q[0];
sx q[0];
rz(-0.94267541) q[0];
sx q[0];
rz(2.511456) q[0];
rz(-0.12763003) q[1];
sx q[1];
rz(-0.6487414) q[1];
sx q[1];
rz(-2.4198467) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5979413) q[0];
sx q[0];
rz(-1.8231892) q[0];
sx q[0];
rz(0.6681722) q[0];
rz(-2.6038405) q[2];
sx q[2];
rz(-0.76965145) q[2];
sx q[2];
rz(0.69307454) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.1303781) q[1];
sx q[1];
rz(-1.5323258) q[1];
sx q[1];
rz(-2.6973666) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0201163) q[3];
sx q[3];
rz(-2.46393) q[3];
sx q[3];
rz(0.47798702) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.3815986) q[2];
sx q[2];
rz(-2.1888032) q[2];
sx q[2];
rz(-0.90908137) q[2];
rz(2.6233853) q[3];
sx q[3];
rz(-1.0639233) q[3];
sx q[3];
rz(-0.28809965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58406126) q[0];
sx q[0];
rz(-2.4052305) q[0];
sx q[0];
rz(-3.0990565) q[0];
rz(-2.361239) q[1];
sx q[1];
rz(-2.6413554) q[1];
sx q[1];
rz(-0.76400486) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4064179) q[0];
sx q[0];
rz(-1.58332) q[0];
sx q[0];
rz(2.2494715) q[0];
x q[1];
rz(0.1673844) q[2];
sx q[2];
rz(-0.65132729) q[2];
sx q[2];
rz(-1.1812047) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.82356794) q[1];
sx q[1];
rz(-1.9293702) q[1];
sx q[1];
rz(-1.093868) q[1];
x q[2];
rz(-2.571991) q[3];
sx q[3];
rz(-0.74865018) q[3];
sx q[3];
rz(0.99018712) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.2970695) q[2];
sx q[2];
rz(-1.917118) q[2];
sx q[2];
rz(2.6085473) q[2];
rz(0.26238966) q[3];
sx q[3];
rz(-0.636594) q[3];
sx q[3];
rz(2.6791402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
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
rz(1.9168636) q[0];
sx q[0];
rz(-0.30039766) q[0];
sx q[0];
rz(-1.2438783) q[0];
rz(2.9149756) q[1];
sx q[1];
rz(-0.77426568) q[1];
sx q[1];
rz(-2.7382543) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4606844) q[0];
sx q[0];
rz(-1.7366689) q[0];
sx q[0];
rz(-0.18780639) q[0];
rz(-1.9344159) q[2];
sx q[2];
rz(-0.91797963) q[2];
sx q[2];
rz(1.2448685) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.4407318) q[1];
sx q[1];
rz(-1.7248389) q[1];
sx q[1];
rz(-0.53149077) q[1];
x q[2];
rz(-2.5087318) q[3];
sx q[3];
rz(-0.78213464) q[3];
sx q[3];
rz(2.6485505) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.32101813) q[2];
sx q[2];
rz(-2.0901188) q[2];
sx q[2];
rz(-0.78249758) q[2];
rz(2.0292422) q[3];
sx q[3];
rz(-2.7045043) q[3];
sx q[3];
rz(-1.0085683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7111506) q[0];
sx q[0];
rz(-0.40200457) q[0];
sx q[0];
rz(-0.45561403) q[0];
rz(3.0420711) q[1];
sx q[1];
rz(-2.0030622) q[1];
sx q[1];
rz(0.0064370357) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59505263) q[0];
sx q[0];
rz(-1.8227302) q[0];
sx q[0];
rz(2.3810054) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9474895) q[2];
sx q[2];
rz(-1.4720535) q[2];
sx q[2];
rz(1.313098) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.8840864) q[1];
sx q[1];
rz(-1.3370561) q[1];
sx q[1];
rz(-1.9952378) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3859343) q[3];
sx q[3];
rz(-2.14058) q[3];
sx q[3];
rz(2.8017686) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.36879888) q[2];
sx q[2];
rz(-1.6364748) q[2];
sx q[2];
rz(2.2951365) q[2];
rz(-2.1438697) q[3];
sx q[3];
rz(-0.80934757) q[3];
sx q[3];
rz(1.6830106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.098175123) q[0];
sx q[0];
rz(-0.88413969) q[0];
sx q[0];
rz(3.085882) q[0];
rz(0.78272351) q[1];
sx q[1];
rz(-1.1306154) q[1];
sx q[1];
rz(1.1605211) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6768778) q[0];
sx q[0];
rz(-2.5589295) q[0];
sx q[0];
rz(-2.23784) q[0];
rz(-0.84476446) q[2];
sx q[2];
rz(-0.65462199) q[2];
sx q[2];
rz(-0.95915937) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.9142368) q[1];
sx q[1];
rz(-2.1319234) q[1];
sx q[1];
rz(-0.44740541) q[1];
x q[2];
rz(-1.9556932) q[3];
sx q[3];
rz(-1.2166096) q[3];
sx q[3];
rz(0.76364309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.065585) q[2];
sx q[2];
rz(-0.92131725) q[2];
sx q[2];
rz(-0.78545061) q[2];
rz(-0.75585946) q[3];
sx q[3];
rz(-1.8463585) q[3];
sx q[3];
rz(-0.41539645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8287559) q[0];
sx q[0];
rz(-0.090933196) q[0];
sx q[0];
rz(-1.998741) q[0];
rz(-1.3061334) q[1];
sx q[1];
rz(-1.1885234) q[1];
sx q[1];
rz(2.725504) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5160822) q[0];
sx q[0];
rz(-1.372638) q[0];
sx q[0];
rz(1.0869736) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1025238) q[2];
sx q[2];
rz(-1.4434012) q[2];
sx q[2];
rz(2.1786736) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.5059698) q[1];
sx q[1];
rz(-1.9700248) q[1];
sx q[1];
rz(-1.3969621) q[1];
x q[2];
rz(-3.0307426) q[3];
sx q[3];
rz(-1.1006315) q[3];
sx q[3];
rz(1.6925616) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1064421) q[2];
sx q[2];
rz(-1.687655) q[2];
sx q[2];
rz(-2.8184334) q[2];
rz(-0.20251003) q[3];
sx q[3];
rz(-1.2812307) q[3];
sx q[3];
rz(-2.5642853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.768196) q[0];
sx q[0];
rz(-2.51238) q[0];
sx q[0];
rz(1.1449822) q[0];
rz(-1.9454983) q[1];
sx q[1];
rz(-0.15592608) q[1];
sx q[1];
rz(-2.6224565) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8828744) q[0];
sx q[0];
rz(-1.9472194) q[0];
sx q[0];
rz(-0.14278485) q[0];
rz(-pi) q[1];
rz(2.3245415) q[2];
sx q[2];
rz(-2.8231986) q[2];
sx q[2];
rz(1.5794093) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.054159315) q[1];
sx q[1];
rz(-0.79257733) q[1];
sx q[1];
rz(-1.4448318) q[1];
rz(2.5095652) q[3];
sx q[3];
rz(-0.29187376) q[3];
sx q[3];
rz(-0.29575086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.8273932) q[2];
sx q[2];
rz(-1.2131571) q[2];
sx q[2];
rz(-0.040977565) q[2];
rz(-0.86769062) q[3];
sx q[3];
rz(-0.65338457) q[3];
sx q[3];
rz(-0.51122558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7460019) q[0];
sx q[0];
rz(-0.87930167) q[0];
sx q[0];
rz(1.5270773) q[0];
rz(-1.7136259) q[1];
sx q[1];
rz(-0.36179301) q[1];
sx q[1];
rz(1.4987) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1536627) q[0];
sx q[0];
rz(-1.6462353) q[0];
sx q[0];
rz(0.020214202) q[0];
rz(-pi) q[1];
rz(0.047949009) q[2];
sx q[2];
rz(-1.4158632) q[2];
sx q[2];
rz(0.63873728) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.26951075) q[1];
sx q[1];
rz(-1.2334358) q[1];
sx q[1];
rz(-0.81706394) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.54159553) q[3];
sx q[3];
rz(-2.9062727) q[3];
sx q[3];
rz(3.1129587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.3830118) q[2];
sx q[2];
rz(-2.0643533) q[2];
sx q[2];
rz(0.14653462) q[2];
rz(2.3274029) q[3];
sx q[3];
rz(-2.345293) q[3];
sx q[3];
rz(-2.7248785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2232589) q[0];
sx q[0];
rz(-1.8948566) q[0];
sx q[0];
rz(-2.1444453) q[0];
rz(-1.8756443) q[1];
sx q[1];
rz(-2.1389778) q[1];
sx q[1];
rz(-1.9139342) q[1];
rz(1.8101495) q[2];
sx q[2];
rz(-1.6806921) q[2];
sx q[2];
rz(0.42721911) q[2];
rz(2.4155865) q[3];
sx q[3];
rz(-0.86961679) q[3];
sx q[3];
rz(1.1179954) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
