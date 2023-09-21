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
rz(2.0547325) q[0];
sx q[0];
rz(7.6261043) q[0];
rz(1.5496594) q[1];
sx q[1];
rz(-0.078443371) q[1];
sx q[1];
rz(2.4989541) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4896048) q[0];
sx q[0];
rz(-1.542924) q[0];
sx q[0];
rz(2.6023988) q[0];
rz(2.6334555) q[2];
sx q[2];
rz(-1.7388441) q[2];
sx q[2];
rz(-2.7009168) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.1145775) q[1];
sx q[1];
rz(-1.6903965) q[1];
sx q[1];
rz(-0.12160614) q[1];
rz(-2.9075165) q[3];
sx q[3];
rz(-1.5353068) q[3];
sx q[3];
rz(1.2897829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.6671483) q[2];
sx q[2];
rz(-0.64626226) q[2];
sx q[2];
rz(1.8623964) q[2];
rz(-0.71875087) q[3];
sx q[3];
rz(-1.5712534) q[3];
sx q[3];
rz(-0.32354245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
sx q[2];
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
rz(-1.4973367) q[1];
sx q[1];
rz(2.352879) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65608998) q[0];
sx q[0];
rz(-2.9733109) q[0];
sx q[0];
rz(-2.3165354) q[0];
rz(0.2561432) q[2];
sx q[2];
rz(-1.5400585) q[2];
sx q[2];
rz(2.5275633) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.1566369) q[1];
sx q[1];
rz(-2.5002694) q[1];
sx q[1];
rz(2.6167469) q[1];
rz(-pi) q[2];
rz(-0.73123587) q[3];
sx q[3];
rz(-0.9160708) q[3];
sx q[3];
rz(-0.54192858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.7754037) q[2];
sx q[2];
rz(-1.0752233) q[2];
sx q[2];
rz(-0.186084) q[2];
rz(-0.6535334) q[3];
sx q[3];
rz(-1.7553522) q[3];
sx q[3];
rz(-0.11793605) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9300951) q[0];
sx q[0];
rz(-0.94267541) q[0];
sx q[0];
rz(0.63013664) q[0];
rz(3.0139626) q[1];
sx q[1];
rz(-0.6487414) q[1];
sx q[1];
rz(-2.4198467) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5436514) q[0];
sx q[0];
rz(-1.8231892) q[0];
sx q[0];
rz(2.4734205) q[0];
rz(-0.53775215) q[2];
sx q[2];
rz(-0.76965145) q[2];
sx q[2];
rz(2.4485181) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.4221146) q[1];
sx q[1];
rz(-1.1269224) q[1];
sx q[1];
rz(-1.5281954) q[1];
rz(-2.0201163) q[3];
sx q[3];
rz(-0.67766261) q[3];
sx q[3];
rz(0.47798702) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.7599941) q[2];
sx q[2];
rz(-2.1888032) q[2];
sx q[2];
rz(2.2325113) q[2];
rz(2.6233853) q[3];
sx q[3];
rz(-1.0639233) q[3];
sx q[3];
rz(-0.28809965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(2.5575314) q[0];
sx q[0];
rz(-0.73636213) q[0];
sx q[0];
rz(0.042536143) q[0];
rz(0.78035367) q[1];
sx q[1];
rz(-2.6413554) q[1];
sx q[1];
rz(2.3775878) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4064179) q[0];
sx q[0];
rz(-1.58332) q[0];
sx q[0];
rz(-0.89212117) q[0];
x q[1];
rz(0.64455428) q[2];
sx q[2];
rz(-1.6719712) q[2];
sx q[2];
rz(2.8855756) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.82356794) q[1];
sx q[1];
rz(-1.9293702) q[1];
sx q[1];
rz(-2.0477247) q[1];
rz(-0.56960168) q[3];
sx q[3];
rz(-0.74865018) q[3];
sx q[3];
rz(2.1514055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.2970695) q[2];
sx q[2];
rz(-1.917118) q[2];
sx q[2];
rz(-2.6085473) q[2];
rz(2.879203) q[3];
sx q[3];
rz(-2.5049987) q[3];
sx q[3];
rz(-0.46245241) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
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
rz(-0.22661701) q[1];
sx q[1];
rz(-0.77426568) q[1];
sx q[1];
rz(0.40333834) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82542244) q[0];
sx q[0];
rz(-0.24992019) q[0];
sx q[0];
rz(0.73096801) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4558805) q[2];
sx q[2];
rz(-1.2843686) q[2];
sx q[2];
rz(-0.098766947) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.7008608) q[1];
sx q[1];
rz(-1.4167538) q[1];
sx q[1];
rz(-0.53149077) q[1];
rz(-pi) q[2];
rz(2.1020528) q[3];
sx q[3];
rz(-0.96635339) q[3];
sx q[3];
rz(-1.8464551) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.32101813) q[2];
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
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(1.7111506) q[0];
sx q[0];
rz(-0.40200457) q[0];
sx q[0];
rz(-0.45561403) q[0];
rz(-0.09952155) q[1];
sx q[1];
rz(-2.0030622) q[1];
sx q[1];
rz(-3.1351556) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59505263) q[0];
sx q[0];
rz(-1.8227302) q[0];
sx q[0];
rz(2.3810054) q[0];
x q[1];
rz(1.1941031) q[2];
sx q[2];
rz(-1.4720535) q[2];
sx q[2];
rz(-1.313098) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.9812614) q[1];
sx q[1];
rz(-0.48109522) q[1];
sx q[1];
rz(2.0950003) q[1];
rz(-pi) q[2];
rz(-0.84900093) q[3];
sx q[3];
rz(-0.95522049) q[3];
sx q[3];
rz(-1.4403696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.7727938) q[2];
sx q[2];
rz(-1.6364748) q[2];
sx q[2];
rz(-0.84645611) q[2];
rz(0.99772292) q[3];
sx q[3];
rz(-2.3322451) q[3];
sx q[3];
rz(1.458582) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0434175) q[0];
sx q[0];
rz(-2.257453) q[0];
sx q[0];
rz(3.085882) q[0];
rz(0.78272351) q[1];
sx q[1];
rz(-2.0109773) q[1];
sx q[1];
rz(-1.1605211) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8502064) q[0];
sx q[0];
rz(-1.123748) q[0];
sx q[0];
rz(0.38711754) q[0];
rz(-pi) q[1];
rz(-2.2968282) q[2];
sx q[2];
rz(-0.65462199) q[2];
sx q[2];
rz(0.95915937) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.6351663) q[1];
sx q[1];
rz(-0.70235683) q[1];
sx q[1];
rz(-2.173645) q[1];
x q[2];
rz(1.9556932) q[3];
sx q[3];
rz(-1.9249831) q[3];
sx q[3];
rz(-2.3779496) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.0760076) q[2];
sx q[2];
rz(-2.2202754) q[2];
sx q[2];
rz(-0.78545061) q[2];
rz(2.3857332) q[3];
sx q[3];
rz(-1.8463585) q[3];
sx q[3];
rz(-0.41539645) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3128368) q[0];
sx q[0];
rz(-3.0506595) q[0];
sx q[0];
rz(1.998741) q[0];
rz(1.3061334) q[1];
sx q[1];
rz(-1.1885234) q[1];
sx q[1];
rz(0.41608861) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.04836719) q[0];
sx q[0];
rz(-2.0443633) q[0];
sx q[0];
rz(-0.22305365) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2942737) q[2];
sx q[2];
rz(-2.6575436) q[2];
sx q[2];
rz(0.36177847) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.9302952) q[1];
sx q[1];
rz(-2.7080309) q[1];
sx q[1];
rz(-2.7525206) q[1];
x q[2];
rz(3.0307426) q[3];
sx q[3];
rz(-2.0409611) q[3];
sx q[3];
rz(-1.449031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.0351506) q[2];
sx q[2];
rz(-1.687655) q[2];
sx q[2];
rz(-2.8184334) q[2];
rz(-0.20251003) q[3];
sx q[3];
rz(-1.860362) q[3];
sx q[3];
rz(2.5642853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.768196) q[0];
sx q[0];
rz(-2.51238) q[0];
sx q[0];
rz(-1.1449822) q[0];
rz(1.9454983) q[1];
sx q[1];
rz(-2.9856666) q[1];
sx q[1];
rz(-2.6224565) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11408344) q[0];
sx q[0];
rz(-2.7402096) q[0];
sx q[0];
rz(1.2252349) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8066508) q[2];
sx q[2];
rz(-1.354885) q[2];
sx q[2];
rz(-0.73667919) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(3.0172826) q[1];
sx q[1];
rz(-0.78622422) q[1];
sx q[1];
rz(3.0148274) q[1];
rz(-pi) q[2];
rz(2.5095652) q[3];
sx q[3];
rz(-2.8497189) q[3];
sx q[3];
rz(-2.8458418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8273932) q[2];
sx q[2];
rz(-1.2131571) q[2];
sx q[2];
rz(-0.040977565) q[2];
rz(-2.273902) q[3];
sx q[3];
rz(-2.4882081) q[3];
sx q[3];
rz(-0.51122558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39559078) q[0];
sx q[0];
rz(-0.87930167) q[0];
sx q[0];
rz(-1.5270773) q[0];
rz(1.7136259) q[1];
sx q[1];
rz(-0.36179301) q[1];
sx q[1];
rz(1.6428927) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1536627) q[0];
sx q[0];
rz(-1.4953574) q[0];
sx q[0];
rz(3.1213785) q[0];
x q[1];
rz(-1.2730359) q[2];
sx q[2];
rz(-2.9794663) q[2];
sx q[2];
rz(-2.8043384) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.6403221) q[1];
sx q[1];
rz(-0.81201279) q[1];
sx q[1];
rz(-1.0971607) q[1];
rz(-1.6937709) q[3];
sx q[3];
rz(-1.3696559) q[3];
sx q[3];
rz(2.5589383) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.3830118) q[2];
sx q[2];
rz(-1.0772394) q[2];
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
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9183337) q[0];
sx q[0];
rz(-1.8948566) q[0];
sx q[0];
rz(-2.1444453) q[0];
rz(1.8756443) q[1];
sx q[1];
rz(-1.0026149) q[1];
sx q[1];
rz(1.2276585) q[1];
rz(-3.0284991) q[2];
sx q[2];
rz(-1.332915) q[2];
sx q[2];
rz(-1.1168196) q[2];
rz(0.90445789) q[3];
sx q[3];
rz(-2.1790128) q[3];
sx q[3];
rz(2.0603767) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];