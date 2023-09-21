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
rz(1.5496594) q[1];
sx q[1];
rz(-0.078443371) q[1];
sx q[1];
rz(-0.6426386) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4896048) q[0];
sx q[0];
rz(-1.542924) q[0];
sx q[0];
rz(0.53919381) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7625916) q[2];
sx q[2];
rz(-2.0711053) q[2];
sx q[2];
rz(-1.9185916) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5832311) q[1];
sx q[1];
rz(-1.4500631) q[1];
sx q[1];
rz(1.4503149) q[1];
rz(0.1518941) q[3];
sx q[3];
rz(-0.23670247) q[3];
sx q[3];
rz(-0.13329245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46368018) q[0];
sx q[0];
rz(-1.3543411) q[0];
sx q[0];
rz(-0.98051488) q[0];
rz(2.9837043) q[1];
sx q[1];
rz(-1.4973367) q[1];
sx q[1];
rz(-2.352879) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0448488) q[0];
sx q[0];
rz(-1.6941438) q[0];
sx q[0];
rz(3.0268273) q[0];
rz(-0.2561432) q[2];
sx q[2];
rz(-1.5400585) q[2];
sx q[2];
rz(-2.5275633) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.84856725) q[1];
sx q[1];
rz(-1.8752521) q[1];
sx q[1];
rz(-0.57363631) q[1];
rz(-pi) q[2];
rz(-0.85487811) q[3];
sx q[3];
rz(-0.93920556) q[3];
sx q[3];
rz(-2.7090493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.366189) q[2];
sx q[2];
rz(-2.0663694) q[2];
sx q[2];
rz(-0.186084) q[2];
rz(-0.6535334) q[3];
sx q[3];
rz(-1.3862405) q[3];
sx q[3];
rz(-3.0236566) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2114975) q[0];
sx q[0];
rz(-0.94267541) q[0];
sx q[0];
rz(-0.63013664) q[0];
rz(-3.0139626) q[1];
sx q[1];
rz(-0.6487414) q[1];
sx q[1];
rz(2.4198467) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5436514) q[0];
sx q[0];
rz(-1.3184034) q[0];
sx q[0];
rz(0.6681722) q[0];
rz(-pi) q[1];
rz(2.0314991) q[2];
sx q[2];
rz(-0.93020541) q[2];
sx q[2];
rz(3.1415423) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7194781) q[1];
sx q[1];
rz(-2.0146703) q[1];
sx q[1];
rz(1.6133973) q[1];
rz(1.1214764) q[3];
sx q[3];
rz(-0.67766261) q[3];
sx q[3];
rz(-2.6636056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
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
x q[1];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58406126) q[0];
sx q[0];
rz(-0.73636213) q[0];
sx q[0];
rz(3.0990565) q[0];
rz(-0.78035367) q[1];
sx q[1];
rz(-2.6413554) q[1];
sx q[1];
rz(0.76400486) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9873136) q[0];
sx q[0];
rz(-2.2494082) q[0];
sx q[0];
rz(-0.016088386) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.64455428) q[2];
sx q[2];
rz(-1.6719712) q[2];
sx q[2];
rz(0.25601706) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3180247) q[1];
sx q[1];
rz(-1.9293702) q[1];
sx q[1];
rz(1.093868) q[1];
rz(-pi) q[2];
rz(-1.10631) q[3];
sx q[3];
rz(-0.96040695) q[3];
sx q[3];
rz(-0.27184091) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.2970695) q[2];
sx q[2];
rz(-1.2244747) q[2];
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
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
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
rz(2.7382543) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4606844) q[0];
sx q[0];
rz(-1.4049238) q[0];
sx q[0];
rz(0.18780639) q[0];
x q[1];
rz(0.43535797) q[2];
sx q[2];
rz(-2.407496) q[2];
sx q[2];
rz(-1.8045319) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4407318) q[1];
sx q[1];
rz(-1.4167538) q[1];
sx q[1];
rz(2.6101019) q[1];
x q[2];
rz(2.4661857) q[3];
sx q[3];
rz(-1.1408148) q[3];
sx q[3];
rz(2.5436385) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.32101813) q[2];
sx q[2];
rz(-1.0514739) q[2];
sx q[2];
rz(-2.3590951) q[2];
rz(-2.0292422) q[3];
sx q[3];
rz(-2.7045043) q[3];
sx q[3];
rz(-2.1330244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
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
rz(-1.1385304) q[1];
sx q[1];
rz(3.1351556) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71953668) q[0];
sx q[0];
rz(-2.3483843) q[0];
sx q[0];
rz(2.7842245) q[0];
rz(-pi) q[1];
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
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.9812614) q[1];
sx q[1];
rz(-2.6604974) q[1];
sx q[1];
rz(-1.0465924) q[1];
rz(-pi) q[2];
rz(-2.3859343) q[3];
sx q[3];
rz(-2.14058) q[3];
sx q[3];
rz(-0.33982402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.7727938) q[2];
sx q[2];
rz(-1.6364748) q[2];
sx q[2];
rz(2.2951365) q[2];
rz(-2.1438697) q[3];
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
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0434175) q[0];
sx q[0];
rz(-0.88413969) q[0];
sx q[0];
rz(0.055710677) q[0];
rz(-2.3588691) q[1];
sx q[1];
rz(-1.1306154) q[1];
sx q[1];
rz(1.1605211) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6877277) q[0];
sx q[0];
rz(-1.9181607) q[0];
sx q[0];
rz(2.0485282) q[0];
rz(2.6703228) q[2];
sx q[2];
rz(-2.0435213) q[2];
sx q[2];
rz(-1.3408692) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.5481422) q[1];
sx q[1];
rz(-1.1957809) q[1];
sx q[1];
rz(2.1795991) q[1];
rz(1.1858995) q[3];
sx q[3];
rz(-1.9249831) q[3];
sx q[3];
rz(-0.76364309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.065585) q[2];
sx q[2];
rz(-2.2202754) q[2];
sx q[2];
rz(-0.78545061) q[2];
rz(-2.3857332) q[3];
sx q[3];
rz(-1.8463585) q[3];
sx q[3];
rz(-2.7261962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
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
rz(-1.8287559) q[0];
sx q[0];
rz(-0.090933196) q[0];
sx q[0];
rz(1.1428517) q[0];
rz(1.8354592) q[1];
sx q[1];
rz(-1.1885234) q[1];
sx q[1];
rz(2.725504) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7284262) q[0];
sx q[0];
rz(-0.51983716) q[0];
sx q[0];
rz(1.1632989) q[0];
rz(-pi) q[1];
rz(-1.1025238) q[2];
sx q[2];
rz(-1.6981914) q[2];
sx q[2];
rz(-0.96291908) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.0033274) q[1];
sx q[1];
rz(-1.4107553) q[1];
sx q[1];
rz(-0.40469594) q[1];
x q[2];
rz(-1.0981406) q[3];
sx q[3];
rz(-1.6695767) q[3];
sx q[3];
rz(0.071382513) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.1064421) q[2];
sx q[2];
rz(-1.4539377) q[2];
sx q[2];
rz(2.8184334) q[2];
rz(-0.20251003) q[3];
sx q[3];
rz(-1.2812307) q[3];
sx q[3];
rz(-2.5642853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.768196) q[0];
sx q[0];
rz(-2.51238) q[0];
sx q[0];
rz(-1.9966104) q[0];
rz(-1.1960944) q[1];
sx q[1];
rz(-0.15592608) q[1];
sx q[1];
rz(2.6224565) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0275092) q[0];
sx q[0];
rz(-0.4013831) q[0];
sx q[0];
rz(1.9163577) q[0];
rz(0.22186188) q[2];
sx q[2];
rz(-1.8010745) q[2];
sx q[2];
rz(0.78267539) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5362894) q[1];
sx q[1];
rz(-1.481206) q[1];
sx q[1];
rz(-0.78219608) q[1];
x q[2];
rz(2.5095652) q[3];
sx q[3];
rz(-0.29187376) q[3];
sx q[3];
rz(-0.29575086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3141994) q[2];
sx q[2];
rz(-1.9284356) q[2];
sx q[2];
rz(-3.1006151) q[2];
rz(-0.86769062) q[3];
sx q[3];
rz(-0.65338457) q[3];
sx q[3];
rz(2.6303671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
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
rz(-1.4279667) q[1];
sx q[1];
rz(-2.7797996) q[1];
sx q[1];
rz(-1.6428927) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1536627) q[0];
sx q[0];
rz(-1.6462353) q[0];
sx q[0];
rz(3.1213785) q[0];
x q[1];
rz(-1.8685568) q[2];
sx q[2];
rz(-2.9794663) q[2];
sx q[2];
rz(-0.33725421) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5012706) q[1];
sx q[1];
rz(-2.3295799) q[1];
sx q[1];
rz(-2.0444319) q[1];
rz(1.6937709) q[3];
sx q[3];
rz(-1.3696559) q[3];
sx q[3];
rz(0.5826544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.75858086) q[2];
sx q[2];
rz(-2.0643533) q[2];
sx q[2];
rz(2.995058) q[2];
rz(-2.3274029) q[3];
sx q[3];
rz(-0.79629961) q[3];
sx q[3];
rz(0.41671419) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
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
rz(1.8756443) q[1];
sx q[1];
rz(-1.0026149) q[1];
sx q[1];
rz(1.2276585) q[1];
rz(-1.8101495) q[2];
sx q[2];
rz(-1.4609006) q[2];
sx q[2];
rz(-2.7143735) q[2];
rz(-0.72487763) q[3];
sx q[3];
rz(-1.0387883) q[3];
sx q[3];
rz(-3.0742857) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
