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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0138577) q[0];
sx q[0];
rz(-2.60175) q[0];
sx q[0];
rz(3.0873469) q[0];
x q[1];
rz(-1.379001) q[2];
sx q[2];
rz(-2.0711053) q[2];
sx q[2];
rz(1.9185916) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.3712284) q[1];
sx q[1];
rz(-2.9712354) q[1];
sx q[1];
rz(0.78070663) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5343127) q[3];
sx q[3];
rz(-1.3368703) q[3];
sx q[3];
rz(-0.28947383) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.6671483) q[2];
sx q[2];
rz(-0.64626226) q[2];
sx q[2];
rz(-1.8623964) q[2];
rz(-0.71875087) q[3];
sx q[3];
rz(-1.5703392) q[3];
sx q[3];
rz(0.32354245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6779125) q[0];
sx q[0];
rz(-1.3543411) q[0];
sx q[0];
rz(-0.98051488) q[0];
rz(0.15788831) q[1];
sx q[1];
rz(-1.4973367) q[1];
sx q[1];
rz(-0.78871361) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.096743874) q[0];
sx q[0];
rz(-1.6941438) q[0];
sx q[0];
rz(0.11476536) q[0];
rz(-pi) q[1];
rz(-2.8854495) q[2];
sx q[2];
rz(-1.6015341) q[2];
sx q[2];
rz(0.61402938) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.84856725) q[1];
sx q[1];
rz(-1.8752521) q[1];
sx q[1];
rz(-2.5679563) q[1];
rz(-0.7699645) q[3];
sx q[3];
rz(-1.0125481) q[3];
sx q[3];
rz(1.528873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.7754037) q[2];
sx q[2];
rz(-1.0752233) q[2];
sx q[2];
rz(0.186084) q[2];
rz(-0.6535334) q[3];
sx q[3];
rz(-1.3862405) q[3];
sx q[3];
rz(0.11793605) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2114975) q[0];
sx q[0];
rz(-2.1989172) q[0];
sx q[0];
rz(-2.511456) q[0];
rz(0.12763003) q[1];
sx q[1];
rz(-0.6487414) q[1];
sx q[1];
rz(2.4198467) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33357027) q[0];
sx q[0];
rz(-0.70735065) q[0];
sx q[0];
rz(2.7471514) q[0];
x q[1];
rz(-2.4475054) q[2];
sx q[2];
rz(-1.9352479) q[2];
sx q[2];
rz(-1.8592161) q[2];
rz(-pi/2) q[3];
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
rz(-3.0522703) q[1];
rz(2.1980522) q[3];
sx q[3];
rz(-1.2949847) q[3];
sx q[3];
rz(0.73345473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.7599941) q[2];
sx q[2];
rz(-0.95278946) q[2];
sx q[2];
rz(2.2325113) q[2];
rz(-2.6233853) q[3];
sx q[3];
rz(-2.0776694) q[3];
sx q[3];
rz(-0.28809965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-2.5575314) q[0];
sx q[0];
rz(-2.4052305) q[0];
sx q[0];
rz(0.042536143) q[0];
rz(0.78035367) q[1];
sx q[1];
rz(-2.6413554) q[1];
sx q[1];
rz(2.3775878) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9873136) q[0];
sx q[0];
rz(-0.89218441) q[0];
sx q[0];
rz(-3.1255043) q[0];
x q[1];
rz(-0.1673844) q[2];
sx q[2];
rz(-0.65132729) q[2];
sx q[2];
rz(-1.960388) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.92661392) q[1];
sx q[1];
rz(-1.1264631) q[1];
sx q[1];
rz(2.742393) q[1];
x q[2];
rz(2.4776811) q[3];
sx q[3];
rz(-1.9466562) q[3];
sx q[3];
rz(2.1223048) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.8445231) q[2];
sx q[2];
rz(-1.917118) q[2];
sx q[2];
rz(-2.6085473) q[2];
rz(-2.879203) q[3];
sx q[3];
rz(-2.5049987) q[3];
sx q[3];
rz(0.46245241) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9168636) q[0];
sx q[0];
rz(-0.30039766) q[0];
sx q[0];
rz(1.2438783) q[0];
rz(-0.22661701) q[1];
sx q[1];
rz(-2.367327) q[1];
sx q[1];
rz(2.7382543) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.078743155) q[0];
sx q[0];
rz(-1.7559949) q[0];
sx q[0];
rz(-1.4020105) q[0];
rz(-0.43535797) q[2];
sx q[2];
rz(-0.73409664) q[2];
sx q[2];
rz(-1.8045319) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.9215645) q[1];
sx q[1];
rz(-2.0953396) q[1];
sx q[1];
rz(-1.3925874) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1020528) q[3];
sx q[3];
rz(-0.96635339) q[3];
sx q[3];
rz(1.8464551) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.8205745) q[2];
sx q[2];
rz(-2.0901188) q[2];
sx q[2];
rz(0.78249758) q[2];
rz(2.0292422) q[3];
sx q[3];
rz(-0.4370884) q[3];
sx q[3];
rz(1.0085683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.430442) q[0];
sx q[0];
rz(-0.40200457) q[0];
sx q[0];
rz(2.6859786) q[0];
rz(-0.09952155) q[1];
sx q[1];
rz(-2.0030622) q[1];
sx q[1];
rz(0.0064370357) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71953668) q[0];
sx q[0];
rz(-2.3483843) q[0];
sx q[0];
rz(-0.35736812) q[0];
x q[1];
rz(-3.0354584) q[2];
sx q[2];
rz(-1.9455633) q[2];
sx q[2];
rz(2.8449164) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.7240119) q[1];
sx q[1];
rz(-1.15861) q[1];
sx q[1];
rz(0.25556232) q[1];
rz(2.2925917) q[3];
sx q[3];
rz(-2.1863722) q[3];
sx q[3];
rz(1.4403696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.7727938) q[2];
sx q[2];
rz(-1.6364748) q[2];
sx q[2];
rz(-0.84645611) q[2];
rz(2.1438697) q[3];
sx q[3];
rz(-2.3322451) q[3];
sx q[3];
rz(-1.458582) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0434175) q[0];
sx q[0];
rz(-2.257453) q[0];
sx q[0];
rz(0.055710677) q[0];
rz(2.3588691) q[1];
sx q[1];
rz(-1.1306154) q[1];
sx q[1];
rz(-1.1605211) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46471483) q[0];
sx q[0];
rz(-2.5589295) q[0];
sx q[0];
rz(-2.23784) q[0];
rz(-pi) q[1];
rz(1.0497401) q[2];
sx q[2];
rz(-1.1546635) q[2];
sx q[2];
rz(0.001948826) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.22735587) q[1];
sx q[1];
rz(-2.1319234) q[1];
sx q[1];
rz(-0.44740541) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1858995) q[3];
sx q[3];
rz(-1.9249831) q[3];
sx q[3];
rz(0.76364309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.065585) q[2];
sx q[2];
rz(-2.2202754) q[2];
sx q[2];
rz(2.356142) q[2];
rz(2.3857332) q[3];
sx q[3];
rz(-1.8463585) q[3];
sx q[3];
rz(2.7261962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-1.3128368) q[0];
sx q[0];
rz(-0.090933196) q[0];
sx q[0];
rz(-1.1428517) q[0];
rz(1.3061334) q[1];
sx q[1];
rz(-1.1885234) q[1];
sx q[1];
rz(0.41608861) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5160822) q[0];
sx q[0];
rz(-1.372638) q[0];
sx q[0];
rz(1.0869736) q[0];
x q[1];
rz(-1.1025238) q[2];
sx q[2];
rz(-1.4434012) q[2];
sx q[2];
rz(-2.1786736) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.63562288) q[1];
sx q[1];
rz(-1.9700248) q[1];
sx q[1];
rz(-1.3969621) q[1];
x q[2];
rz(1.3564524) q[3];
sx q[3];
rz(-0.48210258) q[3];
sx q[3];
rz(1.4516423) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.1064421) q[2];
sx q[2];
rz(-1.4539377) q[2];
sx q[2];
rz(0.32315928) q[2];
rz(-0.20251003) q[3];
sx q[3];
rz(-1.860362) q[3];
sx q[3];
rz(-0.57730738) q[3];
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
rz(pi/2) q[0];
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
rz(-0.37339661) q[0];
sx q[0];
rz(-2.51238) q[0];
sx q[0];
rz(1.9966104) q[0];
rz(1.1960944) q[1];
sx q[1];
rz(-0.15592608) q[1];
sx q[1];
rz(0.51913613) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11408344) q[0];
sx q[0];
rz(-0.4013831) q[0];
sx q[0];
rz(-1.2252349) q[0];
rz(-pi) q[1];
rz(-1.8066508) q[2];
sx q[2];
rz(-1.7867076) q[2];
sx q[2];
rz(-0.73667919) q[2];
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
rz(-1.395123) q[3];
sx q[3];
rz(-1.8050977) q[3];
sx q[3];
rz(-0.3570041) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.3141994) q[2];
sx q[2];
rz(-1.9284356) q[2];
sx q[2];
rz(-3.1006151) q[2];
rz(0.86769062) q[3];
sx q[3];
rz(-0.65338457) q[3];
sx q[3];
rz(-2.6303671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7460019) q[0];
sx q[0];
rz(-0.87930167) q[0];
sx q[0];
rz(1.6145153) q[0];
rz(-1.7136259) q[1];
sx q[1];
rz(-0.36179301) q[1];
sx q[1];
rz(1.4987) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4157383) q[0];
sx q[0];
rz(-3.0634974) q[0];
sx q[0];
rz(-1.3094835) q[0];
x q[1];
rz(-1.8685568) q[2];
sx q[2];
rz(-0.16212633) q[2];
sx q[2];
rz(0.33725421) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.5012706) q[1];
sx q[1];
rz(-2.3295799) q[1];
sx q[1];
rz(2.0444319) q[1];
x q[2];
rz(0.20262952) q[3];
sx q[3];
rz(-1.6912795) q[3];
sx q[3];
rz(-1.0128302) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.75858086) q[2];
sx q[2];
rz(-1.0772394) q[2];
sx q[2];
rz(0.14653462) q[2];
rz(0.81418973) q[3];
sx q[3];
rz(-2.345293) q[3];
sx q[3];
rz(-0.41671419) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
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
rz(-2.2371348) q[3];
sx q[3];
rz(-2.1790128) q[3];
sx q[3];
rz(2.0603767) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
