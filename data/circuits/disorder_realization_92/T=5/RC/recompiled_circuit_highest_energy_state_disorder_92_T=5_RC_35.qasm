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
rz(-1.5795213) q[0];
sx q[0];
rz(6.4223839) q[0];
sx q[0];
rz(9.7425707) q[0];
rz(-2.3843482) q[1];
sx q[1];
rz(-1.550996) q[1];
sx q[1];
rz(-1.6928033) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1549281) q[0];
sx q[0];
rz(-1.9408615) q[0];
sx q[0];
rz(1.4246902) q[0];
x q[1];
rz(-2.8002003) q[2];
sx q[2];
rz(-1.72014) q[2];
sx q[2];
rz(0.53440375) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.536126) q[1];
sx q[1];
rz(-1.2647293) q[1];
sx q[1];
rz(1.1121967) q[1];
rz(-pi) q[2];
rz(1.3643614) q[3];
sx q[3];
rz(-0.41461333) q[3];
sx q[3];
rz(-1.1846836) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.455787) q[2];
sx q[2];
rz(-2.0152342) q[2];
sx q[2];
rz(0.2618928) q[2];
rz(-2.453878) q[3];
sx q[3];
rz(-2.3759638) q[3];
sx q[3];
rz(2.8253637) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19551109) q[0];
sx q[0];
rz(-1.8571778) q[0];
sx q[0];
rz(2.2971357) q[0];
rz(0.63956815) q[1];
sx q[1];
rz(-0.25194672) q[1];
sx q[1];
rz(1.2676988) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15877341) q[0];
sx q[0];
rz(-2.4381659) q[0];
sx q[0];
rz(1.3838468) q[0];
rz(-pi) q[1];
rz(2.4652609) q[2];
sx q[2];
rz(-1.6937635) q[2];
sx q[2];
rz(0.42675297) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.69985327) q[1];
sx q[1];
rz(-1.5754712) q[1];
sx q[1];
rz(-0.75231268) q[1];
rz(-0.92314641) q[3];
sx q[3];
rz(-0.83469838) q[3];
sx q[3];
rz(-1.0506964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.1203479) q[2];
sx q[2];
rz(-0.52560386) q[2];
sx q[2];
rz(-2.8525103) q[2];
rz(-2.8287079) q[3];
sx q[3];
rz(-0.94066921) q[3];
sx q[3];
rz(0.59278929) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4631571) q[0];
sx q[0];
rz(-2.1364697) q[0];
sx q[0];
rz(2.4682755) q[0];
rz(-2.8939269) q[1];
sx q[1];
rz(-2.1238056) q[1];
sx q[1];
rz(2.8376875) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.073868) q[0];
sx q[0];
rz(-1.0751599) q[0];
sx q[0];
rz(-0.93511029) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.7627569) q[2];
sx q[2];
rz(-2.1386563) q[2];
sx q[2];
rz(-3.0936846) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5317638) q[1];
sx q[1];
rz(-1.4312688) q[1];
sx q[1];
rz(-0.05910766) q[1];
x q[2];
rz(0.21411333) q[3];
sx q[3];
rz(-1.1592093) q[3];
sx q[3];
rz(1.0836835) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.5725382) q[2];
sx q[2];
rz(-1.7819541) q[2];
sx q[2];
rz(3.1247826) q[2];
rz(-2.861764) q[3];
sx q[3];
rz(-0.40585104) q[3];
sx q[3];
rz(0.60389891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0458321) q[0];
sx q[0];
rz(-1.1149167) q[0];
sx q[0];
rz(-1.8810077) q[0];
rz(2.6450805) q[1];
sx q[1];
rz(-1.1505726) q[1];
sx q[1];
rz(1.5042492) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.838546) q[0];
sx q[0];
rz(-1.9936512) q[0];
sx q[0];
rz(2.1648876) q[0];
rz(-pi) q[1];
rz(-0.38522933) q[2];
sx q[2];
rz(-2.5653815) q[2];
sx q[2];
rz(0.50265898) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.97442409) q[1];
sx q[1];
rz(-1.7361871) q[1];
sx q[1];
rz(-2.1163634) q[1];
rz(-pi) q[2];
rz(-2.6332985) q[3];
sx q[3];
rz(-2.5807305) q[3];
sx q[3];
rz(-1.9702553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.090195) q[2];
sx q[2];
rz(-2.6020256) q[2];
sx q[2];
rz(-2.9810193) q[2];
rz(-1.4422902) q[3];
sx q[3];
rz(-1.620404) q[3];
sx q[3];
rz(-0.035339385) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8869239) q[0];
sx q[0];
rz(-1.2441607) q[0];
sx q[0];
rz(0.96920335) q[0];
rz(1.8907549) q[1];
sx q[1];
rz(-0.67986095) q[1];
sx q[1];
rz(-0.23076898) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31393181) q[0];
sx q[0];
rz(-1.5761762) q[0];
sx q[0];
rz(1.5350585) q[0];
x q[1];
rz(-1.5838301) q[2];
sx q[2];
rz(-2.763204) q[2];
sx q[2];
rz(-2.3005444) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.240162) q[1];
sx q[1];
rz(-1.104875) q[1];
sx q[1];
rz(0.45216143) q[1];
x q[2];
rz(-0.67471116) q[3];
sx q[3];
rz(-0.58739788) q[3];
sx q[3];
rz(-0.34492465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.1988924) q[2];
sx q[2];
rz(-2.5653699) q[2];
sx q[2];
rz(0.25299859) q[2];
rz(-1.7871855) q[3];
sx q[3];
rz(-1.1043045) q[3];
sx q[3];
rz(1.3033006) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81786466) q[0];
sx q[0];
rz(-2.4619894) q[0];
sx q[0];
rz(-0.060081765) q[0];
rz(-0.59728638) q[1];
sx q[1];
rz(-2.2164454) q[1];
sx q[1];
rz(-3.0937321) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5756719) q[0];
sx q[0];
rz(-1.4788155) q[0];
sx q[0];
rz(-3.1264114) q[0];
rz(0.20345236) q[2];
sx q[2];
rz(-1.5849539) q[2];
sx q[2];
rz(-0.79256533) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.3678494) q[1];
sx q[1];
rz(-1.7715104) q[1];
sx q[1];
rz(-1.9108921) q[1];
rz(-pi) q[2];
rz(2.8196857) q[3];
sx q[3];
rz(-1.6471146) q[3];
sx q[3];
rz(0.25086153) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.73703274) q[2];
sx q[2];
rz(-0.76475778) q[2];
sx q[2];
rz(-0.16727373) q[2];
rz(3.0625999) q[3];
sx q[3];
rz(-1.1827712) q[3];
sx q[3];
rz(0.8305859) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
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
rz(-1.1792184) q[0];
sx q[0];
rz(-3.0308864) q[0];
sx q[0];
rz(-3.0141444) q[0];
rz(2.2178862) q[1];
sx q[1];
rz(-0.5674924) q[1];
sx q[1];
rz(-2.4097402) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29402052) q[0];
sx q[0];
rz(-2.5811727) q[0];
sx q[0];
rz(0.27604485) q[0];
x q[1];
rz(0.62035608) q[2];
sx q[2];
rz(-0.51232238) q[2];
sx q[2];
rz(-2.870395) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.7596563) q[1];
sx q[1];
rz(-2.3664673) q[1];
sx q[1];
rz(-1.0785021) q[1];
x q[2];
rz(2.022677) q[3];
sx q[3];
rz(-2.2377399) q[3];
sx q[3];
rz(-0.42740145) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.38714108) q[2];
sx q[2];
rz(-1.0820505) q[2];
sx q[2];
rz(-2.2608536) q[2];
rz(2.8494917) q[3];
sx q[3];
rz(-2.5160774) q[3];
sx q[3];
rz(-0.97091278) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77325118) q[0];
sx q[0];
rz(-1.0057978) q[0];
sx q[0];
rz(-0.7811501) q[0];
rz(1.7225522) q[1];
sx q[1];
rz(-0.96839372) q[1];
sx q[1];
rz(0.17793812) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3969361) q[0];
sx q[0];
rz(-1.2632003) q[0];
sx q[0];
rz(-1.164695) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0974085) q[2];
sx q[2];
rz(-1.1605613) q[2];
sx q[2];
rz(0.1178478) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.7922349) q[1];
sx q[1];
rz(-0.48949896) q[1];
sx q[1];
rz(1.6952747) q[1];
rz(1.7573886) q[3];
sx q[3];
rz(-1.8301788) q[3];
sx q[3];
rz(0.6849388) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.6302744) q[2];
sx q[2];
rz(-1.7061468) q[2];
sx q[2];
rz(-1.1159631) q[2];
rz(2.337194) q[3];
sx q[3];
rz(-0.65785995) q[3];
sx q[3];
rz(0.72949725) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6804752) q[0];
sx q[0];
rz(-2.456433) q[0];
sx q[0];
rz(-2.6035736) q[0];
rz(-0.33942014) q[1];
sx q[1];
rz(-1.598204) q[1];
sx q[1];
rz(-3.0982049) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26608426) q[0];
sx q[0];
rz(-1.636063) q[0];
sx q[0];
rz(-1.4738183) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.814079) q[2];
sx q[2];
rz(-0.56381153) q[2];
sx q[2];
rz(-0.76432894) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.2895537) q[1];
sx q[1];
rz(-2.3869136) q[1];
sx q[1];
rz(-0.11615495) q[1];
rz(-2.6049117) q[3];
sx q[3];
rz(-2.4701475) q[3];
sx q[3];
rz(-1.0332024) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.033120774) q[2];
sx q[2];
rz(-2.7816009) q[2];
sx q[2];
rz(-1.6610422) q[2];
rz(0.2964274) q[3];
sx q[3];
rz(-1.1400283) q[3];
sx q[3];
rz(-2.3886783) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4794889) q[0];
sx q[0];
rz(-0.39283735) q[0];
sx q[0];
rz(-1.9774849) q[0];
rz(-2.8083943) q[1];
sx q[1];
rz(-1.4771799) q[1];
sx q[1];
rz(2.3936757) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5912022) q[0];
sx q[0];
rz(-1.8699614) q[0];
sx q[0];
rz(-2.5183866) q[0];
rz(-pi) q[1];
rz(-0.21166734) q[2];
sx q[2];
rz(-0.94918409) q[2];
sx q[2];
rz(1.7424316) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.038422) q[1];
sx q[1];
rz(-1.5161247) q[1];
sx q[1];
rz(1.5780539) q[1];
x q[2];
rz(1.6259057) q[3];
sx q[3];
rz(-2.4477473) q[3];
sx q[3];
rz(-1.9195279) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.36772874) q[2];
sx q[2];
rz(-1.7325956) q[2];
sx q[2];
rz(1.9749953) q[2];
rz(1.5424607) q[3];
sx q[3];
rz(-1.5365994) q[3];
sx q[3];
rz(2.0994073) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
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
rz(-1.790697) q[0];
sx q[0];
rz(-2.6804374) q[0];
sx q[0];
rz(2.8275369) q[0];
rz(-0.028388609) q[1];
sx q[1];
rz(-2.8969565) q[1];
sx q[1];
rz(1.2895186) q[1];
rz(0.51043541) q[2];
sx q[2];
rz(-1.0344997) q[2];
sx q[2];
rz(-2.79881) q[2];
rz(2.3343347) q[3];
sx q[3];
rz(-1.1716598) q[3];
sx q[3];
rz(-1.4902478) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
