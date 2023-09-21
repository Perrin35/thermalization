OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.1768567) q[0];
sx q[0];
rz(5.4194874) q[0];
sx q[0];
rz(9.6258862) q[0];
rz(1.7445298) q[1];
sx q[1];
rz(-1.3367329) q[1];
sx q[1];
rz(0.62682682) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7456822) q[0];
sx q[0];
rz(-2.9183309) q[0];
sx q[0];
rz(-1.0049055) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5457382) q[2];
sx q[2];
rz(-2.7239954) q[2];
sx q[2];
rz(-2.8516304) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.37571733) q[1];
sx q[1];
rz(-0.98343508) q[1];
sx q[1];
rz(-0.37382965) q[1];
rz(-pi) q[2];
rz(-2.9653373) q[3];
sx q[3];
rz(-1.3000543) q[3];
sx q[3];
rz(0.83812974) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.8864002) q[2];
sx q[2];
rz(-0.93078405) q[2];
sx q[2];
rz(1.8480776) q[2];
rz(-1.1039929) q[3];
sx q[3];
rz(-2.2938426) q[3];
sx q[3];
rz(-2.3867992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
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
rz(-2.8936546) q[0];
sx q[0];
rz(-1.0590483) q[0];
sx q[0];
rz(-0.98518103) q[0];
rz(2.7424116) q[1];
sx q[1];
rz(-1.8789411) q[1];
sx q[1];
rz(3.0342297) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0221314) q[0];
sx q[0];
rz(-0.71170002) q[0];
sx q[0];
rz(2.2244781) q[0];
rz(0.58789247) q[2];
sx q[2];
rz(-1.8111472) q[2];
sx q[2];
rz(2.8628778) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.0774035) q[1];
sx q[1];
rz(-1.940155) q[1];
sx q[1];
rz(2.2662152) q[1];
x q[2];
rz(0.3932088) q[3];
sx q[3];
rz(-2.7244096) q[3];
sx q[3];
rz(-1.259491) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5169516) q[2];
sx q[2];
rz(-1.7525643) q[2];
sx q[2];
rz(0.66169935) q[2];
rz(0.055756904) q[3];
sx q[3];
rz(-0.35900933) q[3];
sx q[3];
rz(-0.24584809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4645638) q[0];
sx q[0];
rz(-0.55389535) q[0];
sx q[0];
rz(0.04034986) q[0];
rz(2.5616052) q[1];
sx q[1];
rz(-2.2433387) q[1];
sx q[1];
rz(2.0592164) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4989935) q[0];
sx q[0];
rz(-1.9369619) q[0];
sx q[0];
rz(0.55143349) q[0];
rz(-2.9897887) q[2];
sx q[2];
rz(-1.9215343) q[2];
sx q[2];
rz(2.0098067) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.1632417) q[1];
sx q[1];
rz(-0.32838531) q[1];
sx q[1];
rz(2.3204625) q[1];
x q[2];
rz(0.2228959) q[3];
sx q[3];
rz(-1.8867023) q[3];
sx q[3];
rz(-0.99881682) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.0064156) q[2];
sx q[2];
rz(-1.2536896) q[2];
sx q[2];
rz(2.5877) q[2];
rz(1.4931549) q[3];
sx q[3];
rz(-1.0529543) q[3];
sx q[3];
rz(2.5890787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4937113) q[0];
sx q[0];
rz(-2.557502) q[0];
sx q[0];
rz(-0.035382263) q[0];
rz(-2.1454051) q[1];
sx q[1];
rz(-1.532998) q[1];
sx q[1];
rz(-2.6534973) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1195388) q[0];
sx q[0];
rz(-1.4931803) q[0];
sx q[0];
rz(2.0366497) q[0];
rz(-0.98056356) q[2];
sx q[2];
rz(-2.0681941) q[2];
sx q[2];
rz(-0.24888466) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.19255895) q[1];
sx q[1];
rz(-2.9395736) q[1];
sx q[1];
rz(2.2662524) q[1];
rz(-pi) q[2];
rz(1.0933236) q[3];
sx q[3];
rz(-2.100088) q[3];
sx q[3];
rz(0.57100163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.35486832) q[2];
sx q[2];
rz(-1.0151261) q[2];
sx q[2];
rz(-2.6341237) q[2];
rz(1.1821702) q[3];
sx q[3];
rz(-1.2927262) q[3];
sx q[3];
rz(-1.2549843) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5421211) q[0];
sx q[0];
rz(-2.4276955) q[0];
sx q[0];
rz(0.3702634) q[0];
rz(-2.8778991) q[1];
sx q[1];
rz(-1.5779243) q[1];
sx q[1];
rz(-0.19651861) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4465966) q[0];
sx q[0];
rz(-1.2906404) q[0];
sx q[0];
rz(-0.3145991) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0318803) q[2];
sx q[2];
rz(-1.5421252) q[2];
sx q[2];
rz(1.920514) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.6424751) q[1];
sx q[1];
rz(-0.31458464) q[1];
sx q[1];
rz(2.5423074) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5971562) q[3];
sx q[3];
rz(-2.8215373) q[3];
sx q[3];
rz(-0.65359945) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.013441) q[2];
sx q[2];
rz(-0.95169607) q[2];
sx q[2];
rz(-2.2244942) q[2];
rz(2.3069978) q[3];
sx q[3];
rz(-1.3093427) q[3];
sx q[3];
rz(-2.807907) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31345263) q[0];
sx q[0];
rz(-1.4414635) q[0];
sx q[0];
rz(0.20203461) q[0];
rz(-2.0687885) q[1];
sx q[1];
rz(-1.5067116) q[1];
sx q[1];
rz(1.3938168) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4112339) q[0];
sx q[0];
rz(-0.3404091) q[0];
sx q[0];
rz(-2.911724) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9453085) q[2];
sx q[2];
rz(-1.6746582) q[2];
sx q[2];
rz(-1.9749157) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.2653633) q[1];
sx q[1];
rz(-0.47788844) q[1];
sx q[1];
rz(-2.1761314) q[1];
x q[2];
rz(-2.6447547) q[3];
sx q[3];
rz(-2.3416069) q[3];
sx q[3];
rz(-2.5132688) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.94971925) q[2];
sx q[2];
rz(-1.4800025) q[2];
sx q[2];
rz(-0.8423155) q[2];
rz(-1.0874282) q[3];
sx q[3];
rz(-0.73173404) q[3];
sx q[3];
rz(-1.6585763) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52952805) q[0];
sx q[0];
rz(-1.8719215) q[0];
sx q[0];
rz(-2.1202309) q[0];
rz(-1.1313324) q[1];
sx q[1];
rz(-2.1919057) q[1];
sx q[1];
rz(0.21025118) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4410941) q[0];
sx q[0];
rz(-2.9919618) q[0];
sx q[0];
rz(1.11087) q[0];
x q[1];
rz(-1.636593) q[2];
sx q[2];
rz(-1.4074667) q[2];
sx q[2];
rz(-2.3492299) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.6639858) q[1];
sx q[1];
rz(-0.37316445) q[1];
sx q[1];
rz(2.5928156) q[1];
rz(-pi) q[2];
rz(-0.57742124) q[3];
sx q[3];
rz(-0.93772674) q[3];
sx q[3];
rz(-3.0949743) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.82688275) q[2];
sx q[2];
rz(-1.751739) q[2];
sx q[2];
rz(-1.7604527) q[2];
rz(-0.76210493) q[3];
sx q[3];
rz(-2.6657181) q[3];
sx q[3];
rz(2.7281318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
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
rz(-0.49522266) q[0];
sx q[0];
rz(-0.78283993) q[0];
sx q[0];
rz(2.5543509) q[0];
rz(-2.6953221) q[1];
sx q[1];
rz(-1.8621567) q[1];
sx q[1];
rz(-2.1648724) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8551089) q[0];
sx q[0];
rz(-2.6628462) q[0];
sx q[0];
rz(-0.32740645) q[0];
rz(-pi) q[1];
rz(-1.4090528) q[2];
sx q[2];
rz(-1.4913034) q[2];
sx q[2];
rz(3.1331935) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.8969438) q[1];
sx q[1];
rz(-1.9172501) q[1];
sx q[1];
rz(1.1036554) q[1];
rz(-2.1920491) q[3];
sx q[3];
rz(-2.3275259) q[3];
sx q[3];
rz(2.4809588) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.4010767) q[2];
sx q[2];
rz(-2.433321) q[2];
sx q[2];
rz(0.55244279) q[2];
rz(0.46594122) q[3];
sx q[3];
rz(-1.7539141) q[3];
sx q[3];
rz(2.1140816) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66822806) q[0];
sx q[0];
rz(-2.3478274) q[0];
sx q[0];
rz(-2.2209432) q[0];
rz(0.051041516) q[1];
sx q[1];
rz(-1.5565245) q[1];
sx q[1];
rz(2.2424973) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5024922) q[0];
sx q[0];
rz(-0.99196767) q[0];
sx q[0];
rz(2.6274908) q[0];
x q[1];
rz(-1.012072) q[2];
sx q[2];
rz(-1.2860824) q[2];
sx q[2];
rz(-0.65288359) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.2832665) q[1];
sx q[1];
rz(-2.7910633) q[1];
sx q[1];
rz(1.9117029) q[1];
rz(-pi) q[2];
rz(-3.1413792) q[3];
sx q[3];
rz(-1.9409688) q[3];
sx q[3];
rz(2.7406524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.94545323) q[2];
sx q[2];
rz(-0.070572704) q[2];
sx q[2];
rz(3.0370039) q[2];
rz(2.3032522) q[3];
sx q[3];
rz(-0.74931562) q[3];
sx q[3];
rz(-0.88275638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33912441) q[0];
sx q[0];
rz(-2.323928) q[0];
sx q[0];
rz(-2.0237645) q[0];
rz(-0.71240187) q[1];
sx q[1];
rz(-2.5103705) q[1];
sx q[1];
rz(-2.7446279) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2755651) q[0];
sx q[0];
rz(-0.70789982) q[0];
sx q[0];
rz(-0.78535725) q[0];
x q[1];
rz(-2.8357387) q[2];
sx q[2];
rz(-1.2199739) q[2];
sx q[2];
rz(-2.2985814) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.18394477) q[1];
sx q[1];
rz(-0.34788222) q[1];
sx q[1];
rz(1.3089048) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7361717) q[3];
sx q[3];
rz(-1.1076895) q[3];
sx q[3];
rz(-0.4511569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.1931856) q[2];
sx q[2];
rz(-2.2414175) q[2];
sx q[2];
rz(-0.096654264) q[2];
rz(1.4139253) q[3];
sx q[3];
rz(-1.1380514) q[3];
sx q[3];
rz(-2.0146501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5861355) q[0];
sx q[0];
rz(-1.1145034) q[0];
sx q[0];
rz(0.282423) q[0];
rz(-1.2697521) q[1];
sx q[1];
rz(-1.3602263) q[1];
sx q[1];
rz(0.7855986) q[1];
rz(0.10791049) q[2];
sx q[2];
rz(-1.4866598) q[2];
sx q[2];
rz(1.286081) q[2];
rz(1.1453015) q[3];
sx q[3];
rz(-0.28278657) q[3];
sx q[3];
rz(0.84859802) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
