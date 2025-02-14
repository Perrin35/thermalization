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
rz(0.12662521) q[0];
sx q[0];
rz(-1.5669444) q[0];
sx q[0];
rz(2.6259165) q[0];
rz(0.22817837) q[1];
sx q[1];
rz(5.5362267) q[1];
sx q[1];
rz(9.8510392) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92862064) q[0];
sx q[0];
rz(-0.80982319) q[0];
sx q[0];
rz(2.9024603) q[0];
rz(-pi) q[1];
rz(-1.8581763) q[2];
sx q[2];
rz(-1.0848019) q[2];
sx q[2];
rz(3.0241242) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.8378046) q[1];
sx q[1];
rz(-0.96988867) q[1];
sx q[1];
rz(-1.5008885) q[1];
rz(-2.0419042) q[3];
sx q[3];
rz(-2.503241) q[3];
sx q[3];
rz(1.9698576) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.2446186) q[2];
sx q[2];
rz(-1.8316869) q[2];
sx q[2];
rz(2.8986325) q[2];
rz(-0.36863676) q[3];
sx q[3];
rz(-2.536085) q[3];
sx q[3];
rz(-1.2322371) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31545562) q[0];
sx q[0];
rz(-2.9189126) q[0];
sx q[0];
rz(-0.36732236) q[0];
rz(-0.79633725) q[1];
sx q[1];
rz(-1.0581191) q[1];
sx q[1];
rz(-2.499089) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2655661) q[0];
sx q[0];
rz(-2.4867704) q[0];
sx q[0];
rz(-0.73237082) q[0];
rz(-2.5902469) q[2];
sx q[2];
rz(-1.9187821) q[2];
sx q[2];
rz(2.6642193) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.8290303) q[1];
sx q[1];
rz(-1.0110657) q[1];
sx q[1];
rz(2.1055431) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0890555) q[3];
sx q[3];
rz(-1.6828487) q[3];
sx q[3];
rz(2.5555536) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.6392886) q[2];
sx q[2];
rz(-2.2059811) q[2];
sx q[2];
rz(-1.7712234) q[2];
rz(-0.227452) q[3];
sx q[3];
rz(-1.8919614) q[3];
sx q[3];
rz(1.9500218) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74533904) q[0];
sx q[0];
rz(-0.17530137) q[0];
sx q[0];
rz(1.9434209) q[0];
rz(2.0612969) q[1];
sx q[1];
rz(-2.9273169) q[1];
sx q[1];
rz(-0.094873039) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32786518) q[0];
sx q[0];
rz(-1.828232) q[0];
sx q[0];
rz(-1.0125005) q[0];
rz(-pi) q[1];
rz(-2.2302365) q[2];
sx q[2];
rz(-1.9007287) q[2];
sx q[2];
rz(2.099769) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.97040365) q[1];
sx q[1];
rz(-2.6034149) q[1];
sx q[1];
rz(-1.6506877) q[1];
rz(2.8646144) q[3];
sx q[3];
rz(-1.9414895) q[3];
sx q[3];
rz(2.1408368) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0692811) q[2];
sx q[2];
rz(-2.5371964) q[2];
sx q[2];
rz(-2.7351232) q[2];
rz(1.1786002) q[3];
sx q[3];
rz(-1.4402163) q[3];
sx q[3];
rz(1.9788205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-2.5793295) q[0];
sx q[0];
rz(-0.13450204) q[0];
sx q[0];
rz(0.33600268) q[0];
rz(-0.52040368) q[1];
sx q[1];
rz(-2.2774179) q[1];
sx q[1];
rz(3.0373108) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0020427) q[0];
sx q[0];
rz(-0.487953) q[0];
sx q[0];
rz(-0.81172184) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5849131) q[2];
sx q[2];
rz(-1.4143362) q[2];
sx q[2];
rz(-0.9597646) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.6426864) q[1];
sx q[1];
rz(-2.3489174) q[1];
sx q[1];
rz(-1.0426056) q[1];
rz(-0.4533259) q[3];
sx q[3];
rz(-1.4651872) q[3];
sx q[3];
rz(0.80704122) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.61863724) q[2];
sx q[2];
rz(-1.0982265) q[2];
sx q[2];
rz(-1.1445507) q[2];
rz(0.54730225) q[3];
sx q[3];
rz(-1.2283044) q[3];
sx q[3];
rz(-1.087629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-0.8892141) q[0];
sx q[0];
rz(-0.19110282) q[0];
sx q[0];
rz(0.55437535) q[0];
rz(2.2303708) q[1];
sx q[1];
rz(-1.2634042) q[1];
sx q[1];
rz(0.30141452) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3867971) q[0];
sx q[0];
rz(-0.63969958) q[0];
sx q[0];
rz(1.9596837) q[0];
x q[1];
rz(0.38827916) q[2];
sx q[2];
rz(-1.7634321) q[2];
sx q[2];
rz(1.1757869) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5382945) q[1];
sx q[1];
rz(-2.7009101) q[1];
sx q[1];
rz(1.0540109) q[1];
rz(-pi) q[2];
x q[2];
rz(0.31378515) q[3];
sx q[3];
rz(-1.9289534) q[3];
sx q[3];
rz(-1.4729983) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.1201799) q[2];
sx q[2];
rz(-2.9605949) q[2];
sx q[2];
rz(3.1346698) q[2];
rz(-0.93112469) q[3];
sx q[3];
rz(-2.0102863) q[3];
sx q[3];
rz(-1.1091703) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6329426) q[0];
sx q[0];
rz(-0.83704346) q[0];
sx q[0];
rz(-1.5337926) q[0];
rz(-2.4389229) q[1];
sx q[1];
rz(-2.6103795) q[1];
sx q[1];
rz(-0.30219561) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.142573) q[0];
sx q[0];
rz(-2.4732865) q[0];
sx q[0];
rz(1.4060941) q[0];
rz(1.2185368) q[2];
sx q[2];
rz(-0.88636905) q[2];
sx q[2];
rz(2.5619439) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.2875048) q[1];
sx q[1];
rz(-2.4937428) q[1];
sx q[1];
rz(1.8604398) q[1];
rz(-pi) q[2];
x q[2];
rz(0.052107776) q[3];
sx q[3];
rz(-1.4192389) q[3];
sx q[3];
rz(1.0406756) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.2493784) q[2];
sx q[2];
rz(-1.0503146) q[2];
sx q[2];
rz(2.2155679) q[2];
rz(-1.9269491) q[3];
sx q[3];
rz(-0.49321431) q[3];
sx q[3];
rz(-2.8652969) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72325426) q[0];
sx q[0];
rz(-0.76699081) q[0];
sx q[0];
rz(-1.1085229) q[0];
rz(0.33379894) q[1];
sx q[1];
rz(-1.8148986) q[1];
sx q[1];
rz(1.0858067) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.007581) q[0];
sx q[0];
rz(-1.7756878) q[0];
sx q[0];
rz(-1.4896859) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9606107) q[2];
sx q[2];
rz(-1.9662734) q[2];
sx q[2];
rz(2.4534695) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.27583308) q[1];
sx q[1];
rz(-2.1001022) q[1];
sx q[1];
rz(0.60739002) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7560739) q[3];
sx q[3];
rz(-2.4967125) q[3];
sx q[3];
rz(1.4035937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.6256025) q[2];
sx q[2];
rz(-2.7770999) q[2];
sx q[2];
rz(1.1642574) q[2];
rz(0.26816756) q[3];
sx q[3];
rz(-1.2428913) q[3];
sx q[3];
rz(-2.3837762) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5371573) q[0];
sx q[0];
rz(-2.6217961) q[0];
sx q[0];
rz(-1.9926158) q[0];
rz(0.2941429) q[1];
sx q[1];
rz(-1.5831169) q[1];
sx q[1];
rz(-1.0424967) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9669078) q[0];
sx q[0];
rz(-2.3174441) q[0];
sx q[0];
rz(2.0281726) q[0];
rz(2.3518219) q[2];
sx q[2];
rz(-1.4068479) q[2];
sx q[2];
rz(-2.3254834) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.6175872) q[1];
sx q[1];
rz(-1.7442787) q[1];
sx q[1];
rz(-0.37131997) q[1];
rz(-pi) q[2];
rz(-0.028826272) q[3];
sx q[3];
rz(-2.2346063) q[3];
sx q[3];
rz(0.038906038) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.5440172) q[2];
sx q[2];
rz(-0.67108265) q[2];
sx q[2];
rz(1.3647122) q[2];
rz(-2.4890066) q[3];
sx q[3];
rz(-2.3502974) q[3];
sx q[3];
rz(2.4932388) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6259916) q[0];
sx q[0];
rz(-0.37726548) q[0];
sx q[0];
rz(-0.01734497) q[0];
rz(-3.1248202) q[1];
sx q[1];
rz(-2.3544632) q[1];
sx q[1];
rz(2.9877072) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0506929) q[0];
sx q[0];
rz(-2.8150301) q[0];
sx q[0];
rz(-2.9039277) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4520313) q[2];
sx q[2];
rz(-2.29792) q[2];
sx q[2];
rz(-0.52983701) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.9388403) q[1];
sx q[1];
rz(-1.8649264) q[1];
sx q[1];
rz(-2.7619656) q[1];
rz(1.2653874) q[3];
sx q[3];
rz(-3.0620259) q[3];
sx q[3];
rz(2.9935392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.953557) q[2];
sx q[2];
rz(-0.81622684) q[2];
sx q[2];
rz(-2.6970862) q[2];
rz(-0.19017531) q[3];
sx q[3];
rz(-1.5580274) q[3];
sx q[3];
rz(-1.7688513) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(-1.001215) q[0];
sx q[0];
rz(-1.2495406) q[0];
sx q[0];
rz(-2.0709399) q[0];
rz(0.045914687) q[1];
sx q[1];
rz(-1.670198) q[1];
sx q[1];
rz(-2.3505223) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5673253) q[0];
sx q[0];
rz(-1.2358032) q[0];
sx q[0];
rz(-0.092760249) q[0];
rz(-2.7439762) q[2];
sx q[2];
rz(-1.8194345) q[2];
sx q[2];
rz(1.4225117) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.2998878) q[1];
sx q[1];
rz(-1.1763402) q[1];
sx q[1];
rz(1.9576555) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.13274712) q[3];
sx q[3];
rz(-2.0721966) q[3];
sx q[3];
rz(2.8341334) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.644824) q[2];
sx q[2];
rz(-2.9569929) q[2];
sx q[2];
rz(1.4727288) q[2];
rz(-1.5276927) q[3];
sx q[3];
rz(-1.0563285) q[3];
sx q[3];
rz(1.24019) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2033757) q[0];
sx q[0];
rz(-1.4980409) q[0];
sx q[0];
rz(-0.71312755) q[0];
rz(-2.6279502) q[1];
sx q[1];
rz(-1.7084264) q[1];
sx q[1];
rz(1.4485566) q[1];
rz(-1.8088874) q[2];
sx q[2];
rz(-2.0621962) q[2];
sx q[2];
rz(-0.65218492) q[2];
rz(2.0075825) q[3];
sx q[3];
rz(-2.3306049) q[3];
sx q[3];
rz(-2.8931531) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
