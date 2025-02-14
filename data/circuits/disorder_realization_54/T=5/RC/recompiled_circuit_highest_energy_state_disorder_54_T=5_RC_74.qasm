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
rz(0.74066585) q[0];
sx q[0];
rz(-1.292115) q[0];
sx q[0];
rz(2.3722755) q[0];
rz(-1.5462592) q[1];
sx q[1];
rz(-0.21472628) q[1];
sx q[1];
rz(0.62155849) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7249994) q[0];
sx q[0];
rz(-2.160797) q[0];
sx q[0];
rz(-0.9585462) q[0];
x q[1];
rz(1.5311422) q[2];
sx q[2];
rz(-0.41799212) q[2];
sx q[2];
rz(0.57340535) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.5947508) q[1];
sx q[1];
rz(-1.2870645) q[1];
sx q[1];
rz(0.8304412) q[1];
x q[2];
rz(1.4906472) q[3];
sx q[3];
rz(-1.4372131) q[3];
sx q[3];
rz(-2.7770673) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2810716) q[2];
sx q[2];
rz(-0.29025429) q[2];
sx q[2];
rz(2.2104635) q[2];
rz(2.9562601) q[3];
sx q[3];
rz(-1.0328181) q[3];
sx q[3];
rz(1.4144271) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0952048) q[0];
sx q[0];
rz(-0.33400184) q[0];
sx q[0];
rz(-0.97292501) q[0];
rz(-2.3880549) q[1];
sx q[1];
rz(-2.2851508) q[1];
sx q[1];
rz(0.10261745) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.590693) q[0];
sx q[0];
rz(-1.4859938) q[0];
sx q[0];
rz(-1.4728611) q[0];
rz(-pi) q[1];
rz(0.13394103) q[2];
sx q[2];
rz(-1.6689577) q[2];
sx q[2];
rz(2.5636473) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.2778138) q[1];
sx q[1];
rz(-1.4881388) q[1];
sx q[1];
rz(1.7277) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9063977) q[3];
sx q[3];
rz(-1.7017659) q[3];
sx q[3];
rz(-1.5168911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7332581) q[2];
sx q[2];
rz(-1.4380941) q[2];
sx q[2];
rz(-0.56780887) q[2];
rz(-0.37219498) q[3];
sx q[3];
rz(-1.8122383) q[3];
sx q[3];
rz(1.6224434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5837625) q[0];
sx q[0];
rz(-2.3985641) q[0];
sx q[0];
rz(-1.8999735) q[0];
rz(-2.325233) q[1];
sx q[1];
rz(-2.7528449) q[1];
sx q[1];
rz(-2.1519318) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0123204) q[0];
sx q[0];
rz(-0.43827692) q[0];
sx q[0];
rz(2.6196036) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0895592) q[2];
sx q[2];
rz(-0.18480572) q[2];
sx q[2];
rz(1.081274) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.4892927) q[1];
sx q[1];
rz(-1.1649141) q[1];
sx q[1];
rz(-1.7651221) q[1];
x q[2];
rz(-1.196546) q[3];
sx q[3];
rz(-1.6658226) q[3];
sx q[3];
rz(2.414741) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.376754) q[2];
sx q[2];
rz(-1.2412485) q[2];
sx q[2];
rz(1.494701) q[2];
rz(-2.3946848) q[3];
sx q[3];
rz(-0.61498314) q[3];
sx q[3];
rz(-2.1493105) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84585369) q[0];
sx q[0];
rz(-2.5339412) q[0];
sx q[0];
rz(2.1920152) q[0];
rz(-1.9751366) q[1];
sx q[1];
rz(-1.2587073) q[1];
sx q[1];
rz(-0.20225254) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1631781) q[0];
sx q[0];
rz(-1.5401773) q[0];
sx q[0];
rz(-3.091673) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4522547) q[2];
sx q[2];
rz(-1.7007014) q[2];
sx q[2];
rz(2.4914233) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.7864322) q[1];
sx q[1];
rz(-2.0789685) q[1];
sx q[1];
rz(-2.2340005) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1465591) q[3];
sx q[3];
rz(-2.0821794) q[3];
sx q[3];
rz(0.56727876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.64968455) q[2];
sx q[2];
rz(-1.4338355) q[2];
sx q[2];
rz(1.3362308) q[2];
rz(2.6797471) q[3];
sx q[3];
rz(-1.0577842) q[3];
sx q[3];
rz(1.0285146) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.860054) q[0];
sx q[0];
rz(-0.10801948) q[0];
sx q[0];
rz(2.7463013) q[0];
rz(-2.1832502) q[1];
sx q[1];
rz(-1.7610565) q[1];
sx q[1];
rz(-0.39452943) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0392691) q[0];
sx q[0];
rz(-1.8365055) q[0];
sx q[0];
rz(0.018331176) q[0];
x q[1];
rz(0.5455234) q[2];
sx q[2];
rz(-2.6370271) q[2];
sx q[2];
rz(-1.2009753) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.9745512) q[1];
sx q[1];
rz(-2.2675505) q[1];
sx q[1];
rz(2.7147994) q[1];
rz(-pi) q[2];
rz(-2.5234918) q[3];
sx q[3];
rz(-1.7458946) q[3];
sx q[3];
rz(-1.3503981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.40256527) q[2];
sx q[2];
rz(-1.6141011) q[2];
sx q[2];
rz(-1.019545) q[2];
rz(1.3022276) q[3];
sx q[3];
rz(-3.0890833) q[3];
sx q[3];
rz(2.4359865) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5810982) q[0];
sx q[0];
rz(-0.56466931) q[0];
sx q[0];
rz(-1.0127006) q[0];
rz(-0.46214354) q[1];
sx q[1];
rz(-1.4116762) q[1];
sx q[1];
rz(0.046646811) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29829866) q[0];
sx q[0];
rz(-0.12353914) q[0];
sx q[0];
rz(-1.6575097) q[0];
x q[1];
rz(2.88575) q[2];
sx q[2];
rz(-1.8579972) q[2];
sx q[2];
rz(-0.89509667) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.6973027) q[1];
sx q[1];
rz(-1.238916) q[1];
sx q[1];
rz(0.96849558) q[1];
rz(-pi) q[2];
rz(-0.024123419) q[3];
sx q[3];
rz(-1.1937448) q[3];
sx q[3];
rz(-2.0029062) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8320273) q[2];
sx q[2];
rz(-0.088531606) q[2];
sx q[2];
rz(0.95477611) q[2];
rz(-0.66983062) q[3];
sx q[3];
rz(-1.1845183) q[3];
sx q[3];
rz(0.36439103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(1.9949812) q[0];
sx q[0];
rz(-2.8611188) q[0];
sx q[0];
rz(1.9914419) q[0];
rz(-3.0448044) q[1];
sx q[1];
rz(-2.0014747) q[1];
sx q[1];
rz(1.6614301) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1065001) q[0];
sx q[0];
rz(-2.273002) q[0];
sx q[0];
rz(-1.0649135) q[0];
rz(0.63235967) q[2];
sx q[2];
rz(-0.39146921) q[2];
sx q[2];
rz(-2.7923358) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4188618) q[1];
sx q[1];
rz(-1.7081385) q[1];
sx q[1];
rz(1.2773499) q[1];
x q[2];
rz(-1.893793) q[3];
sx q[3];
rz(-0.45489254) q[3];
sx q[3];
rz(0.98857075) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.4074771) q[2];
sx q[2];
rz(-0.87782562) q[2];
sx q[2];
rz(0.62848282) q[2];
rz(2.8001522) q[3];
sx q[3];
rz(-2.1554558) q[3];
sx q[3];
rz(2.3054874) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(1.2899365) q[0];
sx q[0];
rz(-0.66806942) q[0];
sx q[0];
rz(0.063902721) q[0];
rz(2.5996161) q[1];
sx q[1];
rz(-2.0109476) q[1];
sx q[1];
rz(2.7625387) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0444894) q[0];
sx q[0];
rz(-1.7428723) q[0];
sx q[0];
rz(-1.8015566) q[0];
rz(0.13001534) q[2];
sx q[2];
rz(-2.3283544) q[2];
sx q[2];
rz(0.38947916) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.8547092) q[1];
sx q[1];
rz(-1.5552551) q[1];
sx q[1];
rz(2.5980224) q[1];
rz(-pi) q[2];
rz(0.98446357) q[3];
sx q[3];
rz(-1.2098312) q[3];
sx q[3];
rz(2.933055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.9717676) q[2];
sx q[2];
rz(-1.1173893) q[2];
sx q[2];
rz(-2.8525412) q[2];
rz(-2.1257832) q[3];
sx q[3];
rz(-2.2063875) q[3];
sx q[3];
rz(-2.9142006) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0357901) q[0];
sx q[0];
rz(-0.38689026) q[0];
sx q[0];
rz(0.33045688) q[0];
rz(-0.48078787) q[1];
sx q[1];
rz(-1.5590706) q[1];
sx q[1];
rz(2.6343583) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0704267) q[0];
sx q[0];
rz(-2.2615763) q[0];
sx q[0];
rz(-0.69657495) q[0];
rz(-pi) q[1];
rz(2.6823235) q[2];
sx q[2];
rz(-1.9063564) q[2];
sx q[2];
rz(-2.2989049) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5937336) q[1];
sx q[1];
rz(-0.47858176) q[1];
sx q[1];
rz(-1.0793769) q[1];
rz(-0.60023579) q[3];
sx q[3];
rz(-1.5215989) q[3];
sx q[3];
rz(0.95358301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.5759739) q[2];
sx q[2];
rz(-2.2810292) q[2];
sx q[2];
rz(-2.7536075) q[2];
rz(3.0787789) q[3];
sx q[3];
rz(-2.4020436) q[3];
sx q[3];
rz(-1.1886103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(2.0906618) q[0];
sx q[0];
rz(-0.020314038) q[0];
sx q[0];
rz(-0.055572979) q[0];
rz(0.22124258) q[1];
sx q[1];
rz(-1.0098927) q[1];
sx q[1];
rz(1.6411068) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3474172) q[0];
sx q[0];
rz(-0.58484924) q[0];
sx q[0];
rz(1.0561159) q[0];
rz(-pi) q[1];
rz(-1.6204206) q[2];
sx q[2];
rz(-1.7073586) q[2];
sx q[2];
rz(2.905576) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.4909497) q[1];
sx q[1];
rz(-1.2526647) q[1];
sx q[1];
rz(1.9716119) q[1];
x q[2];
rz(2.5440823) q[3];
sx q[3];
rz(-2.198285) q[3];
sx q[3];
rz(-1.5016055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.9556433) q[2];
sx q[2];
rz(-2.9538437) q[2];
sx q[2];
rz(-1.1090247) q[2];
rz(0.016131314) q[3];
sx q[3];
rz(-1.707209) q[3];
sx q[3];
rz(-0.46512887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7598509) q[0];
sx q[0];
rz(-2.1325337) q[0];
sx q[0];
rz(2.7251563) q[0];
rz(2.3460559) q[1];
sx q[1];
rz(-1.7557314) q[1];
sx q[1];
rz(3.0605127) q[1];
rz(-0.14079413) q[2];
sx q[2];
rz(-1.9119605) q[2];
sx q[2];
rz(-0.20902363) q[2];
rz(-1.0352739) q[3];
sx q[3];
rz(-2.2761619) q[3];
sx q[3];
rz(-0.63009562) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
