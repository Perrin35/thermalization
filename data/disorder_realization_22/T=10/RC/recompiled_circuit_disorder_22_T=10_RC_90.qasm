OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.7093631) q[0];
sx q[0];
rz(-2.1837283) q[0];
sx q[0];
rz(-0.14444484) q[0];
rz(0.56675178) q[1];
sx q[1];
rz(-0.52539879) q[1];
sx q[1];
rz(0.97775835) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35740556) q[0];
sx q[0];
rz(-1.4854684) q[0];
sx q[0];
rz(0.25934319) q[0];
rz(-pi) q[1];
rz(-1.058504) q[2];
sx q[2];
rz(-0.57527486) q[2];
sx q[2];
rz(-3.1303867) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.2797151) q[1];
sx q[1];
rz(-1.4289083) q[1];
sx q[1];
rz(-3.0008297) q[1];
x q[2];
rz(0.17840673) q[3];
sx q[3];
rz(-2.3462786) q[3];
sx q[3];
rz(2.2291396) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.67291659) q[2];
sx q[2];
rz(-1.1489931) q[2];
sx q[2];
rz(0.93227512) q[2];
rz(0.19876924) q[3];
sx q[3];
rz(-1.1105744) q[3];
sx q[3];
rz(-0.96536243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7070049) q[0];
sx q[0];
rz(-2.2362237) q[0];
sx q[0];
rz(0.36112753) q[0];
rz(1.7065642) q[1];
sx q[1];
rz(-1.3577434) q[1];
sx q[1];
rz(0.8180058) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6011071) q[0];
sx q[0];
rz(-1.846082) q[0];
sx q[0];
rz(2.7406373) q[0];
rz(-3.0188574) q[2];
sx q[2];
rz(-1.2018179) q[2];
sx q[2];
rz(-0.69532794) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2745143) q[1];
sx q[1];
rz(-2.3550905) q[1];
sx q[1];
rz(2.2134476) q[1];
rz(-pi) q[2];
rz(-1.7198256) q[3];
sx q[3];
rz(-2.7741198) q[3];
sx q[3];
rz(3.0847103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.8890185) q[2];
sx q[2];
rz(-2.694464) q[2];
sx q[2];
rz(-1.4206295) q[2];
rz(1.3160926) q[3];
sx q[3];
rz(-2.3836453) q[3];
sx q[3];
rz(-2.7533598) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4017568) q[0];
sx q[0];
rz(-0.49210423) q[0];
sx q[0];
rz(2.2706568) q[0];
rz(0.31618205) q[1];
sx q[1];
rz(-0.28156391) q[1];
sx q[1];
rz(0.2972163) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7446211) q[0];
sx q[0];
rz(-2.731785) q[0];
sx q[0];
rz(2.9817392) q[0];
rz(-pi) q[1];
rz(-2.259953) q[2];
sx q[2];
rz(-2.214553) q[2];
sx q[2];
rz(-1.2780485) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.19993648) q[1];
sx q[1];
rz(-1.3841277) q[1];
sx q[1];
rz(-0.048528683) q[1];
x q[2];
rz(2.3171595) q[3];
sx q[3];
rz(-1.0472877) q[3];
sx q[3];
rz(1.5200966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3729942) q[2];
sx q[2];
rz(-2.3186389) q[2];
sx q[2];
rz(0.42759839) q[2];
rz(1.1887431) q[3];
sx q[3];
rz(-2.5116428) q[3];
sx q[3];
rz(0.52743131) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7052085) q[0];
sx q[0];
rz(-1.5773062) q[0];
sx q[0];
rz(-2.3676681) q[0];
rz(-2.4286843) q[1];
sx q[1];
rz(-2.1122825) q[1];
sx q[1];
rz(-0.68177044) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43142372) q[0];
sx q[0];
rz(-1.7792601) q[0];
sx q[0];
rz(-1.8375977) q[0];
rz(-pi) q[1];
rz(-2.1749928) q[2];
sx q[2];
rz(-1.3912429) q[2];
sx q[2];
rz(-2.0008848) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.75979739) q[1];
sx q[1];
rz(-2.0054381) q[1];
sx q[1];
rz(-2.1684907) q[1];
x q[2];
rz(-2.3247129) q[3];
sx q[3];
rz(-0.48135346) q[3];
sx q[3];
rz(0.64827418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.008808) q[2];
sx q[2];
rz(-0.71295732) q[2];
sx q[2];
rz(2.183389) q[2];
rz(-1.0664252) q[3];
sx q[3];
rz(-1.8582148) q[3];
sx q[3];
rz(0.3796033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.585007) q[0];
sx q[0];
rz(-1.3469232) q[0];
sx q[0];
rz(0.56030309) q[0];
rz(-2.141748) q[1];
sx q[1];
rz(-0.20345774) q[1];
sx q[1];
rz(-1.5195742) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9902089) q[0];
sx q[0];
rz(-2.4575893) q[0];
sx q[0];
rz(-2.3955976) q[0];
rz(-2.4155596) q[2];
sx q[2];
rz(-2.3880929) q[2];
sx q[2];
rz(1.3940648) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.49415627) q[1];
sx q[1];
rz(-1.372822) q[1];
sx q[1];
rz(-0.76311771) q[1];
x q[2];
rz(1.7468466) q[3];
sx q[3];
rz(-2.5213443) q[3];
sx q[3];
rz(0.32905096) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.4524298) q[2];
sx q[2];
rz(-0.55183691) q[2];
sx q[2];
rz(-0.2229283) q[2];
rz(0.034742268) q[3];
sx q[3];
rz(-1.3870753) q[3];
sx q[3];
rz(-3.070014) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3361622) q[0];
sx q[0];
rz(-2.8283089) q[0];
sx q[0];
rz(-1.0700595) q[0];
rz(-1.3609715) q[1];
sx q[1];
rz(-0.37934163) q[1];
sx q[1];
rz(-1.8575352) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7362471) q[0];
sx q[0];
rz(-2.1655472) q[0];
sx q[0];
rz(-0.88881641) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.96096303) q[2];
sx q[2];
rz(-2.6764538) q[2];
sx q[2];
rz(-2.7679408) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.9575427) q[1];
sx q[1];
rz(-1.8584538) q[1];
sx q[1];
rz(-0.48344739) q[1];
x q[2];
rz(2.0606023) q[3];
sx q[3];
rz(-2.5488857) q[3];
sx q[3];
rz(2.5268775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.6665035) q[2];
sx q[2];
rz(-1.6978426) q[2];
sx q[2];
rz(2.5040023) q[2];
rz(-2.3049138) q[3];
sx q[3];
rz(-2.0357318) q[3];
sx q[3];
rz(-1.3440514) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56629431) q[0];
sx q[0];
rz(-0.42995444) q[0];
sx q[0];
rz(0.56754011) q[0];
rz(2.7138846) q[1];
sx q[1];
rz(-1.5274915) q[1];
sx q[1];
rz(-2.2033851) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.095031247) q[0];
sx q[0];
rz(-1.1633658) q[0];
sx q[0];
rz(0.33336063) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0023408) q[2];
sx q[2];
rz(-2.217514) q[2];
sx q[2];
rz(1.2127753) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.92749121) q[1];
sx q[1];
rz(-0.92988211) q[1];
sx q[1];
rz(-2.6919041) q[1];
rz(0.70721831) q[3];
sx q[3];
rz(-1.8054609) q[3];
sx q[3];
rz(-0.16682391) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.2618835) q[2];
sx q[2];
rz(-1.1738913) q[2];
sx q[2];
rz(1.3860469) q[2];
rz(-1.8188247) q[3];
sx q[3];
rz(-1.1498007) q[3];
sx q[3];
rz(0.02903207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8742074) q[0];
sx q[0];
rz(-2.8081942) q[0];
sx q[0];
rz(-1.4338795) q[0];
rz(-1.8677615) q[1];
sx q[1];
rz(-1.1359943) q[1];
sx q[1];
rz(0.83126718) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6600663) q[0];
sx q[0];
rz(-0.43018451) q[0];
sx q[0];
rz(2.7823886) q[0];
rz(-pi) q[1];
rz(1.3966884) q[2];
sx q[2];
rz(-0.6422407) q[2];
sx q[2];
rz(1.0345392) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.9030994) q[1];
sx q[1];
rz(-1.0817173) q[1];
sx q[1];
rz(0.78519435) q[1];
rz(-pi) q[2];
rz(0.23037489) q[3];
sx q[3];
rz(-1.7883693) q[3];
sx q[3];
rz(2.2090467) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.5780118) q[2];
sx q[2];
rz(-1.3728023) q[2];
sx q[2];
rz(-2.1772299) q[2];
rz(-1.1635121) q[3];
sx q[3];
rz(-0.61146277) q[3];
sx q[3];
rz(-0.65892974) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37373856) q[0];
sx q[0];
rz(-2.4066194) q[0];
sx q[0];
rz(0.97737616) q[0];
rz(1.7550229) q[1];
sx q[1];
rz(-1.8354548) q[1];
sx q[1];
rz(1.1057373) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0375263) q[0];
sx q[0];
rz(-0.46188799) q[0];
sx q[0];
rz(-0.26432963) q[0];
rz(2.2374723) q[2];
sx q[2];
rz(-1.4881926) q[2];
sx q[2];
rz(0.77043515) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.22509174) q[1];
sx q[1];
rz(-1.8019925) q[1];
sx q[1];
rz(2.9514312) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2826142) q[3];
sx q[3];
rz(-1.9087221) q[3];
sx q[3];
rz(-1.6328788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.608312) q[2];
sx q[2];
rz(-2.7567342) q[2];
sx q[2];
rz(0.14979714) q[2];
rz(1.3730565) q[3];
sx q[3];
rz(-1.7356197) q[3];
sx q[3];
rz(1.3214553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49366429) q[0];
sx q[0];
rz(-2.6239008) q[0];
sx q[0];
rz(0.1272442) q[0];
rz(-1.4808902) q[1];
sx q[1];
rz(-2.7791185) q[1];
sx q[1];
rz(0.1677992) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0509335) q[0];
sx q[0];
rz(-1.5516084) q[0];
sx q[0];
rz(3.1113383) q[0];
rz(-pi) q[1];
rz(-2.9211505) q[2];
sx q[2];
rz(-0.61969212) q[2];
sx q[2];
rz(-2.9881791) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.4976616) q[1];
sx q[1];
rz(-1.5776458) q[1];
sx q[1];
rz(2.5531406) q[1];
rz(0.76370244) q[3];
sx q[3];
rz(-0.16994952) q[3];
sx q[3];
rz(-0.22596879) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.026713513) q[2];
sx q[2];
rz(-0.9388887) q[2];
sx q[2];
rz(2.3804469) q[2];
rz(-3.051565) q[3];
sx q[3];
rz(-2.138425) q[3];
sx q[3];
rz(-2.1910523) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54820838) q[0];
sx q[0];
rz(-1.1593288) q[0];
sx q[0];
rz(2.819084) q[0];
rz(0.38800115) q[1];
sx q[1];
rz(-1.7419659) q[1];
sx q[1];
rz(2.3566125) q[1];
rz(2.2446752) q[2];
sx q[2];
rz(-2.2581836) q[2];
sx q[2];
rz(2.7473292) q[2];
rz(-1.0317867) q[3];
sx q[3];
rz(-1.7179334) q[3];
sx q[3];
rz(-2.4739305) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
