OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(3.0290282) q[0];
sx q[0];
rz(-1.5052786) q[0];
sx q[0];
rz(2.3583052) q[0];
rz(1.0822436) q[1];
sx q[1];
rz(-1.487027) q[1];
sx q[1];
rz(-2.3496871) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5740252) q[0];
sx q[0];
rz(-1.3114531) q[0];
sx q[0];
rz(-1.3784842) q[0];
rz(-pi) q[1];
rz(1.3375086) q[2];
sx q[2];
rz(-0.15809862) q[2];
sx q[2];
rz(-2.5776517) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.72630771) q[1];
sx q[1];
rz(-1.0276762) q[1];
sx q[1];
rz(1.5341137) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1551829) q[3];
sx q[3];
rz(-2.274868) q[3];
sx q[3];
rz(2.5147284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5504494) q[2];
sx q[2];
rz(-2.998896) q[2];
sx q[2];
rz(-3.0421416) q[2];
rz(2.8948696) q[3];
sx q[3];
rz(-1.6171425) q[3];
sx q[3];
rz(-2.7439086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0196911) q[0];
sx q[0];
rz(-2.1653403) q[0];
sx q[0];
rz(2.2130527) q[0];
rz(-1.4506725) q[1];
sx q[1];
rz(-1.848315) q[1];
sx q[1];
rz(2.4931152) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.649594) q[0];
sx q[0];
rz(-1.0186983) q[0];
sx q[0];
rz(2.2158438) q[0];
rz(-pi) q[1];
rz(1.8869867) q[2];
sx q[2];
rz(-2.11497) q[2];
sx q[2];
rz(-0.7989102) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.8806421) q[1];
sx q[1];
rz(-1.7127258) q[1];
sx q[1];
rz(1.5550343) q[1];
x q[2];
rz(-1.1901598) q[3];
sx q[3];
rz(-2.8320304) q[3];
sx q[3];
rz(0.42282399) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.4497946) q[2];
sx q[2];
rz(-2.3851676) q[2];
sx q[2];
rz(-0.85255426) q[2];
rz(-0.84613386) q[3];
sx q[3];
rz(-1.6328014) q[3];
sx q[3];
rz(2.1107296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6344675) q[0];
sx q[0];
rz(-1.2256624) q[0];
sx q[0];
rz(-1.0852098) q[0];
rz(2.197544) q[1];
sx q[1];
rz(-1.4986821) q[1];
sx q[1];
rz(1.5012213) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0414934) q[0];
sx q[0];
rz(-1.3736808) q[0];
sx q[0];
rz(1.6105349) q[0];
rz(-0.19135059) q[2];
sx q[2];
rz(-1.8008968) q[2];
sx q[2];
rz(-1.5887345) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.2861917) q[1];
sx q[1];
rz(-0.58127379) q[1];
sx q[1];
rz(-1.3894677) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7872058) q[3];
sx q[3];
rz(-1.4167656) q[3];
sx q[3];
rz(1.1104402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7368855) q[2];
sx q[2];
rz(-0.53386226) q[2];
sx q[2];
rz(-2.286818) q[2];
rz(-0.9084304) q[3];
sx q[3];
rz(-1.5646489) q[3];
sx q[3];
rz(2.0250208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92723769) q[0];
sx q[0];
rz(-2.7404009) q[0];
sx q[0];
rz(-1.9450564) q[0];
rz(-1.6664956) q[1];
sx q[1];
rz(-0.42088446) q[1];
sx q[1];
rz(1.9570785) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6086585) q[0];
sx q[0];
rz(-2.2581165) q[0];
sx q[0];
rz(1.88009) q[0];
rz(-pi) q[1];
x q[1];
rz(0.20303161) q[2];
sx q[2];
rz(-2.3782999) q[2];
sx q[2];
rz(-0.53084669) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.2503916) q[1];
sx q[1];
rz(-1.7970835) q[1];
sx q[1];
rz(2.0865284) q[1];
rz(-1.7538025) q[3];
sx q[3];
rz(-1.840045) q[3];
sx q[3];
rz(1.315801) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.92675942) q[2];
sx q[2];
rz(-2.6439809) q[2];
sx q[2];
rz(-1.2614177) q[2];
rz(-0.43241209) q[3];
sx q[3];
rz(-1.132248) q[3];
sx q[3];
rz(2.6692218) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
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
rz(2.0686907) q[0];
sx q[0];
rz(-1.9111159) q[0];
sx q[0];
rz(-0.029408971) q[0];
rz(0.75621653) q[1];
sx q[1];
rz(-0.58719802) q[1];
sx q[1];
rz(0.75278935) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2761951) q[0];
sx q[0];
rz(-1.5704234) q[0];
sx q[0];
rz(-1.5651817) q[0];
rz(1.8389614) q[2];
sx q[2];
rz(-1.6577621) q[2];
sx q[2];
rz(0.11448569) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.636105) q[1];
sx q[1];
rz(-1.4661745) q[1];
sx q[1];
rz(-2.5701447) q[1];
rz(-2.6008368) q[3];
sx q[3];
rz(-1.1782559) q[3];
sx q[3];
rz(-0.3730216) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0011562) q[2];
sx q[2];
rz(-1.6480646) q[2];
sx q[2];
rz(0.39503869) q[2];
rz(2.6664074) q[3];
sx q[3];
rz(-1.7893712) q[3];
sx q[3];
rz(0.37676677) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3756556) q[0];
sx q[0];
rz(-0.7203311) q[0];
sx q[0];
rz(-2.1790867) q[0];
rz(-2.6267701) q[1];
sx q[1];
rz(-1.5341026) q[1];
sx q[1];
rz(-2.8033676) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.048934919) q[0];
sx q[0];
rz(-0.83507628) q[0];
sx q[0];
rz(-0.2661163) q[0];
rz(-pi) q[1];
rz(-2.9043496) q[2];
sx q[2];
rz(-1.3297992) q[2];
sx q[2];
rz(-0.70912305) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.64778642) q[1];
sx q[1];
rz(-1.5739417) q[1];
sx q[1];
rz(2.6814744e-05) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7231483) q[3];
sx q[3];
rz(-2.7151516) q[3];
sx q[3];
rz(1.3063198) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6414791) q[2];
sx q[2];
rz(-1.1621472) q[2];
sx q[2];
rz(-0.55434736) q[2];
rz(-0.85136271) q[3];
sx q[3];
rz(-2.7832289) q[3];
sx q[3];
rz(2.5135777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3250658) q[0];
sx q[0];
rz(-2.5287703) q[0];
sx q[0];
rz(-2.391173) q[0];
rz(-0.57506192) q[1];
sx q[1];
rz(-1.776418) q[1];
sx q[1];
rz(-0.65779984) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0938134) q[0];
sx q[0];
rz(-1.9792611) q[0];
sx q[0];
rz(-0.7042709) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6148241) q[2];
sx q[2];
rz(-0.86404534) q[2];
sx q[2];
rz(-2.6296774) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.9283674) q[1];
sx q[1];
rz(-1.2847273) q[1];
sx q[1];
rz(0.017315344) q[1];
x q[2];
rz(-0.21605394) q[3];
sx q[3];
rz(-1.9441868) q[3];
sx q[3];
rz(1.4114398) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.49360069) q[2];
sx q[2];
rz(-1.0085663) q[2];
sx q[2];
rz(1.6592525) q[2];
rz(3.0209387) q[3];
sx q[3];
rz(-1.3841265) q[3];
sx q[3];
rz(0.83546662) q[3];
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
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3487314) q[0];
sx q[0];
rz(-0.86935765) q[0];
sx q[0];
rz(-2.8435775) q[0];
rz(-2.1030078) q[1];
sx q[1];
rz(-2.5972001) q[1];
sx q[1];
rz(-2.9878152) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5030432) q[0];
sx q[0];
rz(-0.26493236) q[0];
sx q[0];
rz(0.28980906) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5030382) q[2];
sx q[2];
rz(-2.2306135) q[2];
sx q[2];
rz(-1.7121079) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.41339918) q[1];
sx q[1];
rz(-1.4811885) q[1];
sx q[1];
rz(2.7941568) q[1];
rz(-pi) q[2];
rz(1.3078717) q[3];
sx q[3];
rz(-0.44028966) q[3];
sx q[3];
rz(-1.057098) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.095857233) q[2];
sx q[2];
rz(-2.6726674) q[2];
sx q[2];
rz(-1.1886965) q[2];
rz(-1.8111604) q[3];
sx q[3];
rz(-1.1878139) q[3];
sx q[3];
rz(-0.73330283) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34340149) q[0];
sx q[0];
rz(-2.5372086) q[0];
sx q[0];
rz(-1.6424302) q[0];
rz(-0.54939735) q[1];
sx q[1];
rz(-1.5645942) q[1];
sx q[1];
rz(-0.25064358) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9662387) q[0];
sx q[0];
rz(-2.4099775) q[0];
sx q[0];
rz(2.0425955) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0876686) q[2];
sx q[2];
rz(-0.26130518) q[2];
sx q[2];
rz(-3.1045632) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.10560606) q[1];
sx q[1];
rz(-0.79933724) q[1];
sx q[1];
rz(-1.4738333) q[1];
rz(2.9467877) q[3];
sx q[3];
rz(-1.6945632) q[3];
sx q[3];
rz(0.74461246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.15829076) q[2];
sx q[2];
rz(-3.0088708) q[2];
sx q[2];
rz(-1.0164725) q[2];
rz(-3.0554092) q[3];
sx q[3];
rz(-2.1250171) q[3];
sx q[3];
rz(-2.3763954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8330399) q[0];
sx q[0];
rz(-0.85423952) q[0];
sx q[0];
rz(2.7225851) q[0];
rz(-0.53681701) q[1];
sx q[1];
rz(-2.5173126) q[1];
sx q[1];
rz(2.4057665) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6572172) q[0];
sx q[0];
rz(-1.7630944) q[0];
sx q[0];
rz(2.4853287) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1081829) q[2];
sx q[2];
rz(-2.1436286) q[2];
sx q[2];
rz(-2.9705641) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.6169301) q[1];
sx q[1];
rz(-1.3549651) q[1];
sx q[1];
rz(-2.4991577) q[1];
rz(-0.60634585) q[3];
sx q[3];
rz(-0.56713533) q[3];
sx q[3];
rz(0.66886574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.1824823) q[2];
sx q[2];
rz(-2.3190494) q[2];
sx q[2];
rz(-0.62270069) q[2];
rz(-0.73838082) q[3];
sx q[3];
rz(-1.9564956) q[3];
sx q[3];
rz(-2.8625989) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57711346) q[0];
sx q[0];
rz(-0.27272419) q[0];
sx q[0];
rz(0.96584366) q[0];
rz(0.64361698) q[1];
sx q[1];
rz(-1.8543961) q[1];
sx q[1];
rz(-2.054945) q[1];
rz(2.116133) q[2];
sx q[2];
rz(-2.7164216) q[2];
sx q[2];
rz(2.807775) q[2];
rz(0.24091992) q[3];
sx q[3];
rz(-1.708247) q[3];
sx q[3];
rz(0.74549992) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
