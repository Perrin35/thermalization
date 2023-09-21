OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.74270785) q[0];
sx q[0];
rz(2.3595915) q[0];
sx q[0];
rz(10.695988) q[0];
rz(3.4186163) q[1];
sx q[1];
rz(3.613598) q[1];
sx q[1];
rz(9.4233905) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.815925) q[0];
sx q[0];
rz(-0.40580931) q[0];
sx q[0];
rz(-1.9947467) q[0];
x q[1];
rz(0.087287993) q[2];
sx q[2];
rz(-0.44863551) q[2];
sx q[2];
rz(-2.0729614) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5056155) q[1];
sx q[1];
rz(-2.639289) q[1];
sx q[1];
rz(2.99519) q[1];
rz(-pi) q[2];
rz(-1.9418342) q[3];
sx q[3];
rz(-1.3317809) q[3];
sx q[3];
rz(-2.1109076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.9871621) q[2];
sx q[2];
rz(-0.61750948) q[2];
sx q[2];
rz(0.74938613) q[2];
rz(1.0162214) q[3];
sx q[3];
rz(-1.1775492) q[3];
sx q[3];
rz(-0.40482503) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7063023) q[0];
sx q[0];
rz(-2.3162233) q[0];
sx q[0];
rz(2.170927) q[0];
rz(-1.0372112) q[1];
sx q[1];
rz(-1.4379921) q[1];
sx q[1];
rz(-0.81545365) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6204651) q[0];
sx q[0];
rz(-1.6308115) q[0];
sx q[0];
rz(0.64882664) q[0];
x q[1];
rz(-1.843812) q[2];
sx q[2];
rz(-2.2815939) q[2];
sx q[2];
rz(-3.0836011) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.7533469) q[1];
sx q[1];
rz(-1.8716295) q[1];
sx q[1];
rz(1.5335598) q[1];
rz(0.29173298) q[3];
sx q[3];
rz(-2.4749304) q[3];
sx q[3];
rz(-0.093689703) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.6796391) q[2];
sx q[2];
rz(-1.5660428) q[2];
sx q[2];
rz(-2.5088076) q[2];
rz(1.9880382) q[3];
sx q[3];
rz(-2.3735235) q[3];
sx q[3];
rz(2.8320584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84045029) q[0];
sx q[0];
rz(-1.9294894) q[0];
sx q[0];
rz(-0.87483037) q[0];
rz(1.8114999) q[1];
sx q[1];
rz(-1.7069838) q[1];
sx q[1];
rz(0.99951807) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2886575) q[0];
sx q[0];
rz(-1.7203727) q[0];
sx q[0];
rz(1.8729217) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2314084) q[2];
sx q[2];
rz(-2.0816457) q[2];
sx q[2];
rz(0.78188932) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.42400186) q[1];
sx q[1];
rz(-1.0832936) q[1];
sx q[1];
rz(1.8600149) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0796892) q[3];
sx q[3];
rz(-1.6403927) q[3];
sx q[3];
rz(-0.059046179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.53753608) q[2];
sx q[2];
rz(-0.92210046) q[2];
sx q[2];
rz(-0.58004722) q[2];
rz(2.3245658) q[3];
sx q[3];
rz(-1.3823119) q[3];
sx q[3];
rz(1.9918611) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76628768) q[0];
sx q[0];
rz(-1.5910609) q[0];
sx q[0];
rz(-2.2312009) q[0];
rz(0.45122775) q[1];
sx q[1];
rz(-1.5463566) q[1];
sx q[1];
rz(-2.8667563) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7461473) q[0];
sx q[0];
rz(-1.6593885) q[0];
sx q[0];
rz(-3.0625312) q[0];
x q[1];
rz(2.0312112) q[2];
sx q[2];
rz(-1.5722256) q[2];
sx q[2];
rz(-1.187385) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.53510016) q[1];
sx q[1];
rz(-2.1537158) q[1];
sx q[1];
rz(0.54775723) q[1];
rz(1.7219909) q[3];
sx q[3];
rz(-2.4393743) q[3];
sx q[3];
rz(2.7694626) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.3466907) q[2];
sx q[2];
rz(-1.1179504) q[2];
sx q[2];
rz(-1.6332731) q[2];
rz(1.1446965) q[3];
sx q[3];
rz(-2.4016524) q[3];
sx q[3];
rz(0.16170734) q[3];
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
rz(0.65790025) q[0];
sx q[0];
rz(-1.9265441) q[0];
sx q[0];
rz(2.143798) q[0];
rz(-0.18355852) q[1];
sx q[1];
rz(-1.654637) q[1];
sx q[1];
rz(1.516974) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8109587) q[0];
sx q[0];
rz(-1.9348382) q[0];
sx q[0];
rz(-0.90322687) q[0];
rz(-pi) q[1];
rz(-2.5324608) q[2];
sx q[2];
rz(-2.5172148) q[2];
sx q[2];
rz(2.6513211) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.5114054) q[1];
sx q[1];
rz(-1.1668219) q[1];
sx q[1];
rz(-1.1955111) q[1];
rz(-pi) q[2];
x q[2];
rz(2.153271) q[3];
sx q[3];
rz(-1.054793) q[3];
sx q[3];
rz(0.9048681) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.8020442) q[2];
sx q[2];
rz(-0.99207726) q[2];
sx q[2];
rz(2.8175763) q[2];
rz(1.8185395) q[3];
sx q[3];
rz(-0.75606212) q[3];
sx q[3];
rz(1.6103305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3762387) q[0];
sx q[0];
rz(-2.1026251) q[0];
sx q[0];
rz(-1.2639686) q[0];
rz(2.2309247) q[1];
sx q[1];
rz(-1.2011386) q[1];
sx q[1];
rz(2.8009159) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7621988) q[0];
sx q[0];
rz(-1.6025935) q[0];
sx q[0];
rz(1.6367153) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2773877) q[2];
sx q[2];
rz(-0.83653203) q[2];
sx q[2];
rz(-0.66213995) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.1250455) q[1];
sx q[1];
rz(-0.53495896) q[1];
sx q[1];
rz(-0.46334456) q[1];
rz(-pi) q[2];
rz(2.2716899) q[3];
sx q[3];
rz(-0.60855908) q[3];
sx q[3];
rz(-2.9121947) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1910151) q[2];
sx q[2];
rz(-2.5577929) q[2];
sx q[2];
rz(-0.77159709) q[2];
rz(-0.54780444) q[3];
sx q[3];
rz(-2.1717725) q[3];
sx q[3];
rz(2.5781393) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9724378) q[0];
sx q[0];
rz(-1.4470402) q[0];
sx q[0];
rz(-0.34564885) q[0];
rz(0.06282839) q[1];
sx q[1];
rz(-2.6627916) q[1];
sx q[1];
rz(2.6766434) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9879887) q[0];
sx q[0];
rz(-1.1055595) q[0];
sx q[0];
rz(2.9265755) q[0];
rz(0.49352383) q[2];
sx q[2];
rz(-2.0888622) q[2];
sx q[2];
rz(1.354419) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.3685776) q[1];
sx q[1];
rz(-1.7891208) q[1];
sx q[1];
rz(-1.6792084) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.025860272) q[3];
sx q[3];
rz(-1.6169294) q[3];
sx q[3];
rz(2.9010454) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.1376301) q[2];
sx q[2];
rz(-1.6370862) q[2];
sx q[2];
rz(0.31759343) q[2];
rz(-2.5701304) q[3];
sx q[3];
rz(-1.0390037) q[3];
sx q[3];
rz(2.8542744) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6222318) q[0];
sx q[0];
rz(-1.2943635) q[0];
sx q[0];
rz(-0.28433329) q[0];
rz(0.55150664) q[1];
sx q[1];
rz(-2.9998144) q[1];
sx q[1];
rz(-0.078358738) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6358444) q[0];
sx q[0];
rz(-2.3730179) q[0];
sx q[0];
rz(1.6705253) q[0];
x q[1];
rz(-0.23711726) q[2];
sx q[2];
rz(-1.3318921) q[2];
sx q[2];
rz(0.6616) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.2625418) q[1];
sx q[1];
rz(-1.7444376) q[1];
sx q[1];
rz(0.91211984) q[1];
x q[2];
rz(-3.1073242) q[3];
sx q[3];
rz(-1.7274324) q[3];
sx q[3];
rz(-0.53965118) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.4006965) q[2];
sx q[2];
rz(-0.90831465) q[2];
sx q[2];
rz(0.25137869) q[2];
rz(-0.58319432) q[3];
sx q[3];
rz(-2.0299032) q[3];
sx q[3];
rz(-1.4096227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
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
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3437929) q[0];
sx q[0];
rz(-3.0594337) q[0];
sx q[0];
rz(0.051368512) q[0];
rz(-0.92357606) q[1];
sx q[1];
rz(-0.66134614) q[1];
sx q[1];
rz(0.87402469) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.066814518) q[0];
sx q[0];
rz(-1.5621119) q[0];
sx q[0];
rz(-2.3735223) q[0];
rz(-pi) q[1];
rz(0.36058493) q[2];
sx q[2];
rz(-1.9546095) q[2];
sx q[2];
rz(-1.4521445) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2487138) q[1];
sx q[1];
rz(-1.6422179) q[1];
sx q[1];
rz(-2.0715269) q[1];
rz(-pi) q[2];
rz(-1.2285352) q[3];
sx q[3];
rz(-2.1517793) q[3];
sx q[3];
rz(-0.57623219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.727227) q[2];
sx q[2];
rz(-2.3870769) q[2];
sx q[2];
rz(-2.5218463) q[2];
rz(-1.184458) q[3];
sx q[3];
rz(-1.8871566) q[3];
sx q[3];
rz(-1.363389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1353564) q[0];
sx q[0];
rz(-2.0993435) q[0];
sx q[0];
rz(-2.4172879) q[0];
rz(-2.9528217) q[1];
sx q[1];
rz(-2.962208) q[1];
sx q[1];
rz(1.9627409) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41335426) q[0];
sx q[0];
rz(-1.2984707) q[0];
sx q[0];
rz(-2.9288835) q[0];
rz(-1.6640501) q[2];
sx q[2];
rz(-1.6972731) q[2];
sx q[2];
rz(1.7699514) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.2420173) q[1];
sx q[1];
rz(-1.3841108) q[1];
sx q[1];
rz(2.1980397) q[1];
x q[2];
rz(-2.6188649) q[3];
sx q[3];
rz(-1.3739112) q[3];
sx q[3];
rz(2.6678391) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.05802352) q[2];
sx q[2];
rz(-1.0408137) q[2];
sx q[2];
rz(0.89938346) q[2];
rz(2.2670238) q[3];
sx q[3];
rz(-0.42566291) q[3];
sx q[3];
rz(-1.4609059) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5205004) q[0];
sx q[0];
rz(-2.7324471) q[0];
sx q[0];
rz(0.24656217) q[0];
rz(0.75795603) q[1];
sx q[1];
rz(-1.6592204) q[1];
sx q[1];
rz(1.6827676) q[1];
rz(-0.035049546) q[2];
sx q[2];
rz(-2.0106342) q[2];
sx q[2];
rz(1.0019279) q[2];
rz(-0.27579565) q[3];
sx q[3];
rz(-1.7134568) q[3];
sx q[3];
rz(-1.4207763) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];