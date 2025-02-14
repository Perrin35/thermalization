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
rz(-0.15260881) q[0];
sx q[0];
rz(-1.1996562) q[0];
sx q[0];
rz(-0.3056404) q[0];
rz(-2.8513554) q[1];
sx q[1];
rz(-1.4479535) q[1];
sx q[1];
rz(-1.0040959) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5619156) q[0];
sx q[0];
rz(-1.5439057) q[0];
sx q[0];
rz(3.0843601) q[0];
x q[1];
rz(0.19308108) q[2];
sx q[2];
rz(-1.8154789) q[2];
sx q[2];
rz(0.31665238) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.43875445) q[1];
sx q[1];
rz(-1.251529) q[1];
sx q[1];
rz(1.1478018) q[1];
rz(-pi) q[2];
rz(3.0707487) q[3];
sx q[3];
rz(-2.6016782) q[3];
sx q[3];
rz(1.2666463) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.64757887) q[2];
sx q[2];
rz(-2.5250489) q[2];
sx q[2];
rz(-0.54440633) q[2];
rz(0.38862774) q[3];
sx q[3];
rz(-1.2231239) q[3];
sx q[3];
rz(-0.94038928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87043864) q[0];
sx q[0];
rz(-2.1156023) q[0];
sx q[0];
rz(1.1218659) q[0];
rz(-1.059996) q[1];
sx q[1];
rz(-2.6385939) q[1];
sx q[1];
rz(-0.88567919) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87410418) q[0];
sx q[0];
rz(-2.7120706) q[0];
sx q[0];
rz(-2.3699479) q[0];
rz(-0.76876872) q[2];
sx q[2];
rz(-0.93009206) q[2];
sx q[2];
rz(-0.51962432) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.6180775) q[1];
sx q[1];
rz(-1.6517868) q[1];
sx q[1];
rz(0.063132719) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2996385) q[3];
sx q[3];
rz(-0.15781395) q[3];
sx q[3];
rz(2.006881) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.1765959) q[2];
sx q[2];
rz(-1.4623888) q[2];
sx q[2];
rz(3.1122567) q[2];
rz(0.36597478) q[3];
sx q[3];
rz(-0.47593203) q[3];
sx q[3];
rz(-0.16987814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.733424) q[0];
sx q[0];
rz(-1.8086731) q[0];
sx q[0];
rz(0.15342203) q[0];
rz(0.19639213) q[1];
sx q[1];
rz(-1.6729665) q[1];
sx q[1];
rz(0.82082716) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0401413) q[0];
sx q[0];
rz(-2.2635604) q[0];
sx q[0];
rz(-2.470284) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.94541855) q[2];
sx q[2];
rz(-2.4074005) q[2];
sx q[2];
rz(1.9482833) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.7559636) q[1];
sx q[1];
rz(-1.5308994) q[1];
sx q[1];
rz(-0.41547687) q[1];
rz(-pi) q[2];
rz(3.0140424) q[3];
sx q[3];
rz(-2.7722428) q[3];
sx q[3];
rz(2.7238059) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(3.0241432) q[2];
sx q[2];
rz(-1.3261745) q[2];
sx q[2];
rz(-2.2085025) q[2];
rz(0.34705958) q[3];
sx q[3];
rz(-1.1091899) q[3];
sx q[3];
rz(-0.6412653) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1953122) q[0];
sx q[0];
rz(-2.0144137) q[0];
sx q[0];
rz(-1.0585002) q[0];
rz(-3.0463386) q[1];
sx q[1];
rz(-1.1493827) q[1];
sx q[1];
rz(2.9375295) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6883012) q[0];
sx q[0];
rz(-2.8729962) q[0];
sx q[0];
rz(-2.1902188) q[0];
rz(-pi) q[1];
rz(-1.0957031) q[2];
sx q[2];
rz(-1.9852601) q[2];
sx q[2];
rz(-1.955065) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.82641027) q[1];
sx q[1];
rz(-2.1433804) q[1];
sx q[1];
rz(-1.8099643) q[1];
x q[2];
rz(-2.4434632) q[3];
sx q[3];
rz(-1.7295726) q[3];
sx q[3];
rz(2.862249) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-3.0846587) q[2];
sx q[2];
rz(-1.6997507) q[2];
sx q[2];
rz(1.320768) q[2];
rz(-1.002958) q[3];
sx q[3];
rz(-1.2489677) q[3];
sx q[3];
rz(0.5654208) q[3];
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
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.716662) q[0];
sx q[0];
rz(-1.2053763) q[0];
sx q[0];
rz(2.5004814) q[0];
rz(2.5504316) q[1];
sx q[1];
rz(-1.7008275) q[1];
sx q[1];
rz(-1.6966746) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.145744) q[0];
sx q[0];
rz(-1.543049) q[0];
sx q[0];
rz(-2.7553245) q[0];
rz(-pi) q[1];
x q[1];
rz(0.88052184) q[2];
sx q[2];
rz(-1.4052999) q[2];
sx q[2];
rz(1.6146162) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.29825764) q[1];
sx q[1];
rz(-1.3239256) q[1];
sx q[1];
rz(2.6554537) q[1];
rz(-pi) q[2];
rz(2.0956482) q[3];
sx q[3];
rz(-1.305745) q[3];
sx q[3];
rz(1.2983398) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.558202) q[2];
sx q[2];
rz(-2.148197) q[2];
sx q[2];
rz(-2.2885585) q[2];
rz(-0.17523003) q[3];
sx q[3];
rz(-1.8195189) q[3];
sx q[3];
rz(2.5186553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(2.919642) q[0];
sx q[0];
rz(-2.3355244) q[0];
sx q[0];
rz(-2.6958418) q[0];
rz(-2.1469965) q[1];
sx q[1];
rz(-2.1371195) q[1];
sx q[1];
rz(-1.4882784) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19916473) q[0];
sx q[0];
rz(-0.86948778) q[0];
sx q[0];
rz(1.2851932) q[0];
rz(-pi) q[1];
rz(-1.4966665) q[2];
sx q[2];
rz(-1.6152116) q[2];
sx q[2];
rz(-1.6629459) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1057652) q[1];
sx q[1];
rz(-1.2700541) q[1];
sx q[1];
rz(-0.2161628) q[1];
rz(1.3429232) q[3];
sx q[3];
rz(-0.3349492) q[3];
sx q[3];
rz(-3.0216274) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.970924) q[2];
sx q[2];
rz(-1.8985775) q[2];
sx q[2];
rz(-2.1155913) q[2];
rz(0.79801997) q[3];
sx q[3];
rz(-2.3012216) q[3];
sx q[3];
rz(0.26871267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84665027) q[0];
sx q[0];
rz(-1.3947399) q[0];
sx q[0];
rz(0.41644874) q[0];
rz(1.396748) q[1];
sx q[1];
rz(-1.6114707) q[1];
sx q[1];
rz(-1.6646615) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0234982) q[0];
sx q[0];
rz(-2.5865285) q[0];
sx q[0];
rz(-1.5007625) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.32768048) q[2];
sx q[2];
rz(-1.9679356) q[2];
sx q[2];
rz(-0.59573345) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.0020963) q[1];
sx q[1];
rz(-1.1285892) q[1];
sx q[1];
rz(-0.79331974) q[1];
x q[2];
rz(-2.8028071) q[3];
sx q[3];
rz(-1.4926406) q[3];
sx q[3];
rz(1.6878267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.2485409) q[2];
sx q[2];
rz(-2.7613642) q[2];
sx q[2];
rz(-0.60379544) q[2];
rz(0.90886146) q[3];
sx q[3];
rz(-1.7085608) q[3];
sx q[3];
rz(-1.0994256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
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
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4029082) q[0];
sx q[0];
rz(-1.7576907) q[0];
sx q[0];
rz(-2.4201194) q[0];
rz(-1.8062704) q[1];
sx q[1];
rz(-1.1386917) q[1];
sx q[1];
rz(-1.8849751) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54859358) q[0];
sx q[0];
rz(-1.8701435) q[0];
sx q[0];
rz(-2.3936098) q[0];
x q[1];
rz(-1.0851791) q[2];
sx q[2];
rz(-2.0943421) q[2];
sx q[2];
rz(2.8741037) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.9049038) q[1];
sx q[1];
rz(-1.3373378) q[1];
sx q[1];
rz(-1.1462565) q[1];
rz(-1.6001892) q[3];
sx q[3];
rz(-1.6163858) q[3];
sx q[3];
rz(-2.8025093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.1544372) q[2];
sx q[2];
rz(-0.64266959) q[2];
sx q[2];
rz(-1.2745693) q[2];
rz(0.68862033) q[3];
sx q[3];
rz(-1.6504811) q[3];
sx q[3];
rz(3.137818) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1275198) q[0];
sx q[0];
rz(-1.8919683) q[0];
sx q[0];
rz(2.4697812) q[0];
rz(-1.5348684) q[1];
sx q[1];
rz(-0.22767362) q[1];
sx q[1];
rz(0.52275503) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1296165) q[0];
sx q[0];
rz(-2.2253941) q[0];
sx q[0];
rz(-1.2415166) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9143064) q[2];
sx q[2];
rz(-1.6186423) q[2];
sx q[2];
rz(-2.7005418) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8485496) q[1];
sx q[1];
rz(-2.5575752) q[1];
sx q[1];
rz(-0.73941458) q[1];
x q[2];
rz(-1.7789442) q[3];
sx q[3];
rz(-1.8878332) q[3];
sx q[3];
rz(3.1274743) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.25040024) q[2];
sx q[2];
rz(-2.4947391) q[2];
sx q[2];
rz(0.35987535) q[2];
rz(-1.322809) q[3];
sx q[3];
rz(-1.4578994) q[3];
sx q[3];
rz(-0.047164269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5322402) q[0];
sx q[0];
rz(-0.9541963) q[0];
sx q[0];
rz(2.1954319) q[0];
rz(-0.040945176) q[1];
sx q[1];
rz(-1.2766726) q[1];
sx q[1];
rz(1.2650222) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0505943) q[0];
sx q[0];
rz(-2.6532905) q[0];
sx q[0];
rz(2.9191892) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6870747) q[2];
sx q[2];
rz(-1.3474166) q[2];
sx q[2];
rz(-0.73385161) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.62502669) q[1];
sx q[1];
rz(-1.0426765) q[1];
sx q[1];
rz(0.81104802) q[1];
rz(-pi) q[2];
x q[2];
rz(0.048790292) q[3];
sx q[3];
rz(-1.1975761) q[3];
sx q[3];
rz(2.9874955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.1653183) q[2];
sx q[2];
rz(-1.1555305) q[2];
sx q[2];
rz(-2.9602642) q[2];
rz(2.2881962) q[3];
sx q[3];
rz(-0.80297628) q[3];
sx q[3];
rz(1.8027423) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4637909) q[0];
sx q[0];
rz(-1.3297357) q[0];
sx q[0];
rz(1.7672675) q[0];
rz(1.455066) q[1];
sx q[1];
rz(-1.4604026) q[1];
sx q[1];
rz(-0.30053465) q[1];
rz(-0.64114943) q[2];
sx q[2];
rz(-1.9998301) q[2];
sx q[2];
rz(2.395973) q[2];
rz(0.72814124) q[3];
sx q[3];
rz(-0.37708382) q[3];
sx q[3];
rz(2.0711475) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
