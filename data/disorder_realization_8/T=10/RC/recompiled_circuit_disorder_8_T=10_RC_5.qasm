OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.33558694) q[0];
sx q[0];
rz(-2.196329) q[0];
sx q[0];
rz(0.52559108) q[0];
rz(-2.8984012) q[1];
sx q[1];
rz(-1.2326198) q[1];
sx q[1];
rz(-0.90484172) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5072767) q[0];
sx q[0];
rz(-2.0024558) q[0];
sx q[0];
rz(1.0213724) q[0];
rz(-pi) q[1];
rz(-2.4279847) q[2];
sx q[2];
rz(-0.2954233) q[2];
sx q[2];
rz(1.7507391) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.3416268) q[1];
sx q[1];
rz(-2.6744665) q[1];
sx q[1];
rz(-2.1145691) q[1];
x q[2];
rz(1.8815133) q[3];
sx q[3];
rz(-1.7508535) q[3];
sx q[3];
rz(-2.3678399) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.9378172) q[2];
sx q[2];
rz(-1.4062466) q[2];
sx q[2];
rz(-0.096244372) q[2];
rz(-1.0359267) q[3];
sx q[3];
rz(-2.7544498) q[3];
sx q[3];
rz(2.9878785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44089833) q[0];
sx q[0];
rz(-0.39114025) q[0];
sx q[0];
rz(0.76517117) q[0];
rz(-1.8493429) q[1];
sx q[1];
rz(-2.6563829) q[1];
sx q[1];
rz(-0.66295019) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11766079) q[0];
sx q[0];
rz(-1.5045325) q[0];
sx q[0];
rz(0.06225417) q[0];
rz(-pi) q[1];
x q[1];
rz(0.40132482) q[2];
sx q[2];
rz(-2.1839645) q[2];
sx q[2];
rz(-2.3344628) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7141124) q[1];
sx q[1];
rz(-0.44891) q[1];
sx q[1];
rz(-0.19101363) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4088267) q[3];
sx q[3];
rz(-1.0893981) q[3];
sx q[3];
rz(3.0961406) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-3.092676) q[2];
sx q[2];
rz(-1.0485704) q[2];
sx q[2];
rz(-0.32901397) q[2];
rz(-0.66550231) q[3];
sx q[3];
rz(-0.21829675) q[3];
sx q[3];
rz(1.8238508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5488141) q[0];
sx q[0];
rz(-2.9187262) q[0];
sx q[0];
rz(-2.9192525) q[0];
rz(2.1242583) q[1];
sx q[1];
rz(-2.4203114) q[1];
sx q[1];
rz(2.6229048) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1919353) q[0];
sx q[0];
rz(-1.2149095) q[0];
sx q[0];
rz(2.3200672) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8027595) q[2];
sx q[2];
rz(-0.28140861) q[2];
sx q[2];
rz(2.0237405) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.7171214) q[1];
sx q[1];
rz(-0.79454225) q[1];
sx q[1];
rz(0.24309991) q[1];
rz(-pi) q[2];
x q[2];
rz(0.014702602) q[3];
sx q[3];
rz(-3.0869752) q[3];
sx q[3];
rz(0.86044776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.8609994) q[2];
sx q[2];
rz(-2.7797647) q[2];
sx q[2];
rz(-0.068543531) q[2];
rz(-2.5391501) q[3];
sx q[3];
rz(-2.3790363) q[3];
sx q[3];
rz(0.13901916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
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
rz(-1.4908726) q[0];
sx q[0];
rz(-2.3644709) q[0];
sx q[0];
rz(-2.9673476) q[0];
rz(0.53025591) q[1];
sx q[1];
rz(-1.6590051) q[1];
sx q[1];
rz(0.51309103) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4671191) q[0];
sx q[0];
rz(-1.6498483) q[0];
sx q[0];
rz(-0.62069686) q[0];
x q[1];
rz(2.3344343) q[2];
sx q[2];
rz(-2.0370738) q[2];
sx q[2];
rz(1.5392787) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.7375813) q[1];
sx q[1];
rz(-1.1176795) q[1];
sx q[1];
rz(1.416942) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4245093) q[3];
sx q[3];
rz(-0.20855599) q[3];
sx q[3];
rz(-1.8640765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.6667368) q[2];
sx q[2];
rz(-2.7754144) q[2];
sx q[2];
rz(0.22988698) q[2];
rz(0.41904467) q[3];
sx q[3];
rz(-1.34904) q[3];
sx q[3];
rz(0.45927799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6722365) q[0];
sx q[0];
rz(-2.3903963) q[0];
sx q[0];
rz(0.67681926) q[0];
rz(0.49304402) q[1];
sx q[1];
rz(-2.1926011) q[1];
sx q[1];
rz(0.61606032) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13000935) q[0];
sx q[0];
rz(-1.6066215) q[0];
sx q[0];
rz(1.623276) q[0];
rz(-pi) q[1];
rz(1.4511756) q[2];
sx q[2];
rz(-1.9336485) q[2];
sx q[2];
rz(0.56064831) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.39561135) q[1];
sx q[1];
rz(-0.99318722) q[1];
sx q[1];
rz(-1.5773768) q[1];
rz(-pi) q[2];
rz(1.1957937) q[3];
sx q[3];
rz(-0.71205322) q[3];
sx q[3];
rz(-0.28111162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.73413509) q[2];
sx q[2];
rz(-0.55441684) q[2];
sx q[2];
rz(-0.25536728) q[2];
rz(1.5363961) q[3];
sx q[3];
rz(-0.95241773) q[3];
sx q[3];
rz(0.77409625) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42246321) q[0];
sx q[0];
rz(-0.96619773) q[0];
sx q[0];
rz(-0.47250026) q[0];
rz(0.52608144) q[1];
sx q[1];
rz(-0.20985797) q[1];
sx q[1];
rz(0.88476673) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9677744) q[0];
sx q[0];
rz(-2.7533555) q[0];
sx q[0];
rz(-1.3740747) q[0];
rz(1.5849437) q[2];
sx q[2];
rz(-1.3607585) q[2];
sx q[2];
rz(-1.1426136) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2192167) q[1];
sx q[1];
rz(-1.1616542) q[1];
sx q[1];
rz(0.31422305) q[1];
x q[2];
rz(-2.1045477) q[3];
sx q[3];
rz(-1.3579206) q[3];
sx q[3];
rz(0.45149976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.3970967) q[2];
sx q[2];
rz(-2.3366191) q[2];
sx q[2];
rz(2.8302144) q[2];
rz(1.3686251) q[3];
sx q[3];
rz(-0.45752782) q[3];
sx q[3];
rz(0.5293203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7235274) q[0];
sx q[0];
rz(-2.8630246) q[0];
sx q[0];
rz(0.061070651) q[0];
rz(-3.1014077) q[1];
sx q[1];
rz(-1.9804852) q[1];
sx q[1];
rz(0.73289245) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7948579) q[0];
sx q[0];
rz(-1.96694) q[0];
sx q[0];
rz(2.5651155) q[0];
rz(-pi) q[1];
rz(1.0211208) q[2];
sx q[2];
rz(-2.3602544) q[2];
sx q[2];
rz(0.48175016) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.43436189) q[1];
sx q[1];
rz(-1.5773298) q[1];
sx q[1];
rz(-2.1965501) q[1];
rz(-pi) q[2];
rz(-0.041386889) q[3];
sx q[3];
rz(-1.4635411) q[3];
sx q[3];
rz(-1.0556575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.0733033) q[2];
sx q[2];
rz(-1.1969593) q[2];
sx q[2];
rz(0.51458365) q[2];
rz(1.9040646) q[3];
sx q[3];
rz(-2.5585744) q[3];
sx q[3];
rz(0.54491836) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.056203689) q[0];
sx q[0];
rz(-1.6563002) q[0];
sx q[0];
rz(-0.7094267) q[0];
rz(1.5052694) q[1];
sx q[1];
rz(-2.067833) q[1];
sx q[1];
rz(0.27871305) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7473135) q[0];
sx q[0];
rz(-2.8907667) q[0];
sx q[0];
rz(0.5745116) q[0];
rz(-pi) q[1];
x q[1];
rz(0.3785554) q[2];
sx q[2];
rz(-2.5494529) q[2];
sx q[2];
rz(-2.2373667) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7396001) q[1];
sx q[1];
rz(-2.7584689) q[1];
sx q[1];
rz(-0.95462228) q[1];
rz(-1.7940984) q[3];
sx q[3];
rz(-1.6468862) q[3];
sx q[3];
rz(-3.1267816) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.30148208) q[2];
sx q[2];
rz(-1.8436517) q[2];
sx q[2];
rz(3.0855132) q[2];
rz(-0.85514832) q[3];
sx q[3];
rz(-2.6931098) q[3];
sx q[3];
rz(2.7364031) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99496019) q[0];
sx q[0];
rz(-2.9652847) q[0];
sx q[0];
rz(2.1726998) q[0];
rz(0.47337198) q[1];
sx q[1];
rz(-2.3469766) q[1];
sx q[1];
rz(-1.1425346) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0363077) q[0];
sx q[0];
rz(-1.6616016) q[0];
sx q[0];
rz(-1.275195) q[0];
rz(-pi) q[1];
x q[1];
rz(0.70334401) q[2];
sx q[2];
rz(-1.568927) q[2];
sx q[2];
rz(1.6714931) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.7361703) q[1];
sx q[1];
rz(-1.0338963) q[1];
sx q[1];
rz(-1.7169397) q[1];
x q[2];
rz(1.8248796) q[3];
sx q[3];
rz(-2.4554606) q[3];
sx q[3];
rz(0.38476598) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2450976) q[2];
sx q[2];
rz(-1.2197887) q[2];
sx q[2];
rz(-1.0207821) q[2];
rz(-2.8178689) q[3];
sx q[3];
rz(-2.3886069) q[3];
sx q[3];
rz(-3.1304205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97994119) q[0];
sx q[0];
rz(-0.027898235) q[0];
sx q[0];
rz(2.4401869) q[0];
rz(2.2258863) q[1];
sx q[1];
rz(-2.1332824) q[1];
sx q[1];
rz(-1.2385626) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4558444) q[0];
sx q[0];
rz(-1.0554753) q[0];
sx q[0];
rz(2.8932163) q[0];
rz(-pi) q[1];
rz(-0.752991) q[2];
sx q[2];
rz(-1.5333311) q[2];
sx q[2];
rz(2.3245036) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.30675948) q[1];
sx q[1];
rz(-0.66654897) q[1];
sx q[1];
rz(-1.3798316) q[1];
x q[2];
rz(0.978312) q[3];
sx q[3];
rz(-1.6654135) q[3];
sx q[3];
rz(2.832151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.4548268) q[2];
sx q[2];
rz(-2.0451615) q[2];
sx q[2];
rz(-3.0977541) q[2];
rz(-1.94058) q[3];
sx q[3];
rz(-0.73533708) q[3];
sx q[3];
rz(2.1380077) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4713521) q[0];
sx q[0];
rz(-0.72605194) q[0];
sx q[0];
rz(-1.3656021) q[0];
rz(-3.1162221) q[1];
sx q[1];
rz(-1.8352958) q[1];
sx q[1];
rz(-1.8713554) q[1];
rz(-2.2123443) q[2];
sx q[2];
rz(-2.6570448) q[2];
sx q[2];
rz(-1.6397283) q[2];
rz(1.745789) q[3];
sx q[3];
rz(-0.77944118) q[3];
sx q[3];
rz(-0.13164095) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
