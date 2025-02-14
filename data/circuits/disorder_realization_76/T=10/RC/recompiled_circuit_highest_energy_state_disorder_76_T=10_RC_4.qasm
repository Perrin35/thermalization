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
rz(-0.79987502) q[0];
sx q[0];
rz(-0.81698155) q[0];
sx q[0];
rz(-0.45720994) q[0];
rz(0.41181052) q[1];
sx q[1];
rz(4.7635912) q[1];
sx q[1];
rz(6.8430321) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.354411) q[0];
sx q[0];
rz(-2.5776064) q[0];
sx q[0];
rz(-0.63040479) q[0];
rz(-pi) q[1];
x q[1];
rz(1.067808) q[2];
sx q[2];
rz(-1.3185698) q[2];
sx q[2];
rz(-2.5258738) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.7840883) q[1];
sx q[1];
rz(-2.5604446) q[1];
sx q[1];
rz(-0.52304348) q[1];
rz(-pi) q[2];
rz(1.4234957) q[3];
sx q[3];
rz(-1.788874) q[3];
sx q[3];
rz(0.17373057) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.37602279) q[2];
sx q[2];
rz(-2.1799808) q[2];
sx q[2];
rz(-1.8454856) q[2];
rz(-0.62093312) q[3];
sx q[3];
rz(-2.4758078) q[3];
sx q[3];
rz(-0.92500979) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45944443) q[0];
sx q[0];
rz(-0.58732533) q[0];
sx q[0];
rz(-2.4272954) q[0];
rz(2.4192877) q[1];
sx q[1];
rz(-2.0562833) q[1];
sx q[1];
rz(1.3246271) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6178169) q[0];
sx q[0];
rz(-1.1259698) q[0];
sx q[0];
rz(-2.8033048) q[0];
rz(2.4087853) q[2];
sx q[2];
rz(-1.5298973) q[2];
sx q[2];
rz(0.78860229) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.87323071) q[1];
sx q[1];
rz(-1.9385846) q[1];
sx q[1];
rz(1.7413543) q[1];
rz(0.12217317) q[3];
sx q[3];
rz(-2.2524912) q[3];
sx q[3];
rz(0.091191779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.033919949) q[2];
sx q[2];
rz(-1.718797) q[2];
sx q[2];
rz(0.99204341) q[2];
rz(0.77573675) q[3];
sx q[3];
rz(-2.3596767) q[3];
sx q[3];
rz(1.0824664) q[3];
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
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42763212) q[0];
sx q[0];
rz(-1.7949224) q[0];
sx q[0];
rz(0.91208518) q[0];
rz(2.7569547) q[1];
sx q[1];
rz(-2.1802528) q[1];
sx q[1];
rz(-0.49547637) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3289483) q[0];
sx q[0];
rz(-1.616713) q[0];
sx q[0];
rz(-1.4981235) q[0];
rz(-pi) q[1];
x q[1];
rz(1.762359) q[2];
sx q[2];
rz(-2.3108916) q[2];
sx q[2];
rz(0.57410115) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.66049536) q[1];
sx q[1];
rz(-0.7259136) q[1];
sx q[1];
rz(-3.0654415) q[1];
x q[2];
rz(2.1927102) q[3];
sx q[3];
rz(-1.5601741) q[3];
sx q[3];
rz(1.1158021) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.89874011) q[2];
sx q[2];
rz(-1.1246559) q[2];
sx q[2];
rz(0.43453547) q[2];
rz(3.0430326) q[3];
sx q[3];
rz(-1.7193272) q[3];
sx q[3];
rz(-2.3677473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74715215) q[0];
sx q[0];
rz(-2.9065865) q[0];
sx q[0];
rz(2.8727942) q[0];
rz(2.0647743) q[1];
sx q[1];
rz(-1.7348758) q[1];
sx q[1];
rz(-1.9487618) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4736611) q[0];
sx q[0];
rz(-2.8369378) q[0];
sx q[0];
rz(-0.61078914) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0054686) q[2];
sx q[2];
rz(-1.4413712) q[2];
sx q[2];
rz(-0.65107513) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.0620246) q[1];
sx q[1];
rz(-2.8753198) q[1];
sx q[1];
rz(-2.5063199) q[1];
rz(-2.4402839) q[3];
sx q[3];
rz(-1.4150672) q[3];
sx q[3];
rz(2.6913672) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.25527915) q[2];
sx q[2];
rz(-0.89526075) q[2];
sx q[2];
rz(-2.8423584) q[2];
rz(-0.29087654) q[3];
sx q[3];
rz(-1.8894922) q[3];
sx q[3];
rz(2.4860184) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3000779) q[0];
sx q[0];
rz(-0.97989196) q[0];
sx q[0];
rz(-3.063524) q[0];
rz(0.95371753) q[1];
sx q[1];
rz(-2.6225312) q[1];
sx q[1];
rz(-1.5740707) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6758976) q[0];
sx q[0];
rz(-2.2392096) q[0];
sx q[0];
rz(-2.9107679) q[0];
rz(-pi) q[1];
rz(2.7360953) q[2];
sx q[2];
rz(-2.4031696) q[2];
sx q[2];
rz(-0.15106311) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.003338699) q[1];
sx q[1];
rz(-1.6032156) q[1];
sx q[1];
rz(2.1829672) q[1];
rz(1.6408345) q[3];
sx q[3];
rz(-0.87890714) q[3];
sx q[3];
rz(0.50229154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.8629525) q[2];
sx q[2];
rz(-2.7095257) q[2];
sx q[2];
rz(-0.26818177) q[2];
rz(-0.48926085) q[3];
sx q[3];
rz(-1.093995) q[3];
sx q[3];
rz(-1.6930273) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0443403) q[0];
sx q[0];
rz(-0.02359979) q[0];
sx q[0];
rz(-0.46446717) q[0];
rz(-2.6181472) q[1];
sx q[1];
rz(-2.372066) q[1];
sx q[1];
rz(0.73062599) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7775041) q[0];
sx q[0];
rz(-2.3100762) q[0];
sx q[0];
rz(2.9502908) q[0];
rz(0.67209216) q[2];
sx q[2];
rz(-1.8037829) q[2];
sx q[2];
rz(-1.1868049) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.043482232) q[1];
sx q[1];
rz(-0.93096369) q[1];
sx q[1];
rz(2.1049064) q[1];
x q[2];
rz(0.2517638) q[3];
sx q[3];
rz(-1.6986966) q[3];
sx q[3];
rz(-2.3826016) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.0940493) q[2];
sx q[2];
rz(-1.0796248) q[2];
sx q[2];
rz(-0.13761061) q[2];
rz(1.1329457) q[3];
sx q[3];
rz(-1.9652941) q[3];
sx q[3];
rz(-1.8784116) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15261821) q[0];
sx q[0];
rz(-1.6292097) q[0];
sx q[0];
rz(3.1076987) q[0];
rz(2.7821817) q[1];
sx q[1];
rz(-1.5243328) q[1];
sx q[1];
rz(2.3072037) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9890404) q[0];
sx q[0];
rz(-1.9473528) q[0];
sx q[0];
rz(1.2370308) q[0];
x q[1];
rz(-0.57247803) q[2];
sx q[2];
rz(-0.77478638) q[2];
sx q[2];
rz(3.0511193) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.5339116) q[1];
sx q[1];
rz(-1.4754731) q[1];
sx q[1];
rz(1.3168112) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5079751) q[3];
sx q[3];
rz(-2.2414226) q[3];
sx q[3];
rz(0.84298979) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.72839165) q[2];
sx q[2];
rz(-0.30343702) q[2];
sx q[2];
rz(2.0835853) q[2];
rz(0.47438619) q[3];
sx q[3];
rz(-0.79550231) q[3];
sx q[3];
rz(2.0856196) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0742842) q[0];
sx q[0];
rz(-1.5165167) q[0];
sx q[0];
rz(-1.7838595) q[0];
rz(-0.43074295) q[1];
sx q[1];
rz(-1.6669225) q[1];
sx q[1];
rz(-0.46708435) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3245832) q[0];
sx q[0];
rz(-1.5276434) q[0];
sx q[0];
rz(3.140121) q[0];
x q[1];
rz(0.22612382) q[2];
sx q[2];
rz(-0.71561049) q[2];
sx q[2];
rz(-2.0307045) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.704485) q[1];
sx q[1];
rz(-2.6117755) q[1];
sx q[1];
rz(-0.036661224) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.93293087) q[3];
sx q[3];
rz(-1.5278421) q[3];
sx q[3];
rz(1.4900946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.0507386) q[2];
sx q[2];
rz(-2.723912) q[2];
sx q[2];
rz(1.0806855) q[2];
rz(-2.4596227) q[3];
sx q[3];
rz(-2.3365648) q[3];
sx q[3];
rz(-0.69847703) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4526116) q[0];
sx q[0];
rz(-0.51849759) q[0];
sx q[0];
rz(0.80775753) q[0];
rz(-1.0001596) q[1];
sx q[1];
rz(-2.2554485) q[1];
sx q[1];
rz(0.79661405) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2444435) q[0];
sx q[0];
rz(-1.5287491) q[0];
sx q[0];
rz(-1.6680662) q[0];
rz(-pi) q[1];
rz(1.3299204) q[2];
sx q[2];
rz(-0.95930525) q[2];
sx q[2];
rz(-0.73630737) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.1353242) q[1];
sx q[1];
rz(-2.709143) q[1];
sx q[1];
rz(1.2342374) q[1];
rz(-pi) q[2];
rz(-2.8878941) q[3];
sx q[3];
rz(-1.2658409) q[3];
sx q[3];
rz(2.0010468) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.54958582) q[2];
sx q[2];
rz(-1.0909811) q[2];
sx q[2];
rz(2.8263212) q[2];
rz(1.5677876) q[3];
sx q[3];
rz(-1.1476293) q[3];
sx q[3];
rz(-1.1228336) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10298097) q[0];
sx q[0];
rz(-0.6655612) q[0];
sx q[0];
rz(0.82157201) q[0];
rz(0.483825) q[1];
sx q[1];
rz(-1.4204104) q[1];
sx q[1];
rz(0.22463591) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.076775) q[0];
sx q[0];
rz(-0.358264) q[0];
sx q[0];
rz(-2.0117156) q[0];
rz(-2.2251525) q[2];
sx q[2];
rz(-1.0447557) q[2];
sx q[2];
rz(1.9182084) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.0269751) q[1];
sx q[1];
rz(-1.7791505) q[1];
sx q[1];
rz(2.2592553) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3566689) q[3];
sx q[3];
rz(-1.4487293) q[3];
sx q[3];
rz(-1.7904953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3839174) q[2];
sx q[2];
rz(-3.0249247) q[2];
sx q[2];
rz(3.0465916) q[2];
rz(1.4875686) q[3];
sx q[3];
rz(-0.58949685) q[3];
sx q[3];
rz(0.76475638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8944396) q[0];
sx q[0];
rz(-2.7338487) q[0];
sx q[0];
rz(-0.26640531) q[0];
rz(-2.5447625) q[1];
sx q[1];
rz(-1.6886371) q[1];
sx q[1];
rz(1.4773038) q[1];
rz(-0.3565557) q[2];
sx q[2];
rz(-2.2168163) q[2];
sx q[2];
rz(1.4067895) q[2];
rz(0.37551914) q[3];
sx q[3];
rz(-1.569034) q[3];
sx q[3];
rz(1.9636596) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
