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
rz(-2.6744106) q[0];
sx q[0];
rz(-1.8888357) q[0];
sx q[0];
rz(-0.80656615) q[0];
rz(-0.84713495) q[1];
sx q[1];
rz(-0.48589125) q[1];
sx q[1];
rz(-2.2810305) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8019077) q[0];
sx q[0];
rz(-0.97754495) q[0];
sx q[0];
rz(2.9723333) q[0];
rz(-pi) q[1];
rz(1.8612618) q[2];
sx q[2];
rz(-0.72735255) q[2];
sx q[2];
rz(1.0199821) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.5570506) q[1];
sx q[1];
rz(-1.2963386) q[1];
sx q[1];
rz(1.518979) q[1];
x q[2];
rz(-2.2673554) q[3];
sx q[3];
rz(-2.4421066) q[3];
sx q[3];
rz(1.4892088) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.2412771) q[2];
sx q[2];
rz(-2.6130455) q[2];
sx q[2];
rz(0.47145525) q[2];
rz(2.6491162) q[3];
sx q[3];
rz(-1.9086647) q[3];
sx q[3];
rz(-1.8878149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(-2.8660698) q[0];
sx q[0];
rz(-0.38551426) q[0];
sx q[0];
rz(0.27819124) q[0];
rz(-1.1909852) q[1];
sx q[1];
rz(-1.3326125) q[1];
sx q[1];
rz(1.2954378) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71466509) q[0];
sx q[0];
rz(-2.1708827) q[0];
sx q[0];
rz(2.1467208) q[0];
x q[1];
rz(3.0619925) q[2];
sx q[2];
rz(-1.6222553) q[2];
sx q[2];
rz(-2.1934794) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.56486121) q[1];
sx q[1];
rz(-0.38034359) q[1];
sx q[1];
rz(1.5647792) q[1];
rz(-pi) q[2];
rz(2.1987887) q[3];
sx q[3];
rz(-2.0217388) q[3];
sx q[3];
rz(2.1357) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.94464716) q[2];
sx q[2];
rz(-1.4713919) q[2];
sx q[2];
rz(0.50986457) q[2];
rz(1.4465796) q[3];
sx q[3];
rz(-2.6810985) q[3];
sx q[3];
rz(-2.6367326) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6079717) q[0];
sx q[0];
rz(-1.5190268) q[0];
sx q[0];
rz(-0.98440379) q[0];
rz(-3.1254752) q[1];
sx q[1];
rz(-2.0033483) q[1];
sx q[1];
rz(1.3498397) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70671073) q[0];
sx q[0];
rz(-0.87315403) q[0];
sx q[0];
rz(1.7876704) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0060293) q[2];
sx q[2];
rz(-1.2899219) q[2];
sx q[2];
rz(1.8233521) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.3193839) q[1];
sx q[1];
rz(-2.2693386) q[1];
sx q[1];
rz(-1.6234267) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4593342) q[3];
sx q[3];
rz(-0.80815274) q[3];
sx q[3];
rz(-2.7956642) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.4237889) q[2];
sx q[2];
rz(-0.70222792) q[2];
sx q[2];
rz(-0.88199893) q[2];
rz(-1.1841904) q[3];
sx q[3];
rz(-0.94112527) q[3];
sx q[3];
rz(-0.73630303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0940359) q[0];
sx q[0];
rz(-1.2268257) q[0];
sx q[0];
rz(-0.077022821) q[0];
rz(2.9432964) q[1];
sx q[1];
rz(-1.6402596) q[1];
sx q[1];
rz(0.1870627) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27320751) q[0];
sx q[0];
rz(-1.4131695) q[0];
sx q[0];
rz(0.34021838) q[0];
rz(-pi) q[1];
rz(2.0346626) q[2];
sx q[2];
rz(-2.7527546) q[2];
sx q[2];
rz(2.1269682) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.9976207) q[1];
sx q[1];
rz(-1.629843) q[1];
sx q[1];
rz(0.70934341) q[1];
x q[2];
rz(1.2310394) q[3];
sx q[3];
rz(-1.8700784) q[3];
sx q[3];
rz(-1.2359985) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.5932811) q[2];
sx q[2];
rz(-0.45479861) q[2];
sx q[2];
rz(-1.456267) q[2];
rz(2.4233387) q[3];
sx q[3];
rz(-1.7536438) q[3];
sx q[3];
rz(2.3072402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89183557) q[0];
sx q[0];
rz(-1.3351853) q[0];
sx q[0];
rz(-0.43858132) q[0];
rz(0.35722411) q[1];
sx q[1];
rz(-1.2780739) q[1];
sx q[1];
rz(1.8908148) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6092458) q[0];
sx q[0];
rz(-0.43918434) q[0];
sx q[0];
rz(1.8076375) q[0];
x q[1];
rz(1.3322796) q[2];
sx q[2];
rz(-0.3198238) q[2];
sx q[2];
rz(2.8826098) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.4692396) q[1];
sx q[1];
rz(-1.1732983) q[1];
sx q[1];
rz(-0.40594493) q[1];
rz(0.26024466) q[3];
sx q[3];
rz(-1.6992927) q[3];
sx q[3];
rz(-0.63333007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.232051) q[2];
sx q[2];
rz(-2.6614058) q[2];
sx q[2];
rz(-1.6298182) q[2];
rz(0.016228598) q[3];
sx q[3];
rz(-2.0461693) q[3];
sx q[3];
rz(-0.79286638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4474354) q[0];
sx q[0];
rz(-1.2801535) q[0];
sx q[0];
rz(-1.9919027) q[0];
rz(2.7206874) q[1];
sx q[1];
rz(-1.4525388) q[1];
sx q[1];
rz(3.0973869) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51254184) q[0];
sx q[0];
rz(-2.3568332) q[0];
sx q[0];
rz(2.6757338) q[0];
rz(1.030613) q[2];
sx q[2];
rz(-1.1655131) q[2];
sx q[2];
rz(-1.1073529) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6730576) q[1];
sx q[1];
rz(-1.4899163) q[1];
sx q[1];
rz(1.1799501) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6930674) q[3];
sx q[3];
rz(-2.5618346) q[3];
sx q[3];
rz(0.019252456) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.9629024) q[2];
sx q[2];
rz(-0.89683214) q[2];
sx q[2];
rz(0.94513354) q[2];
rz(1.6804228) q[3];
sx q[3];
rz(-1.1472568) q[3];
sx q[3];
rz(0.51819658) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1951676) q[0];
sx q[0];
rz(-2.8312046) q[0];
sx q[0];
rz(-2.2921966) q[0];
rz(1.7646344) q[1];
sx q[1];
rz(-2.5701249) q[1];
sx q[1];
rz(-1.2845385) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97235128) q[0];
sx q[0];
rz(-1.961876) q[0];
sx q[0];
rz(-2.611586) q[0];
rz(-pi) q[1];
rz(-2.8340802) q[2];
sx q[2];
rz(-1.4685681) q[2];
sx q[2];
rz(-0.44164613) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.7276638) q[1];
sx q[1];
rz(-1.5352121) q[1];
sx q[1];
rz(-2.6160081) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6553978) q[3];
sx q[3];
rz(-1.0665585) q[3];
sx q[3];
rz(2.274226) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.66594243) q[2];
sx q[2];
rz(-1.7721756) q[2];
sx q[2];
rz(1.5671889) q[2];
rz(0.33268467) q[3];
sx q[3];
rz(-2.0761469) q[3];
sx q[3];
rz(0.98193297) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57629267) q[0];
sx q[0];
rz(-1.1197634) q[0];
sx q[0];
rz(-2.0530307) q[0];
rz(2.2125878) q[1];
sx q[1];
rz(-1.6477511) q[1];
sx q[1];
rz(-0.30379024) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.740363) q[0];
sx q[0];
rz(-2.1308194) q[0];
sx q[0];
rz(1.1267046) q[0];
rz(-pi) q[1];
rz(-0.11907507) q[2];
sx q[2];
rz(-1.4085253) q[2];
sx q[2];
rz(-1.281126) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.5227844) q[1];
sx q[1];
rz(-1.4780103) q[1];
sx q[1];
rz(-1.1397847) q[1];
x q[2];
rz(1.0904543) q[3];
sx q[3];
rz(-2.1340202) q[3];
sx q[3];
rz(-2.8728812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.8888237) q[2];
sx q[2];
rz(-1.6656275) q[2];
sx q[2];
rz(2.8775173) q[2];
rz(-2.0090328) q[3];
sx q[3];
rz(-2.8045636) q[3];
sx q[3];
rz(1.0549217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7842512) q[0];
sx q[0];
rz(-2.7989474) q[0];
sx q[0];
rz(-2.1270879) q[0];
rz(-2.2134589) q[1];
sx q[1];
rz(-1.4200297) q[1];
sx q[1];
rz(0.49044213) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3532012) q[0];
sx q[0];
rz(-1.7176796) q[0];
sx q[0];
rz(-3.0753972) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1855256) q[2];
sx q[2];
rz(-1.2495572) q[2];
sx q[2];
rz(2.3209077) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.52940581) q[1];
sx q[1];
rz(-1.2703964) q[1];
sx q[1];
rz(-0.043073853) q[1];
x q[2];
rz(1.0317627) q[3];
sx q[3];
rz(-1.7787053) q[3];
sx q[3];
rz(-2.5202005) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.32533112) q[2];
sx q[2];
rz(-2.0897431) q[2];
sx q[2];
rz(-0.10247792) q[2];
rz(-1.8898194) q[3];
sx q[3];
rz(-1.3636369) q[3];
sx q[3];
rz(-1.6583091) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3935299) q[0];
sx q[0];
rz(-3.0952125) q[0];
sx q[0];
rz(-1.8512132) q[0];
rz(-0.45404592) q[1];
sx q[1];
rz(-1.3244649) q[1];
sx q[1];
rz(-2.9581199) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.030122) q[0];
sx q[0];
rz(-2.9577177) q[0];
sx q[0];
rz(-1.4669815) q[0];
rz(2.543599) q[2];
sx q[2];
rz(-1.5584297) q[2];
sx q[2];
rz(-1.7747161) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.390229) q[1];
sx q[1];
rz(-1.5469264) q[1];
sx q[1];
rz(-0.49225251) q[1];
x q[2];
rz(-2.1366664) q[3];
sx q[3];
rz(-1.6477895) q[3];
sx q[3];
rz(-1.0371829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.4219249) q[2];
sx q[2];
rz(-2.0823961) q[2];
sx q[2];
rz(3.1165519) q[2];
rz(-0.54343623) q[3];
sx q[3];
rz(-2.4380324) q[3];
sx q[3];
rz(2.8652625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80171361) q[0];
sx q[0];
rz(-2.1262953) q[0];
sx q[0];
rz(-2.3332818) q[0];
rz(0.33357757) q[1];
sx q[1];
rz(-1.7671276) q[1];
sx q[1];
rz(-0.50552013) q[1];
rz(0.12814604) q[2];
sx q[2];
rz(-1.2432877) q[2];
sx q[2];
rz(0.23040085) q[2];
rz(-0.89743817) q[3];
sx q[3];
rz(-1.031395) q[3];
sx q[3];
rz(-0.54678834) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
