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
rz(2.3350265) q[0];
rz(2.2944577) q[1];
sx q[1];
rz(-2.6557014) q[1];
sx q[1];
rz(2.2810305) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8019077) q[0];
sx q[0];
rz(-0.97754495) q[0];
sx q[0];
rz(0.16925933) q[0];
rz(-pi) q[1];
rz(1.2803308) q[2];
sx q[2];
rz(-2.4142401) q[2];
sx q[2];
rz(1.0199821) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.7736176) q[1];
sx q[1];
rz(-0.27918511) q[1];
sx q[1];
rz(-0.18191819) q[1];
rz(-pi) q[2];
rz(-2.2673554) q[3];
sx q[3];
rz(-2.4421066) q[3];
sx q[3];
rz(1.4892088) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.2412771) q[2];
sx q[2];
rz(-2.6130455) q[2];
sx q[2];
rz(2.6701374) q[2];
rz(-0.49247646) q[3];
sx q[3];
rz(-1.9086647) q[3];
sx q[3];
rz(-1.8878149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27552283) q[0];
sx q[0];
rz(-2.7560784) q[0];
sx q[0];
rz(-2.8634014) q[0];
rz(1.9506075) q[1];
sx q[1];
rz(-1.3326125) q[1];
sx q[1];
rz(-1.8461548) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9339855) q[0];
sx q[0];
rz(-1.1046358) q[0];
sx q[0];
rz(-2.4572479) q[0];
rz(-1.5191742) q[2];
sx q[2];
rz(-1.4913017) q[2];
sx q[2];
rz(0.62678601) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.5767314) q[1];
sx q[1];
rz(-0.38034359) q[1];
sx q[1];
rz(1.5647792) q[1];
rz(0.53923082) q[3];
sx q[3];
rz(-2.1279716) q[3];
sx q[3];
rz(2.2702366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.94464716) q[2];
sx q[2];
rz(-1.6702007) q[2];
sx q[2];
rz(2.6317281) q[2];
rz(-1.4465796) q[3];
sx q[3];
rz(-0.46049419) q[3];
sx q[3];
rz(-2.6367326) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6079717) q[0];
sx q[0];
rz(-1.6225659) q[0];
sx q[0];
rz(2.1571889) q[0];
rz(-0.016117485) q[1];
sx q[1];
rz(-1.1382444) q[1];
sx q[1];
rz(-1.791753) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70671073) q[0];
sx q[0];
rz(-0.87315403) q[0];
sx q[0];
rz(-1.3539223) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0651814) q[2];
sx q[2];
rz(-0.62389031) q[2];
sx q[2];
rz(-2.9817493) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.74047725) q[1];
sx q[1];
rz(-0.7001895) q[1];
sx q[1];
rz(-3.0790331) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0257019) q[3];
sx q[3];
rz(-2.3724764) q[3];
sx q[3];
rz(0.18530857) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.71780378) q[2];
sx q[2];
rz(-0.70222792) q[2];
sx q[2];
rz(0.88199893) q[2];
rz(-1.1841904) q[3];
sx q[3];
rz(-2.2004674) q[3];
sx q[3];
rz(-2.4052896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0475567) q[0];
sx q[0];
rz(-1.9147669) q[0];
sx q[0];
rz(-0.077022821) q[0];
rz(2.9432964) q[1];
sx q[1];
rz(-1.501333) q[1];
sx q[1];
rz(2.95453) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27320751) q[0];
sx q[0];
rz(-1.4131695) q[0];
sx q[0];
rz(-0.34021838) q[0];
rz(1.2195829) q[2];
sx q[2];
rz(-1.7412392) q[2];
sx q[2];
rz(2.1518681) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.6641561) q[1];
sx q[1];
rz(-0.86294791) q[1];
sx q[1];
rz(1.6485467) q[1];
rz(-pi) q[2];
rz(-1.2310394) q[3];
sx q[3];
rz(-1.2715142) q[3];
sx q[3];
rz(1.9055942) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.5932811) q[2];
sx q[2];
rz(-2.686794) q[2];
sx q[2];
rz(1.456267) q[2];
rz(2.4233387) q[3];
sx q[3];
rz(-1.3879489) q[3];
sx q[3];
rz(0.83435241) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2497571) q[0];
sx q[0];
rz(-1.3351853) q[0];
sx q[0];
rz(2.7030113) q[0];
rz(0.35722411) q[1];
sx q[1];
rz(-1.2780739) q[1];
sx q[1];
rz(-1.2507778) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8698602) q[0];
sx q[0];
rz(-1.1446867) q[0];
sx q[0];
rz(-0.10978384) q[0];
x q[1];
rz(1.8093131) q[2];
sx q[2];
rz(-0.3198238) q[2];
sx q[2];
rz(-2.8826098) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.631613) q[1];
sx q[1];
rz(-0.56029472) q[1];
sx q[1];
rz(-2.3255682) q[1];
x q[2];
rz(0.46535551) q[3];
sx q[3];
rz(-0.28959238) q[3];
sx q[3];
rz(-0.48894879) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.232051) q[2];
sx q[2];
rz(-2.6614058) q[2];
sx q[2];
rz(-1.5117744) q[2];
rz(-0.016228598) q[3];
sx q[3];
rz(-1.0954233) q[3];
sx q[3];
rz(-0.79286638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4474354) q[0];
sx q[0];
rz(-1.8614391) q[0];
sx q[0];
rz(-1.1496899) q[0];
rz(-2.7206874) q[1];
sx q[1];
rz(-1.4525388) q[1];
sx q[1];
rz(0.044205753) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51254184) q[0];
sx q[0];
rz(-0.78475941) q[0];
sx q[0];
rz(-2.6757338) q[0];
rz(-2.1109796) q[2];
sx q[2];
rz(-1.1655131) q[2];
sx q[2];
rz(2.0342397) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.468535) q[1];
sx q[1];
rz(-1.6516764) q[1];
sx q[1];
rz(1.1799501) q[1];
rz(-0.53307311) q[3];
sx q[3];
rz(-1.8106451) q[3];
sx q[3];
rz(1.1687973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.9629024) q[2];
sx q[2];
rz(-2.2447605) q[2];
sx q[2];
rz(2.1964591) q[2];
rz(-1.6804228) q[3];
sx q[3];
rz(-1.1472568) q[3];
sx q[3];
rz(2.6233961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1951676) q[0];
sx q[0];
rz(-2.8312046) q[0];
sx q[0];
rz(0.84939605) q[0];
rz(-1.3769582) q[1];
sx q[1];
rz(-2.5701249) q[1];
sx q[1];
rz(-1.2845385) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3234069) q[0];
sx q[0];
rz(-2.0571097) q[0];
sx q[0];
rz(1.125) q[0];
rz(-pi) q[1];
rz(1.6780186) q[2];
sx q[2];
rz(-1.8766512) q[2];
sx q[2];
rz(-2.0448409) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(3.0459874) q[1];
sx q[1];
rz(-0.52667499) q[1];
sx q[1];
rz(0.070835872) q[1];
rz(-pi) q[2];
rz(-2.128781) q[3];
sx q[3];
rz(-1.9922755) q[3];
sx q[3];
rz(0.95343219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.66594243) q[2];
sx q[2];
rz(-1.3694171) q[2];
sx q[2];
rz(1.5671889) q[2];
rz(-2.808908) q[3];
sx q[3];
rz(-1.0654457) q[3];
sx q[3];
rz(-0.98193297) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57629267) q[0];
sx q[0];
rz(-1.1197634) q[0];
sx q[0];
rz(1.0885619) q[0];
rz(2.2125878) q[1];
sx q[1];
rz(-1.6477511) q[1];
sx q[1];
rz(2.8378024) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.740363) q[0];
sx q[0];
rz(-1.0107733) q[0];
sx q[0];
rz(-1.1267046) q[0];
rz(-pi) q[1];
x q[1];
rz(0.94307282) q[2];
sx q[2];
rz(-2.94063) q[2];
sx q[2];
rz(-0.64370868) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.090592314) q[1];
sx q[1];
rz(-1.9998308) q[1];
sx q[1];
rz(-0.10206435) q[1];
x q[2];
rz(1.0904543) q[3];
sx q[3];
rz(-1.0075724) q[3];
sx q[3];
rz(-0.26871142) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.2527689) q[2];
sx q[2];
rz(-1.6656275) q[2];
sx q[2];
rz(0.26407537) q[2];
rz(-1.1325599) q[3];
sx q[3];
rz(-0.3370291) q[3];
sx q[3];
rz(-2.0866709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7842512) q[0];
sx q[0];
rz(-0.34264523) q[0];
sx q[0];
rz(-1.0145048) q[0];
rz(-0.92813379) q[1];
sx q[1];
rz(-1.4200297) q[1];
sx q[1];
rz(2.6511505) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2078932) q[0];
sx q[0];
rz(-1.636278) q[0];
sx q[0];
rz(1.4235953) q[0];
rz(0.84607203) q[2];
sx q[2];
rz(-2.6451689) q[2];
sx q[2];
rz(1.4116532) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.6121868) q[1];
sx q[1];
rz(-1.8711963) q[1];
sx q[1];
rz(-3.0985188) q[1];
x q[2];
rz(2.9005592) q[3];
sx q[3];
rz(-1.0445945) q[3];
sx q[3];
rz(-1.0722425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.8162615) q[2];
sx q[2];
rz(-1.0518495) q[2];
sx q[2];
rz(-3.0391147) q[2];
rz(1.2517733) q[3];
sx q[3];
rz(-1.7779558) q[3];
sx q[3];
rz(1.6583091) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74806279) q[0];
sx q[0];
rz(-0.04638014) q[0];
sx q[0];
rz(-1.8512132) q[0];
rz(-2.6875467) q[1];
sx q[1];
rz(-1.8171277) q[1];
sx q[1];
rz(0.18347278) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35724872) q[0];
sx q[0];
rz(-1.5518477) q[0];
sx q[0];
rz(-1.7537033) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5857592) q[2];
sx q[2];
rz(-2.1687379) q[2];
sx q[2];
rz(-0.19549616) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.390229) q[1];
sx q[1];
rz(-1.5946663) q[1];
sx q[1];
rz(2.6493401) q[1];
x q[2];
rz(3.0504543) q[3];
sx q[3];
rz(-2.1347858) q[3];
sx q[3];
rz(-0.48479776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.4219249) q[2];
sx q[2];
rz(-2.0823961) q[2];
sx q[2];
rz(-3.1165519) q[2];
rz(-2.5981564) q[3];
sx q[3];
rz(-0.70356026) q[3];
sx q[3];
rz(2.8652625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80171361) q[0];
sx q[0];
rz(-1.0152974) q[0];
sx q[0];
rz(0.80831084) q[0];
rz(2.8080151) q[1];
sx q[1];
rz(-1.3744651) q[1];
sx q[1];
rz(2.6360725) q[1];
rz(1.9008209) q[2];
sx q[2];
rz(-1.6920964) q[2];
sx q[2];
rz(1.7597711) q[2];
rz(2.2441545) q[3];
sx q[3];
rz(-1.031395) q[3];
sx q[3];
rz(-0.54678834) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
