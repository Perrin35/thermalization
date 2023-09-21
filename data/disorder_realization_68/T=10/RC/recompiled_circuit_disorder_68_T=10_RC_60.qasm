OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.39419898) q[0];
sx q[0];
rz(-0.49180254) q[0];
sx q[0];
rz(0.1879745) q[0];
rz(2.0239053) q[1];
sx q[1];
rz(4.6586577) q[1];
sx q[1];
rz(12.933856) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.039283218) q[0];
sx q[0];
rz(-1.6741721) q[0];
sx q[0];
rz(-2.1084059) q[0];
rz(-pi) q[1];
rz(3.0066178) q[2];
sx q[2];
rz(-2.0581323) q[2];
sx q[2];
rz(-0.48263532) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.623466) q[1];
sx q[1];
rz(-1.7382442) q[1];
sx q[1];
rz(1.3256339) q[1];
rz(-pi) q[2];
rz(0.46842694) q[3];
sx q[3];
rz(-2.763063) q[3];
sx q[3];
rz(-0.98640501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.1774896) q[2];
sx q[2];
rz(-0.51012817) q[2];
sx q[2];
rz(-0.5509848) q[2];
rz(-1.8356813) q[3];
sx q[3];
rz(-1.6492313) q[3];
sx q[3];
rz(-1.3163542) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
rz(2.6630163) q[0];
sx q[0];
rz(-1.0000279) q[0];
sx q[0];
rz(0.4719032) q[0];
rz(0.42981237) q[1];
sx q[1];
rz(-1.2496354) q[1];
sx q[1];
rz(2.205251) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9443003) q[0];
sx q[0];
rz(-2.8951277) q[0];
sx q[0];
rz(0.36578567) q[0];
x q[1];
rz(0.41994862) q[2];
sx q[2];
rz(-0.83050767) q[2];
sx q[2];
rz(2.3967957) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.1658926) q[1];
sx q[1];
rz(-2.4411538) q[1];
sx q[1];
rz(-2.9794934) q[1];
rz(1.068088) q[3];
sx q[3];
rz(-1.5912676) q[3];
sx q[3];
rz(-0.81077829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.3669746) q[2];
sx q[2];
rz(-2.8145511) q[2];
sx q[2];
rz(2.7152087) q[2];
rz(1.9042227) q[3];
sx q[3];
rz(-2.5137413) q[3];
sx q[3];
rz(-0.0330851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8957829) q[0];
sx q[0];
rz(-1.3170467) q[0];
sx q[0];
rz(-2.202503) q[0];
rz(0.89871961) q[1];
sx q[1];
rz(-2.6627314) q[1];
sx q[1];
rz(-0.59392196) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1547326) q[0];
sx q[0];
rz(-1.8911456) q[0];
sx q[0];
rz(-2.8264168) q[0];
rz(-0.8823231) q[2];
sx q[2];
rz(-1.6154628) q[2];
sx q[2];
rz(-2.4633212) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.18332874) q[1];
sx q[1];
rz(-1.7839583) q[1];
sx q[1];
rz(-2.0002685) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2265615) q[3];
sx q[3];
rz(-1.1262745) q[3];
sx q[3];
rz(-1.3821186) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.5014191) q[2];
sx q[2];
rz(-2.4121425) q[2];
sx q[2];
rz(-1.7017986) q[2];
rz(2.7539608) q[3];
sx q[3];
rz(-1.6250316) q[3];
sx q[3];
rz(0.38813996) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7664292) q[0];
sx q[0];
rz(-1.5901934) q[0];
sx q[0];
rz(-0.50278062) q[0];
rz(0.76820961) q[1];
sx q[1];
rz(-2.6380824) q[1];
sx q[1];
rz(-0.75685135) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2198974) q[0];
sx q[0];
rz(-1.9031525) q[0];
sx q[0];
rz(2.7601348) q[0];
x q[1];
rz(3.1286131) q[2];
sx q[2];
rz(-2.0364967) q[2];
sx q[2];
rz(-2.4109858) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.6824324) q[1];
sx q[1];
rz(-1.9285413) q[1];
sx q[1];
rz(2.1898502) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4079011) q[3];
sx q[3];
rz(-1.9568223) q[3];
sx q[3];
rz(1.2804077) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.7148774) q[2];
sx q[2];
rz(-1.2299512) q[2];
sx q[2];
rz(-1.654401) q[2];
rz(2.5590844) q[3];
sx q[3];
rz(-2.0472066) q[3];
sx q[3];
rz(0.55707651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29397598) q[0];
sx q[0];
rz(-1.0852381) q[0];
sx q[0];
rz(2.3838682) q[0];
rz(1.2879397) q[1];
sx q[1];
rz(-2.2133591) q[1];
sx q[1];
rz(-2.0910738) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3462853) q[0];
sx q[0];
rz(-1.2588358) q[0];
sx q[0];
rz(-0.26766582) q[0];
x q[1];
rz(2.6987223) q[2];
sx q[2];
rz(-1.0204698) q[2];
sx q[2];
rz(-2.0572822) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.037462385) q[1];
sx q[1];
rz(-2.3298652) q[1];
sx q[1];
rz(-0.8123668) q[1];
x q[2];
rz(0.80661185) q[3];
sx q[3];
rz(-1.9730554) q[3];
sx q[3];
rz(0.81392399) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.22333764) q[2];
sx q[2];
rz(-0.35327062) q[2];
sx q[2];
rz(-2.5081432) q[2];
rz(1.9472306) q[3];
sx q[3];
rz(-1.4712237) q[3];
sx q[3];
rz(2.4244394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.441992) q[0];
sx q[0];
rz(-2.6514335) q[0];
sx q[0];
rz(2.8884086) q[0];
rz(1.5340012) q[1];
sx q[1];
rz(-1.4350767) q[1];
sx q[1];
rz(-1.4621428) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40186858) q[0];
sx q[0];
rz(-1.3603633) q[0];
sx q[0];
rz(-1.6519288) q[0];
x q[1];
rz(-0.24121933) q[2];
sx q[2];
rz(-2.2058645) q[2];
sx q[2];
rz(-1.6068174) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.51965442) q[1];
sx q[1];
rz(-0.42236537) q[1];
sx q[1];
rz(0.42610355) q[1];
rz(1.2268279) q[3];
sx q[3];
rz(-1.868639) q[3];
sx q[3];
rz(-2.1732268) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.9399461) q[2];
sx q[2];
rz(-2.3942409) q[2];
sx q[2];
rz(2.3366826) q[2];
rz(-1.1770052) q[3];
sx q[3];
rz(-1.0369119) q[3];
sx q[3];
rz(-2.0578407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8686304) q[0];
sx q[0];
rz(-1.0694163) q[0];
sx q[0];
rz(2.2139363) q[0];
rz(1.0246798) q[1];
sx q[1];
rz(-1.506348) q[1];
sx q[1];
rz(-2.129508) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1597848) q[0];
sx q[0];
rz(-2.4111528) q[0];
sx q[0];
rz(0.83321379) q[0];
rz(-pi) q[1];
rz(1.7094678) q[2];
sx q[2];
rz(-1.9562634) q[2];
sx q[2];
rz(0.0027545714) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.5092897) q[1];
sx q[1];
rz(-1.0273233) q[1];
sx q[1];
rz(1.5322881) q[1];
rz(-pi) q[2];
rz(1.4885694) q[3];
sx q[3];
rz(-2.1279018) q[3];
sx q[3];
rz(-1.132387) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.6107789) q[2];
sx q[2];
rz(-1.6616219) q[2];
sx q[2];
rz(-0.42993316) q[2];
rz(-1.0144462) q[3];
sx q[3];
rz(-2.7323664) q[3];
sx q[3];
rz(2.608192) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6124509) q[0];
sx q[0];
rz(-0.93944678) q[0];
sx q[0];
rz(-2.9274143) q[0];
rz(1.051349) q[1];
sx q[1];
rz(-2.9290757) q[1];
sx q[1];
rz(-0.28373757) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3389694) q[0];
sx q[0];
rz(-1.6049275) q[0];
sx q[0];
rz(2.8404833) q[0];
x q[1];
rz(2.9494638) q[2];
sx q[2];
rz(-1.2461975) q[2];
sx q[2];
rz(-2.0786376) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.6730496) q[1];
sx q[1];
rz(-2.0646411) q[1];
sx q[1];
rz(-1.7935351) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.92054263) q[3];
sx q[3];
rz(-2.1170108) q[3];
sx q[3];
rz(-1.9073245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.6909137) q[2];
sx q[2];
rz(-2.6288855) q[2];
sx q[2];
rz(-1.696375) q[2];
rz(-1.5444267) q[3];
sx q[3];
rz(-1.701136) q[3];
sx q[3];
rz(-2.8022695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63672367) q[0];
sx q[0];
rz(-3.0915785) q[0];
sx q[0];
rz(-3.0723363) q[0];
rz(-1.6537369) q[1];
sx q[1];
rz(-1.2830877) q[1];
sx q[1];
rz(-1.5690631) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5208961) q[0];
sx q[0];
rz(-0.81917742) q[0];
sx q[0];
rz(2.2197414) q[0];
x q[1];
rz(-1.4885159) q[2];
sx q[2];
rz(-1.8956208) q[2];
sx q[2];
rz(-0.61330739) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.931817) q[1];
sx q[1];
rz(-1.4869542) q[1];
sx q[1];
rz(-0.2340338) q[1];
rz(-pi) q[2];
rz(-0.2089573) q[3];
sx q[3];
rz(-2.2088802) q[3];
sx q[3];
rz(1.0554505) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.9562324) q[2];
sx q[2];
rz(-2.9113443) q[2];
sx q[2];
rz(-0.13988477) q[2];
rz(0.36758962) q[3];
sx q[3];
rz(-1.9544173) q[3];
sx q[3];
rz(-0.99115133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96520987) q[0];
sx q[0];
rz(-2.7503224) q[0];
sx q[0];
rz(-0.64176732) q[0];
rz(-1.9104674) q[1];
sx q[1];
rz(-1.1522013) q[1];
sx q[1];
rz(-0.26783255) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56775996) q[0];
sx q[0];
rz(-1.1288252) q[0];
sx q[0];
rz(2.3964336) q[0];
rz(-3.015976) q[2];
sx q[2];
rz(-1.3889424) q[2];
sx q[2];
rz(-0.4609209) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.2512868) q[1];
sx q[1];
rz(-0.44154134) q[1];
sx q[1];
rz(0.73015405) q[1];
rz(-1.0288826) q[3];
sx q[3];
rz(-1.8174603) q[3];
sx q[3];
rz(2.820462) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.3315167) q[2];
sx q[2];
rz(-1.5248359) q[2];
sx q[2];
rz(1.4769185) q[2];
rz(2.8752575) q[3];
sx q[3];
rz(-0.24644066) q[3];
sx q[3];
rz(2.5951071) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1289566) q[0];
sx q[0];
rz(-2.2205882) q[0];
sx q[0];
rz(0.90482774) q[0];
rz(0.77990445) q[1];
sx q[1];
rz(-2.654568) q[1];
sx q[1];
rz(1.754896) q[1];
rz(-2.4308464) q[2];
sx q[2];
rz(-1.1087316) q[2];
sx q[2];
rz(1.0089594) q[2];
rz(1.7846617) q[3];
sx q[3];
rz(-2.2371117) q[3];
sx q[3];
rz(-0.37123751) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];