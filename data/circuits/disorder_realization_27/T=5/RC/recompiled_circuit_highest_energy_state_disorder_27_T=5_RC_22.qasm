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
rz(1.2394387) q[0];
sx q[0];
rz(-1.8128938) q[0];
sx q[0];
rz(-0.21188307) q[0];
rz(-1.4172685) q[1];
sx q[1];
rz(-2.6091726) q[1];
sx q[1];
rz(0.37791696) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0402474) q[0];
sx q[0];
rz(-1.5703526) q[0];
sx q[0];
rz(-1.5612649) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9340408) q[2];
sx q[2];
rz(-0.9245199) q[2];
sx q[2];
rz(-2.5763047) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.35740543) q[1];
sx q[1];
rz(-2.2012246) q[1];
sx q[1];
rz(-1.9131843) q[1];
x q[2];
rz(2.3842281) q[3];
sx q[3];
rz(-2.2038282) q[3];
sx q[3];
rz(2.4773134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.8523031) q[2];
sx q[2];
rz(-2.0630344) q[2];
sx q[2];
rz(-0.079785384) q[2];
rz(2.1763109) q[3];
sx q[3];
rz(-1.4502757) q[3];
sx q[3];
rz(-2.4108346) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6642283) q[0];
sx q[0];
rz(-2.1381162) q[0];
sx q[0];
rz(3.0526414) q[0];
rz(-1.3720007) q[1];
sx q[1];
rz(-1.4274495) q[1];
sx q[1];
rz(1.1044097) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8744132) q[0];
sx q[0];
rz(-0.66775087) q[0];
sx q[0];
rz(-2.7655168) q[0];
rz(-pi) q[1];
rz(-2.3588793) q[2];
sx q[2];
rz(-1.9316976) q[2];
sx q[2];
rz(2.1825143) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.6767533) q[1];
sx q[1];
rz(-2.0130806) q[1];
sx q[1];
rz(-0.16525903) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4711963) q[3];
sx q[3];
rz(-1.0336813) q[3];
sx q[3];
rz(2.296534) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.264512) q[2];
sx q[2];
rz(-2.0161714) q[2];
sx q[2];
rz(-1.6748927) q[2];
rz(0.029021164) q[3];
sx q[3];
rz(-1.0163739) q[3];
sx q[3];
rz(-0.69028729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2324209) q[0];
sx q[0];
rz(-0.066815289) q[0];
sx q[0];
rz(2.8636041) q[0];
rz(1.7104644) q[1];
sx q[1];
rz(-2.2093096) q[1];
sx q[1];
rz(-0.40036449) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58996449) q[0];
sx q[0];
rz(-1.1271521) q[0];
sx q[0];
rz(-0.69728627) q[0];
x q[1];
rz(-0.68636559) q[2];
sx q[2];
rz(-1.7698506) q[2];
sx q[2];
rz(-0.65716568) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.4440205) q[1];
sx q[1];
rz(-2.6409147) q[1];
sx q[1];
rz(2.8050007) q[1];
rz(-pi) q[2];
x q[2];
rz(0.72680803) q[3];
sx q[3];
rz(-0.44379674) q[3];
sx q[3];
rz(-0.82647317) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.6843188) q[2];
sx q[2];
rz(-1.6088586) q[2];
sx q[2];
rz(-0.45480967) q[2];
rz(1.8999892) q[3];
sx q[3];
rz(-2.1865215) q[3];
sx q[3];
rz(-2.4212867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1860745) q[0];
sx q[0];
rz(-1.3294514) q[0];
sx q[0];
rz(-0.00051001471) q[0];
rz(-2.5406802) q[1];
sx q[1];
rz(-0.83952236) q[1];
sx q[1];
rz(2.9972163) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91286425) q[0];
sx q[0];
rz(-2.1179541) q[0];
sx q[0];
rz(-2.8990101) q[0];
rz(-pi) q[1];
rz(1.7208485) q[2];
sx q[2];
rz(-1.320082) q[2];
sx q[2];
rz(1.5992407) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.2541009) q[1];
sx q[1];
rz(-2.2057167) q[1];
sx q[1];
rz(1.8545521) q[1];
rz(0.5471121) q[3];
sx q[3];
rz(-1.662385) q[3];
sx q[3];
rz(1.509059) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.1986177) q[2];
sx q[2];
rz(-1.696442) q[2];
sx q[2];
rz(2.1790738) q[2];
rz(1.6019542) q[3];
sx q[3];
rz(-1.736015) q[3];
sx q[3];
rz(-0.34077728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2365504) q[0];
sx q[0];
rz(-1.6860697) q[0];
sx q[0];
rz(-1.1849674) q[0];
rz(0.22625893) q[1];
sx q[1];
rz(-0.87892756) q[1];
sx q[1];
rz(-2.102898) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5796367) q[0];
sx q[0];
rz(-0.59369864) q[0];
sx q[0];
rz(-0.14621347) q[0];
rz(-pi) q[1];
rz(2.4657927) q[2];
sx q[2];
rz(-2.2797109) q[2];
sx q[2];
rz(2.076864) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.86335582) q[1];
sx q[1];
rz(-2.0788631) q[1];
sx q[1];
rz(2.3418576) q[1];
x q[2];
rz(1.1306954) q[3];
sx q[3];
rz(-1.6301042) q[3];
sx q[3];
rz(1.9594994) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.7096536) q[2];
sx q[2];
rz(-0.69810549) q[2];
sx q[2];
rz(0.31614885) q[2];
rz(-1.6759253) q[3];
sx q[3];
rz(-2.0114055) q[3];
sx q[3];
rz(-1.3795615) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50452152) q[0];
sx q[0];
rz(-2.2383454) q[0];
sx q[0];
rz(-1.1035408) q[0];
rz(0.48031131) q[1];
sx q[1];
rz(-0.66245285) q[1];
sx q[1];
rz(0.59741098) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7561853) q[0];
sx q[0];
rz(-0.28235897) q[0];
sx q[0];
rz(2.307111) q[0];
rz(0.2603724) q[2];
sx q[2];
rz(-0.7625167) q[2];
sx q[2];
rz(1.1077566) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.5525517) q[1];
sx q[1];
rz(-1.3197761) q[1];
sx q[1];
rz(-1.4000721) q[1];
x q[2];
rz(-1.1998981) q[3];
sx q[3];
rz(-1.6565595) q[3];
sx q[3];
rz(-2.896912) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.9395113) q[2];
sx q[2];
rz(-1.5152405) q[2];
sx q[2];
rz(-1.0763947) q[2];
rz(-1.1427897) q[3];
sx q[3];
rz(-0.79311526) q[3];
sx q[3];
rz(0.49016652) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48653212) q[0];
sx q[0];
rz(-1.158411) q[0];
sx q[0];
rz(-0.038473815) q[0];
rz(-0.064727457) q[1];
sx q[1];
rz(-1.3836626) q[1];
sx q[1];
rz(0.23385349) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42002871) q[0];
sx q[0];
rz(-1.7917243) q[0];
sx q[0];
rz(1.3807644) q[0];
rz(-2.8296956) q[2];
sx q[2];
rz(-1.2205178) q[2];
sx q[2];
rz(2.6147982) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.95215248) q[1];
sx q[1];
rz(-0.92467151) q[1];
sx q[1];
rz(1.4586071) q[1];
rz(-pi) q[2];
rz(-0.12988083) q[3];
sx q[3];
rz(-1.6280773) q[3];
sx q[3];
rz(-2.2735325) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.6528299) q[2];
sx q[2];
rz(-1.5583928) q[2];
sx q[2];
rz(-2.5433507) q[2];
rz(3.0277142) q[3];
sx q[3];
rz(-1.3860393) q[3];
sx q[3];
rz(-2.2940476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7365731) q[0];
sx q[0];
rz(-0.84252715) q[0];
sx q[0];
rz(2.8952428) q[0];
rz(-1.8440638) q[1];
sx q[1];
rz(-2.0228701) q[1];
sx q[1];
rz(0.49682239) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6778567) q[0];
sx q[0];
rz(-1.1564768) q[0];
sx q[0];
rz(0.68092771) q[0];
rz(-pi) q[1];
rz(-1.4250579) q[2];
sx q[2];
rz(-2.2706804) q[2];
sx q[2];
rz(0.041797195) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.5649181) q[1];
sx q[1];
rz(-1.0905301) q[1];
sx q[1];
rz(2.8757921) q[1];
rz(-1.854135) q[3];
sx q[3];
rz(-0.50931286) q[3];
sx q[3];
rz(2.4278502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.42314998) q[2];
sx q[2];
rz(-2.6079874) q[2];
sx q[2];
rz(-1.9971087) q[2];
rz(-1.9874969) q[3];
sx q[3];
rz(-1.7172979) q[3];
sx q[3];
rz(2.9437039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6941187) q[0];
sx q[0];
rz(-0.23451528) q[0];
sx q[0];
rz(-0.06614729) q[0];
rz(-1.9937531) q[1];
sx q[1];
rz(-1.7639953) q[1];
sx q[1];
rz(-2.537421) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19949958) q[0];
sx q[0];
rz(-1.2954933) q[0];
sx q[0];
rz(1.3574187) q[0];
rz(-pi) q[1];
rz(-1.4140903) q[2];
sx q[2];
rz(-2.092166) q[2];
sx q[2];
rz(-0.40796134) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.79163247) q[1];
sx q[1];
rz(-0.59795982) q[1];
sx q[1];
rz(-2.7954742) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.82138942) q[3];
sx q[3];
rz(-2.2940972) q[3];
sx q[3];
rz(2.5625474) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.8784647) q[2];
sx q[2];
rz(-0.85004127) q[2];
sx q[2];
rz(-2.7086332) q[2];
rz(1.912502) q[3];
sx q[3];
rz(-1.9357598) q[3];
sx q[3];
rz(-1.3273201) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1935254) q[0];
sx q[0];
rz(-2.0663517) q[0];
sx q[0];
rz(1.2731592) q[0];
rz(-2.676447) q[1];
sx q[1];
rz(-1.7756614) q[1];
sx q[1];
rz(-0.21496162) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0419768) q[0];
sx q[0];
rz(-0.93659217) q[0];
sx q[0];
rz(-1.8862104) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6097492) q[2];
sx q[2];
rz(-1.9491247) q[2];
sx q[2];
rz(-1.7500306) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7219639) q[1];
sx q[1];
rz(-0.89744324) q[1];
sx q[1];
rz(-0.90760214) q[1];
rz(-pi) q[2];
rz(-0.56193476) q[3];
sx q[3];
rz(-1.3519545) q[3];
sx q[3];
rz(-0.17499017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.8421858) q[2];
sx q[2];
rz(-0.81373787) q[2];
sx q[2];
rz(-0.95883933) q[2];
rz(1.1622608) q[3];
sx q[3];
rz(-1.6543038) q[3];
sx q[3];
rz(1.6963262) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
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
rz(1.5398298) q[0];
sx q[0];
rz(-1.5400664) q[0];
sx q[0];
rz(-1.6590317) q[0];
rz(2.4330347) q[1];
sx q[1];
rz(-0.19150145) q[1];
sx q[1];
rz(2.3932744) q[1];
rz(-0.49025771) q[2];
sx q[2];
rz(-0.62050642) q[2];
sx q[2];
rz(-1.4668589) q[2];
rz(1.7583354) q[3];
sx q[3];
rz(-2.2839727) q[3];
sx q[3];
rz(-1.4061389) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
