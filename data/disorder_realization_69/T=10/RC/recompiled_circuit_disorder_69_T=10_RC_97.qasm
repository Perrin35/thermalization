OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.7944613) q[0];
sx q[0];
rz(-2.1262655) q[0];
sx q[0];
rz(-0.46749687) q[0];
rz(-0.52019083) q[1];
sx q[1];
rz(1.7953035) q[1];
sx q[1];
rz(10.105159) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2517393) q[0];
sx q[0];
rz(-1.7008874) q[0];
sx q[0];
rz(2.2962909) q[0];
x q[1];
rz(0.99256398) q[2];
sx q[2];
rz(-2.7928068) q[2];
sx q[2];
rz(1.055583) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.2441794) q[1];
sx q[1];
rz(-1.3379828) q[1];
sx q[1];
rz(-2.0884872) q[1];
rz(-pi) q[2];
rz(1.3939788) q[3];
sx q[3];
rz(-1.8999945) q[3];
sx q[3];
rz(-0.65229177) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.15158571) q[2];
sx q[2];
rz(-2.1534584) q[2];
sx q[2];
rz(0.087466784) q[2];
rz(-0.72922373) q[3];
sx q[3];
rz(-0.37344033) q[3];
sx q[3];
rz(-0.27174404) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7049578) q[0];
sx q[0];
rz(-0.74902642) q[0];
sx q[0];
rz(2.1402284) q[0];
rz(0.17240605) q[1];
sx q[1];
rz(-1.1162076) q[1];
sx q[1];
rz(-0.52406812) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.24633) q[0];
sx q[0];
rz(-1.2753914) q[0];
sx q[0];
rz(-2.2569879) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7153154) q[2];
sx q[2];
rz(-0.91765109) q[2];
sx q[2];
rz(-0.041235812) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.9429051) q[1];
sx q[1];
rz(-2.6932979) q[1];
sx q[1];
rz(-0.56221902) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4161795) q[3];
sx q[3];
rz(-1.039045) q[3];
sx q[3];
rz(-2.1173409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.15741631) q[2];
sx q[2];
rz(-2.6817862) q[2];
sx q[2];
rz(-1.8015507) q[2];
rz(-0.79483461) q[3];
sx q[3];
rz(-2.0017616) q[3];
sx q[3];
rz(3.1047344) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79725093) q[0];
sx q[0];
rz(-2.4448555) q[0];
sx q[0];
rz(2.5168193) q[0];
rz(0.96427381) q[1];
sx q[1];
rz(-0.48502973) q[1];
sx q[1];
rz(2.952081) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7150869) q[0];
sx q[0];
rz(-0.77930342) q[0];
sx q[0];
rz(-2.7515609) q[0];
x q[1];
rz(2.7663642) q[2];
sx q[2];
rz(-1.669075) q[2];
sx q[2];
rz(0.74429846) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.5513969) q[1];
sx q[1];
rz(-1.8101748) q[1];
sx q[1];
rz(1.6281284) q[1];
rz(-2.1720445) q[3];
sx q[3];
rz(-1.7161955) q[3];
sx q[3];
rz(1.5090514) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.0345962) q[2];
sx q[2];
rz(-0.84356374) q[2];
sx q[2];
rz(-1.754388) q[2];
rz(-0.43131367) q[3];
sx q[3];
rz(-1.8547736) q[3];
sx q[3];
rz(-1.0881902) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65961924) q[0];
sx q[0];
rz(-0.94835931) q[0];
sx q[0];
rz(1.5959928) q[0];
rz(-2.003147) q[1];
sx q[1];
rz(-2.3661416) q[1];
sx q[1];
rz(-2.0514354) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2290303) q[0];
sx q[0];
rz(-1.9711442) q[0];
sx q[0];
rz(-1.1926665) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1074739) q[2];
sx q[2];
rz(-0.82674971) q[2];
sx q[2];
rz(0.85988753) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.480099) q[1];
sx q[1];
rz(-1.9839994) q[1];
sx q[1];
rz(2.5109629) q[1];
x q[2];
rz(0.4269883) q[3];
sx q[3];
rz(-1.7592351) q[3];
sx q[3];
rz(-0.45005902) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1426992) q[2];
sx q[2];
rz(-0.45140758) q[2];
sx q[2];
rz(-0.85582716) q[2];
rz(-1.9479729) q[3];
sx q[3];
rz(-1.6222298) q[3];
sx q[3];
rz(0.91317552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6511433) q[0];
sx q[0];
rz(-2.1676846) q[0];
sx q[0];
rz(-0.14973101) q[0];
rz(0.99114746) q[1];
sx q[1];
rz(-1.2065572) q[1];
sx q[1];
rz(-1.1700464) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4867602) q[0];
sx q[0];
rz(-2.1954384) q[0];
sx q[0];
rz(-2.7420068) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9636376) q[2];
sx q[2];
rz(-0.71937865) q[2];
sx q[2];
rz(-2.9658085) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3586732) q[1];
sx q[1];
rz(-1.1766608) q[1];
sx q[1];
rz(2.6637117) q[1];
rz(-pi) q[2];
rz(1.7647469) q[3];
sx q[3];
rz(-2.100482) q[3];
sx q[3];
rz(2.291631) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.871792) q[2];
sx q[2];
rz(-1.5912) q[2];
sx q[2];
rz(-0.13723792) q[2];
rz(1.8042701) q[3];
sx q[3];
rz(-0.58115712) q[3];
sx q[3];
rz(-1.2924682) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18734922) q[0];
sx q[0];
rz(-1.7631148) q[0];
sx q[0];
rz(1.7911918) q[0];
rz(2.2987135) q[1];
sx q[1];
rz(-0.73892361) q[1];
sx q[1];
rz(0.67289105) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85149375) q[0];
sx q[0];
rz(-2.5585472) q[0];
sx q[0];
rz(1.7593988) q[0];
x q[1];
rz(3.1168675) q[2];
sx q[2];
rz(-0.85660663) q[2];
sx q[2];
rz(-1.6753472) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.8012961) q[1];
sx q[1];
rz(-0.9749037) q[1];
sx q[1];
rz(-1.3188386) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4937431) q[3];
sx q[3];
rz(-1.9814081) q[3];
sx q[3];
rz(1.3239469) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.3488397) q[2];
sx q[2];
rz(-1.9323843) q[2];
sx q[2];
rz(2.7065281) q[2];
rz(-1.3600291) q[3];
sx q[3];
rz(-0.74917787) q[3];
sx q[3];
rz(0.24766651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.172794) q[0];
sx q[0];
rz(-2.2491169) q[0];
sx q[0];
rz(1.0621747) q[0];
rz(2.0299714) q[1];
sx q[1];
rz(-1.9042791) q[1];
sx q[1];
rz(-1.7395082) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3081449) q[0];
sx q[0];
rz(-1.2045367) q[0];
sx q[0];
rz(1.2783947) q[0];
rz(-1.8066605) q[2];
sx q[2];
rz(-2.3216322) q[2];
sx q[2];
rz(1.5136528) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.2981616) q[1];
sx q[1];
rz(-1.2097675) q[1];
sx q[1];
rz(1.4069188) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3659031) q[3];
sx q[3];
rz(-2.5435102) q[3];
sx q[3];
rz(0.3860592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.7193675) q[2];
sx q[2];
rz(-2.5740467) q[2];
sx q[2];
rz(2.4613703) q[2];
rz(-2.7133572) q[3];
sx q[3];
rz(-1.8869583) q[3];
sx q[3];
rz(-1.359882) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5087886) q[0];
sx q[0];
rz(-1.7254242) q[0];
sx q[0];
rz(1.5493786) q[0];
rz(0.24958615) q[1];
sx q[1];
rz(-1.9892178) q[1];
sx q[1];
rz(0.54135281) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4562948) q[0];
sx q[0];
rz(-1.0066427) q[0];
sx q[0];
rz(1.5575404) q[0];
rz(-1.2301684) q[2];
sx q[2];
rz(-1.1414141) q[2];
sx q[2];
rz(2.9609749) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.036877) q[1];
sx q[1];
rz(-2.2762205) q[1];
sx q[1];
rz(0.55286644) q[1];
rz(1.5726611) q[3];
sx q[3];
rz(-1.1904753) q[3];
sx q[3];
rz(-0.62455458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(3.0970739) q[2];
sx q[2];
rz(-0.35733435) q[2];
sx q[2];
rz(0.77735916) q[2];
rz(0.85123953) q[3];
sx q[3];
rz(-2.080353) q[3];
sx q[3];
rz(-1.8204934) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7580496) q[0];
sx q[0];
rz(-1.8109011) q[0];
sx q[0];
rz(0.0099649075) q[0];
rz(-1.0154356) q[1];
sx q[1];
rz(-0.76428691) q[1];
sx q[1];
rz(-1.4452971) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24628595) q[0];
sx q[0];
rz(-2.200897) q[0];
sx q[0];
rz(-2.2340328) q[0];
rz(-2.2465835) q[2];
sx q[2];
rz(-1.0969321) q[2];
sx q[2];
rz(2.5788139) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.6362793) q[1];
sx q[1];
rz(-1.5552102) q[1];
sx q[1];
rz(-2.511335) q[1];
x q[2];
rz(2.5310998) q[3];
sx q[3];
rz(-0.99184147) q[3];
sx q[3];
rz(1.0510774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.4426667) q[2];
sx q[2];
rz(-0.93124229) q[2];
sx q[2];
rz(-2.5342069) q[2];
rz(-1.7025042) q[3];
sx q[3];
rz(-1.3834229) q[3];
sx q[3];
rz(0.79090345) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33567515) q[0];
sx q[0];
rz(-1.6573925) q[0];
sx q[0];
rz(0.6138531) q[0];
rz(-1.0461668) q[1];
sx q[1];
rz(-0.26509735) q[1];
sx q[1];
rz(0.39224958) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0454355) q[0];
sx q[0];
rz(-0.7924315) q[0];
sx q[0];
rz(-2.5655454) q[0];
x q[1];
rz(0.73108436) q[2];
sx q[2];
rz(-0.77027551) q[2];
sx q[2];
rz(2.0172271) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.85234648) q[1];
sx q[1];
rz(-2.1834071) q[1];
sx q[1];
rz(-2.5958519) q[1];
rz(-pi) q[2];
rz(-0.13855374) q[3];
sx q[3];
rz(-0.20692736) q[3];
sx q[3];
rz(-3.0272527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.067651) q[2];
sx q[2];
rz(-1.3839046) q[2];
sx q[2];
rz(-0.60047853) q[2];
rz(1.0673808) q[3];
sx q[3];
rz(-1.3200656) q[3];
sx q[3];
rz(-2.0480806) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7077211) q[0];
sx q[0];
rz(-2.7086471) q[0];
sx q[0];
rz(1.4824296) q[0];
rz(2.1451163) q[1];
sx q[1];
rz(-1.6468208) q[1];
sx q[1];
rz(-1.5368808) q[1];
rz(2.375013) q[2];
sx q[2];
rz(-1.8324413) q[2];
sx q[2];
rz(1.2698297) q[2];
rz(2.2157833) q[3];
sx q[3];
rz(-2.3498597) q[3];
sx q[3];
rz(1.3510977) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];