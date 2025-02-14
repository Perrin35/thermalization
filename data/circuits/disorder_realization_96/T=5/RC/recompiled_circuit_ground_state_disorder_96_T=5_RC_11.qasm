OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.191303) q[0];
sx q[0];
rz(-2.8714955) q[0];
sx q[0];
rz(-0.88859963) q[0];
rz(-1.3130045) q[1];
sx q[1];
rz(-1.5994025) q[1];
sx q[1];
rz(-1.3808274) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.439405) q[0];
sx q[0];
rz(-0.26964615) q[0];
sx q[0];
rz(0.7650956) q[0];
rz(0.51989748) q[2];
sx q[2];
rz(-0.87018229) q[2];
sx q[2];
rz(1.1288647) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.0128768) q[1];
sx q[1];
rz(-2.8059462) q[1];
sx q[1];
rz(3.1124093) q[1];
rz(0.22486658) q[3];
sx q[3];
rz(-1.9266085) q[3];
sx q[3];
rz(2.0840621) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.8218653) q[2];
sx q[2];
rz(-1.3250019) q[2];
sx q[2];
rz(-0.74903178) q[2];
rz(2.9876515) q[3];
sx q[3];
rz(-1.0356244) q[3];
sx q[3];
rz(0.099893959) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.446949) q[0];
sx q[0];
rz(-1.4780937) q[0];
sx q[0];
rz(-1.2167759) q[0];
rz(-2.079839) q[1];
sx q[1];
rz(-0.80445015) q[1];
sx q[1];
rz(-0.42713508) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0674853) q[0];
sx q[0];
rz(-2.4541313) q[0];
sx q[0];
rz(0.56337728) q[0];
x q[1];
rz(0.92556503) q[2];
sx q[2];
rz(-1.3750018) q[2];
sx q[2];
rz(0.93709556) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4915575) q[1];
sx q[1];
rz(-1.2833529) q[1];
sx q[1];
rz(-2.9950055) q[1];
rz(-pi) q[2];
rz(-1.3668381) q[3];
sx q[3];
rz(-1.5195623) q[3];
sx q[3];
rz(2.8683087) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.079166807) q[2];
sx q[2];
rz(-1.7288952) q[2];
sx q[2];
rz(-1.3986577) q[2];
rz(-2.9178197) q[3];
sx q[3];
rz(-0.90884915) q[3];
sx q[3];
rz(1.5435425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-0.021521213) q[0];
sx q[0];
rz(-1.3503617) q[0];
sx q[0];
rz(0.045510005) q[0];
rz(1.1687763) q[1];
sx q[1];
rz(-1.903542) q[1];
sx q[1];
rz(2.5659836) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2164477) q[0];
sx q[0];
rz(-2.5024104) q[0];
sx q[0];
rz(0.018763824) q[0];
rz(-pi) q[1];
x q[1];
rz(0.91012886) q[2];
sx q[2];
rz(-1.5572539) q[2];
sx q[2];
rz(-0.34156583) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.2896881) q[1];
sx q[1];
rz(-0.59483268) q[1];
sx q[1];
rz(-1.965305) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.1166385) q[3];
sx q[3];
rz(-0.91812274) q[3];
sx q[3];
rz(-3.0252473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.0098972926) q[2];
sx q[2];
rz(-0.74247777) q[2];
sx q[2];
rz(0.95727813) q[2];
rz(2.5668528) q[3];
sx q[3];
rz(-1.5301907) q[3];
sx q[3];
rz(-1.8360229) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15605536) q[0];
sx q[0];
rz(-2.1888581) q[0];
sx q[0];
rz(-1.8804469) q[0];
rz(0.082322923) q[1];
sx q[1];
rz(-2.0539093) q[1];
sx q[1];
rz(1.7877158) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60786965) q[0];
sx q[0];
rz(-1.3114616) q[0];
sx q[0];
rz(-1.3277256) q[0];
rz(-pi) q[1];
rz(1.7548082) q[2];
sx q[2];
rz(-2.0495601) q[2];
sx q[2];
rz(1.1637853) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.9685843) q[1];
sx q[1];
rz(-0.49793303) q[1];
sx q[1];
rz(-3.0281316) q[1];
rz(-pi) q[2];
rz(-0.4889826) q[3];
sx q[3];
rz(-2.1900898) q[3];
sx q[3];
rz(-1.796738) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7202619) q[2];
sx q[2];
rz(-2.223184) q[2];
sx q[2];
rz(1.2565695) q[2];
rz(2.4300857) q[3];
sx q[3];
rz(-1.3308176) q[3];
sx q[3];
rz(2.8594657) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8307777) q[0];
sx q[0];
rz(-2.2956235) q[0];
sx q[0];
rz(0.34314439) q[0];
rz(-3.0798196) q[1];
sx q[1];
rz(-0.97007483) q[1];
sx q[1];
rz(-1.7247346) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2767216) q[0];
sx q[0];
rz(-1.7407932) q[0];
sx q[0];
rz(-1.2890588) q[0];
rz(-pi) q[1];
rz(2.4409962) q[2];
sx q[2];
rz(-1.3178409) q[2];
sx q[2];
rz(-0.41159901) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7570863) q[1];
sx q[1];
rz(-1.0913186) q[1];
sx q[1];
rz(-2.9887385) q[1];
x q[2];
rz(2.0569818) q[3];
sx q[3];
rz(-0.47322464) q[3];
sx q[3];
rz(-1.5718232) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.5938277) q[2];
sx q[2];
rz(-1.6739028) q[2];
sx q[2];
rz(1.0859547) q[2];
rz(-0.91935277) q[3];
sx q[3];
rz(-1.7707526) q[3];
sx q[3];
rz(0.61387387) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3980961) q[0];
sx q[0];
rz(-1.9323823) q[0];
sx q[0];
rz(-2.3486775) q[0];
rz(2.4010557) q[1];
sx q[1];
rz(-2.1426327) q[1];
sx q[1];
rz(0.83121306) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3073458) q[0];
sx q[0];
rz(-1.9458658) q[0];
sx q[0];
rz(-0.22386472) q[0];
rz(-pi) q[1];
rz(1.5446051) q[2];
sx q[2];
rz(-0.27715836) q[2];
sx q[2];
rz(-1.8089393) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.97396353) q[1];
sx q[1];
rz(-2.2669753) q[1];
sx q[1];
rz(0.40145282) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6115723) q[3];
sx q[3];
rz(-0.47104657) q[3];
sx q[3];
rz(-0.99179964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.071216019) q[2];
sx q[2];
rz(-1.2717609) q[2];
sx q[2];
rz(2.0224723) q[2];
rz(-1.2707155) q[3];
sx q[3];
rz(-2.0344574) q[3];
sx q[3];
rz(1.3814111) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1374461) q[0];
sx q[0];
rz(-2.9747712) q[0];
sx q[0];
rz(-1.5995837) q[0];
rz(-1.0150602) q[1];
sx q[1];
rz(-1.5559745) q[1];
sx q[1];
rz(-2.9687845) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8625582) q[0];
sx q[0];
rz(-2.1365215) q[0];
sx q[0];
rz(0.59086694) q[0];
rz(-1.9394933) q[2];
sx q[2];
rz(-2.2208344) q[2];
sx q[2];
rz(-1.4024613) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.037464868) q[1];
sx q[1];
rz(-1.3224673) q[1];
sx q[1];
rz(1.1315956) q[1];
x q[2];
rz(0.76877131) q[3];
sx q[3];
rz(-1.988058) q[3];
sx q[3];
rz(2.4989243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.66199866) q[2];
sx q[2];
rz(-2.2981503) q[2];
sx q[2];
rz(2.1638828) q[2];
rz(-2.5665723) q[3];
sx q[3];
rz(-1.2049371) q[3];
sx q[3];
rz(-1.2303111) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32790312) q[0];
sx q[0];
rz(-2.8459025) q[0];
sx q[0];
rz(-3.1258702) q[0];
rz(-1.9533336) q[1];
sx q[1];
rz(-0.63242042) q[1];
sx q[1];
rz(-1.1994919) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9415814) q[0];
sx q[0];
rz(-1.1304378) q[0];
sx q[0];
rz(0.26047996) q[0];
rz(-pi) q[1];
rz(-2.0616777) q[2];
sx q[2];
rz(-1.5245617) q[2];
sx q[2];
rz(-0.85109988) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.6024692) q[1];
sx q[1];
rz(-1.1014465) q[1];
sx q[1];
rz(-2.5825809) q[1];
rz(0.99908546) q[3];
sx q[3];
rz(-1.3093108) q[3];
sx q[3];
rz(-0.7657683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.1477995) q[2];
sx q[2];
rz(-1.2387929) q[2];
sx q[2];
rz(1.7564868) q[2];
rz(-2.5177453) q[3];
sx q[3];
rz(-1.8892989) q[3];
sx q[3];
rz(2.0417716) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1856336) q[0];
sx q[0];
rz(-2.3210242) q[0];
sx q[0];
rz(-2.4483335) q[0];
rz(-0.26501003) q[1];
sx q[1];
rz(-0.82844228) q[1];
sx q[1];
rz(1.5230491) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.023848195) q[0];
sx q[0];
rz(-1.5818412) q[0];
sx q[0];
rz(-2.3680229) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.62317836) q[2];
sx q[2];
rz(-1.0747391) q[2];
sx q[2];
rz(3.0219699) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.48623006) q[1];
sx q[1];
rz(-1.2312447) q[1];
sx q[1];
rz(1.5454253) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.070793666) q[3];
sx q[3];
rz(-0.99081836) q[3];
sx q[3];
rz(-1.9656612) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.8388464) q[2];
sx q[2];
rz(-1.6281717) q[2];
sx q[2];
rz(1.4576853) q[2];
rz(2.1863106) q[3];
sx q[3];
rz(-2.6421319) q[3];
sx q[3];
rz(-2.3063851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38462287) q[0];
sx q[0];
rz(-2.5769825) q[0];
sx q[0];
rz(2.6589174) q[0];
rz(-0.92974281) q[1];
sx q[1];
rz(-2.1371806) q[1];
sx q[1];
rz(-2.7453056) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38448725) q[0];
sx q[0];
rz(-1.1570017) q[0];
sx q[0];
rz(1.1767469) q[0];
x q[1];
rz(-1.7223139) q[2];
sx q[2];
rz(-1.7477525) q[2];
sx q[2];
rz(2.3529476) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.63782802) q[1];
sx q[1];
rz(-0.85496584) q[1];
sx q[1];
rz(-2.9010335) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.40376264) q[3];
sx q[3];
rz(-1.2376533) q[3];
sx q[3];
rz(-0.008283188) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.79260176) q[2];
sx q[2];
rz(-1.4395809) q[2];
sx q[2];
rz(-0.3271884) q[2];
rz(0.78091019) q[3];
sx q[3];
rz(-1.9802997) q[3];
sx q[3];
rz(1.4769295) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
sx q[3];
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
rz(1.1031716) q[0];
sx q[0];
rz(-2.1506943) q[0];
sx q[0];
rz(-2.8839169) q[0];
rz(-0.36956638) q[1];
sx q[1];
rz(-0.8820487) q[1];
sx q[1];
rz(-0.65912156) q[1];
rz(0.054389537) q[2];
sx q[2];
rz(-2.4475606) q[2];
sx q[2];
rz(-2.9771752) q[2];
rz(2.0965602) q[3];
sx q[3];
rz(-2.4637194) q[3];
sx q[3];
rz(1.2135492) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
