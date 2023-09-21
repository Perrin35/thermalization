OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.45286173) q[0];
sx q[0];
rz(-0.17012574) q[0];
sx q[0];
rz(-0.78599077) q[0];
rz(-2.535948) q[1];
sx q[1];
rz(-0.65094596) q[1];
sx q[1];
rz(-2.5075066) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85672985) q[0];
sx q[0];
rz(-0.75548178) q[0];
sx q[0];
rz(1.612624) q[0];
rz(-pi) q[1];
rz(-1.6526821) q[2];
sx q[2];
rz(-0.6559283) q[2];
sx q[2];
rz(1.8413078) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7635599) q[1];
sx q[1];
rz(-1.130192) q[1];
sx q[1];
rz(-0.9159169) q[1];
rz(-3.0862023) q[3];
sx q[3];
rz(-2.7408544) q[3];
sx q[3];
rz(1.3415568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.2259851) q[2];
sx q[2];
rz(-2.6245124) q[2];
sx q[2];
rz(1.8784286) q[2];
rz(1.2849215) q[3];
sx q[3];
rz(-1.6804755) q[3];
sx q[3];
rz(2.8485956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15853515) q[0];
sx q[0];
rz(-1.8479269) q[0];
sx q[0];
rz(2.7040226) q[0];
rz(-0.63105398) q[1];
sx q[1];
rz(-2.8129306) q[1];
sx q[1];
rz(-2.9204869) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76925812) q[0];
sx q[0];
rz(-0.49494574) q[0];
sx q[0];
rz(1.9702205) q[0];
rz(-pi) q[1];
rz(1.2731304) q[2];
sx q[2];
rz(-2.7781099) q[2];
sx q[2];
rz(3.1306981) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.088614956) q[1];
sx q[1];
rz(-0.47514519) q[1];
sx q[1];
rz(0.26762025) q[1];
rz(-pi) q[2];
rz(-2.9651871) q[3];
sx q[3];
rz(-0.50900148) q[3];
sx q[3];
rz(0.41364663) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.9235886) q[2];
sx q[2];
rz(-1.4322586) q[2];
sx q[2];
rz(-0.58829266) q[2];
rz(-0.44899392) q[3];
sx q[3];
rz(-0.42789999) q[3];
sx q[3];
rz(-1.0480405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56851971) q[0];
sx q[0];
rz(-0.097671106) q[0];
sx q[0];
rz(-1.0004689) q[0];
rz(-2.4160642) q[1];
sx q[1];
rz(-1.0745445) q[1];
sx q[1];
rz(-0.75769889) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.706447) q[0];
sx q[0];
rz(-0.62951127) q[0];
sx q[0];
rz(-2.5802617) q[0];
rz(-pi) q[1];
rz(-0.20747848) q[2];
sx q[2];
rz(-1.7277272) q[2];
sx q[2];
rz(-2.4776138) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.1670907) q[1];
sx q[1];
rz(-2.3579512) q[1];
sx q[1];
rz(2.2553315) q[1];
x q[2];
rz(-2.9919639) q[3];
sx q[3];
rz(-0.64390874) q[3];
sx q[3];
rz(0.097188918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.2086601) q[2];
sx q[2];
rz(-0.22298403) q[2];
sx q[2];
rz(2.3051252) q[2];
rz(-1.6992735) q[3];
sx q[3];
rz(-1.6675555) q[3];
sx q[3];
rz(0.20382717) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1214685) q[0];
sx q[0];
rz(-0.080105372) q[0];
sx q[0];
rz(0.51112038) q[0];
rz(0.077443667) q[1];
sx q[1];
rz(-0.42235342) q[1];
sx q[1];
rz(1.3285332) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1140095) q[0];
sx q[0];
rz(-2.8767577) q[0];
sx q[0];
rz(0.94075216) q[0];
rz(-2.0650234) q[2];
sx q[2];
rz(-2.2708714) q[2];
sx q[2];
rz(-1.1717403) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.7698344) q[1];
sx q[1];
rz(-1.8375988) q[1];
sx q[1];
rz(3.1410602) q[1];
x q[2];
rz(2.02416) q[3];
sx q[3];
rz(-0.24586596) q[3];
sx q[3];
rz(2.9256431) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-3.0478583) q[2];
sx q[2];
rz(-1.2183943) q[2];
sx q[2];
rz(2.034534) q[2];
rz(2.6691061) q[3];
sx q[3];
rz(-1.6146336) q[3];
sx q[3];
rz(0.85737491) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
x q[3];
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
rz(0.040314019) q[0];
sx q[0];
rz(-0.89656985) q[0];
sx q[0];
rz(-0.3381981) q[0];
rz(1.2942554) q[1];
sx q[1];
rz(-1.5816403) q[1];
sx q[1];
rz(-0.50450528) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.550068) q[0];
sx q[0];
rz(-1.3378157) q[0];
sx q[0];
rz(-1.0250807) q[0];
rz(-2.4733876) q[2];
sx q[2];
rz(-0.56958157) q[2];
sx q[2];
rz(-0.44797209) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.5096163) q[1];
sx q[1];
rz(-2.7512433) q[1];
sx q[1];
rz(2.2401287) q[1];
rz(-pi) q[2];
rz(-2.3584064) q[3];
sx q[3];
rz(-0.83055701) q[3];
sx q[3];
rz(0.98379788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.0304886) q[2];
sx q[2];
rz(-0.8384853) q[2];
sx q[2];
rz(1.4939235) q[2];
rz(-1.6882287) q[3];
sx q[3];
rz(-2.2036392) q[3];
sx q[3];
rz(-1.3822314) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21022739) q[0];
sx q[0];
rz(-0.87451044) q[0];
sx q[0];
rz(1.0908303) q[0];
rz(2.6018654) q[1];
sx q[1];
rz(-1.9878186) q[1];
sx q[1];
rz(0.18879034) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81297368) q[0];
sx q[0];
rz(-2.6972174) q[0];
sx q[0];
rz(-0.61291738) q[0];
rz(-3.0374755) q[2];
sx q[2];
rz(-1.8407485) q[2];
sx q[2];
rz(1.5013258) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.69953883) q[1];
sx q[1];
rz(-2.8901254) q[1];
sx q[1];
rz(2.8738408) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.24174989) q[3];
sx q[3];
rz(-1.7484192) q[3];
sx q[3];
rz(-1.0279946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.9699036) q[2];
sx q[2];
rz(-2.1128555) q[2];
sx q[2];
rz(1.8264654) q[2];
rz(-1.5054437) q[3];
sx q[3];
rz(-0.35741487) q[3];
sx q[3];
rz(-1.2833387) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0514907) q[0];
sx q[0];
rz(-0.83054709) q[0];
sx q[0];
rz(-2.8175957) q[0];
rz(-1.8404768) q[1];
sx q[1];
rz(-2.3062861) q[1];
sx q[1];
rz(-1.3791929) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24810476) q[0];
sx q[0];
rz(-1.4890492) q[0];
sx q[0];
rz(1.3260613) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4263779) q[2];
sx q[2];
rz(-2.5589802) q[2];
sx q[2];
rz(1.6974534) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.0566237) q[1];
sx q[1];
rz(-1.8693722) q[1];
sx q[1];
rz(-1.741239) q[1];
rz(-pi) q[2];
rz(2.561065) q[3];
sx q[3];
rz(-2.1248098) q[3];
sx q[3];
rz(-2.011727) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(3.1088915) q[2];
sx q[2];
rz(-1.6872493) q[2];
sx q[2];
rz(0.94318715) q[2];
rz(-2.8105248) q[3];
sx q[3];
rz(-1.3676684) q[3];
sx q[3];
rz(2.846472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.023495) q[0];
sx q[0];
rz(-2.5546615) q[0];
sx q[0];
rz(2.3773637) q[0];
rz(-0.14097342) q[1];
sx q[1];
rz(-0.39607513) q[1];
sx q[1];
rz(1.8483298) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9613567) q[0];
sx q[0];
rz(-2.0016252) q[0];
sx q[0];
rz(-0.21813099) q[0];
rz(-2.5404732) q[2];
sx q[2];
rz(-1.6496611) q[2];
sx q[2];
rz(2.2476946) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1845491) q[1];
sx q[1];
rz(-2.7126985) q[1];
sx q[1];
rz(-1.9098319) q[1];
x q[2];
rz(-1.771365) q[3];
sx q[3];
rz(-1.1784394) q[3];
sx q[3];
rz(-1.6605103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.76688898) q[2];
sx q[2];
rz(-2.1024487) q[2];
sx q[2];
rz(-0.58132201) q[2];
rz(0.86822048) q[3];
sx q[3];
rz(-0.41728443) q[3];
sx q[3];
rz(-1.8301615) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54810369) q[0];
sx q[0];
rz(-2.5618401) q[0];
sx q[0];
rz(-1.2506437) q[0];
rz(-0.99682322) q[1];
sx q[1];
rz(-2.2579028) q[1];
sx q[1];
rz(1.2876127) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8087013) q[0];
sx q[0];
rz(-1.5455855) q[0];
sx q[0];
rz(1.0813792) q[0];
x q[1];
rz(-2.655517) q[2];
sx q[2];
rz(-1.4176148) q[2];
sx q[2];
rz(-1.7024405) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.77546706) q[1];
sx q[1];
rz(-0.93535103) q[1];
sx q[1];
rz(2.8833564) q[1];
rz(3.0887293) q[3];
sx q[3];
rz(-0.82517805) q[3];
sx q[3];
rz(-2.4406274) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8210956) q[2];
sx q[2];
rz(-2.7536776) q[2];
sx q[2];
rz(-1.256475) q[2];
rz(0.16658941) q[3];
sx q[3];
rz(-1.5766671) q[3];
sx q[3];
rz(2.0690209) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
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
rz(-0.88084108) q[0];
sx q[0];
rz(-0.80417019) q[0];
sx q[0];
rz(-2.8163731) q[0];
rz(-1.1351769) q[1];
sx q[1];
rz(-2.4644641) q[1];
sx q[1];
rz(-0.36718711) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67883915) q[0];
sx q[0];
rz(-0.05719962) q[0];
sx q[0];
rz(-2.5662533) q[0];
rz(-2.514421) q[2];
sx q[2];
rz(-1.3417202) q[2];
sx q[2];
rz(2.0723745) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.6348833) q[1];
sx q[1];
rz(-2.2193529) q[1];
sx q[1];
rz(-1.8222005) q[1];
rz(-2.2121067) q[3];
sx q[3];
rz(-1.1747642) q[3];
sx q[3];
rz(-3.0505153) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.8964768) q[2];
sx q[2];
rz(-0.8090691) q[2];
sx q[2];
rz(-2.0754576) q[2];
rz(0.079244763) q[3];
sx q[3];
rz(-2.5388122) q[3];
sx q[3];
rz(1.665303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8463678) q[0];
sx q[0];
rz(-1.5867148) q[0];
sx q[0];
rz(-1.5750194) q[0];
rz(-0.29905839) q[1];
sx q[1];
rz(-0.60332861) q[1];
sx q[1];
rz(0.52287846) q[1];
rz(-1.4916186) q[2];
sx q[2];
rz(-2.3186602) q[2];
sx q[2];
rz(3.0849948) q[2];
rz(-2.9885837) q[3];
sx q[3];
rz(-0.85481337) q[3];
sx q[3];
rz(-0.92683642) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
