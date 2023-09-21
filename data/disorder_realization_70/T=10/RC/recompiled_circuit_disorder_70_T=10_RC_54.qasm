OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.6887309) q[0];
sx q[0];
rz(-2.9714669) q[0];
sx q[0];
rz(-2.3556019) q[0];
rz(-2.535948) q[1];
sx q[1];
rz(-0.65094596) q[1];
sx q[1];
rz(0.63408607) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2274322) q[0];
sx q[0];
rz(-0.81613805) q[0];
sx q[0];
rz(-3.1022275) q[0];
rz(-2.2251031) q[2];
sx q[2];
rz(-1.5208897) q[2];
sx q[2];
rz(-2.8061342) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3780327) q[1];
sx q[1];
rz(-1.130192) q[1];
sx q[1];
rz(-0.9159169) q[1];
x q[2];
rz(2.7414054) q[3];
sx q[3];
rz(-1.5923946) q[3];
sx q[3];
rz(-2.9633629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.2259851) q[2];
sx q[2];
rz(-0.51708022) q[2];
sx q[2];
rz(-1.8784286) q[2];
rz(1.8566711) q[3];
sx q[3];
rz(-1.6804755) q[3];
sx q[3];
rz(-2.8485956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-2.9830575) q[0];
sx q[0];
rz(-1.8479269) q[0];
sx q[0];
rz(-0.43757004) q[0];
rz(-0.63105398) q[1];
sx q[1];
rz(-0.32866207) q[1];
sx q[1];
rz(2.9204869) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76925812) q[0];
sx q[0];
rz(-0.49494574) q[0];
sx q[0];
rz(1.1713722) q[0];
x q[1];
rz(-1.2220076) q[2];
sx q[2];
rz(-1.6752599) q[2];
sx q[2];
rz(-1.2806569) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.7213388) q[1];
sx q[1];
rz(-1.6920648) q[1];
sx q[1];
rz(-2.6810357) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-1.709334) q[2];
sx q[2];
rz(0.58829266) q[2];
rz(-0.44899392) q[3];
sx q[3];
rz(-2.7136927) q[3];
sx q[3];
rz(-2.0935521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5730729) q[0];
sx q[0];
rz(-0.097671106) q[0];
sx q[0];
rz(1.0004689) q[0];
rz(2.4160642) q[1];
sx q[1];
rz(-2.0670481) q[1];
sx q[1];
rz(2.3838938) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0962778) q[0];
sx q[0];
rz(-1.0490388) q[0];
sx q[0];
rz(1.9406712) q[0];
x q[1];
rz(1.7311086) q[2];
sx q[2];
rz(-1.3659039) q[2];
sx q[2];
rz(-2.2018873) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.1670907) q[1];
sx q[1];
rz(-2.3579512) q[1];
sx q[1];
rz(-0.88626119) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.63852255) q[3];
sx q[3];
rz(-1.4811852) q[3];
sx q[3];
rz(-1.5479969) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.9329325) q[2];
sx q[2];
rz(-0.22298403) q[2];
sx q[2];
rz(2.3051252) q[2];
rz(1.6992735) q[3];
sx q[3];
rz(-1.6675555) q[3];
sx q[3];
rz(2.9377655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.02012415) q[0];
sx q[0];
rz(-0.080105372) q[0];
sx q[0];
rz(-0.51112038) q[0];
rz(3.064149) q[1];
sx q[1];
rz(-0.42235342) q[1];
sx q[1];
rz(1.8130594) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9851345) q[0];
sx q[0];
rz(-1.4159604) q[0];
sx q[0];
rz(-1.3550718) q[0];
rz(-pi) q[1];
rz(1.0765692) q[2];
sx q[2];
rz(-2.2708714) q[2];
sx q[2];
rz(-1.1717403) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.36973876) q[1];
sx q[1];
rz(-0.26680294) q[1];
sx q[1];
rz(-1.5688483) q[1];
rz(-pi) q[2];
rz(3.0321211) q[3];
sx q[3];
rz(-1.3502035) q[3];
sx q[3];
rz(-0.68144875) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.093734309) q[2];
sx q[2];
rz(-1.9231984) q[2];
sx q[2];
rz(-1.1070586) q[2];
rz(-0.47248653) q[3];
sx q[3];
rz(-1.6146336) q[3];
sx q[3];
rz(0.85737491) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.040314019) q[0];
sx q[0];
rz(-0.89656985) q[0];
sx q[0];
rz(-0.3381981) q[0];
rz(-1.2942554) q[1];
sx q[1];
rz(-1.5816403) q[1];
sx q[1];
rz(0.50450528) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0230334) q[0];
sx q[0];
rz(-2.1001864) q[0];
sx q[0];
rz(-0.27079196) q[0];
rz(-0.46576969) q[2];
sx q[2];
rz(-1.911474) q[2];
sx q[2];
rz(-2.6054232) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5704755) q[1];
sx q[1];
rz(-1.332453) q[1];
sx q[1];
rz(1.8829324) q[1];
rz(-pi) q[2];
rz(-2.3584064) q[3];
sx q[3];
rz(-0.83055701) q[3];
sx q[3];
rz(0.98379788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.111104) q[2];
sx q[2];
rz(-2.3031074) q[2];
sx q[2];
rz(-1.4939235) q[2];
rz(-1.6882287) q[3];
sx q[3];
rz(-2.2036392) q[3];
sx q[3];
rz(-1.3822314) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9313653) q[0];
sx q[0];
rz(-2.2670822) q[0];
sx q[0];
rz(-2.0507623) q[0];
rz(0.53972721) q[1];
sx q[1];
rz(-1.153774) q[1];
sx q[1];
rz(-2.9528023) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.328619) q[0];
sx q[0];
rz(-0.44437528) q[0];
sx q[0];
rz(-2.5286753) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2115057) q[2];
sx q[2];
rz(-0.28887666) q[2];
sx q[2];
rz(1.8747683) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.4420538) q[1];
sx q[1];
rz(-0.25146723) q[1];
sx q[1];
rz(-0.26775189) q[1];
x q[2];
rz(-0.24174989) q[3];
sx q[3];
rz(-1.3931735) q[3];
sx q[3];
rz(1.0279946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.9699036) q[2];
sx q[2];
rz(-1.0287372) q[2];
sx q[2];
rz(-1.3151273) q[2];
rz(-1.6361489) q[3];
sx q[3];
rz(-0.35741487) q[3];
sx q[3];
rz(-1.858254) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.090102) q[0];
sx q[0];
rz(-0.83054709) q[0];
sx q[0];
rz(-2.8175957) q[0];
rz(-1.8404768) q[1];
sx q[1];
rz(-0.83530656) q[1];
sx q[1];
rz(1.3791929) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8934879) q[0];
sx q[0];
rz(-1.4890492) q[0];
sx q[0];
rz(-1.3260613) q[0];
rz(-pi) q[1];
rz(-1.9786644) q[2];
sx q[2];
rz(-1.1424354) q[2];
sx q[2];
rz(0.63901627) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.56475793) q[1];
sx q[1];
rz(-1.4079637) q[1];
sx q[1];
rz(-2.8388883) q[1];
x q[2];
rz(-2.2961388) q[3];
sx q[3];
rz(-0.77973706) q[3];
sx q[3];
rz(0.23507915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.032701187) q[2];
sx q[2];
rz(-1.6872493) q[2];
sx q[2];
rz(0.94318715) q[2];
rz(0.33106783) q[3];
sx q[3];
rz(-1.7739242) q[3];
sx q[3];
rz(0.29512063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.023495) q[0];
sx q[0];
rz(-0.58693111) q[0];
sx q[0];
rz(2.3773637) q[0];
rz(-0.14097342) q[1];
sx q[1];
rz(-2.7455175) q[1];
sx q[1];
rz(1.2932628) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8338776) q[0];
sx q[0];
rz(-0.47979646) q[0];
sx q[0];
rz(2.0108372) q[0];
x q[1];
rz(-2.5404732) q[2];
sx q[2];
rz(-1.6496611) q[2];
sx q[2];
rz(2.2476946) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.075899374) q[1];
sx q[1];
rz(-1.432044) q[1];
sx q[1];
rz(-1.1636415) q[1];
rz(-pi) q[2];
rz(-2.7420298) q[3];
sx q[3];
rz(-1.385653) q[3];
sx q[3];
rz(-0.012133908) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.76688898) q[2];
sx q[2];
rz(-2.1024487) q[2];
sx q[2];
rz(-2.5602706) q[2];
rz(0.86822048) q[3];
sx q[3];
rz(-2.7243082) q[3];
sx q[3];
rz(1.8301615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
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
rz(-0.54810369) q[0];
sx q[0];
rz(-2.5618401) q[0];
sx q[0];
rz(1.2506437) q[0];
rz(2.1447694) q[1];
sx q[1];
rz(-2.2579028) q[1];
sx q[1];
rz(1.2876127) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9171137) q[0];
sx q[0];
rz(-2.0600442) q[0];
sx q[0];
rz(-3.1130303) q[0];
rz(-pi) q[1];
x q[1];
rz(0.31918819) q[2];
sx q[2];
rz(-0.50779283) q[2];
sx q[2];
rz(-2.9920981) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.3567644) q[1];
sx q[1];
rz(-2.462466) q[1];
sx q[1];
rz(-1.9041512) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5136396) q[3];
sx q[3];
rz(-0.74712979) q[3];
sx q[3];
rz(0.77880083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.8210956) q[2];
sx q[2];
rz(-0.38791502) q[2];
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
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(2.2607516) q[0];
sx q[0];
rz(-2.3374225) q[0];
sx q[0];
rz(2.8163731) q[0];
rz(1.1351769) q[1];
sx q[1];
rz(-0.67712855) q[1];
sx q[1];
rz(2.7744055) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8866667) q[0];
sx q[0];
rz(-1.6187795) q[0];
sx q[0];
rz(-1.6019437) q[0];
rz(-0.62717168) q[2];
sx q[2];
rz(-1.7998724) q[2];
sx q[2];
rz(2.0723745) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.6348833) q[1];
sx q[1];
rz(-2.2193529) q[1];
sx q[1];
rz(-1.8222005) q[1];
x q[2];
rz(-0.96079798) q[3];
sx q[3];
rz(-0.73878091) q[3];
sx q[3];
rz(-1.184954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8964768) q[2];
sx q[2];
rz(-2.3325236) q[2];
sx q[2];
rz(-1.0661351) q[2];
rz(3.0623479) q[3];
sx q[3];
rz(-2.5388122) q[3];
sx q[3];
rz(-1.665303) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
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
rz(0.29905839) q[1];
sx q[1];
rz(-2.538264) q[1];
sx q[1];
rz(-2.6187142) q[1];
rz(1.6499741) q[2];
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
