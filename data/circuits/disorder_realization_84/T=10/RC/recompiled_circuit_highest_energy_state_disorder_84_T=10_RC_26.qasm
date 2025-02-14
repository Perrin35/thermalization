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
rz(-2.6211205) q[0];
sx q[0];
rz(-1.0360798) q[0];
sx q[0];
rz(-2.5007024) q[0];
rz(2.9042397) q[1];
sx q[1];
rz(-0.64164716) q[1];
sx q[1];
rz(0.084903804) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9176506) q[0];
sx q[0];
rz(-0.71086796) q[0];
sx q[0];
rz(-1.4213597) q[0];
rz(3.0831111) q[2];
sx q[2];
rz(-1.6019434) q[2];
sx q[2];
rz(0.29171388) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.1168659) q[1];
sx q[1];
rz(-1.4566531) q[1];
sx q[1];
rz(2.8495339) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0155979) q[3];
sx q[3];
rz(-2.4732504) q[3];
sx q[3];
rz(-2.4618142) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.25195965) q[2];
sx q[2];
rz(-1.2026938) q[2];
sx q[2];
rz(1.5084051) q[2];
rz(2.3097322) q[3];
sx q[3];
rz(-2.8155477) q[3];
sx q[3];
rz(1.7970201) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93083301) q[0];
sx q[0];
rz(-0.63527125) q[0];
sx q[0];
rz(0.26500901) q[0];
rz(-0.79568544) q[1];
sx q[1];
rz(-2.7964451) q[1];
sx q[1];
rz(1.4211242) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4583295) q[0];
sx q[0];
rz(-0.70752549) q[0];
sx q[0];
rz(-1.9842771) q[0];
rz(-0.95395835) q[2];
sx q[2];
rz(-1.1797311) q[2];
sx q[2];
rz(2.3352438) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.4231734) q[1];
sx q[1];
rz(-1.3191901) q[1];
sx q[1];
rz(0.045028506) q[1];
x q[2];
rz(-0.029985241) q[3];
sx q[3];
rz(-1.7383778) q[3];
sx q[3];
rz(-0.94887892) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.77218324) q[2];
sx q[2];
rz(-0.99516827) q[2];
sx q[2];
rz(-2.4824202) q[2];
rz(-2.2199421) q[3];
sx q[3];
rz(-2.0989336) q[3];
sx q[3];
rz(-1.5825745) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42940656) q[0];
sx q[0];
rz(-2.5262008) q[0];
sx q[0];
rz(-1.2834826) q[0];
rz(0.42690024) q[1];
sx q[1];
rz(-0.59395298) q[1];
sx q[1];
rz(1.1675507) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3385462) q[0];
sx q[0];
rz(-1.9340705) q[0];
sx q[0];
rz(0.35420322) q[0];
rz(-pi) q[1];
x q[1];
rz(0.53041544) q[2];
sx q[2];
rz(-1.1258719) q[2];
sx q[2];
rz(-2.4702768) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.1192273) q[1];
sx q[1];
rz(-1.4767547) q[1];
sx q[1];
rz(2.8507807) q[1];
x q[2];
rz(-2.5646788) q[3];
sx q[3];
rz(-0.31539279) q[3];
sx q[3];
rz(-1.9384991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.20649642) q[2];
sx q[2];
rz(-1.6556211) q[2];
sx q[2];
rz(3.0237954) q[2];
rz(-1.3055034) q[3];
sx q[3];
rz(-0.85956231) q[3];
sx q[3];
rz(-0.99087244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-0.30462581) q[0];
sx q[0];
rz(-1.0517629) q[0];
sx q[0];
rz(0.032531746) q[0];
rz(-2.1578728) q[1];
sx q[1];
rz(-2.3301221) q[1];
sx q[1];
rz(0.54289114) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5158096) q[0];
sx q[0];
rz(-0.90321945) q[0];
sx q[0];
rz(1.9208917) q[0];
rz(-pi) q[1];
rz(-0.10323287) q[2];
sx q[2];
rz(-1.2599753) q[2];
sx q[2];
rz(0.23676591) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.11527744) q[1];
sx q[1];
rz(-1.3168854) q[1];
sx q[1];
rz(2.3350969) q[1];
rz(2.9916688) q[3];
sx q[3];
rz(-2.3036851) q[3];
sx q[3];
rz(2.2265707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.5618318) q[2];
sx q[2];
rz(-1.1299955) q[2];
sx q[2];
rz(-0.96060166) q[2];
rz(-0.97551712) q[3];
sx q[3];
rz(-0.87149182) q[3];
sx q[3];
rz(0.52198854) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8303216) q[0];
sx q[0];
rz(-2.4595342) q[0];
sx q[0];
rz(2.8090546) q[0];
rz(-0.63811103) q[1];
sx q[1];
rz(-2.5159914) q[1];
sx q[1];
rz(2.6854551) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77959427) q[0];
sx q[0];
rz(-1.5080305) q[0];
sx q[0];
rz(-1.4383352) q[0];
rz(2.664723) q[2];
sx q[2];
rz(-0.84706351) q[2];
sx q[2];
rz(2.4349231) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3553473) q[1];
sx q[1];
rz(-1.659871) q[1];
sx q[1];
rz(2.1554355) q[1];
rz(-pi) q[2];
rz(0.36231492) q[3];
sx q[3];
rz(-2.1862595) q[3];
sx q[3];
rz(-0.42847732) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.0629603) q[2];
sx q[2];
rz(-0.76138622) q[2];
sx q[2];
rz(-2.6893993) q[2];
rz(-2.8771628) q[3];
sx q[3];
rz(-2.9954438) q[3];
sx q[3];
rz(-1.2368081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
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
rz(2.0753157) q[0];
sx q[0];
rz(-0.83106581) q[0];
sx q[0];
rz(0.7178632) q[0];
rz(-1.5526937) q[1];
sx q[1];
rz(-1.540686) q[1];
sx q[1];
rz(-2.9343361) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7826876) q[0];
sx q[0];
rz(-2.0126632) q[0];
sx q[0];
rz(1.5405517) q[0];
rz(-pi) q[1];
rz(1.147993) q[2];
sx q[2];
rz(-1.5059274) q[2];
sx q[2];
rz(-2.1366304) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.26970999) q[1];
sx q[1];
rz(-1.9970532) q[1];
sx q[1];
rz(1.0428509) q[1];
rz(-pi) q[2];
rz(3.0318694) q[3];
sx q[3];
rz(-1.8302813) q[3];
sx q[3];
rz(-0.55486996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.3805716) q[2];
sx q[2];
rz(-2.0863159) q[2];
sx q[2];
rz(2.5804248) q[2];
rz(-2.0821954) q[3];
sx q[3];
rz(-1.1166409) q[3];
sx q[3];
rz(1.4580956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92784944) q[0];
sx q[0];
rz(-0.97074592) q[0];
sx q[0];
rz(-2.2947776) q[0];
rz(0.98833409) q[1];
sx q[1];
rz(-1.7653932) q[1];
sx q[1];
rz(-2.3086937) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4109548) q[0];
sx q[0];
rz(-0.22408501) q[0];
sx q[0];
rz(-0.48416324) q[0];
x q[1];
rz(-2.6360699) q[2];
sx q[2];
rz(-0.16211432) q[2];
sx q[2];
rz(-0.41451987) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0497694) q[1];
sx q[1];
rz(-0.97477978) q[1];
sx q[1];
rz(2.6906162) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0924721) q[3];
sx q[3];
rz(-1.6843154) q[3];
sx q[3];
rz(0.037179557) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.46574584) q[2];
sx q[2];
rz(-0.70285672) q[2];
sx q[2];
rz(0.77009002) q[2];
rz(-0.13104023) q[3];
sx q[3];
rz(-1.659487) q[3];
sx q[3];
rz(-1.2469163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7939664) q[0];
sx q[0];
rz(-2.2063231) q[0];
sx q[0];
rz(-0.4916218) q[0];
rz(-0.31279102) q[1];
sx q[1];
rz(-0.28952315) q[1];
sx q[1];
rz(-3.0591931) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0534957) q[0];
sx q[0];
rz(-1.4812569) q[0];
sx q[0];
rz(-3.128818) q[0];
rz(-pi) q[1];
rz(-0.92026842) q[2];
sx q[2];
rz(-0.81020497) q[2];
sx q[2];
rz(-1.6595464) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.5746153) q[1];
sx q[1];
rz(-1.7321371) q[1];
sx q[1];
rz(-1.3938422) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2639323) q[3];
sx q[3];
rz(-0.32684775) q[3];
sx q[3];
rz(0.89034789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.4440492) q[2];
sx q[2];
rz(-1.9335582) q[2];
sx q[2];
rz(3.1235798) q[2];
rz(0.87120122) q[3];
sx q[3];
rz(-0.33659354) q[3];
sx q[3];
rz(0.66756162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1910601) q[0];
sx q[0];
rz(-0.55452269) q[0];
sx q[0];
rz(-2.6771123) q[0];
rz(2.65061) q[1];
sx q[1];
rz(-1.4717439) q[1];
sx q[1];
rz(-1.6741265) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53297645) q[0];
sx q[0];
rz(-2.5166844) q[0];
sx q[0];
rz(3.141538) q[0];
x q[1];
rz(1.8764795) q[2];
sx q[2];
rz(-1.0846429) q[2];
sx q[2];
rz(0.4403688) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.3882093) q[1];
sx q[1];
rz(-0.25553726) q[1];
sx q[1];
rz(-1.8123883) q[1];
x q[2];
rz(-0.76045658) q[3];
sx q[3];
rz(-1.0257313) q[3];
sx q[3];
rz(2.2732002) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(3.1269425) q[2];
sx q[2];
rz(-1.4913968) q[2];
sx q[2];
rz(-1.3313741) q[2];
rz(-0.93975449) q[3];
sx q[3];
rz(-2.0846114) q[3];
sx q[3];
rz(-1.1872928) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2347655) q[0];
sx q[0];
rz(-0.48521388) q[0];
sx q[0];
rz(-2.1038798) q[0];
rz(2.0872133) q[1];
sx q[1];
rz(-1.8260006) q[1];
sx q[1];
rz(-2.0920928) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.755167) q[0];
sx q[0];
rz(-1.5706129) q[0];
sx q[0];
rz(-1.6205233) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0825572) q[2];
sx q[2];
rz(-2.0768713) q[2];
sx q[2];
rz(2.807694) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.29925525) q[1];
sx q[1];
rz(-1.5399057) q[1];
sx q[1];
rz(0.96986356) q[1];
rz(-pi) q[2];
x q[2];
rz(0.58934642) q[3];
sx q[3];
rz(-2.0799619) q[3];
sx q[3];
rz(1.2065534) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.99027571) q[2];
sx q[2];
rz(-1.8345202) q[2];
sx q[2];
rz(-1.6952197) q[2];
rz(2.9826048) q[3];
sx q[3];
rz(-1.1587326) q[3];
sx q[3];
rz(-0.75001636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0010407) q[0];
sx q[0];
rz(-0.97698553) q[0];
sx q[0];
rz(-2.3401596) q[0];
rz(1.4871545) q[1];
sx q[1];
rz(-2.9947037) q[1];
sx q[1];
rz(-0.53302232) q[1];
rz(-2.2764789) q[2];
sx q[2];
rz(-0.8556753) q[2];
sx q[2];
rz(-0.52134261) q[2];
rz(-2.4132397) q[3];
sx q[3];
rz(-0.77941685) q[3];
sx q[3];
rz(0.70913915) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
