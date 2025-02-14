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
rz(2.4616315) q[0];
sx q[0];
rz(-2.2357219) q[0];
sx q[0];
rz(1.0933956) q[0];
rz(1.5054585) q[1];
sx q[1];
rz(-1.9688164) q[1];
sx q[1];
rz(-2.7170031) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.070804) q[0];
sx q[0];
rz(-1.6089526) q[0];
sx q[0];
rz(1.7724994) q[0];
rz(3.0908747) q[2];
sx q[2];
rz(-1.7246304) q[2];
sx q[2];
rz(2.5209629) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.72480768) q[1];
sx q[1];
rz(-2.8832286) q[1];
sx q[1];
rz(-1.8838253) q[1];
rz(-0.86058094) q[3];
sx q[3];
rz(-1.672555) q[3];
sx q[3];
rz(1.4836131) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.7666185) q[2];
sx q[2];
rz(-1.5593636) q[2];
sx q[2];
rz(1.1589104) q[2];
rz(2.9186987) q[3];
sx q[3];
rz(-1.736172) q[3];
sx q[3];
rz(0.29002732) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3463773) q[0];
sx q[0];
rz(-1.5257436) q[0];
sx q[0];
rz(1.079153) q[0];
rz(1.4214628) q[1];
sx q[1];
rz(-0.68669569) q[1];
sx q[1];
rz(3.0770643) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20708974) q[0];
sx q[0];
rz(-0.10182589) q[0];
sx q[0];
rz(-1.0290716) q[0];
rz(-pi) q[1];
rz(0.47625292) q[2];
sx q[2];
rz(-1.5956399) q[2];
sx q[2];
rz(0.78643878) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.77212438) q[1];
sx q[1];
rz(-1.8126948) q[1];
sx q[1];
rz(0.78472225) q[1];
rz(-0.82134473) q[3];
sx q[3];
rz(-1.708416) q[3];
sx q[3];
rz(2.2087165) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.1530389) q[2];
sx q[2];
rz(-1.3554074) q[2];
sx q[2];
rz(1.1174196) q[2];
rz(2.2828263) q[3];
sx q[3];
rz(-2.358181) q[3];
sx q[3];
rz(2.7045238) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3808463) q[0];
sx q[0];
rz(-0.28457156) q[0];
sx q[0];
rz(-0.17307702) q[0];
rz(-0.95678798) q[1];
sx q[1];
rz(-1.0280949) q[1];
sx q[1];
rz(2.079336) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0883368) q[0];
sx q[0];
rz(-2.2940141) q[0];
sx q[0];
rz(0.0084612554) q[0];
x q[1];
rz(0.64086242) q[2];
sx q[2];
rz(-1.0932361) q[2];
sx q[2];
rz(-1.8030082) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.77327496) q[1];
sx q[1];
rz(-0.89491208) q[1];
sx q[1];
rz(-0.85974093) q[1];
rz(-pi) q[2];
x q[2];
rz(0.62543243) q[3];
sx q[3];
rz(-2.1533826) q[3];
sx q[3];
rz(0.39804493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.411285) q[2];
sx q[2];
rz(-1.2016808) q[2];
sx q[2];
rz(-0.74472204) q[2];
rz(0.64106411) q[3];
sx q[3];
rz(-1.3499667) q[3];
sx q[3];
rz(2.6509269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4500126) q[0];
sx q[0];
rz(-0.55679655) q[0];
sx q[0];
rz(-0.84521729) q[0];
rz(-2.7382964) q[1];
sx q[1];
rz(-1.6019628) q[1];
sx q[1];
rz(0.32803112) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83932861) q[0];
sx q[0];
rz(-1.7516881) q[0];
sx q[0];
rz(2.9995287) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5045067) q[2];
sx q[2];
rz(-1.8083329) q[2];
sx q[2];
rz(1.4927166) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.5535132) q[1];
sx q[1];
rz(-0.57826406) q[1];
sx q[1];
rz(0.2167313) q[1];
rz(-2.8680236) q[3];
sx q[3];
rz(-2.0988587) q[3];
sx q[3];
rz(0.4116962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.8959117) q[2];
sx q[2];
rz(-0.74378219) q[2];
sx q[2];
rz(-0.15920676) q[2];
rz(0.016247449) q[3];
sx q[3];
rz(-1.0280321) q[3];
sx q[3];
rz(-2.4102559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3943587) q[0];
sx q[0];
rz(-1.9488229) q[0];
sx q[0];
rz(1.0900981) q[0];
rz(3.0116426) q[1];
sx q[1];
rz(-2.3268685) q[1];
sx q[1];
rz(-1.6158339) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65924683) q[0];
sx q[0];
rz(-1.0291983) q[0];
sx q[0];
rz(2.0603176) q[0];
x q[1];
rz(0.075399727) q[2];
sx q[2];
rz(-1.7345718) q[2];
sx q[2];
rz(-0.57367708) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.9435038) q[1];
sx q[1];
rz(-0.63423613) q[1];
sx q[1];
rz(3.0953118) q[1];
rz(-pi) q[2];
rz(-1.3899743) q[3];
sx q[3];
rz(-1.6055593) q[3];
sx q[3];
rz(-1.8145479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.7502363) q[2];
sx q[2];
rz(-2.3199234) q[2];
sx q[2];
rz(2.6643122) q[2];
rz(2.9376049) q[3];
sx q[3];
rz(-2.9450649) q[3];
sx q[3];
rz(0.6240713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9452962) q[0];
sx q[0];
rz(-0.81552234) q[0];
sx q[0];
rz(-0.77769172) q[0];
rz(2.5241191) q[1];
sx q[1];
rz(-1.6505046) q[1];
sx q[1];
rz(-2.2106574) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6361178) q[0];
sx q[0];
rz(-1.1606998) q[0];
sx q[0];
rz(-2.0652886) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.41895509) q[2];
sx q[2];
rz(-2.0134996) q[2];
sx q[2];
rz(-1.0087763) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.8691784) q[1];
sx q[1];
rz(-2.8003516) q[1];
sx q[1];
rz(-1.7029087) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.191675) q[3];
sx q[3];
rz(-1.2395879) q[3];
sx q[3];
rz(-3.1049867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4882539) q[2];
sx q[2];
rz(-1.808017) q[2];
sx q[2];
rz(-0.2674357) q[2];
rz(-2.2087162) q[3];
sx q[3];
rz(-1.2837912) q[3];
sx q[3];
rz(2.647184) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47144181) q[0];
sx q[0];
rz(-1.9981367) q[0];
sx q[0];
rz(-2.4757521) q[0];
rz(-0.2039856) q[1];
sx q[1];
rz(-0.98122707) q[1];
sx q[1];
rz(-1.1873672) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68722938) q[0];
sx q[0];
rz(-0.47894127) q[0];
sx q[0];
rz(2.9402551) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.97839956) q[2];
sx q[2];
rz(-0.80139388) q[2];
sx q[2];
rz(-1.3260384) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.59661181) q[1];
sx q[1];
rz(-2.2447733) q[1];
sx q[1];
rz(1.9121714) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3755619) q[3];
sx q[3];
rz(-0.25849202) q[3];
sx q[3];
rz(-3.0743161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.5953956) q[2];
sx q[2];
rz(-0.51837817) q[2];
sx q[2];
rz(-0.67016822) q[2];
rz(2.8403122) q[3];
sx q[3];
rz(-1.2944841) q[3];
sx q[3];
rz(2.7231351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30348521) q[0];
sx q[0];
rz(-2.1391588) q[0];
sx q[0];
rz(2.0759034) q[0];
rz(0.65713716) q[1];
sx q[1];
rz(-0.70825759) q[1];
sx q[1];
rz(1.4297952) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2524714) q[0];
sx q[0];
rz(-1.9074797) q[0];
sx q[0];
rz(-2.1648429) q[0];
rz(-pi) q[1];
rz(-0.33139511) q[2];
sx q[2];
rz(-1.7426531) q[2];
sx q[2];
rz(-0.21769014) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.2473879) q[1];
sx q[1];
rz(-2.6726843) q[1];
sx q[1];
rz(-0.50007485) q[1];
rz(-2.9518806) q[3];
sx q[3];
rz(-1.4290775) q[3];
sx q[3];
rz(2.313368) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.8810205) q[2];
sx q[2];
rz(-1.2994956) q[2];
sx q[2];
rz(1.0183498) q[2];
rz(1.5492505) q[3];
sx q[3];
rz(-1.6377623) q[3];
sx q[3];
rz(0.75756592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-2.8757979) q[0];
sx q[0];
rz(-0.51545155) q[0];
sx q[0];
rz(2.1739668) q[0];
rz(-0.34010092) q[1];
sx q[1];
rz(-2.3022771) q[1];
sx q[1];
rz(-1.0338773) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5115625) q[0];
sx q[0];
rz(-2.1242737) q[0];
sx q[0];
rz(0.16279499) q[0];
x q[1];
rz(-2.5754355) q[2];
sx q[2];
rz(-0.70933178) q[2];
sx q[2];
rz(-0.69401238) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.441923) q[1];
sx q[1];
rz(-1.960454) q[1];
sx q[1];
rz(0.81894919) q[1];
x q[2];
rz(2.9176209) q[3];
sx q[3];
rz(-1.4272318) q[3];
sx q[3];
rz(-0.75907133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.10244441) q[2];
sx q[2];
rz(-1.6951963) q[2];
sx q[2];
rz(-0.038912494) q[2];
rz(0.14032042) q[3];
sx q[3];
rz(-0.084429927) q[3];
sx q[3];
rz(-1.3154359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4843531) q[0];
sx q[0];
rz(-1.0056647) q[0];
sx q[0];
rz(-1.2580309) q[0];
rz(0.21615061) q[1];
sx q[1];
rz(-2.436147) q[1];
sx q[1];
rz(0.36453882) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.011019) q[0];
sx q[0];
rz(-2.6219061) q[0];
sx q[0];
rz(1.0566061) q[0];
x q[1];
rz(-0.44224085) q[2];
sx q[2];
rz(-0.85390515) q[2];
sx q[2];
rz(2.2539162) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.6072104) q[1];
sx q[1];
rz(-3.0845007) q[1];
sx q[1];
rz(-2.9109138) q[1];
rz(-pi) q[2];
rz(-1.5384501) q[3];
sx q[3];
rz(-1.5129287) q[3];
sx q[3];
rz(-0.75266121) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.2898499) q[2];
sx q[2];
rz(-1.936329) q[2];
sx q[2];
rz(0.64129788) q[2];
rz(1.6614206) q[3];
sx q[3];
rz(-1.8931754) q[3];
sx q[3];
rz(-0.67354584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
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
rz(0.73364532) q[0];
sx q[0];
rz(-0.49325627) q[0];
sx q[0];
rz(2.7706964) q[0];
rz(0.98450487) q[1];
sx q[1];
rz(-1.6592818) q[1];
sx q[1];
rz(2.0936113) q[1];
rz(-1.7683239) q[2];
sx q[2];
rz(-1.6850204) q[2];
sx q[2];
rz(-1.9584283) q[2];
rz(0.52624191) q[3];
sx q[3];
rz(-1.146823) q[3];
sx q[3];
rz(-1.5592195) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
