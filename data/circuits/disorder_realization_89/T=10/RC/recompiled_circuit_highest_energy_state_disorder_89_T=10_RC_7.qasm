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
rz(3.2193174) q[0];
sx q[0];
rz(5.2585703) q[0];
sx q[0];
rz(10.724714) q[0];
rz(-0.44252244) q[1];
sx q[1];
rz(-1.1392925) q[1];
sx q[1];
rz(-2.7055969) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0511507) q[0];
sx q[0];
rz(-1.465475) q[0];
sx q[0];
rz(1.9645343) q[0];
rz(-pi) q[1];
rz(0.73968621) q[2];
sx q[2];
rz(-1.7285182) q[2];
sx q[2];
rz(-0.24631234) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.2548542) q[1];
sx q[1];
rz(-1.5553932) q[1];
sx q[1];
rz(-0.40285691) q[1];
rz(-pi) q[2];
rz(-0.95753448) q[3];
sx q[3];
rz(-1.6744348) q[3];
sx q[3];
rz(0.14305015) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.062833) q[2];
sx q[2];
rz(-1.6597513) q[2];
sx q[2];
rz(0.01595846) q[2];
rz(1.6503085) q[3];
sx q[3];
rz(-1.242638) q[3];
sx q[3];
rz(-2.7119467) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4780739) q[0];
sx q[0];
rz(-2.0500545) q[0];
sx q[0];
rz(1.3462521) q[0];
rz(-1.1675872) q[1];
sx q[1];
rz(-2.0474032) q[1];
sx q[1];
rz(1.8928554) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7715147) q[0];
sx q[0];
rz(-1.9569719) q[0];
sx q[0];
rz(-2.2313928) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6925519) q[2];
sx q[2];
rz(-1.7261862) q[2];
sx q[2];
rz(2.4933585) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.84184835) q[1];
sx q[1];
rz(-1.9711291) q[1];
sx q[1];
rz(3.0818066) q[1];
rz(-pi) q[2];
rz(1.0661032) q[3];
sx q[3];
rz(-2.4923996) q[3];
sx q[3];
rz(2.0673942) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.0060129082) q[2];
sx q[2];
rz(-2.4485782) q[2];
sx q[2];
rz(2.6683624) q[2];
rz(3.0114975) q[3];
sx q[3];
rz(-1.3869163) q[3];
sx q[3];
rz(-0.010802833) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3267645) q[0];
sx q[0];
rz(-0.15693754) q[0];
sx q[0];
rz(-2.9275295) q[0];
rz(1.3940943) q[1];
sx q[1];
rz(-0.97624818) q[1];
sx q[1];
rz(-0.086437978) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6753767) q[0];
sx q[0];
rz(-1.6965742) q[0];
sx q[0];
rz(-1.2745538) q[0];
x q[1];
rz(-2.222074) q[2];
sx q[2];
rz(-0.79139564) q[2];
sx q[2];
rz(0.00080303116) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.4519388) q[1];
sx q[1];
rz(-1.1517593) q[1];
sx q[1];
rz(-1.3801912) q[1];
rz(0.92126927) q[3];
sx q[3];
rz(-2.0100397) q[3];
sx q[3];
rz(0.64542239) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.67132407) q[2];
sx q[2];
rz(-0.9321804) q[2];
sx q[2];
rz(1.1364802) q[2];
rz(1.350435) q[3];
sx q[3];
rz(-1.8108436) q[3];
sx q[3];
rz(2.7854846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7478624) q[0];
sx q[0];
rz(-0.50734729) q[0];
sx q[0];
rz(0.90079975) q[0];
rz(-1.2334709) q[1];
sx q[1];
rz(-0.8546468) q[1];
sx q[1];
rz(2.5305117) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.184) q[0];
sx q[0];
rz(-1.5779571) q[0];
sx q[0];
rz(1.0953127) q[0];
x q[1];
rz(0.22902352) q[2];
sx q[2];
rz(-1.5364858) q[2];
sx q[2];
rz(-1.2364309) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.722695) q[1];
sx q[1];
rz(-1.1686106) q[1];
sx q[1];
rz(-2.6821892) q[1];
rz(-1.6875721) q[3];
sx q[3];
rz(-1.2216785) q[3];
sx q[3];
rz(-1.4064096) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.5078807) q[2];
sx q[2];
rz(-0.66358006) q[2];
sx q[2];
rz(-3.0207685) q[2];
rz(-0.027033022) q[3];
sx q[3];
rz(-2.9398672) q[3];
sx q[3];
rz(-0.87593186) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19288572) q[0];
sx q[0];
rz(-0.8256194) q[0];
sx q[0];
rz(-1.8008308) q[0];
rz(1.6138529) q[1];
sx q[1];
rz(-1.8283045) q[1];
sx q[1];
rz(2.5921879) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8698688) q[0];
sx q[0];
rz(-2.3921674) q[0];
sx q[0];
rz(-3.0282227) q[0];
rz(-pi) q[1];
x q[1];
rz(3.1099186) q[2];
sx q[2];
rz(-2.1018545) q[2];
sx q[2];
rz(2.2344294) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.0548965) q[1];
sx q[1];
rz(-0.48829309) q[1];
sx q[1];
rz(1.6773305) q[1];
rz(-1.109455) q[3];
sx q[3];
rz(-1.6008988) q[3];
sx q[3];
rz(0.96772675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.01650979) q[2];
sx q[2];
rz(-1.5089401) q[2];
sx q[2];
rz(-2.5588918) q[2];
rz(-1.6216283) q[3];
sx q[3];
rz(-2.426332) q[3];
sx q[3];
rz(-1.3165855) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3358066) q[0];
sx q[0];
rz(-1.0842706) q[0];
sx q[0];
rz(1.8871319) q[0];
rz(-1.1955903) q[1];
sx q[1];
rz(-1.7343727) q[1];
sx q[1];
rz(1.5171299) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74745893) q[0];
sx q[0];
rz(-1.582495) q[0];
sx q[0];
rz(1.6375084) q[0];
rz(-pi) q[1];
rz(-0.56052358) q[2];
sx q[2];
rz(-0.88489489) q[2];
sx q[2];
rz(-2.1328515) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.9395917) q[1];
sx q[1];
rz(-2.4701405) q[1];
sx q[1];
rz(-0.47452773) q[1];
rz(-3.1213753) q[3];
sx q[3];
rz(-2.6264694) q[3];
sx q[3];
rz(2.8216854) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.967041) q[2];
sx q[2];
rz(-1.5341362) q[2];
sx q[2];
rz(-3.0957481) q[2];
rz(2.6144419) q[3];
sx q[3];
rz(-1.0322626) q[3];
sx q[3];
rz(2.3003858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6531649) q[0];
sx q[0];
rz(-2.4808352) q[0];
sx q[0];
rz(2.7556457) q[0];
rz(1.7653607) q[1];
sx q[1];
rz(-2.0538581) q[1];
sx q[1];
rz(-1.2765346) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4285884) q[0];
sx q[0];
rz(-2.3516293) q[0];
sx q[0];
rz(-2.6053564) q[0];
rz(0.047983147) q[2];
sx q[2];
rz(-1.1986102) q[2];
sx q[2];
rz(-1.3822777) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.3500665) q[1];
sx q[1];
rz(-3.1050395) q[1];
sx q[1];
rz(1.3514446) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5064729) q[3];
sx q[3];
rz(-1.5350047) q[3];
sx q[3];
rz(2.2944642) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4253) q[2];
sx q[2];
rz(-1.425681) q[2];
sx q[2];
rz(-1.7745793) q[2];
rz(-3.0573209) q[3];
sx q[3];
rz(-1.1140946) q[3];
sx q[3];
rz(-3.058694) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.080634236) q[0];
sx q[0];
rz(-3.0301889) q[0];
sx q[0];
rz(-1.8633307) q[0];
rz(0.06079611) q[1];
sx q[1];
rz(-2.1863329) q[1];
sx q[1];
rz(2.0416868) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5583145) q[0];
sx q[0];
rz(-0.41039) q[0];
sx q[0];
rz(-1.9463842) q[0];
rz(2.3448337) q[2];
sx q[2];
rz(-2.1646059) q[2];
sx q[2];
rz(-1.6006921) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.5775313) q[1];
sx q[1];
rz(-1.8558335) q[1];
sx q[1];
rz(-0.37071812) q[1];
rz(-pi) q[2];
rz(2.527929) q[3];
sx q[3];
rz(-1.1432768) q[3];
sx q[3];
rz(0.18117426) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.78201571) q[2];
sx q[2];
rz(-0.64543739) q[2];
sx q[2];
rz(0.47425708) q[2];
rz(-0.54840243) q[3];
sx q[3];
rz(-2.4689398) q[3];
sx q[3];
rz(1.3801581) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7062374) q[0];
sx q[0];
rz(-2.7147003) q[0];
sx q[0];
rz(-1.2109582) q[0];
rz(-1.8636761) q[1];
sx q[1];
rz(-1.5958818) q[1];
sx q[1];
rz(-0.011215297) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29182157) q[0];
sx q[0];
rz(-1.0106083) q[0];
sx q[0];
rz(1.4683506) q[0];
x q[1];
rz(0.36478784) q[2];
sx q[2];
rz(-1.786149) q[2];
sx q[2];
rz(-0.47206934) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-3.1234731) q[1];
sx q[1];
rz(-1.5168861) q[1];
sx q[1];
rz(2.3192295) q[1];
x q[2];
rz(-1.8770921) q[3];
sx q[3];
rz(-2.838475) q[3];
sx q[3];
rz(2.1949286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.32030216) q[2];
sx q[2];
rz(-1.236329) q[2];
sx q[2];
rz(-0.75766364) q[2];
rz(-0.6997987) q[3];
sx q[3];
rz(-2.564513) q[3];
sx q[3];
rz(0.9052161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59117544) q[0];
sx q[0];
rz(-1.2400405) q[0];
sx q[0];
rz(-2.8072939) q[0];
rz(-0.39407691) q[1];
sx q[1];
rz(-2.1912992) q[1];
sx q[1];
rz(-1.9091512) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52215965) q[0];
sx q[0];
rz(-1.30511) q[0];
sx q[0];
rz(-2.9152318) q[0];
x q[1];
rz(-2.903995) q[2];
sx q[2];
rz(-2.4294649) q[2];
sx q[2];
rz(-2.1472296) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.10127549) q[1];
sx q[1];
rz(-0.35666944) q[1];
sx q[1];
rz(0.89445313) q[1];
x q[2];
rz(-1.0717137) q[3];
sx q[3];
rz(-1.8129375) q[3];
sx q[3];
rz(-1.4361962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.95840994) q[2];
sx q[2];
rz(-0.56768688) q[2];
sx q[2];
rz(-3.0779823) q[2];
rz(2.375864) q[3];
sx q[3];
rz(-2.016341) q[3];
sx q[3];
rz(0.061802797) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3381989) q[0];
sx q[0];
rz(-0.82427187) q[0];
sx q[0];
rz(-1.6765539) q[0];
rz(1.4631396) q[1];
sx q[1];
rz(-1.2668162) q[1];
sx q[1];
rz(-0.75513671) q[1];
rz(-2.9207567) q[2];
sx q[2];
rz(-2.4855843) q[2];
sx q[2];
rz(1.9334855) q[2];
rz(-2.3018671) q[3];
sx q[3];
rz(-1.4473549) q[3];
sx q[3];
rz(2.2458129) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
