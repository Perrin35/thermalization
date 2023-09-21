OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.4322296) q[0];
sx q[0];
rz(-0.95786434) q[0];
sx q[0];
rz(0.14444484) q[0];
rz(0.56675178) q[1];
sx q[1];
rz(-0.52539879) q[1];
sx q[1];
rz(0.97775835) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5241961) q[0];
sx q[0];
rz(-2.8688736) q[0];
sx q[0];
rz(-2.8196536) q[0];
x q[1];
rz(-1.058504) q[2];
sx q[2];
rz(-0.57527486) q[2];
sx q[2];
rz(0.011205999) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.2797151) q[1];
sx q[1];
rz(-1.4289083) q[1];
sx q[1];
rz(-0.14076294) q[1];
x q[2];
rz(-2.3542777) q[3];
sx q[3];
rz(-1.4437321) q[3];
sx q[3];
rz(-0.78391677) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.67291659) q[2];
sx q[2];
rz(-1.9925995) q[2];
sx q[2];
rz(0.93227512) q[2];
rz(2.9428234) q[3];
sx q[3];
rz(-2.0310183) q[3];
sx q[3];
rz(2.1762302) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4345877) q[0];
sx q[0];
rz(-0.90536896) q[0];
sx q[0];
rz(-0.36112753) q[0];
rz(1.7065642) q[1];
sx q[1];
rz(-1.7838493) q[1];
sx q[1];
rz(2.3235869) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0571787) q[0];
sx q[0];
rz(-1.1857496) q[0];
sx q[0];
rz(1.2731228) q[0];
rz(-pi) q[1];
rz(1.8774162) q[2];
sx q[2];
rz(-2.7536256) q[2];
sx q[2];
rz(-1.0248794) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6815731) q[1];
sx q[1];
rz(-0.96833723) q[1];
sx q[1];
rz(-0.54089344) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7198256) q[3];
sx q[3];
rz(-0.36747284) q[3];
sx q[3];
rz(-0.056882337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.25257418) q[2];
sx q[2];
rz(-2.694464) q[2];
sx q[2];
rz(1.4206295) q[2];
rz(1.3160926) q[3];
sx q[3];
rz(-2.3836453) q[3];
sx q[3];
rz(0.38823286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7398359) q[0];
sx q[0];
rz(-0.49210423) q[0];
sx q[0];
rz(-2.2706568) q[0];
rz(2.8254106) q[1];
sx q[1];
rz(-0.28156391) q[1];
sx q[1];
rz(2.8443764) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0270099) q[0];
sx q[0];
rz(-1.5073338) q[0];
sx q[0];
rz(0.40513904) q[0];
x q[1];
rz(-0.88163968) q[2];
sx q[2];
rz(-2.214553) q[2];
sx q[2];
rz(1.2780485) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.19993648) q[1];
sx q[1];
rz(-1.7574649) q[1];
sx q[1];
rz(3.093064) q[1];
rz(2.3171595) q[3];
sx q[3];
rz(-2.094305) q[3];
sx q[3];
rz(-1.5200966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.7685984) q[2];
sx q[2];
rz(-0.82295376) q[2];
sx q[2];
rz(-2.7139943) q[2];
rz(1.9528495) q[3];
sx q[3];
rz(-0.62994981) q[3];
sx q[3];
rz(0.52743131) q[3];
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
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7052085) q[0];
sx q[0];
rz(-1.5642865) q[0];
sx q[0];
rz(2.3676681) q[0];
rz(0.71290839) q[1];
sx q[1];
rz(-2.1122825) q[1];
sx q[1];
rz(2.4598222) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1958774) q[0];
sx q[0];
rz(-1.3099075) q[0];
sx q[0];
rz(-2.9257141) q[0];
x q[1];
rz(2.9245124) q[2];
sx q[2];
rz(-2.1639369) q[2];
sx q[2];
rz(-0.3074239) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.75979739) q[1];
sx q[1];
rz(-1.1361546) q[1];
sx q[1];
rz(-2.1684907) q[1];
rz(-pi) q[2];
rz(2.3247129) q[3];
sx q[3];
rz(-2.6602392) q[3];
sx q[3];
rz(0.64827418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.1327847) q[2];
sx q[2];
rz(-0.71295732) q[2];
sx q[2];
rz(-2.183389) q[2];
rz(-1.0664252) q[3];
sx q[3];
rz(-1.2833779) q[3];
sx q[3];
rz(-0.3796033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.585007) q[0];
sx q[0];
rz(-1.3469232) q[0];
sx q[0];
rz(-0.56030309) q[0];
rz(0.99984461) q[1];
sx q[1];
rz(-2.9381349) q[1];
sx q[1];
rz(1.5195742) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.863127) q[0];
sx q[0];
rz(-2.0533877) q[0];
sx q[0];
rz(-2.0762073) q[0];
rz(-pi) q[1];
rz(-0.72603307) q[2];
sx q[2];
rz(-2.3880929) q[2];
sx q[2];
rz(1.7475278) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.49415627) q[1];
sx q[1];
rz(-1.372822) q[1];
sx q[1];
rz(-0.76311771) q[1];
x q[2];
rz(0.12445478) q[3];
sx q[3];
rz(-0.96154562) q[3];
sx q[3];
rz(-2.5973158) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.6891629) q[2];
sx q[2];
rz(-0.55183691) q[2];
sx q[2];
rz(0.2229283) q[2];
rz(-3.1068504) q[3];
sx q[3];
rz(-1.3870753) q[3];
sx q[3];
rz(0.071578659) q[3];
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
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3361622) q[0];
sx q[0];
rz(-0.31328377) q[0];
sx q[0];
rz(-1.0700595) q[0];
rz(1.7806212) q[1];
sx q[1];
rz(-2.762251) q[1];
sx q[1];
rz(-1.2840575) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5491935) q[0];
sx q[0];
rz(-2.1201029) q[0];
sx q[0];
rz(-2.4248289) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1806296) q[2];
sx q[2];
rz(-2.6764538) q[2];
sx q[2];
rz(-0.37365183) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.18404993) q[1];
sx q[1];
rz(-1.2831389) q[1];
sx q[1];
rz(-2.6581453) q[1];
rz(-1.0345801) q[3];
sx q[3];
rz(-1.3048733) q[3];
sx q[3];
rz(1.769161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.47508919) q[2];
sx q[2];
rz(-1.6978426) q[2];
sx q[2];
rz(0.63759032) q[2];
rz(2.3049138) q[3];
sx q[3];
rz(-1.1058608) q[3];
sx q[3];
rz(1.7975413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5752983) q[0];
sx q[0];
rz(-2.7116382) q[0];
sx q[0];
rz(-0.56754011) q[0];
rz(0.42770806) q[1];
sx q[1];
rz(-1.6141012) q[1];
sx q[1];
rz(-2.2033851) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.095031247) q[0];
sx q[0];
rz(-1.9782269) q[0];
sx q[0];
rz(-2.808232) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0023408) q[2];
sx q[2];
rz(-0.92407862) q[2];
sx q[2];
rz(1.2127753) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7792556) q[1];
sx q[1];
rz(-1.926683) q[1];
sx q[1];
rz(2.2625655) q[1];
rz(1.8754962) q[3];
sx q[3];
rz(-0.88677553) q[3];
sx q[3];
rz(-1.5414433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.2618835) q[2];
sx q[2];
rz(-1.1738913) q[2];
sx q[2];
rz(-1.7555457) q[2];
rz(1.322768) q[3];
sx q[3];
rz(-1.991792) q[3];
sx q[3];
rz(-0.02903207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26738527) q[0];
sx q[0];
rz(-2.8081942) q[0];
sx q[0];
rz(-1.7077131) q[0];
rz(1.2738312) q[1];
sx q[1];
rz(-2.0055983) q[1];
sx q[1];
rz(2.3103255) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6600663) q[0];
sx q[0];
rz(-0.43018451) q[0];
sx q[0];
rz(-0.35920401) q[0];
rz(-pi) q[1];
rz(2.2057461) q[2];
sx q[2];
rz(-1.4668462) q[2];
sx q[2];
rz(-2.7452591) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.9030994) q[1];
sx q[1];
rz(-2.0598754) q[1];
sx q[1];
rz(2.3563983) q[1];
rz(-1.3475111) q[3];
sx q[3];
rz(-1.3459473) q[3];
sx q[3];
rz(0.68883483) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.5635809) q[2];
sx q[2];
rz(-1.7687904) q[2];
sx q[2];
rz(0.9643628) q[2];
rz(-1.1635121) q[3];
sx q[3];
rz(-2.5301299) q[3];
sx q[3];
rz(0.65892974) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37373856) q[0];
sx q[0];
rz(-2.4066194) q[0];
sx q[0];
rz(-2.1642165) q[0];
rz(1.7550229) q[1];
sx q[1];
rz(-1.3061378) q[1];
sx q[1];
rz(-1.1057373) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39767299) q[0];
sx q[0];
rz(-1.1261252) q[0];
sx q[0];
rz(-1.7001274) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0366304) q[2];
sx q[2];
rz(-2.2347921) q[2];
sx q[2];
rz(-0.86519372) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.2178206) q[1];
sx q[1];
rz(-2.8433228) q[1];
sx q[1];
rz(2.24733) q[1];
rz(-pi) q[2];
rz(1.2826142) q[3];
sx q[3];
rz(-1.2328706) q[3];
sx q[3];
rz(1.6328788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.5332807) q[2];
sx q[2];
rz(-2.7567342) q[2];
sx q[2];
rz(-2.9917955) q[2];
rz(-1.3730565) q[3];
sx q[3];
rz(-1.7356197) q[3];
sx q[3];
rz(-1.3214553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49366429) q[0];
sx q[0];
rz(-2.6239008) q[0];
sx q[0];
rz(0.1272442) q[0];
rz(1.6607025) q[1];
sx q[1];
rz(-2.7791185) q[1];
sx q[1];
rz(-2.9737934) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6608749) q[0];
sx q[0];
rz(-1.5405476) q[0];
sx q[0];
rz(1.589993) q[0];
x q[1];
rz(2.9211505) q[2];
sx q[2];
rz(-2.5219005) q[2];
sx q[2];
rz(0.1534136) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.92229453) q[1];
sx q[1];
rz(-0.98235992) q[1];
sx q[1];
rz(-1.5625619) q[1];
rz(-pi) q[2];
rz(0.76370244) q[3];
sx q[3];
rz(-2.9716431) q[3];
sx q[3];
rz(-2.9156239) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(3.1148791) q[2];
sx q[2];
rz(-2.202704) q[2];
sx q[2];
rz(2.3804469) q[2];
rz(3.051565) q[3];
sx q[3];
rz(-1.0031676) q[3];
sx q[3];
rz(-2.1910523) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54820838) q[0];
sx q[0];
rz(-1.9822639) q[0];
sx q[0];
rz(-0.32250861) q[0];
rz(2.7535915) q[1];
sx q[1];
rz(-1.3996268) q[1];
sx q[1];
rz(-0.7849801) q[1];
rz(0.81007304) q[2];
sx q[2];
rz(-1.0675061) q[2];
sx q[2];
rz(1.6455417) q[2];
rz(1.0317867) q[3];
sx q[3];
rz(-1.4236593) q[3];
sx q[3];
rz(0.66766213) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];