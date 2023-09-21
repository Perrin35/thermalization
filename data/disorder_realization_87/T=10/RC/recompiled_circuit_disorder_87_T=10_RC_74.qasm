OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.23624578) q[0];
sx q[0];
rz(-2.4155004) q[0];
sx q[0];
rz(0.2015764) q[0];
rz(0.4959313) q[1];
sx q[1];
rz(-0.5402686) q[1];
sx q[1];
rz(2.2044866) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0721604) q[0];
sx q[0];
rz(-1.1496815) q[0];
sx q[0];
rz(-0.62299563) q[0];
rz(-pi) q[1];
rz(1.0785157) q[2];
sx q[2];
rz(-1.0168076) q[2];
sx q[2];
rz(2.2433777) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.5697437) q[1];
sx q[1];
rz(-1.7535926) q[1];
sx q[1];
rz(3.0621431) q[1];
rz(-pi) q[2];
rz(1.3017544) q[3];
sx q[3];
rz(-1.4770513) q[3];
sx q[3];
rz(0.71414381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.15930882) q[2];
sx q[2];
rz(-0.26580492) q[2];
sx q[2];
rz(0.086159555) q[2];
rz(-2.384095) q[3];
sx q[3];
rz(-2.3870654) q[3];
sx q[3];
rz(-1.0936201) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57698292) q[0];
sx q[0];
rz(-1.4819205) q[0];
sx q[0];
rz(-1.3264867) q[0];
rz(-1.2558698) q[1];
sx q[1];
rz(-1.5763667) q[1];
sx q[1];
rz(0.27145162) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2651228) q[0];
sx q[0];
rz(-2.070283) q[0];
sx q[0];
rz(2.6861517) q[0];
x q[1];
rz(-2.6255252) q[2];
sx q[2];
rz(-1.8415673) q[2];
sx q[2];
rz(-2.5004636) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.506362) q[1];
sx q[1];
rz(-2.2681384) q[1];
sx q[1];
rz(0.49763775) q[1];
x q[2];
rz(2.7852887) q[3];
sx q[3];
rz(-2.2998527) q[3];
sx q[3];
rz(2.9475398) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0945956) q[2];
sx q[2];
rz(-1.7589898) q[2];
sx q[2];
rz(0.57717741) q[2];
rz(-0.92352891) q[3];
sx q[3];
rz(-0.90463343) q[3];
sx q[3];
rz(1.7318168) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8975824) q[0];
sx q[0];
rz(-2.0799461) q[0];
sx q[0];
rz(0.14257167) q[0];
rz(1.3525195) q[1];
sx q[1];
rz(-2.1058154) q[1];
sx q[1];
rz(-2.9325063) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3554879) q[0];
sx q[0];
rz(-2.3007563) q[0];
sx q[0];
rz(-2.3441322) q[0];
x q[1];
rz(-0.44720165) q[2];
sx q[2];
rz(-0.7913835) q[2];
sx q[2];
rz(-2.7176822) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.28029728) q[1];
sx q[1];
rz(-1.1251083) q[1];
sx q[1];
rz(-0.46511005) q[1];
rz(1.3554058) q[3];
sx q[3];
rz(-1.8672018) q[3];
sx q[3];
rz(-2.4733558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.7769988) q[2];
sx q[2];
rz(-2.2475188) q[2];
sx q[2];
rz(-0.40412942) q[2];
rz(-1.8557619) q[3];
sx q[3];
rz(-1.1288246) q[3];
sx q[3];
rz(-2.9742441) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7053112) q[0];
sx q[0];
rz(-3.0349338) q[0];
sx q[0];
rz(-0.96281111) q[0];
rz(2.6722233) q[1];
sx q[1];
rz(-0.58987394) q[1];
sx q[1];
rz(-3.1406291) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9199333) q[0];
sx q[0];
rz(-2.1533305) q[0];
sx q[0];
rz(-0.68908738) q[0];
rz(-pi) q[1];
rz(1.5713111) q[2];
sx q[2];
rz(-1.6948023) q[2];
sx q[2];
rz(2.4027783) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2767267) q[1];
sx q[1];
rz(-2.556567) q[1];
sx q[1];
rz(1.3328972) q[1];
rz(-0.13929316) q[3];
sx q[3];
rz(-2.5677498) q[3];
sx q[3];
rz(-1.7593256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.0756388) q[2];
sx q[2];
rz(-1.6299738) q[2];
sx q[2];
rz(0.11165079) q[2];
rz(2.3305437) q[3];
sx q[3];
rz(-2.6871197) q[3];
sx q[3];
rz(0.013899175) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.951293) q[0];
sx q[0];
rz(-0.90068978) q[0];
sx q[0];
rz(-2.7850889) q[0];
rz(0.50645343) q[1];
sx q[1];
rz(-1.3566147) q[1];
sx q[1];
rz(-2.8809663) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5379932) q[0];
sx q[0];
rz(-1.4289745) q[0];
sx q[0];
rz(1.5285368) q[0];
x q[1];
rz(0.25116253) q[2];
sx q[2];
rz(-1.3248982) q[2];
sx q[2];
rz(0.22566251) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(3.0780328) q[1];
sx q[1];
rz(-1.5884001) q[1];
sx q[1];
rz(-0.5743919) q[1];
rz(-pi) q[2];
rz(-0.77633206) q[3];
sx q[3];
rz(-0.20271248) q[3];
sx q[3];
rz(2.4622038) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.9115209) q[2];
sx q[2];
rz(-1.6610961) q[2];
sx q[2];
rz(-0.22932209) q[2];
rz(2.5991332) q[3];
sx q[3];
rz(-2.8312603) q[3];
sx q[3];
rz(2.1054161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77263537) q[0];
sx q[0];
rz(-2.2326523) q[0];
sx q[0];
rz(1.4916346) q[0];
rz(-2.1024599) q[1];
sx q[1];
rz(-1.2960478) q[1];
sx q[1];
rz(-1.414149) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4947858) q[0];
sx q[0];
rz(-1.8456869) q[0];
sx q[0];
rz(-1.6199714) q[0];
x q[1];
rz(-0.51775198) q[2];
sx q[2];
rz(-1.8038097) q[2];
sx q[2];
rz(-1.053996) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.8687739) q[1];
sx q[1];
rz(-1.0827853) q[1];
sx q[1];
rz(0.52656071) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.066091) q[3];
sx q[3];
rz(-1.7864979) q[3];
sx q[3];
rz(-0.84738934) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.002939) q[2];
sx q[2];
rz(-2.5230375) q[2];
sx q[2];
rz(3.1138528) q[2];
rz(-2.6489143) q[3];
sx q[3];
rz(-1.490482) q[3];
sx q[3];
rz(-1.7475351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3843) q[0];
sx q[0];
rz(-1.5889656) q[0];
sx q[0];
rz(-3.0199155) q[0];
rz(1.1514459) q[1];
sx q[1];
rz(-0.45184389) q[1];
sx q[1];
rz(2.7391403) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50169045) q[0];
sx q[0];
rz(-2.2417289) q[0];
sx q[0];
rz(-2.7946266) q[0];
rz(1.5791113) q[2];
sx q[2];
rz(-1.6168211) q[2];
sx q[2];
rz(0.41551057) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.788159) q[1];
sx q[1];
rz(-0.33553365) q[1];
sx q[1];
rz(-2.5787756) q[1];
x q[2];
rz(-1.4943069) q[3];
sx q[3];
rz(-2.7923931) q[3];
sx q[3];
rz(2.4793712) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.014331269) q[2];
sx q[2];
rz(-1.971259) q[2];
sx q[2];
rz(-2.2793615) q[2];
rz(2.6640653) q[3];
sx q[3];
rz(-1.9600441) q[3];
sx q[3];
rz(2.2235218) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7858793) q[0];
sx q[0];
rz(-2.9861351) q[0];
sx q[0];
rz(-0.73295897) q[0];
rz(-3.0015302) q[1];
sx q[1];
rz(-2.1439794) q[1];
sx q[1];
rz(2.1070811) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2037172) q[0];
sx q[0];
rz(-2.8563742) q[0];
sx q[0];
rz(0.9271778) q[0];
rz(0.10266281) q[2];
sx q[2];
rz(-1.5005939) q[2];
sx q[2];
rz(1.7662802) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.95755277) q[1];
sx q[1];
rz(-1.5337481) q[1];
sx q[1];
rz(1.0047822) q[1];
rz(-1.6789867) q[3];
sx q[3];
rz(-1.1547609) q[3];
sx q[3];
rz(-3.0707404) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.90116477) q[2];
sx q[2];
rz(-1.1952091) q[2];
sx q[2];
rz(-0.77587664) q[2];
rz(-0.72426978) q[3];
sx q[3];
rz(-2.3359719) q[3];
sx q[3];
rz(1.7075214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4416606) q[0];
sx q[0];
rz(-2.3775546) q[0];
sx q[0];
rz(-0.016816703) q[0];
rz(3.1230714) q[1];
sx q[1];
rz(-1.8073945) q[1];
sx q[1];
rz(0.7787849) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9176365) q[0];
sx q[0];
rz(-0.42376873) q[0];
sx q[0];
rz(0.64026041) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0987894) q[2];
sx q[2];
rz(-1.325377) q[2];
sx q[2];
rz(2.5185042) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.89214954) q[1];
sx q[1];
rz(-1.8037233) q[1];
sx q[1];
rz(1.1942785) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3498464) q[3];
sx q[3];
rz(-1.6337992) q[3];
sx q[3];
rz(2.8642879) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.6918216) q[2];
sx q[2];
rz(-1.8372953) q[2];
sx q[2];
rz(0.71643913) q[2];
rz(1.4572432) q[3];
sx q[3];
rz(-1.4102035) q[3];
sx q[3];
rz(1.8523857) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0338106) q[0];
sx q[0];
rz(-0.53492117) q[0];
sx q[0];
rz(2.9901436) q[0];
rz(2.9653213) q[1];
sx q[1];
rz(-1.9468032) q[1];
sx q[1];
rz(-2.418628) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82397126) q[0];
sx q[0];
rz(-2.9569607) q[0];
sx q[0];
rz(-2.1872107) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1638219) q[2];
sx q[2];
rz(-0.40728912) q[2];
sx q[2];
rz(-2.2182857) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.1180229) q[1];
sx q[1];
rz(-1.5155063) q[1];
sx q[1];
rz(-1.0667332) q[1];
rz(-pi) q[2];
rz(-1.1687752) q[3];
sx q[3];
rz(-2.3178181) q[3];
sx q[3];
rz(1.0812425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.4621949) q[2];
sx q[2];
rz(-1.2827736) q[2];
sx q[2];
rz(0.52552137) q[2];
rz(0.28371352) q[3];
sx q[3];
rz(-2.0339537) q[3];
sx q[3];
rz(2.7450558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59239607) q[0];
sx q[0];
rz(-2.017673) q[0];
sx q[0];
rz(2.8155433) q[0];
rz(2.0422968) q[1];
sx q[1];
rz(-2.9922843) q[1];
sx q[1];
rz(-0.86984787) q[1];
rz(-0.90080558) q[2];
sx q[2];
rz(-1.1849891) q[2];
sx q[2];
rz(1.6580788) q[2];
rz(2.5476534) q[3];
sx q[3];
rz(-0.38729061) q[3];
sx q[3];
rz(-0.95872986) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];