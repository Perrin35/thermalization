OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.52842927) q[0];
sx q[0];
rz(2.0818721) q[0];
sx q[0];
rz(11.835397) q[0];
rz(1.641474) q[1];
sx q[1];
rz(5.2483622) q[1];
sx q[1];
rz(11.622826) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4306513) q[0];
sx q[0];
rz(-1.5358155) q[0];
sx q[0];
rz(3.0745688) q[0];
x q[1];
rz(-0.11124723) q[2];
sx q[2];
rz(-2.4876378) q[2];
sx q[2];
rz(2.7788869) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2343443) q[1];
sx q[1];
rz(-2.611479) q[1];
sx q[1];
rz(2.2905473) q[1];
rz(-1.7581851) q[3];
sx q[3];
rz(-1.2710147) q[3];
sx q[3];
rz(-0.0084458394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8865108) q[2];
sx q[2];
rz(-1.7604897) q[2];
sx q[2];
rz(-1.250766) q[2];
rz(-1.7154153) q[3];
sx q[3];
rz(-2.2255247) q[3];
sx q[3];
rz(-0.9799408) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.132906) q[0];
sx q[0];
rz(-2.0693021) q[0];
sx q[0];
rz(0.59894484) q[0];
rz(-1.8006181) q[1];
sx q[1];
rz(-0.95021617) q[1];
sx q[1];
rz(-0.96639955) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4616529) q[0];
sx q[0];
rz(-1.1302117) q[0];
sx q[0];
rz(-0.64042129) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0722926) q[2];
sx q[2];
rz(-0.6233223) q[2];
sx q[2];
rz(2.9181883) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.6858789) q[1];
sx q[1];
rz(-3.0027632) q[1];
sx q[1];
rz(0.60771897) q[1];
rz(2.6838449) q[3];
sx q[3];
rz(-2.3380937) q[3];
sx q[3];
rz(0.97704923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0559343) q[2];
sx q[2];
rz(-2.3264383) q[2];
sx q[2];
rz(-0.30109626) q[2];
rz(1.1931233) q[3];
sx q[3];
rz(-1.5914702) q[3];
sx q[3];
rz(1.1863856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
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
rz(0.10467228) q[0];
sx q[0];
rz(-1.5043229) q[0];
sx q[0];
rz(1.1874636) q[0];
rz(1.2359515) q[1];
sx q[1];
rz(-1.0373479) q[1];
sx q[1];
rz(1.3175861) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7147489) q[0];
sx q[0];
rz(-2.2509529) q[0];
sx q[0];
rz(1.8450518) q[0];
x q[1];
rz(2.3972829) q[2];
sx q[2];
rz(-2.8544606) q[2];
sx q[2];
rz(-0.14936514) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.6685851) q[1];
sx q[1];
rz(-0.43273941) q[1];
sx q[1];
rz(0.24344484) q[1];
rz(-pi) q[2];
rz(-0.1173238) q[3];
sx q[3];
rz(-2.085272) q[3];
sx q[3];
rz(-2.8770212) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.0111771) q[2];
sx q[2];
rz(-1.3935564) q[2];
sx q[2];
rz(2.9349566) q[2];
rz(2.4335499) q[3];
sx q[3];
rz(-2.9338624) q[3];
sx q[3];
rz(-2.2487683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3669423) q[0];
sx q[0];
rz(-2.9270524) q[0];
sx q[0];
rz(-2.2553717) q[0];
rz(1.0097424) q[1];
sx q[1];
rz(-2.2354398) q[1];
sx q[1];
rz(1.9151691) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0013989) q[0];
sx q[0];
rz(-1.6370019) q[0];
sx q[0];
rz(2.2982236) q[0];
x q[1];
rz(1.0272155) q[2];
sx q[2];
rz(-2.377844) q[2];
sx q[2];
rz(-0.56994146) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.36818477) q[1];
sx q[1];
rz(-0.28973026) q[1];
sx q[1];
rz(-2.3875931) q[1];
x q[2];
rz(-2.6552116) q[3];
sx q[3];
rz(-1.2380621) q[3];
sx q[3];
rz(-2.2330724) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.4975138) q[2];
sx q[2];
rz(-1.3138694) q[2];
sx q[2];
rz(-2.148596) q[2];
rz(1.8289061) q[3];
sx q[3];
rz(-2.1195181) q[3];
sx q[3];
rz(2.657857) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
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
rz(3.118367) q[0];
sx q[0];
rz(-2.826773) q[0];
sx q[0];
rz(-1.9556048) q[0];
rz(-1.7182619) q[1];
sx q[1];
rz(-1.6512197) q[1];
sx q[1];
rz(-0.58247724) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7459481) q[0];
sx q[0];
rz(-2.1714604) q[0];
sx q[0];
rz(-3.0546741) q[0];
rz(0.14898665) q[2];
sx q[2];
rz(-2.2431231) q[2];
sx q[2];
rz(2.7186) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.7354048) q[1];
sx q[1];
rz(-1.5516073) q[1];
sx q[1];
rz(-2.4270942) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.1547847) q[3];
sx q[3];
rz(-0.41066658) q[3];
sx q[3];
rz(-1.1759024) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.4867268) q[2];
sx q[2];
rz(-1.606769) q[2];
sx q[2];
rz(3.0598818) q[2];
rz(-2.667526) q[3];
sx q[3];
rz(-1.877955) q[3];
sx q[3];
rz(-1.6430395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.896647) q[0];
sx q[0];
rz(-1.1279673) q[0];
sx q[0];
rz(3.1337877) q[0];
rz(-1.4004978) q[1];
sx q[1];
rz(-2.2875319) q[1];
sx q[1];
rz(-1.1046462) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35690755) q[0];
sx q[0];
rz(-2.4570358) q[0];
sx q[0];
rz(0.25951578) q[0];
x q[1];
rz(-0.64013021) q[2];
sx q[2];
rz(-1.374561) q[2];
sx q[2];
rz(3.047608) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.1659516) q[1];
sx q[1];
rz(-2.2659321) q[1];
sx q[1];
rz(-2.402311) q[1];
x q[2];
rz(1.2475345) q[3];
sx q[3];
rz(-2.6544016) q[3];
sx q[3];
rz(-0.67684735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.2509987) q[2];
sx q[2];
rz(-2.4119174) q[2];
sx q[2];
rz(-0.55523038) q[2];
rz(-2.9688719) q[3];
sx q[3];
rz(-1.828086) q[3];
sx q[3];
rz(1.593332) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-1.5884488) q[0];
sx q[0];
rz(-0.48184904) q[0];
sx q[0];
rz(3.0798262) q[0];
rz(-2.899509) q[1];
sx q[1];
rz(-2.7663019) q[1];
sx q[1];
rz(2.0297091) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2568946) q[0];
sx q[0];
rz(-1.1567133) q[0];
sx q[0];
rz(-2.3143682) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7467473) q[2];
sx q[2];
rz(-0.07882747) q[2];
sx q[2];
rz(1.6579171) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5824077) q[1];
sx q[1];
rz(-1.294194) q[1];
sx q[1];
rz(-0.72699593) q[1];
x q[2];
rz(-0.63515969) q[3];
sx q[3];
rz(-0.40024647) q[3];
sx q[3];
rz(-1.1349585) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.3322488) q[2];
sx q[2];
rz(-2.377254) q[2];
sx q[2];
rz(-0.81364441) q[2];
rz(1.404473) q[3];
sx q[3];
rz(-0.22870326) q[3];
sx q[3];
rz(0.61541921) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5161045) q[0];
sx q[0];
rz(-1.7913211) q[0];
sx q[0];
rz(1.425449) q[0];
rz(-1.5215993) q[1];
sx q[1];
rz(-2.3896673) q[1];
sx q[1];
rz(2.5040748) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91598788) q[0];
sx q[0];
rz(-2.2921352) q[0];
sx q[0];
rz(-2.3774873) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2136739) q[2];
sx q[2];
rz(-2.2631858) q[2];
sx q[2];
rz(-1.1536319) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.058640826) q[1];
sx q[1];
rz(-0.73459) q[1];
sx q[1];
rz(-2.525108) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4333126) q[3];
sx q[3];
rz(-1.004389) q[3];
sx q[3];
rz(2.8979104) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.1311538) q[2];
sx q[2];
rz(-2.4413979) q[2];
sx q[2];
rz(1.1361702) q[2];
rz(1.4853959) q[3];
sx q[3];
rz(-2.6172726) q[3];
sx q[3];
rz(-1.2095399) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-1.6256325) q[0];
sx q[0];
rz(-2.3598598) q[0];
sx q[0];
rz(-0.67614722) q[0];
rz(0.82538429) q[1];
sx q[1];
rz(-0.4239347) q[1];
sx q[1];
rz(1.0151781) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0494173) q[0];
sx q[0];
rz(-1.1560688) q[0];
sx q[0];
rz(1.7134922) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9416111) q[2];
sx q[2];
rz(-0.66255424) q[2];
sx q[2];
rz(2.8560864) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.3489192) q[1];
sx q[1];
rz(-1.9209849) q[1];
sx q[1];
rz(0.37383553) q[1];
x q[2];
rz(1.2328524) q[3];
sx q[3];
rz(-0.83331185) q[3];
sx q[3];
rz(-2.1144707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.72145808) q[2];
sx q[2];
rz(-0.2162424) q[2];
sx q[2];
rz(2.1742415) q[2];
rz(-1.597065) q[3];
sx q[3];
rz(-1.8287851) q[3];
sx q[3];
rz(2.7887204) q[3];
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
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5230781) q[0];
sx q[0];
rz(-1.5788364) q[0];
sx q[0];
rz(-3.0850947) q[0];
rz(2.0691195) q[1];
sx q[1];
rz(-1.2697376) q[1];
sx q[1];
rz(1.7369695) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72337379) q[0];
sx q[0];
rz(-1.9310111) q[0];
sx q[0];
rz(-0.13803137) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6925473) q[2];
sx q[2];
rz(-1.8728313) q[2];
sx q[2];
rz(-1.9027325) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.9060668) q[1];
sx q[1];
rz(-1.1065673) q[1];
sx q[1];
rz(-0.81495754) q[1];
rz(-1.5795454) q[3];
sx q[3];
rz(-1.7213271) q[3];
sx q[3];
rz(2.6237486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.29356062) q[2];
sx q[2];
rz(-0.75389391) q[2];
sx q[2];
rz(0.30612293) q[2];
rz(0.39811578) q[3];
sx q[3];
rz(-1.476215) q[3];
sx q[3];
rz(-1.0242296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33481471) q[0];
sx q[0];
rz(-1.0659185) q[0];
sx q[0];
rz(0.48068) q[0];
rz(-1.9819992) q[1];
sx q[1];
rz(-2.1203142) q[1];
sx q[1];
rz(0.28657985) q[1];
rz(-1.2724614) q[2];
sx q[2];
rz(-1.8229501) q[2];
sx q[2];
rz(2.0291871) q[2];
rz(2.1383193) q[3];
sx q[3];
rz(-0.56837396) q[3];
sx q[3];
rz(3.0336998) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];