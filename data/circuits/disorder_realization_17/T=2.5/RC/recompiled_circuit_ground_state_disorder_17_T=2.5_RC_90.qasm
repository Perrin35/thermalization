OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.390653) q[0];
sx q[0];
rz(-0.57354623) q[0];
sx q[0];
rz(3.1107922) q[0];
rz(-2.7544694) q[1];
sx q[1];
rz(-0.20063278) q[1];
sx q[1];
rz(0.90955847) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5873588) q[0];
sx q[0];
rz(-2.2116714) q[0];
sx q[0];
rz(-2.4184722) q[0];
x q[1];
rz(-2.4417905) q[2];
sx q[2];
rz(-1.095311) q[2];
sx q[2];
rz(0.83719992) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.8821669) q[1];
sx q[1];
rz(-1.7710238) q[1];
sx q[1];
rz(3.0012896) q[1];
rz(-1.8424787) q[3];
sx q[3];
rz(-1.6713023) q[3];
sx q[3];
rz(2.1372014) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.51547852) q[2];
sx q[2];
rz(-1.4592183) q[2];
sx q[2];
rz(2.1558351) q[2];
rz(1.23729) q[3];
sx q[3];
rz(-2.3111549) q[3];
sx q[3];
rz(-1.5498836) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80376959) q[0];
sx q[0];
rz(-1.060312) q[0];
sx q[0];
rz(2.7983303) q[0];
rz(-0.23974165) q[1];
sx q[1];
rz(-2.7570351) q[1];
sx q[1];
rz(-0.029031001) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.507258) q[0];
sx q[0];
rz(-2.0478978) q[0];
sx q[0];
rz(1.3647337) q[0];
rz(-pi) q[1];
rz(2.4137561) q[2];
sx q[2];
rz(-1.5994047) q[2];
sx q[2];
rz(-1.4403573) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(3.1085384) q[1];
sx q[1];
rz(-1.5863451) q[1];
sx q[1];
rz(-2.5485746) q[1];
rz(2.4655627) q[3];
sx q[3];
rz(-1.8799729) q[3];
sx q[3];
rz(-1.3622854) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.79641882) q[2];
sx q[2];
rz(-2.5721305) q[2];
sx q[2];
rz(0.84863895) q[2];
rz(0.38550115) q[3];
sx q[3];
rz(-1.6698839) q[3];
sx q[3];
rz(-0.29335415) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.327453) q[0];
sx q[0];
rz(-1.7610981) q[0];
sx q[0];
rz(0.13963786) q[0];
rz(2.8756554) q[1];
sx q[1];
rz(-1.9640924) q[1];
sx q[1];
rz(-1.0136484) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3861198) q[0];
sx q[0];
rz(-1.1869703) q[0];
sx q[0];
rz(-0.95037937) q[0];
rz(-pi) q[1];
x q[1];
rz(2.707194) q[2];
sx q[2];
rz(-1.1936099) q[2];
sx q[2];
rz(2.3247256) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.2964281) q[1];
sx q[1];
rz(-1.0966874) q[1];
sx q[1];
rz(-1.7557322) q[1];
rz(3.072282) q[3];
sx q[3];
rz(-1.6737752) q[3];
sx q[3];
rz(1.5147095) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.2024112) q[2];
sx q[2];
rz(-2.1911669) q[2];
sx q[2];
rz(0.44516304) q[2];
rz(-2.4116481) q[3];
sx q[3];
rz(-1.688262) q[3];
sx q[3];
rz(1.9385519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3438123) q[0];
sx q[0];
rz(-2.7773363) q[0];
sx q[0];
rz(1.9601747) q[0];
rz(1.6254609) q[1];
sx q[1];
rz(-1.5189891) q[1];
sx q[1];
rz(2.6536062) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0868091) q[0];
sx q[0];
rz(-1.8300608) q[0];
sx q[0];
rz(1.4075341) q[0];
x q[1];
rz(0.67537157) q[2];
sx q[2];
rz(-2.3006735) q[2];
sx q[2];
rz(-2.5414064) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5023971) q[1];
sx q[1];
rz(-2.4885414) q[1];
sx q[1];
rz(-2.0060397) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9259629) q[3];
sx q[3];
rz(-1.9653335) q[3];
sx q[3];
rz(-1.153314) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.2130012) q[2];
sx q[2];
rz(-1.7276305) q[2];
sx q[2];
rz(-0.50814) q[2];
rz(-2.7022434) q[3];
sx q[3];
rz(-2.5033247) q[3];
sx q[3];
rz(0.58258575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4947434) q[0];
sx q[0];
rz(-1.9288394) q[0];
sx q[0];
rz(-0.61432046) q[0];
rz(0.99614227) q[1];
sx q[1];
rz(-2.3042945) q[1];
sx q[1];
rz(3.0706629) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.030132564) q[0];
sx q[0];
rz(-1.2637402) q[0];
sx q[0];
rz(-2.2024406) q[0];
rz(-pi) q[1];
rz(-1.0164169) q[2];
sx q[2];
rz(-0.51018894) q[2];
sx q[2];
rz(-1.73908) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6560737) q[1];
sx q[1];
rz(-0.89886256) q[1];
sx q[1];
rz(-2.333332) q[1];
rz(-pi) q[2];
rz(-0.89226858) q[3];
sx q[3];
rz(-2.1315977) q[3];
sx q[3];
rz(-2.4099007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.4225685) q[2];
sx q[2];
rz(-1.4113798) q[2];
sx q[2];
rz(0.076449797) q[2];
rz(0.092175305) q[3];
sx q[3];
rz(-1.1276827) q[3];
sx q[3];
rz(2.5341212) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3344264) q[0];
sx q[0];
rz(-0.18944117) q[0];
sx q[0];
rz(-0.46519753) q[0];
rz(3.0472962) q[1];
sx q[1];
rz(-1.0153208) q[1];
sx q[1];
rz(-1.2634855) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.87288) q[0];
sx q[0];
rz(-1.2194192) q[0];
sx q[0];
rz(-0.12704487) q[0];
rz(-0.63226065) q[2];
sx q[2];
rz(-1.5692729) q[2];
sx q[2];
rz(-1.8830118) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.81663075) q[1];
sx q[1];
rz(-2.0010608) q[1];
sx q[1];
rz(-2.5041786) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6776836) q[3];
sx q[3];
rz(-1.2281443) q[3];
sx q[3];
rz(-3.1009348) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.59244001) q[2];
sx q[2];
rz(-2.0024039) q[2];
sx q[2];
rz(0.69063866) q[2];
rz(1.45951) q[3];
sx q[3];
rz(-3.1091318) q[3];
sx q[3];
rz(-0.261511) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9948147) q[0];
sx q[0];
rz(-2.1187466) q[0];
sx q[0];
rz(-0.47635517) q[0];
rz(-1.9626544) q[1];
sx q[1];
rz(-2.4504688) q[1];
sx q[1];
rz(-1.3271416) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1803197) q[0];
sx q[0];
rz(-1.1256582) q[0];
sx q[0];
rz(2.0044786) q[0];
x q[1];
rz(1.0977237) q[2];
sx q[2];
rz(-1.5247727) q[2];
sx q[2];
rz(-1.0180576) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.33846291) q[1];
sx q[1];
rz(-1.5005365) q[1];
sx q[1];
rz(-0.52609013) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5487212) q[3];
sx q[3];
rz(-2.1855003) q[3];
sx q[3];
rz(0.90329492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.8644766) q[2];
sx q[2];
rz(-1.6552552) q[2];
sx q[2];
rz(0.76249301) q[2];
rz(1.7755194) q[3];
sx q[3];
rz(-1.6906747) q[3];
sx q[3];
rz(-0.33179992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(0.45258006) q[0];
sx q[0];
rz(-3.0735425) q[0];
sx q[0];
rz(-2.3633603) q[0];
rz(-1.4100086) q[1];
sx q[1];
rz(-2.2739613) q[1];
sx q[1];
rz(2.5708139) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2222969) q[0];
sx q[0];
rz(-1.9864829) q[0];
sx q[0];
rz(-2.2908874) q[0];
rz(-pi) q[1];
x q[1];
rz(0.052841803) q[2];
sx q[2];
rz(-1.853301) q[2];
sx q[2];
rz(1.5898012) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.60940424) q[1];
sx q[1];
rz(-1.6488607) q[1];
sx q[1];
rz(1.0482657) q[1];
rz(-pi) q[2];
x q[2];
rz(3.1166534) q[3];
sx q[3];
rz(-0.93526269) q[3];
sx q[3];
rz(2.3221515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0290231) q[2];
sx q[2];
rz(-1.7971635) q[2];
sx q[2];
rz(0.53978187) q[2];
rz(-2.6953186) q[3];
sx q[3];
rz(-2.4036784) q[3];
sx q[3];
rz(-2.3257183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(3.0453294) q[0];
sx q[0];
rz(-0.49880609) q[0];
sx q[0];
rz(0.86736429) q[0];
rz(1.7034886) q[1];
sx q[1];
rz(-2.4080364) q[1];
sx q[1];
rz(2.6677456) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60199262) q[0];
sx q[0];
rz(-1.8162382) q[0];
sx q[0];
rz(0.28781366) q[0];
rz(-pi) q[1];
rz(-2.1146853) q[2];
sx q[2];
rz(-2.4615114) q[2];
sx q[2];
rz(0.56267738) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.074939097) q[1];
sx q[1];
rz(-2.3015296) q[1];
sx q[1];
rz(-1.7615307) q[1];
rz(-0.881265) q[3];
sx q[3];
rz(-1.1598969) q[3];
sx q[3];
rz(-1.9104065) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.69502407) q[2];
sx q[2];
rz(-0.5646255) q[2];
sx q[2];
rz(-2.0107644) q[2];
rz(1.9628149) q[3];
sx q[3];
rz(-1.7480353) q[3];
sx q[3];
rz(2.668837) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9216565) q[0];
sx q[0];
rz(-2.642785) q[0];
sx q[0];
rz(2.1529799) q[0];
rz(-2.0693208) q[1];
sx q[1];
rz(-1.6055454) q[1];
sx q[1];
rz(-1.706121) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50668119) q[0];
sx q[0];
rz(-1.5959863) q[0];
sx q[0];
rz(0.9947302) q[0];
rz(-pi) q[1];
x q[1];
rz(3.1085787) q[2];
sx q[2];
rz(-0.22840127) q[2];
sx q[2];
rz(-1.8637125) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.4354108) q[1];
sx q[1];
rz(-1.5110043) q[1];
sx q[1];
rz(-0.00018187025) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7434968) q[3];
sx q[3];
rz(-0.73258607) q[3];
sx q[3];
rz(-1.7375515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2995305) q[2];
sx q[2];
rz(-0.64075035) q[2];
sx q[2];
rz(-2.2068742) q[2];
rz(-1.0457906) q[3];
sx q[3];
rz(-0.17242923) q[3];
sx q[3];
rz(-2.893253) q[3];
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
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0228148) q[0];
sx q[0];
rz(-1.7775443) q[0];
sx q[0];
rz(-1.4015629) q[0];
rz(-0.10764311) q[1];
sx q[1];
rz(-1.5168774) q[1];
sx q[1];
rz(0.37227896) q[1];
rz(-0.27009985) q[2];
sx q[2];
rz(-2.3348689) q[2];
sx q[2];
rz(0.058090799) q[2];
rz(2.6519898) q[3];
sx q[3];
rz(-2.1402749) q[3];
sx q[3];
rz(-0.28771852) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
