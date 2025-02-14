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
rz(-0.14577785) q[0];
sx q[0];
rz(-1.3480027) q[0];
sx q[0];
rz(-2.7733127) q[0];
rz(-0.085973099) q[1];
sx q[1];
rz(0.84288016) q[1];
sx q[1];
rz(9.1755484) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1220142) q[0];
sx q[0];
rz(-1.459157) q[0];
sx q[0];
rz(-0.48099244) q[0];
rz(-pi) q[1];
rz(-1.0722798) q[2];
sx q[2];
rz(-0.59293737) q[2];
sx q[2];
rz(-0.007615533) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.18663861) q[1];
sx q[1];
rz(-0.74184275) q[1];
sx q[1];
rz(-0.35278503) q[1];
rz(-2.7830368) q[3];
sx q[3];
rz(-1.6772146) q[3];
sx q[3];
rz(-2.2582683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.2117846) q[2];
sx q[2];
rz(-2.4346508) q[2];
sx q[2];
rz(2.9014273) q[2];
rz(0.54720488) q[3];
sx q[3];
rz(-1.6271084) q[3];
sx q[3];
rz(-2.1521294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7996247) q[0];
sx q[0];
rz(-0.37536055) q[0];
sx q[0];
rz(2.4976835) q[0];
rz(-1.0470692) q[1];
sx q[1];
rz(-1.964485) q[1];
sx q[1];
rz(-2.8783096) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9851259) q[0];
sx q[0];
rz(-1.5432285) q[0];
sx q[0];
rz(0.03027244) q[0];
rz(-3.075454) q[2];
sx q[2];
rz(-1.6355343) q[2];
sx q[2];
rz(-1.3548702) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.4519516) q[1];
sx q[1];
rz(-0.76162377) q[1];
sx q[1];
rz(-0.93859251) q[1];
x q[2];
rz(1.0399148) q[3];
sx q[3];
rz(-2.4270319) q[3];
sx q[3];
rz(1.2604063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.6836267) q[2];
sx q[2];
rz(-1.3764952) q[2];
sx q[2];
rz(-2.4158884) q[2];
rz(-2.8343685) q[3];
sx q[3];
rz(-1.3064462) q[3];
sx q[3];
rz(2.0871625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3211408) q[0];
sx q[0];
rz(-1.4844002) q[0];
sx q[0];
rz(-2.3725574) q[0];
rz(0.10737315) q[1];
sx q[1];
rz(-1.4391876) q[1];
sx q[1];
rz(-0.69951397) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25722593) q[0];
sx q[0];
rz(-2.6623137) q[0];
sx q[0];
rz(0.81247654) q[0];
x q[1];
rz(0.6147521) q[2];
sx q[2];
rz(-1.5745133) q[2];
sx q[2];
rz(0.89064179) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.87772885) q[1];
sx q[1];
rz(-0.45719621) q[1];
sx q[1];
rz(1.3844107) q[1];
rz(-pi) q[2];
rz(2.3885257) q[3];
sx q[3];
rz(-2.1478466) q[3];
sx q[3];
rz(-2.6073539) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7677782) q[2];
sx q[2];
rz(-2.5641597) q[2];
sx q[2];
rz(2.2334297) q[2];
rz(2.5670037) q[3];
sx q[3];
rz(-1.8879954) q[3];
sx q[3];
rz(2.7558034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3278367) q[0];
sx q[0];
rz(-1.3264553) q[0];
sx q[0];
rz(-0.5089708) q[0];
rz(-3.0834037) q[1];
sx q[1];
rz(-1.6845082) q[1];
sx q[1];
rz(1.0675272) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.120507) q[0];
sx q[0];
rz(-0.82553804) q[0];
sx q[0];
rz(-0.50936459) q[0];
x q[1];
rz(1.3994965) q[2];
sx q[2];
rz(-0.41738415) q[2];
sx q[2];
rz(2.4462552) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.4219027) q[1];
sx q[1];
rz(-0.1830782) q[1];
sx q[1];
rz(0.94193926) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7948143) q[3];
sx q[3];
rz(-0.30034143) q[3];
sx q[3];
rz(-3.0413922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.031614583) q[2];
sx q[2];
rz(-1.5441394) q[2];
sx q[2];
rz(0.99008375) q[2];
rz(-1.7456985) q[3];
sx q[3];
rz(-3.0315704) q[3];
sx q[3];
rz(-1.8949738) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9083967) q[0];
sx q[0];
rz(-1.7447423) q[0];
sx q[0];
rz(-2.0821849) q[0];
rz(1.1147095) q[1];
sx q[1];
rz(-1.6672986) q[1];
sx q[1];
rz(1.4310736) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3774893) q[0];
sx q[0];
rz(-1.4015348) q[0];
sx q[0];
rz(-0.15248032) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.25974187) q[2];
sx q[2];
rz(-2.0613823) q[2];
sx q[2];
rz(-2.7502613) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.614735) q[1];
sx q[1];
rz(-2.4522374) q[1];
sx q[1];
rz(2.7562642) q[1];
x q[2];
rz(1.7201172) q[3];
sx q[3];
rz(-1.3589371) q[3];
sx q[3];
rz(-2.6184591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.7913671) q[2];
sx q[2];
rz(-0.75608772) q[2];
sx q[2];
rz(-1.4678601) q[2];
rz(-1.0427467) q[3];
sx q[3];
rz(-2.0839033) q[3];
sx q[3];
rz(2.3500672) q[3];
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
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70822155) q[0];
sx q[0];
rz(-1.2935761) q[0];
sx q[0];
rz(-0.21155393) q[0];
rz(2.4987706) q[1];
sx q[1];
rz(-1.1245518) q[1];
sx q[1];
rz(-1.0483673) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.12473561) q[0];
sx q[0];
rz(-2.2679459) q[0];
sx q[0];
rz(-1.7378896) q[0];
rz(0.17906923) q[2];
sx q[2];
rz(-1.6229651) q[2];
sx q[2];
rz(1.6889926) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.6395289) q[1];
sx q[1];
rz(-2.2090753) q[1];
sx q[1];
rz(1.5997574) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6392751) q[3];
sx q[3];
rz(-2.5279495) q[3];
sx q[3];
rz(0.63921038) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.50531498) q[2];
sx q[2];
rz(-1.9376398) q[2];
sx q[2];
rz(0.25406507) q[2];
rz(2.9680179) q[3];
sx q[3];
rz(-2.5901399) q[3];
sx q[3];
rz(1.180163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59488615) q[0];
sx q[0];
rz(-2.7017024) q[0];
sx q[0];
rz(0.79676262) q[0];
rz(-2.2548389) q[1];
sx q[1];
rz(-1.795307) q[1];
sx q[1];
rz(-0.60060445) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89413801) q[0];
sx q[0];
rz(-2.8464212) q[0];
sx q[0];
rz(2.9516793) q[0];
rz(-pi) q[1];
x q[1];
rz(2.372416) q[2];
sx q[2];
rz(-1.666579) q[2];
sx q[2];
rz(-0.2792007) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.010041865) q[1];
sx q[1];
rz(-2.2497592) q[1];
sx q[1];
rz(0.20874397) q[1];
rz(2.7891385) q[3];
sx q[3];
rz(-1.0648994) q[3];
sx q[3];
rz(-0.91010362) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0299006) q[2];
sx q[2];
rz(-1.4375968) q[2];
sx q[2];
rz(-2.0491484) q[2];
rz(-1.368329) q[3];
sx q[3];
rz(-2.5320801) q[3];
sx q[3];
rz(0.4121367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4261037) q[0];
sx q[0];
rz(-1.1142718) q[0];
sx q[0];
rz(-0.45147595) q[0];
rz(-3.0637947) q[1];
sx q[1];
rz(-2.1092238) q[1];
sx q[1];
rz(0.66158867) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84216994) q[0];
sx q[0];
rz(-1.3515858) q[0];
sx q[0];
rz(1.1482394) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.069839283) q[2];
sx q[2];
rz(-2.0681212) q[2];
sx q[2];
rz(2.2004009) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.80967951) q[1];
sx q[1];
rz(-0.47985489) q[1];
sx q[1];
rz(-1.859568) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.39644661) q[3];
sx q[3];
rz(-1.8155193) q[3];
sx q[3];
rz(2.1978767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0563858) q[2];
sx q[2];
rz(-1.2678009) q[2];
sx q[2];
rz(-0.80643225) q[2];
rz(1.8993529) q[3];
sx q[3];
rz(-0.16568383) q[3];
sx q[3];
rz(0.14036673) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(-0.30271444) q[0];
sx q[0];
rz(-1.1352204) q[0];
sx q[0];
rz(-0.43933991) q[0];
rz(1.9413403) q[1];
sx q[1];
rz(-2.3019583) q[1];
sx q[1];
rz(0.95019788) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2586177) q[0];
sx q[0];
rz(-1.1990917) q[0];
sx q[0];
rz(0.88754334) q[0];
x q[1];
rz(-0.15301159) q[2];
sx q[2];
rz(-1.8540314) q[2];
sx q[2];
rz(2.6161414) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.62225759) q[1];
sx q[1];
rz(-0.37168113) q[1];
sx q[1];
rz(1.8628013) q[1];
rz(-pi) q[2];
rz(2.1860649) q[3];
sx q[3];
rz(-0.63378171) q[3];
sx q[3];
rz(-1.2891226) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.3085559) q[2];
sx q[2];
rz(-2.6335282) q[2];
sx q[2];
rz(1.9527831) q[2];
rz(-3.0360119) q[3];
sx q[3];
rz(-2.2002386) q[3];
sx q[3];
rz(-2.1281435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.062716) q[0];
sx q[0];
rz(-2.833241) q[0];
sx q[0];
rz(2.4421332) q[0];
rz(1.7309011) q[1];
sx q[1];
rz(-0.78838333) q[1];
sx q[1];
rz(2.9404822) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7459262) q[0];
sx q[0];
rz(-2.2494083) q[0];
sx q[0];
rz(0.9890438) q[0];
x q[1];
rz(2.6568065) q[2];
sx q[2];
rz(-1.03063) q[2];
sx q[2];
rz(-0.44819427) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6545668) q[1];
sx q[1];
rz(-0.89201614) q[1];
sx q[1];
rz(2.8553633) q[1];
x q[2];
rz(-1.4803545) q[3];
sx q[3];
rz(-1.3916257) q[3];
sx q[3];
rz(2.1018525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.9451399) q[2];
sx q[2];
rz(-1.2597193) q[2];
sx q[2];
rz(-3.073976) q[2];
rz(-1.128528) q[3];
sx q[3];
rz(-2.0566548) q[3];
sx q[3];
rz(1.5170001) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45831281) q[0];
sx q[0];
rz(-1.0673609) q[0];
sx q[0];
rz(1.3491032) q[0];
rz(1.5051399) q[1];
sx q[1];
rz(-2.9121193) q[1];
sx q[1];
rz(3.0519003) q[1];
rz(0.60293052) q[2];
sx q[2];
rz(-1.7431734) q[2];
sx q[2];
rz(-1.0760104) q[2];
rz(-2.454461) q[3];
sx q[3];
rz(-2.7239068) q[3];
sx q[3];
rz(2.1129029) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
