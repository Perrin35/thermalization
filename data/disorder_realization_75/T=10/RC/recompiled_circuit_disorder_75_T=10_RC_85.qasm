OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.5443213) q[0];
sx q[0];
rz(-2.7379524) q[0];
sx q[0];
rz(0.37024745) q[0];
rz(0.20180841) q[1];
sx q[1];
rz(-1.2528074) q[1];
sx q[1];
rz(-1.4189781) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50492935) q[0];
sx q[0];
rz(-1.814827) q[0];
sx q[0];
rz(-1.3958449) q[0];
x q[1];
rz(-2.6909157) q[2];
sx q[2];
rz(-2.5833231) q[2];
sx q[2];
rz(-0.4980586) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.287392) q[1];
sx q[1];
rz(-0.55254793) q[1];
sx q[1];
rz(1.6709177) q[1];
rz(1.8889514) q[3];
sx q[3];
rz(-0.83359026) q[3];
sx q[3];
rz(-2.6300501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.7449164) q[2];
sx q[2];
rz(-1.6884721) q[2];
sx q[2];
rz(-2.7837226) q[2];
rz(2.9499124) q[3];
sx q[3];
rz(-2.7087757) q[3];
sx q[3];
rz(0.55309692) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1382004) q[0];
sx q[0];
rz(-1.9906185) q[0];
sx q[0];
rz(-3.0733118) q[0];
rz(-1.0066907) q[1];
sx q[1];
rz(-0.12193646) q[1];
sx q[1];
rz(-2.4904747) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2383645) q[0];
sx q[0];
rz(-2.4698386) q[0];
sx q[0];
rz(-2.8147459) q[0];
rz(-3.100349) q[2];
sx q[2];
rz(-2.6154499) q[2];
sx q[2];
rz(-1.0893084) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.880289) q[1];
sx q[1];
rz(-1.1297047) q[1];
sx q[1];
rz(1.0198358) q[1];
rz(-pi) q[2];
rz(2.6504374) q[3];
sx q[3];
rz(-1.3782116) q[3];
sx q[3];
rz(1.1112978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.45289257) q[2];
sx q[2];
rz(-2.9850027) q[2];
sx q[2];
rz(-2.2382656) q[2];
rz(-0.7545169) q[3];
sx q[3];
rz(-1.6468331) q[3];
sx q[3];
rz(0.27404684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
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
rz(2.2925401) q[0];
sx q[0];
rz(-0.64340574) q[0];
sx q[0];
rz(0.49355155) q[0];
rz(1.6429398) q[1];
sx q[1];
rz(-0.40619266) q[1];
sx q[1];
rz(1.0292056) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8619336) q[0];
sx q[0];
rz(-2.5760981) q[0];
sx q[0];
rz(-2.4740567) q[0];
x q[1];
rz(2.7846787) q[2];
sx q[2];
rz(-1.1733574) q[2];
sx q[2];
rz(-2.7473161) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5268847) q[1];
sx q[1];
rz(-1.7803368) q[1];
sx q[1];
rz(-0.31063147) q[1];
rz(-2.1152705) q[3];
sx q[3];
rz(-1.4635651) q[3];
sx q[3];
rz(-2.2162007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.9323953) q[2];
sx q[2];
rz(-2.5019427) q[2];
sx q[2];
rz(1.2970682) q[2];
rz(0.57859892) q[3];
sx q[3];
rz(-1.9208627) q[3];
sx q[3];
rz(-2.4659618) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(1.9349174) q[0];
sx q[0];
rz(-2.4949555) q[0];
sx q[0];
rz(0.40859616) q[0];
rz(1.4273377) q[1];
sx q[1];
rz(-1.4988377) q[1];
sx q[1];
rz(-0.88159195) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5170167) q[0];
sx q[0];
rz(-1.3991742) q[0];
sx q[0];
rz(-2.4044328) q[0];
x q[1];
rz(0.66742113) q[2];
sx q[2];
rz(-0.9657269) q[2];
sx q[2];
rz(-2.8990926) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.2243005) q[1];
sx q[1];
rz(-0.91887337) q[1];
sx q[1];
rz(0.93445458) q[1];
rz(-pi) q[2];
x q[2];
rz(0.41739695) q[3];
sx q[3];
rz(-1.6317123) q[3];
sx q[3];
rz(0.062373769) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0980229) q[2];
sx q[2];
rz(-1.6514401) q[2];
sx q[2];
rz(-1.1126474) q[2];
rz(-2.5860795) q[3];
sx q[3];
rz(-1.3534618) q[3];
sx q[3];
rz(-1.6120733) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1495789) q[0];
sx q[0];
rz(-1.1554138) q[0];
sx q[0];
rz(1.7768815) q[0];
rz(-0.31750202) q[1];
sx q[1];
rz(-0.96173871) q[1];
sx q[1];
rz(-0.11725765) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.011615959) q[0];
sx q[0];
rz(-2.5419309) q[0];
sx q[0];
rz(2.9090803) q[0];
rz(1.7848894) q[2];
sx q[2];
rz(-2.1630686) q[2];
sx q[2];
rz(-1.585373) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.6834516) q[1];
sx q[1];
rz(-1.7675753) q[1];
sx q[1];
rz(1.7276006) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4555143) q[3];
sx q[3];
rz(-0.59024631) q[3];
sx q[3];
rz(0.78199996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.34510288) q[2];
sx q[2];
rz(-1.3619276) q[2];
sx q[2];
rz(-1.3797181) q[2];
rz(-1.1896677) q[3];
sx q[3];
rz(-0.15933557) q[3];
sx q[3];
rz(3.0670847) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18774408) q[0];
sx q[0];
rz(-1.6402316) q[0];
sx q[0];
rz(0.70621079) q[0];
rz(-2.0300991) q[1];
sx q[1];
rz(-2.5128384) q[1];
sx q[1];
rz(0.10791735) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4554493) q[0];
sx q[0];
rz(-2.3685072) q[0];
sx q[0];
rz(-2.7129052) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.040060476) q[2];
sx q[2];
rz(-1.5821579) q[2];
sx q[2];
rz(0.60919112) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.38938746) q[1];
sx q[1];
rz(-2.2332193) q[1];
sx q[1];
rz(-1.2090769) q[1];
rz(-pi) q[2];
rz(0.48072731) q[3];
sx q[3];
rz(-2.2016964) q[3];
sx q[3];
rz(2.7709099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.32249054) q[2];
sx q[2];
rz(-2.1671447) q[2];
sx q[2];
rz(-1.1090013) q[2];
rz(1.8479944) q[3];
sx q[3];
rz(-1.785991) q[3];
sx q[3];
rz(-3.0343645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7727707) q[0];
sx q[0];
rz(-1.4819551) q[0];
sx q[0];
rz(-3.1266881) q[0];
rz(0.42124721) q[1];
sx q[1];
rz(-1.0528456) q[1];
sx q[1];
rz(-0.79963911) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5478202) q[0];
sx q[0];
rz(-0.9043588) q[0];
sx q[0];
rz(-1.244546) q[0];
rz(-pi) q[1];
rz(1.854935) q[2];
sx q[2];
rz(-1.26789) q[2];
sx q[2];
rz(-0.45380935) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.15496847) q[1];
sx q[1];
rz(-1.5151549) q[1];
sx q[1];
rz(0.14543047) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3560489) q[3];
sx q[3];
rz(-2.5657006) q[3];
sx q[3];
rz(-1.4240595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.556276) q[2];
sx q[2];
rz(-0.63627807) q[2];
sx q[2];
rz(0.9220534) q[2];
rz(-1.3098035) q[3];
sx q[3];
rz(-1.9366879) q[3];
sx q[3];
rz(0.95782763) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17469445) q[0];
sx q[0];
rz(-2.0860724) q[0];
sx q[0];
rz(2.877537) q[0];
rz(1.7407725) q[1];
sx q[1];
rz(-1.4512647) q[1];
sx q[1];
rz(0.31731269) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6415629) q[0];
sx q[0];
rz(-0.39425685) q[0];
sx q[0];
rz(1.3803598) q[0];
x q[1];
rz(-1.9192341) q[2];
sx q[2];
rz(-1.4197883) q[2];
sx q[2];
rz(2.3210971) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.4006153) q[1];
sx q[1];
rz(-1.4794011) q[1];
sx q[1];
rz(-0.11465794) q[1];
rz(-0.21934261) q[3];
sx q[3];
rz(-0.6570802) q[3];
sx q[3];
rz(0.59857644) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.43549609) q[2];
sx q[2];
rz(-2.2075682) q[2];
sx q[2];
rz(1.6395817) q[2];
rz(2.8912985) q[3];
sx q[3];
rz(-1.7264265) q[3];
sx q[3];
rz(-1.9992453) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89649993) q[0];
sx q[0];
rz(-2.8332906) q[0];
sx q[0];
rz(1.7171575) q[0];
rz(-2.6622488) q[1];
sx q[1];
rz(-1.4762676) q[1];
sx q[1];
rz(0.11553484) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.696366) q[0];
sx q[0];
rz(-2.1292902) q[0];
sx q[0];
rz(0.88748705) q[0];
rz(-0.71474448) q[2];
sx q[2];
rz(-1.7822767) q[2];
sx q[2];
rz(-1.5874869) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.86347307) q[1];
sx q[1];
rz(-1.1520471) q[1];
sx q[1];
rz(-2.8249885) q[1];
rz(-1.6910016) q[3];
sx q[3];
rz(-2.2634441) q[3];
sx q[3];
rz(1.8409178) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.4775548) q[2];
sx q[2];
rz(-1.1096191) q[2];
sx q[2];
rz(1.7158562) q[2];
rz(-1.7410295) q[3];
sx q[3];
rz(-1.0481342) q[3];
sx q[3];
rz(-0.28276309) q[3];
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
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7857159) q[0];
sx q[0];
rz(-0.73369217) q[0];
sx q[0];
rz(2.8724331) q[0];
rz(-2.0458938) q[1];
sx q[1];
rz(-2.2311189) q[1];
sx q[1];
rz(1.3795308) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3740765) q[0];
sx q[0];
rz(-0.19664581) q[0];
sx q[0];
rz(2.1392512) q[0];
rz(-pi) q[1];
rz(-1.1534259) q[2];
sx q[2];
rz(-1.6438612) q[2];
sx q[2];
rz(0.086581143) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1684432) q[1];
sx q[1];
rz(-1.7440376) q[1];
sx q[1];
rz(3.1107305) q[1];
rz(-2.8096301) q[3];
sx q[3];
rz(-1.8649857) q[3];
sx q[3];
rz(2.4320137) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.049008869) q[2];
sx q[2];
rz(-0.2139341) q[2];
sx q[2];
rz(1.6258378) q[2];
rz(1.1670636) q[3];
sx q[3];
rz(-1.556282) q[3];
sx q[3];
rz(-1.0313755) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9243069) q[0];
sx q[0];
rz(-1.8368245) q[0];
sx q[0];
rz(0.40689847) q[0];
rz(-0.45463195) q[1];
sx q[1];
rz(-1.1063207) q[1];
sx q[1];
rz(2.8938821) q[1];
rz(-2.0438647) q[2];
sx q[2];
rz(-1.4164783) q[2];
sx q[2];
rz(2.867792) q[2];
rz(-1.8299673) q[3];
sx q[3];
rz(-1.5117241) q[3];
sx q[3];
rz(-1.1321887) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
