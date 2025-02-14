OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.8425771) q[0];
sx q[0];
rz(-0.17248532) q[0];
sx q[0];
rz(-3.1245998) q[0];
rz(3.6530082) q[1];
sx q[1];
rz(3.6379171) q[1];
sx q[1];
rz(12.080893) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4489789) q[0];
sx q[0];
rz(-2.9831121) q[0];
sx q[0];
rz(-2.3628209) q[0];
x q[1];
rz(-2.4713712) q[2];
sx q[2];
rz(-2.1677758) q[2];
sx q[2];
rz(2.1917997) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.37066165) q[1];
sx q[1];
rz(-1.071613) q[1];
sx q[1];
rz(0.73727495) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9157682) q[3];
sx q[3];
rz(-0.32805035) q[3];
sx q[3];
rz(-1.6102745) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.39712507) q[2];
sx q[2];
rz(-3.0391389) q[2];
sx q[2];
rz(2.3490119) q[2];
rz(2.1553195) q[3];
sx q[3];
rz(-1.4873742) q[3];
sx q[3];
rz(0.13315323) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2546286) q[0];
sx q[0];
rz(-1.5381085) q[0];
sx q[0];
rz(-3.0358553) q[0];
rz(0.75791439) q[1];
sx q[1];
rz(-0.71280232) q[1];
sx q[1];
rz(0.47592083) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71005586) q[0];
sx q[0];
rz(-1.8049631) q[0];
sx q[0];
rz(-1.4708999) q[0];
rz(1.7143102) q[2];
sx q[2];
rz(-2.0594308) q[2];
sx q[2];
rz(-1.2554864) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.8393499) q[1];
sx q[1];
rz(-2.1363597) q[1];
sx q[1];
rz(2.9981705) q[1];
x q[2];
rz(-2.845221) q[3];
sx q[3];
rz(-0.94816531) q[3];
sx q[3];
rz(-0.53106703) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.9953352) q[2];
sx q[2];
rz(-2.6671851) q[2];
sx q[2];
rz(1.8246626) q[2];
rz(2.3138192) q[3];
sx q[3];
rz(-1.0931949) q[3];
sx q[3];
rz(0.83724666) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4227609) q[0];
sx q[0];
rz(-1.7976924) q[0];
sx q[0];
rz(-1.2368917) q[0];
rz(-1.749136) q[1];
sx q[1];
rz(-1.720865) q[1];
sx q[1];
rz(-2.4512591) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8906517) q[0];
sx q[0];
rz(-1.1162045) q[0];
sx q[0];
rz(2.5174052) q[0];
x q[1];
rz(3.1375935) q[2];
sx q[2];
rz(-1.1074813) q[2];
sx q[2];
rz(0.057553854) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.54385932) q[1];
sx q[1];
rz(-1.3784474) q[1];
sx q[1];
rz(-1.400977) q[1];
x q[2];
rz(-0.66392939) q[3];
sx q[3];
rz(-2.2888043) q[3];
sx q[3];
rz(1.688886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.51848015) q[2];
sx q[2];
rz(-1.6192351) q[2];
sx q[2];
rz(-3.0677262) q[2];
rz(1.9833924) q[3];
sx q[3];
rz(-1.9246212) q[3];
sx q[3];
rz(1.7346252) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2288007) q[0];
sx q[0];
rz(-3.085629) q[0];
sx q[0];
rz(-0.19293109) q[0];
rz(1.8979161) q[1];
sx q[1];
rz(-1.7453777) q[1];
sx q[1];
rz(0.23304932) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4503964) q[0];
sx q[0];
rz(-1.47627) q[0];
sx q[0];
rz(2.9942542) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2795882) q[2];
sx q[2];
rz(-1.5420015) q[2];
sx q[2];
rz(0.51143247) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.081542) q[1];
sx q[1];
rz(-2.8257408) q[1];
sx q[1];
rz(0.18209034) q[1];
rz(-pi) q[2];
x q[2];
rz(0.47073029) q[3];
sx q[3];
rz(-1.004537) q[3];
sx q[3];
rz(1.2486718) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.9243246) q[2];
sx q[2];
rz(-1.3793719) q[2];
sx q[2];
rz(3.0775089) q[2];
rz(-0.25121769) q[3];
sx q[3];
rz(-2.678674) q[3];
sx q[3];
rz(-0.29138756) q[3];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0578882) q[0];
sx q[0];
rz(-1.2936445) q[0];
sx q[0];
rz(1.8573014) q[0];
rz(-3.1341556) q[1];
sx q[1];
rz(-1.1859272) q[1];
sx q[1];
rz(0.95710212) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83643267) q[0];
sx q[0];
rz(-2.0070939) q[0];
sx q[0];
rz(-2.7899105) q[0];
rz(-pi) q[1];
rz(0.77218036) q[2];
sx q[2];
rz(-0.77754687) q[2];
sx q[2];
rz(-2.5302056) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.19647287) q[1];
sx q[1];
rz(-0.34356657) q[1];
sx q[1];
rz(0.29702227) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.383833) q[3];
sx q[3];
rz(-2.2174944) q[3];
sx q[3];
rz(1.9892429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.1288422) q[2];
sx q[2];
rz(-0.70657554) q[2];
sx q[2];
rz(2.1112704) q[2];
rz(2.4230867) q[3];
sx q[3];
rz(-2.0593144) q[3];
sx q[3];
rz(-2.4350186) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4265823) q[0];
sx q[0];
rz(-2.3937245) q[0];
sx q[0];
rz(2.7744875) q[0];
rz(-0.35588613) q[1];
sx q[1];
rz(-2.0242736) q[1];
sx q[1];
rz(1.5497367) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9420351) q[0];
sx q[0];
rz(-0.74239391) q[0];
sx q[0];
rz(-3.0846157) q[0];
rz(1.3328959) q[2];
sx q[2];
rz(-0.44778433) q[2];
sx q[2];
rz(-0.16668038) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.8425297) q[1];
sx q[1];
rz(-2.6767414) q[1];
sx q[1];
rz(0.47375394) q[1];
x q[2];
rz(-0.14629062) q[3];
sx q[3];
rz(-0.47463575) q[3];
sx q[3];
rz(2.4865347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.1648569) q[2];
sx q[2];
rz(-1.2931436) q[2];
sx q[2];
rz(3.1381651) q[2];
rz(0.88611832) q[3];
sx q[3];
rz(-2.0485853) q[3];
sx q[3];
rz(1.6974983) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43585983) q[0];
sx q[0];
rz(-1.6678565) q[0];
sx q[0];
rz(0.71520299) q[0];
rz(-3.0192979) q[1];
sx q[1];
rz(-1.0194174) q[1];
sx q[1];
rz(-0.14762793) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1081165) q[0];
sx q[0];
rz(-1.2434992) q[0];
sx q[0];
rz(1.2121588) q[0];
x q[1];
rz(-1.1806025) q[2];
sx q[2];
rz(-1.9762594) q[2];
sx q[2];
rz(0.46145876) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.25179204) q[1];
sx q[1];
rz(-1.7528868) q[1];
sx q[1];
rz(0.86343335) q[1];
rz(-3.0711898) q[3];
sx q[3];
rz(-1.8569267) q[3];
sx q[3];
rz(0.64976684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.3132396) q[2];
sx q[2];
rz(-0.30892631) q[2];
sx q[2];
rz(-3.0301136) q[2];
rz(2.3765423) q[3];
sx q[3];
rz(-1.4234411) q[3];
sx q[3];
rz(-0.033871977) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-2.2400804) q[0];
sx q[0];
rz(-0.96918786) q[0];
sx q[0];
rz(-1.0028268) q[0];
rz(-2.3560933) q[1];
sx q[1];
rz(-1.5248884) q[1];
sx q[1];
rz(1.921382) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0058380143) q[0];
sx q[0];
rz(-1.5501115) q[0];
sx q[0];
rz(-0.51405859) q[0];
rz(-pi) q[1];
x q[1];
rz(0.062510291) q[2];
sx q[2];
rz(-1.4803807) q[2];
sx q[2];
rz(0.16363283) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.8250019) q[1];
sx q[1];
rz(-1.8866393) q[1];
sx q[1];
rz(-2.1027673) q[1];
rz(-pi) q[2];
x q[2];
rz(0.93525993) q[3];
sx q[3];
rz(-1.4339841) q[3];
sx q[3];
rz(-1.2350596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9057374) q[2];
sx q[2];
rz(-1.3964272) q[2];
sx q[2];
rz(3.1077969) q[2];
rz(2.3679521) q[3];
sx q[3];
rz(-1.6126817) q[3];
sx q[3];
rz(-2.0788367) q[3];
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
x q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77699023) q[0];
sx q[0];
rz(-2.5211625) q[0];
sx q[0];
rz(2.799209) q[0];
rz(1.7204334) q[1];
sx q[1];
rz(-2.3183289) q[1];
sx q[1];
rz(-1.8470496) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7086805) q[0];
sx q[0];
rz(-1.4231235) q[0];
sx q[0];
rz(-2.7765034) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2850288) q[2];
sx q[2];
rz(-0.98443778) q[2];
sx q[2];
rz(2.7991653) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.61664904) q[1];
sx q[1];
rz(-1.6934062) q[1];
sx q[1];
rz(0.7864237) q[1];
x q[2];
rz(1.7429084) q[3];
sx q[3];
rz(-1.9034608) q[3];
sx q[3];
rz(1.9486959) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.8483868) q[2];
sx q[2];
rz(-1.3853955) q[2];
sx q[2];
rz(2.6050341) q[2];
rz(2.7287591) q[3];
sx q[3];
rz(-0.53240132) q[3];
sx q[3];
rz(2.0280973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(0.31529108) q[0];
sx q[0];
rz(-1.2696126) q[0];
sx q[0];
rz(1.6444561) q[0];
rz(0.81028384) q[1];
sx q[1];
rz(-1.5850001) q[1];
sx q[1];
rz(2.2946045) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4522471) q[0];
sx q[0];
rz(-2.9072177) q[0];
sx q[0];
rz(2.330392) q[0];
x q[1];
rz(-2.2566206) q[2];
sx q[2];
rz(-0.77407265) q[2];
sx q[2];
rz(-2.1669801) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.21928043) q[1];
sx q[1];
rz(-2.4488291) q[1];
sx q[1];
rz(-2.7119066) q[1];
rz(-pi) q[2];
rz(-0.096135898) q[3];
sx q[3];
rz(-1.8435974) q[3];
sx q[3];
rz(0.48855272) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.102313) q[2];
sx q[2];
rz(-1.3156834) q[2];
sx q[2];
rz(2.5073591) q[2];
rz(-2.5041194) q[3];
sx q[3];
rz(-2.0135148) q[3];
sx q[3];
rz(-0.2171966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(2.1524326) q[0];
sx q[0];
rz(-1.9518873) q[0];
sx q[0];
rz(1.963203) q[0];
rz(0.70869008) q[1];
sx q[1];
rz(-1.3460174) q[1];
sx q[1];
rz(-1.8585471) q[1];
rz(1.218956) q[2];
sx q[2];
rz(-1.2165804) q[2];
sx q[2];
rz(2.658398) q[2];
rz(-0.49571464) q[3];
sx q[3];
rz(-2.442822) q[3];
sx q[3];
rz(-1.9030489) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
