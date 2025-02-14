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
rz(-2.0464309) q[0];
sx q[0];
rz(-0.14552966) q[0];
sx q[0];
rz(-2.6382883) q[0];
rz(-1.9637928) q[1];
sx q[1];
rz(-1.5468532) q[1];
sx q[1];
rz(-1.4454747) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58848042) q[0];
sx q[0];
rz(-2.7403826) q[0];
sx q[0];
rz(-0.11102872) q[0];
x q[1];
rz(-2.4769449) q[2];
sx q[2];
rz(-2.2713813) q[2];
sx q[2];
rz(1.0657016) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.1481944) q[1];
sx q[1];
rz(-2.1752417) q[1];
sx q[1];
rz(1.1682214) q[1];
x q[2];
rz(-1.5863933) q[3];
sx q[3];
rz(-1.5216148) q[3];
sx q[3];
rz(1.5893857) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.4989) q[2];
sx q[2];
rz(-1.2314726) q[2];
sx q[2];
rz(1.6373681) q[2];
rz(-0.73689342) q[3];
sx q[3];
rz(-1.6156018) q[3];
sx q[3];
rz(-2.1604497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40642834) q[0];
sx q[0];
rz(-0.46877113) q[0];
sx q[0];
rz(0.51097393) q[0];
rz(-1.3276395) q[1];
sx q[1];
rz(-1.7402486) q[1];
sx q[1];
rz(-0.32646349) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5559306) q[0];
sx q[0];
rz(-1.3071951) q[0];
sx q[0];
rz(1.2946412) q[0];
rz(-0.20997491) q[2];
sx q[2];
rz(-2.2433503) q[2];
sx q[2];
rz(0.58174101) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0357246) q[1];
sx q[1];
rz(-0.99109036) q[1];
sx q[1];
rz(2.7766262) q[1];
x q[2];
rz(-1.9580919) q[3];
sx q[3];
rz(-1.3069469) q[3];
sx q[3];
rz(1.6843759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.1702801) q[2];
sx q[2];
rz(-2.9176517) q[2];
sx q[2];
rz(-1.8453321) q[2];
rz(-0.41804677) q[3];
sx q[3];
rz(-2.1635735) q[3];
sx q[3];
rz(0.70438284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
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
rz(-2.068785) q[0];
sx q[0];
rz(-2.7588221) q[0];
sx q[0];
rz(2.5991154) q[0];
rz(1.3756649) q[1];
sx q[1];
rz(-2.723697) q[1];
sx q[1];
rz(0.03104041) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2707996) q[0];
sx q[0];
rz(-0.36211553) q[0];
sx q[0];
rz(0.55678456) q[0];
rz(-pi) q[1];
rz(1.8607742) q[2];
sx q[2];
rz(-0.6781247) q[2];
sx q[2];
rz(1.6999619) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.7621029) q[1];
sx q[1];
rz(-0.65564686) q[1];
sx q[1];
rz(1.1608221) q[1];
x q[2];
rz(0.55501819) q[3];
sx q[3];
rz(-1.5848098) q[3];
sx q[3];
rz(-1.1323347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.5151908) q[2];
sx q[2];
rz(-0.73651892) q[2];
sx q[2];
rz(-2.314563) q[2];
rz(-0.51860297) q[3];
sx q[3];
rz(-1.1122455) q[3];
sx q[3];
rz(-1.6270858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9943635) q[0];
sx q[0];
rz(-1.8356859) q[0];
sx q[0];
rz(1.9299141) q[0];
rz(-2.0934824) q[1];
sx q[1];
rz(-2.4217889) q[1];
sx q[1];
rz(1.8009708) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1667092) q[0];
sx q[0];
rz(-2.1926342) q[0];
sx q[0];
rz(1.1896672) q[0];
rz(-pi) q[1];
rz(2.9049314) q[2];
sx q[2];
rz(-1.7156895) q[2];
sx q[2];
rz(-1.3587111) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.27599469) q[1];
sx q[1];
rz(-2.137724) q[1];
sx q[1];
rz(-1.2204942) q[1];
rz(-1.1598489) q[3];
sx q[3];
rz(-1.2839497) q[3];
sx q[3];
rz(-0.2528068) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.94053215) q[2];
sx q[2];
rz(-2.9643855) q[2];
sx q[2];
rz(-2.000467) q[2];
rz(-1.6107669) q[3];
sx q[3];
rz(-2.174236) q[3];
sx q[3];
rz(-2.358986) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.57946) q[0];
sx q[0];
rz(-1.1701522) q[0];
sx q[0];
rz(-2.539047) q[0];
rz(2.5613979) q[1];
sx q[1];
rz(-1.5126901) q[1];
sx q[1];
rz(2.3604732) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3042219) q[0];
sx q[0];
rz(-1.2835555) q[0];
sx q[0];
rz(0.17673136) q[0];
rz(-pi) q[1];
rz(-2.2011312) q[2];
sx q[2];
rz(-2.2302719) q[2];
sx q[2];
rz(2.3465921) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.8153861) q[1];
sx q[1];
rz(-1.657173) q[1];
sx q[1];
rz(-2.2934329) q[1];
rz(-0.72569287) q[3];
sx q[3];
rz(-1.7039434) q[3];
sx q[3];
rz(0.6729047) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.96131229) q[2];
sx q[2];
rz(-2.0013516) q[2];
sx q[2];
rz(-2.9126634) q[2];
rz(1.8217575) q[3];
sx q[3];
rz(-1.6596551) q[3];
sx q[3];
rz(0.84407097) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(-2.9127386) q[0];
sx q[0];
rz(-1.5874533) q[0];
sx q[0];
rz(1.1786906) q[0];
rz(1.9121869) q[1];
sx q[1];
rz(-2.7472159) q[1];
sx q[1];
rz(-1.8454525) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.926149) q[0];
sx q[0];
rz(-2.8279732) q[0];
sx q[0];
rz(-2.0263894) q[0];
rz(-pi) q[1];
x q[1];
rz(0.8528233) q[2];
sx q[2];
rz(-0.97495026) q[2];
sx q[2];
rz(0.15534523) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.7414822) q[1];
sx q[1];
rz(-2.0026836) q[1];
sx q[1];
rz(2.6071965) q[1];
rz(-pi) q[2];
rz(-0.48719897) q[3];
sx q[3];
rz(-1.6604843) q[3];
sx q[3];
rz(-1.7023757) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.65427762) q[2];
sx q[2];
rz(-2.2389905) q[2];
sx q[2];
rz(-0.26228341) q[2];
rz(0.30019635) q[3];
sx q[3];
rz(-1.2223949) q[3];
sx q[3];
rz(2.9330971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7009785) q[0];
sx q[0];
rz(-2.746026) q[0];
sx q[0];
rz(-1.3025008) q[0];
rz(-0.66849661) q[1];
sx q[1];
rz(-1.786247) q[1];
sx q[1];
rz(-2.1065333) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9205017) q[0];
sx q[0];
rz(-1.1701487) q[0];
sx q[0];
rz(-2.9095838) q[0];
rz(2.9547353) q[2];
sx q[2];
rz(-1.2070939) q[2];
sx q[2];
rz(-2.0314558) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.7986732) q[1];
sx q[1];
rz(-1.914115) q[1];
sx q[1];
rz(0.080282057) q[1];
rz(-0.97888246) q[3];
sx q[3];
rz(-2.22284) q[3];
sx q[3];
rz(-1.0197717) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-3.0237026) q[2];
sx q[2];
rz(-1.002545) q[2];
sx q[2];
rz(1.4855851) q[2];
rz(-0.28247908) q[3];
sx q[3];
rz(-1.3586724) q[3];
sx q[3];
rz(-0.43454596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3561803) q[0];
sx q[0];
rz(-1.0419351) q[0];
sx q[0];
rz(2.8670512) q[0];
rz(1.2046332) q[1];
sx q[1];
rz(-2.5614673) q[1];
sx q[1];
rz(-2.5801632) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0208541) q[0];
sx q[0];
rz(-2.0971813) q[0];
sx q[0];
rz(0.83197439) q[0];
rz(1.8738633) q[2];
sx q[2];
rz(-0.91721877) q[2];
sx q[2];
rz(2.4619964) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.8703979) q[1];
sx q[1];
rz(-1.2673089) q[1];
sx q[1];
rz(1.9998026) q[1];
x q[2];
rz(-1.2227603) q[3];
sx q[3];
rz(-1.819918) q[3];
sx q[3];
rz(3.0580229) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.56208912) q[2];
sx q[2];
rz(-2.0359437) q[2];
sx q[2];
rz(-1.4037464) q[2];
rz(1.7956519) q[3];
sx q[3];
rz(-0.74508777) q[3];
sx q[3];
rz(-0.48817202) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76029921) q[0];
sx q[0];
rz(-0.60307044) q[0];
sx q[0];
rz(-0.90989939) q[0];
rz(-1.9445885) q[1];
sx q[1];
rz(-1.0531813) q[1];
sx q[1];
rz(-0.88814703) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.040374856) q[0];
sx q[0];
rz(-1.4618317) q[0];
sx q[0];
rz(-0.094281406) q[0];
x q[1];
rz(3.0641563) q[2];
sx q[2];
rz(-2.4803023) q[2];
sx q[2];
rz(-1.8183625) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.8348114) q[1];
sx q[1];
rz(-1.957292) q[1];
sx q[1];
rz(0.67910925) q[1];
x q[2];
rz(-1.8664594) q[3];
sx q[3];
rz(-1.5348892) q[3];
sx q[3];
rz(1.1334238) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0137279) q[2];
sx q[2];
rz(-1.700054) q[2];
sx q[2];
rz(1.0825276) q[2];
rz(0.23076375) q[3];
sx q[3];
rz(-1.328238) q[3];
sx q[3];
rz(-0.48399353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3103264) q[0];
sx q[0];
rz(-0.75208298) q[0];
sx q[0];
rz(0.58151522) q[0];
rz(-0.61559081) q[1];
sx q[1];
rz(-2.1938727) q[1];
sx q[1];
rz(-1.2967348) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6433577) q[0];
sx q[0];
rz(-1.5603934) q[0];
sx q[0];
rz(-1.5525981) q[0];
x q[1];
rz(-3.1148071) q[2];
sx q[2];
rz(-2.7973632) q[2];
sx q[2];
rz(0.29634288) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.91528581) q[1];
sx q[1];
rz(-1.9939594) q[1];
sx q[1];
rz(-2.8598815) q[1];
rz(-pi) q[2];
rz(2.4695685) q[3];
sx q[3];
rz(-1.8906396) q[3];
sx q[3];
rz(2.8120188) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.41946188) q[2];
sx q[2];
rz(-0.1408793) q[2];
sx q[2];
rz(-0.44931832) q[2];
rz(0.32014534) q[3];
sx q[3];
rz(-1.82205) q[3];
sx q[3];
rz(1.2464574) q[3];
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
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.768059) q[0];
sx q[0];
rz(-0.67714416) q[0];
sx q[0];
rz(-1.9376391) q[0];
rz(-1.3846579) q[1];
sx q[1];
rz(-2.90381) q[1];
sx q[1];
rz(-0.79868383) q[1];
rz(-1.2133219) q[2];
sx q[2];
rz(-0.64897474) q[2];
sx q[2];
rz(-2.8233768) q[2];
rz(1.0998691) q[3];
sx q[3];
rz(-1.9178598) q[3];
sx q[3];
rz(-1.9334067) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
