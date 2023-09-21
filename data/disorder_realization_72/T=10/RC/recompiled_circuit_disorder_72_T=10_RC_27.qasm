OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.5207198) q[0];
sx q[0];
rz(-1.7680661) q[0];
sx q[0];
rz(-1.5078478) q[0];
rz(0.047343407) q[1];
sx q[1];
rz(3.9197796) q[1];
sx q[1];
rz(9.9240886) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0190174) q[0];
sx q[0];
rz(-2.8589773) q[0];
sx q[0];
rz(0.22960381) q[0];
rz(-0.84988611) q[2];
sx q[2];
rz(-1.4404785) q[2];
sx q[2];
rz(0.40638129) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.62568134) q[1];
sx q[1];
rz(-1.7209007) q[1];
sx q[1];
rz(2.4019269) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3189089) q[3];
sx q[3];
rz(-0.10676521) q[3];
sx q[3];
rz(1.8330542) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.9177861) q[2];
sx q[2];
rz(-2.1710158) q[2];
sx q[2];
rz(2.1271558) q[2];
rz(0.23400083) q[3];
sx q[3];
rz(-2.6205385) q[3];
sx q[3];
rz(-0.28449374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
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
rz(-1.6681799) q[0];
sx q[0];
rz(-1.4665335) q[0];
sx q[0];
rz(-1.0043253) q[0];
rz(-1.5197808) q[1];
sx q[1];
rz(-2.2147949) q[1];
sx q[1];
rz(-2.1388334) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8589477) q[0];
sx q[0];
rz(-0.98416057) q[0];
sx q[0];
rz(1.0413917) q[0];
x q[1];
rz(-0.1728671) q[2];
sx q[2];
rz(-1.7555408) q[2];
sx q[2];
rz(2.6181521) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.037307449) q[1];
sx q[1];
rz(-0.98854317) q[1];
sx q[1];
rz(-0.18584713) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6456804) q[3];
sx q[3];
rz(-2.4818588) q[3];
sx q[3];
rz(2.0663313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.26943794) q[2];
sx q[2];
rz(-2.1472011) q[2];
sx q[2];
rz(-0.15094748) q[2];
rz(2.7271467) q[3];
sx q[3];
rz(-0.60025418) q[3];
sx q[3];
rz(3.0533561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8931483) q[0];
sx q[0];
rz(-1.8284766) q[0];
sx q[0];
rz(0.77899581) q[0];
rz(2.7422854) q[1];
sx q[1];
rz(-1.8931959) q[1];
sx q[1];
rz(-0.88358203) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1012264) q[0];
sx q[0];
rz(-1.7352312) q[0];
sx q[0];
rz(0.40584392) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9687037) q[2];
sx q[2];
rz(-1.298561) q[2];
sx q[2];
rz(2.9583601) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.98595847) q[1];
sx q[1];
rz(-1.8556719) q[1];
sx q[1];
rz(-0.023521544) q[1];
rz(1.6885353) q[3];
sx q[3];
rz(-2.2265834) q[3];
sx q[3];
rz(-1.1586787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8399923) q[2];
sx q[2];
rz(-1.2568544) q[2];
sx q[2];
rz(2.7123614) q[2];
rz(0.99003506) q[3];
sx q[3];
rz(-1.0777377) q[3];
sx q[3];
rz(-0.86301962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4317959) q[0];
sx q[0];
rz(-1.3732095) q[0];
sx q[0];
rz(-2.5298932) q[0];
rz(-2.0344095) q[1];
sx q[1];
rz(-2.2829843) q[1];
sx q[1];
rz(0.5501737) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4470091) q[0];
sx q[0];
rz(-2.3059079) q[0];
sx q[0];
rz(1.4674594) q[0];
rz(-0.29454622) q[2];
sx q[2];
rz(-1.6606332) q[2];
sx q[2];
rz(-1.7812658) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.95851129) q[1];
sx q[1];
rz(-1.4554879) q[1];
sx q[1];
rz(2.1469841) q[1];
rz(-1.5189734) q[3];
sx q[3];
rz(-1.7997777) q[3];
sx q[3];
rz(2.1932972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.6198373) q[2];
sx q[2];
rz(-2.6553314) q[2];
sx q[2];
rz(2.8660529) q[2];
rz(3.0299305) q[3];
sx q[3];
rz(-1.2005946) q[3];
sx q[3];
rz(-0.47079852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6977285) q[0];
sx q[0];
rz(-1.0794909) q[0];
sx q[0];
rz(2.2648947) q[0];
rz(-0.69119167) q[1];
sx q[1];
rz(-0.87379876) q[1];
sx q[1];
rz(0.91526389) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1268069) q[0];
sx q[0];
rz(-0.94313398) q[0];
sx q[0];
rz(-0.058237596) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6164242) q[2];
sx q[2];
rz(-1.5319954) q[2];
sx q[2];
rz(0.9135439) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.60587864) q[1];
sx q[1];
rz(-1.5389581) q[1];
sx q[1];
rz(1.4153626) q[1];
x q[2];
rz(-2.3441633) q[3];
sx q[3];
rz(-2.1505822) q[3];
sx q[3];
rz(-0.42339719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.22121945) q[2];
sx q[2];
rz(-2.129014) q[2];
sx q[2];
rz(-0.4635703) q[2];
rz(2.5772337) q[3];
sx q[3];
rz(-2.1488991) q[3];
sx q[3];
rz(0.92818964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0267462) q[0];
sx q[0];
rz(-0.62018728) q[0];
sx q[0];
rz(-1.0790496) q[0];
rz(0.59533978) q[1];
sx q[1];
rz(-0.7535615) q[1];
sx q[1];
rz(-0.39658305) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2561589) q[0];
sx q[0];
rz(-1.116917) q[0];
sx q[0];
rz(-0.28215779) q[0];
rz(-pi) q[1];
rz(1.8814538) q[2];
sx q[2];
rz(-1.1144181) q[2];
sx q[2];
rz(0.059547101) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5571594) q[1];
sx q[1];
rz(-1.4769308) q[1];
sx q[1];
rz(-2.2319016) q[1];
rz(-pi) q[2];
rz(0.35950426) q[3];
sx q[3];
rz(-0.61509575) q[3];
sx q[3];
rz(-2.2744327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.98809272) q[2];
sx q[2];
rz(-0.92553878) q[2];
sx q[2];
rz(1.5552103) q[2];
rz(1.68613) q[3];
sx q[3];
rz(-2.5374135) q[3];
sx q[3];
rz(0.82715183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39847386) q[0];
sx q[0];
rz(-1.159659) q[0];
sx q[0];
rz(0.78480762) q[0];
rz(-1.2706884) q[1];
sx q[1];
rz(-1.7653468) q[1];
sx q[1];
rz(2.0152337) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9046017) q[0];
sx q[0];
rz(-2.8303574) q[0];
sx q[0];
rz(2.8471332) q[0];
rz(-pi) q[1];
rz(-0.14234219) q[2];
sx q[2];
rz(-1.0668313) q[2];
sx q[2];
rz(-0.0056643639) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.2396444) q[1];
sx q[1];
rz(-0.48929729) q[1];
sx q[1];
rz(-0.40892618) q[1];
rz(-pi) q[2];
rz(0.095110006) q[3];
sx q[3];
rz(-1.9169807) q[3];
sx q[3];
rz(2.9444875) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.209098) q[2];
sx q[2];
rz(-1.1065437) q[2];
sx q[2];
rz(0.15360019) q[2];
rz(-2.8364654) q[3];
sx q[3];
rz(-2.025445) q[3];
sx q[3];
rz(-1.7677527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7711733) q[0];
sx q[0];
rz(-2.6035247) q[0];
sx q[0];
rz(-0.11288189) q[0];
rz(1.0007292) q[1];
sx q[1];
rz(-0.74644867) q[1];
sx q[1];
rz(-0.032756068) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5927133) q[0];
sx q[0];
rz(-1.2476876) q[0];
sx q[0];
rz(-2.8723641) q[0];
x q[1];
rz(2.6997386) q[2];
sx q[2];
rz(-2.0896857) q[2];
sx q[2];
rz(1.9372802) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.9084839) q[1];
sx q[1];
rz(-1.1540124) q[1];
sx q[1];
rz(-3.0770739) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7912346) q[3];
sx q[3];
rz(-0.38808295) q[3];
sx q[3];
rz(-0.9006075) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.96413606) q[2];
sx q[2];
rz(-2.5033341) q[2];
sx q[2];
rz(0.62210554) q[2];
rz(-1.1671676) q[3];
sx q[3];
rz(-1.111235) q[3];
sx q[3];
rz(2.7511403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0416097) q[0];
sx q[0];
rz(-2.6384625) q[0];
sx q[0];
rz(1.6149678) q[0];
rz(0.73293066) q[1];
sx q[1];
rz(-0.65299487) q[1];
sx q[1];
rz(-0.48535767) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55886666) q[0];
sx q[0];
rz(-2.9968046) q[0];
sx q[0];
rz(-0.38082122) q[0];
rz(-pi) q[1];
rz(-0.23156667) q[2];
sx q[2];
rz(-0.61083691) q[2];
sx q[2];
rz(-1.4996741) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.8977412) q[1];
sx q[1];
rz(-0.85695367) q[1];
sx q[1];
rz(0.41390151) q[1];
rz(-pi) q[2];
rz(-2.5076809) q[3];
sx q[3];
rz(-0.68448193) q[3];
sx q[3];
rz(-2.4661635) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.0316524) q[2];
sx q[2];
rz(-1.3039219) q[2];
sx q[2];
rz(2.9821441) q[2];
rz(1.4032646) q[3];
sx q[3];
rz(-0.38968971) q[3];
sx q[3];
rz(0.015550912) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6195246) q[0];
sx q[0];
rz(-0.63729006) q[0];
sx q[0];
rz(-1.7074701) q[0];
rz(-1.2592978) q[1];
sx q[1];
rz(-2.1914296) q[1];
sx q[1];
rz(0.5982582) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1372676) q[0];
sx q[0];
rz(-2.8398872) q[0];
sx q[0];
rz(-2.3803677) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5355849) q[2];
sx q[2];
rz(-2.8851644) q[2];
sx q[2];
rz(2.3026349) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.1640537) q[1];
sx q[1];
rz(-2.7643449) q[1];
sx q[1];
rz(1.3110728) q[1];
rz(-1.4062905) q[3];
sx q[3];
rz(-1.7269772) q[3];
sx q[3];
rz(-2.3567049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(3.0593807) q[2];
sx q[2];
rz(-1.8965315) q[2];
sx q[2];
rz(-0.65199488) q[2];
rz(-0.59371289) q[3];
sx q[3];
rz(-1.9605325) q[3];
sx q[3];
rz(1.1317071) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.839529) q[0];
sx q[0];
rz(-2.4840214) q[0];
sx q[0];
rz(-2.1485463) q[0];
rz(1.2364173) q[1];
sx q[1];
rz(-2.1724783) q[1];
sx q[1];
rz(1.9289) q[1];
rz(-0.17157208) q[2];
sx q[2];
rz(-0.62660672) q[2];
sx q[2];
rz(-1.3852711) q[2];
rz(0.89647722) q[3];
sx q[3];
rz(-1.8835246) q[3];
sx q[3];
rz(-2.6261331) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];