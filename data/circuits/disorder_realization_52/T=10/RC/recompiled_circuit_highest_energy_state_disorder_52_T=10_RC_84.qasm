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
rz(-2.4177457) q[0];
sx q[0];
rz(-1.9404193) q[0];
sx q[0];
rz(2.6475651) q[0];
rz(1.5552893) q[1];
sx q[1];
rz(2.2145693) q[1];
sx q[1];
rz(8.5742843) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14334003) q[0];
sx q[0];
rz(-1.7368642) q[0];
sx q[0];
rz(2.2448426) q[0];
x q[1];
rz(1.1095959) q[2];
sx q[2];
rz(-0.96994931) q[2];
sx q[2];
rz(-2.8043562) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.96877484) q[1];
sx q[1];
rz(-2.7451395) q[1];
sx q[1];
rz(-2.5345567) q[1];
rz(-pi) q[2];
rz(-2.7755819) q[3];
sx q[3];
rz(-1.6386392) q[3];
sx q[3];
rz(0.74857601) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.6583307) q[2];
sx q[2];
rz(-0.61507812) q[2];
sx q[2];
rz(-0.94493803) q[2];
rz(0.0065217892) q[3];
sx q[3];
rz(-2.3801453) q[3];
sx q[3];
rz(-2.884088) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1061123) q[0];
sx q[0];
rz(-2.1529614) q[0];
sx q[0];
rz(-2.535787) q[0];
rz(0.93217355) q[1];
sx q[1];
rz(-1.4235539) q[1];
sx q[1];
rz(-0.1056284) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86565844) q[0];
sx q[0];
rz(-1.4167042) q[0];
sx q[0];
rz(2.4529414) q[0];
x q[1];
rz(0.95870734) q[2];
sx q[2];
rz(-2.5664133) q[2];
sx q[2];
rz(-1.5782331) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.82984207) q[1];
sx q[1];
rz(-2.0220322) q[1];
sx q[1];
rz(-1.898502) q[1];
x q[2];
rz(-2.3096931) q[3];
sx q[3];
rz(-2.6100592) q[3];
sx q[3];
rz(0.040415045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.85670829) q[2];
sx q[2];
rz(-1.363089) q[2];
sx q[2];
rz(1.523783) q[2];
rz(2.3273322) q[3];
sx q[3];
rz(-1.7819449) q[3];
sx q[3];
rz(0.8849357) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.043269) q[0];
sx q[0];
rz(-2.3461778) q[0];
sx q[0];
rz(-1.5765618) q[0];
rz(-0.99524975) q[1];
sx q[1];
rz(-2.1627656) q[1];
sx q[1];
rz(2.3547122) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9564215) q[0];
sx q[0];
rz(-1.4354683) q[0];
sx q[0];
rz(1.772373) q[0];
rz(-pi) q[1];
rz(0.22898211) q[2];
sx q[2];
rz(-2.7858284) q[2];
sx q[2];
rz(2.7965429) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7461024) q[1];
sx q[1];
rz(-1.8999892) q[1];
sx q[1];
rz(1.2389061) q[1];
rz(-pi) q[2];
rz(1.6305109) q[3];
sx q[3];
rz(-1.51126) q[3];
sx q[3];
rz(-1.8465896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.83871049) q[2];
sx q[2];
rz(-1.6532712) q[2];
sx q[2];
rz(-0.6380471) q[2];
rz(-0.4979411) q[3];
sx q[3];
rz(-2.0697856) q[3];
sx q[3];
rz(-1.0897442) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8968673) q[0];
sx q[0];
rz(-1.3571955) q[0];
sx q[0];
rz(-3.0778399) q[0];
rz(1.6925192) q[1];
sx q[1];
rz(-1.8359102) q[1];
sx q[1];
rz(1.3099028) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1242715) q[0];
sx q[0];
rz(-1.6505662) q[0];
sx q[0];
rz(-1.5333453) q[0];
rz(-pi) q[1];
rz(-2.9430423) q[2];
sx q[2];
rz(-1.7524002) q[2];
sx q[2];
rz(-1.9924194) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.0367239) q[1];
sx q[1];
rz(-0.83170676) q[1];
sx q[1];
rz(1.402852) q[1];
rz(-2.2579262) q[3];
sx q[3];
rz(-1.851463) q[3];
sx q[3];
rz(1.3786045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.0706851) q[2];
sx q[2];
rz(-2.0347774) q[2];
sx q[2];
rz(3.0359388) q[2];
rz(-1.9581155) q[3];
sx q[3];
rz(-2.2453997) q[3];
sx q[3];
rz(-2.4699874) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79528177) q[0];
sx q[0];
rz(-0.91049894) q[0];
sx q[0];
rz(1.7972535) q[0];
rz(-0.67289871) q[1];
sx q[1];
rz(-1.8310603) q[1];
sx q[1];
rz(3.1281298) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5199842) q[0];
sx q[0];
rz(-2.2604687) q[0];
sx q[0];
rz(-3.1075928) q[0];
rz(-0.014217397) q[2];
sx q[2];
rz(-1.9849615) q[2];
sx q[2];
rz(-2.3384475) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.8720189) q[1];
sx q[1];
rz(-2.7260352) q[1];
sx q[1];
rz(-3.1088105) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7605288) q[3];
sx q[3];
rz(-2.237202) q[3];
sx q[3];
rz(1.2418509) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.54231918) q[2];
sx q[2];
rz(-0.66212526) q[2];
sx q[2];
rz(1.8947961) q[2];
rz(-3.0689734) q[3];
sx q[3];
rz(-2.9473372) q[3];
sx q[3];
rz(2.8323925) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4319864) q[0];
sx q[0];
rz(-2.015634) q[0];
sx q[0];
rz(-3.0288938) q[0];
rz(-0.20507774) q[1];
sx q[1];
rz(-1.3321184) q[1];
sx q[1];
rz(2.233706) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4604707) q[0];
sx q[0];
rz(-2.2325071) q[0];
sx q[0];
rz(2.3813093) q[0];
x q[1];
rz(-1.9992746) q[2];
sx q[2];
rz(-2.1399763) q[2];
sx q[2];
rz(-2.6250397) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.4389613) q[1];
sx q[1];
rz(-1.6621894) q[1];
sx q[1];
rz(0.63321873) q[1];
x q[2];
rz(-0.27067462) q[3];
sx q[3];
rz(-0.68289103) q[3];
sx q[3];
rz(2.0840933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1418566) q[2];
sx q[2];
rz(-2.8047968) q[2];
sx q[2];
rz(-3.0736308) q[2];
rz(-1.9939907) q[3];
sx q[3];
rz(-1.2679029) q[3];
sx q[3];
rz(1.2609153) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1934018) q[0];
sx q[0];
rz(-1.8302487) q[0];
sx q[0];
rz(2.6780658) q[0];
rz(2.9947128) q[1];
sx q[1];
rz(-1.905922) q[1];
sx q[1];
rz(2.1489977) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3459839) q[0];
sx q[0];
rz(-2.3844686) q[0];
sx q[0];
rz(-2.7753995) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2847177) q[2];
sx q[2];
rz(-1.5062638) q[2];
sx q[2];
rz(2.4308824) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.9001683) q[1];
sx q[1];
rz(-2.0219775) q[1];
sx q[1];
rz(-2.3996949) q[1];
rz(2.4086508) q[3];
sx q[3];
rz(-1.9187201) q[3];
sx q[3];
rz(2.3331353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.85840449) q[2];
sx q[2];
rz(-1.7121199) q[2];
sx q[2];
rz(2.6962213) q[2];
rz(-1.3238268) q[3];
sx q[3];
rz(-1.0704853) q[3];
sx q[3];
rz(-1.0858735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0367947) q[0];
sx q[0];
rz(-3.1322271) q[0];
sx q[0];
rz(0.15583663) q[0];
rz(-1.8644631) q[1];
sx q[1];
rz(-1.3469478) q[1];
sx q[1];
rz(0.015930463) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4919419) q[0];
sx q[0];
rz(-1.8109522) q[0];
sx q[0];
rz(0.7911327) q[0];
rz(-pi) q[1];
x q[1];
rz(0.35959776) q[2];
sx q[2];
rz(-2.3162875) q[2];
sx q[2];
rz(1.2957089) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.58305743) q[1];
sx q[1];
rz(-1.0638405) q[1];
sx q[1];
rz(0.81113775) q[1];
rz(-pi) q[2];
rz(1.995581) q[3];
sx q[3];
rz(-2.1652555) q[3];
sx q[3];
rz(-0.068258523) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.5126123) q[2];
sx q[2];
rz(-3.0292558) q[2];
sx q[2];
rz(0.12574276) q[2];
rz(-2.1851152) q[3];
sx q[3];
rz(-1.4230909) q[3];
sx q[3];
rz(2.8023348) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9823031) q[0];
sx q[0];
rz(-0.24022261) q[0];
sx q[0];
rz(-0.45641986) q[0];
rz(-1.5688815) q[1];
sx q[1];
rz(-1.0036889) q[1];
sx q[1];
rz(-0.68663866) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4149041) q[0];
sx q[0];
rz(-2.3285638) q[0];
sx q[0];
rz(-2.4343879) q[0];
rz(-pi) q[1];
rz(-1.3977658) q[2];
sx q[2];
rz(-1.695172) q[2];
sx q[2];
rz(1.5501407) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.2575779) q[1];
sx q[1];
rz(-0.19022372) q[1];
sx q[1];
rz(2.4567) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6410646) q[3];
sx q[3];
rz(-1.3965551) q[3];
sx q[3];
rz(-1.4605923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.75866428) q[2];
sx q[2];
rz(-1.1229346) q[2];
sx q[2];
rz(-1.135896) q[2];
rz(-2.1618333) q[3];
sx q[3];
rz(-1.3330678) q[3];
sx q[3];
rz(-0.69825828) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6887688) q[0];
sx q[0];
rz(-1.8166421) q[0];
sx q[0];
rz(0.45595566) q[0];
rz(-0.042757209) q[1];
sx q[1];
rz(-1.9330934) q[1];
sx q[1];
rz(-2.4483689) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9589577) q[0];
sx q[0];
rz(-1.6696249) q[0];
sx q[0];
rz(1.1846428) q[0];
x q[1];
rz(-1.9749538) q[2];
sx q[2];
rz(-1.9459138) q[2];
sx q[2];
rz(2.1585495) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.0574368) q[1];
sx q[1];
rz(-3.0433214) q[1];
sx q[1];
rz(1.4986503) q[1];
x q[2];
rz(2.6900378) q[3];
sx q[3];
rz(-1.8579036) q[3];
sx q[3];
rz(2.5352458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.2263055) q[2];
sx q[2];
rz(-1.2682356) q[2];
sx q[2];
rz(0.66429663) q[2];
rz(1.213446) q[3];
sx q[3];
rz(-0.72187859) q[3];
sx q[3];
rz(2.2693995) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.450347) q[0];
sx q[0];
rz(-2.3007614) q[0];
sx q[0];
rz(1.7578516) q[0];
rz(1.6830403) q[1];
sx q[1];
rz(-2.8589307) q[1];
sx q[1];
rz(0.9160441) q[1];
rz(-0.073804755) q[2];
sx q[2];
rz(-2.7701785) q[2];
sx q[2];
rz(-2.2899173) q[2];
rz(-0.63587832) q[3];
sx q[3];
rz(-0.81428953) q[3];
sx q[3];
rz(-2.3220111) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
