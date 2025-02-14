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
rz(-0.92702335) q[1];
sx q[1];
rz(0.85049367) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9181831) q[0];
sx q[0];
rz(-0.69108686) q[0];
sx q[0];
rz(1.8331493) q[0];
x q[1];
rz(0.57588864) q[2];
sx q[2];
rz(-2.401899) q[2];
sx q[2];
rz(2.0832555) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.1728178) q[1];
sx q[1];
rz(-0.39645312) q[1];
sx q[1];
rz(-0.60703599) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9539724) q[3];
sx q[3];
rz(-2.7696262) q[3];
sx q[3];
rz(2.1443129) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.4832619) q[2];
sx q[2];
rz(-0.61507812) q[2];
sx q[2];
rz(-0.94493803) q[2];
rz(0.0065217892) q[3];
sx q[3];
rz(-2.3801453) q[3];
sx q[3];
rz(0.25750461) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1061123) q[0];
sx q[0];
rz(-0.98863125) q[0];
sx q[0];
rz(-2.535787) q[0];
rz(0.93217355) q[1];
sx q[1];
rz(-1.7180387) q[1];
sx q[1];
rz(0.1056284) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6208134) q[0];
sx q[0];
rz(-0.70292369) q[0];
sx q[0];
rz(-2.9018794) q[0];
rz(-pi) q[1];
rz(-1.0829686) q[2];
sx q[2];
rz(-1.8886856) q[2];
sx q[2];
rz(-2.6018104) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.16766891) q[1];
sx q[1];
rz(-2.5906256) q[1];
sx q[1];
rz(2.5552804) q[1];
rz(1.1607882) q[3];
sx q[3];
rz(-1.2224397) q[3];
sx q[3];
rz(-0.94545555) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.2848844) q[2];
sx q[2];
rz(-1.7785037) q[2];
sx q[2];
rz(-1.6178097) q[2];
rz(-2.3273322) q[3];
sx q[3];
rz(-1.3596478) q[3];
sx q[3];
rz(0.8849357) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.043269) q[0];
sx q[0];
rz(-2.3461778) q[0];
sx q[0];
rz(-1.5650308) q[0];
rz(-2.1463429) q[1];
sx q[1];
rz(-0.97882706) q[1];
sx q[1];
rz(-0.78688041) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3580619) q[0];
sx q[0];
rz(-1.770505) q[0];
sx q[0];
rz(3.003503) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9126105) q[2];
sx q[2];
rz(-2.7858284) q[2];
sx q[2];
rz(-0.34504978) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.2128558) q[1];
sx q[1];
rz(-0.46312216) q[1];
sx q[1];
rz(-2.379851) q[1];
rz(-pi) q[2];
rz(-3.0819503) q[3];
sx q[3];
rz(-1.5111877) q[3];
sx q[3];
rz(0.27935057) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.3028822) q[2];
sx q[2];
rz(-1.4883214) q[2];
sx q[2];
rz(2.5035456) q[2];
rz(0.4979411) q[3];
sx q[3];
rz(-2.0697856) q[3];
sx q[3];
rz(1.0897442) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2447253) q[0];
sx q[0];
rz(-1.7843972) q[0];
sx q[0];
rz(3.0778399) q[0];
rz(1.4490734) q[1];
sx q[1];
rz(-1.3056825) q[1];
sx q[1];
rz(1.3099028) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44353911) q[0];
sx q[0];
rz(-1.6081282) q[0];
sx q[0];
rz(-3.061767) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3919427) q[2];
sx q[2];
rz(-0.26826619) q[2];
sx q[2];
rz(1.1531354) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.79015649) q[1];
sx q[1];
rz(-0.75439851) q[1];
sx q[1];
rz(-2.960207) q[1];
rz(1.1442361) q[3];
sx q[3];
rz(-2.4080695) q[3];
sx q[3];
rz(-3.008212) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.070907585) q[2];
sx q[2];
rz(-2.0347774) q[2];
sx q[2];
rz(0.10565383) q[2];
rz(-1.9581155) q[3];
sx q[3];
rz(-2.2453997) q[3];
sx q[3];
rz(0.67160523) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79528177) q[0];
sx q[0];
rz(-2.2310937) q[0];
sx q[0];
rz(-1.3443391) q[0];
rz(0.67289871) q[1];
sx q[1];
rz(-1.8310603) q[1];
sx q[1];
rz(-3.1281298) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1707662) q[0];
sx q[0];
rz(-1.544569) q[0];
sx q[0];
rz(-2.2607525) q[0];
rz(-1.6031262) q[2];
sx q[2];
rz(-0.41439498) q[2];
sx q[2];
rz(0.76782819) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.2337444) q[1];
sx q[1];
rz(-1.9861167) q[1];
sx q[1];
rz(1.585258) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1290523) q[3];
sx q[3];
rz(-0.75296445) q[3];
sx q[3];
rz(-0.66679614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.54231918) q[2];
sx q[2];
rz(-2.4794674) q[2];
sx q[2];
rz(1.8947961) q[2];
rz(-0.072619297) q[3];
sx q[3];
rz(-0.19425546) q[3];
sx q[3];
rz(2.8323925) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70960629) q[0];
sx q[0];
rz(-1.1259587) q[0];
sx q[0];
rz(0.11269888) q[0];
rz(2.9365149) q[1];
sx q[1];
rz(-1.3321184) q[1];
sx q[1];
rz(-0.90788666) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7230835) q[0];
sx q[0];
rz(-2.1456133) q[0];
sx q[0];
rz(2.3922343) q[0];
rz(-pi) q[1];
rz(1.9992746) q[2];
sx q[2];
rz(-2.1399763) q[2];
sx q[2];
rz(2.6250397) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.0081387719) q[1];
sx q[1];
rz(-0.63888237) q[1];
sx q[1];
rz(-0.15365803) q[1];
rz(2.870918) q[3];
sx q[3];
rz(-0.68289103) q[3];
sx q[3];
rz(2.0840933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1418566) q[2];
sx q[2];
rz(-0.33679589) q[2];
sx q[2];
rz(3.0736308) q[2];
rz(-1.9939907) q[3];
sx q[3];
rz(-1.8736898) q[3];
sx q[3];
rz(1.8806774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1934018) q[0];
sx q[0];
rz(-1.8302487) q[0];
sx q[0];
rz(2.6780658) q[0];
rz(-2.9947128) q[1];
sx q[1];
rz(-1.905922) q[1];
sx q[1];
rz(-2.1489977) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6382321) q[0];
sx q[0];
rz(-1.8192756) q[0];
sx q[0];
rz(2.4186224) q[0];
rz(-pi) q[1];
rz(3.0743344) q[2];
sx q[2];
rz(-1.8562627) q[2];
sx q[2];
rz(0.87905264) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.0508355) q[1];
sx q[1];
rz(-0.91716424) q[1];
sx q[1];
rz(-0.98934503) q[1];
x q[2];
rz(-2.6448574) q[3];
sx q[3];
rz(-0.79724803) q[3];
sx q[3];
rz(1.1243785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.85840449) q[2];
sx q[2];
rz(-1.7121199) q[2];
sx q[2];
rz(-2.6962213) q[2];
rz(-1.3238268) q[3];
sx q[3];
rz(-2.0711074) q[3];
sx q[3];
rz(-2.0557192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.104798) q[0];
sx q[0];
rz(-3.1322271) q[0];
sx q[0];
rz(2.985756) q[0];
rz(1.8644631) q[1];
sx q[1];
rz(-1.3469478) q[1];
sx q[1];
rz(3.1256622) q[1];
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
rz(-2.35046) q[0];
rz(-pi) q[1];
rz(0.35959776) q[2];
sx q[2];
rz(-2.3162875) q[2];
sx q[2];
rz(1.2957089) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.585938) q[1];
sx q[1];
rz(-0.92472142) q[1];
sx q[1];
rz(-2.487961) q[1];
rz(-pi) q[2];
rz(-2.5032477) q[3];
sx q[3];
rz(-1.2223635) q[3];
sx q[3];
rz(1.750647) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.62898034) q[2];
sx q[2];
rz(-0.11233687) q[2];
sx q[2];
rz(-0.12574276) q[2];
rz(-2.1851152) q[3];
sx q[3];
rz(-1.7185017) q[3];
sx q[3];
rz(-2.8023348) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15928957) q[0];
sx q[0];
rz(-0.24022261) q[0];
sx q[0];
rz(-0.45641986) q[0];
rz(1.5688815) q[1];
sx q[1];
rz(-1.0036889) q[1];
sx q[1];
rz(-2.454954) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68701247) q[0];
sx q[0];
rz(-1.079315) q[0];
sx q[0];
rz(0.67680741) q[0];
rz(-pi) q[1];
rz(-3.0153515) q[2];
sx q[2];
rz(-1.3991157) q[2];
sx q[2];
rz(-3.1405695) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.19018826) q[1];
sx q[1];
rz(-1.7177637) q[1];
sx q[1];
rz(-1.4495871) q[1];
rz(1.500528) q[3];
sx q[3];
rz(-1.3965551) q[3];
sx q[3];
rz(-1.6810004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.75866428) q[2];
sx q[2];
rz(-1.1229346) q[2];
sx q[2];
rz(-1.135896) q[2];
rz(-0.9797594) q[3];
sx q[3];
rz(-1.3330678) q[3];
sx q[3];
rz(-2.4433344) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45282388) q[0];
sx q[0];
rz(-1.3249506) q[0];
sx q[0];
rz(-0.45595566) q[0];
rz(3.0988354) q[1];
sx q[1];
rz(-1.9330934) q[1];
sx q[1];
rz(-2.4483689) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15007818) q[0];
sx q[0];
rz(-0.39798361) q[0];
sx q[0];
rz(-1.8282169) q[0];
x q[1];
rz(-0.78473994) q[2];
sx q[2];
rz(-2.5973136) q[2];
sx q[2];
rz(1.8455102) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.4148414) q[1];
sx q[1];
rz(-1.5637239) q[1];
sx q[1];
rz(-1.6688136) q[1];
rz(-pi) q[2];
rz(-0.59488036) q[3];
sx q[3];
rz(-0.52973807) q[3];
sx q[3];
rz(1.6486419) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.2263055) q[2];
sx q[2];
rz(-1.8733571) q[2];
sx q[2];
rz(-2.477296) q[2];
rz(-1.213446) q[3];
sx q[3];
rz(-2.4197141) q[3];
sx q[3];
rz(-0.87219316) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6912457) q[0];
sx q[0];
rz(-0.84083122) q[0];
sx q[0];
rz(-1.3837411) q[0];
rz(-1.6830403) q[1];
sx q[1];
rz(-0.28266193) q[1];
sx q[1];
rz(-2.2255486) q[1];
rz(2.7710995) q[2];
sx q[2];
rz(-1.5975614) q[2];
sx q[2];
rz(2.3536828) q[2];
rz(2.5057143) q[3];
sx q[3];
rz(-0.81428953) q[3];
sx q[3];
rz(-2.3220111) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
