OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.62087286) q[0];
sx q[0];
rz(1.7680661) q[0];
sx q[0];
rz(11.058523) q[0];
rz(-3.0942492) q[1];
sx q[1];
rz(-0.77818692) q[1];
sx q[1];
rz(-0.49931061) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7802785) q[0];
sx q[0];
rz(-1.2957934) q[0];
sx q[0];
rz(-1.6367903) q[0];
rz(-pi) q[1];
x q[1];
rz(0.17272858) q[2];
sx q[2];
rz(-2.2842801) q[2];
sx q[2];
rz(-1.2781065) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.0608622) q[1];
sx q[1];
rz(-2.3002491) q[1];
sx q[1];
rz(-1.3688449) q[1];
rz(-pi) q[2];
x q[2];
rz(3.1148881) q[3];
sx q[3];
rz(-1.4674125) q[3];
sx q[3];
rz(1.0552693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.22380655) q[2];
sx q[2];
rz(-2.1710158) q[2];
sx q[2];
rz(1.0144368) q[2];
rz(-0.23400083) q[3];
sx q[3];
rz(-2.6205385) q[3];
sx q[3];
rz(-2.8570989) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6681799) q[0];
sx q[0];
rz(-1.6750591) q[0];
sx q[0];
rz(1.0043253) q[0];
rz(1.5197808) q[1];
sx q[1];
rz(-0.92679778) q[1];
sx q[1];
rz(1.0027592) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53045814) q[0];
sx q[0];
rz(-2.3728752) q[0];
sx q[0];
rz(-2.4918633) q[0];
rz(-2.3147896) q[2];
sx q[2];
rz(-0.25233341) q[2];
sx q[2];
rz(-1.2834872) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.36675378) q[1];
sx q[1];
rz(-0.60791053) q[1];
sx q[1];
rz(-1.8444091) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5428883) q[3];
sx q[3];
rz(-1.2748534) q[3];
sx q[3];
rz(-0.89950365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.26943794) q[2];
sx q[2];
rz(-2.1472011) q[2];
sx q[2];
rz(-2.9906452) q[2];
rz(-2.7271467) q[3];
sx q[3];
rz(-2.5413385) q[3];
sx q[3];
rz(3.0533561) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2484444) q[0];
sx q[0];
rz(-1.3131161) q[0];
sx q[0];
rz(-0.77899581) q[0];
rz(-0.39930725) q[1];
sx q[1];
rz(-1.8931959) q[1];
sx q[1];
rz(2.2580106) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1012264) q[0];
sx q[0];
rz(-1.4063615) q[0];
sx q[0];
rz(-2.7357487) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.29404624) q[2];
sx q[2];
rz(-1.9532734) q[2];
sx q[2];
rz(-1.2750212) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.98595847) q[1];
sx q[1];
rz(-1.8556719) q[1];
sx q[1];
rz(-3.1180711) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6885353) q[3];
sx q[3];
rz(-0.91500926) q[3];
sx q[3];
rz(1.9829139) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8399923) q[2];
sx q[2];
rz(-1.2568544) q[2];
sx q[2];
rz(2.7123614) q[2];
rz(-2.1515576) q[3];
sx q[3];
rz(-2.0638549) q[3];
sx q[3];
rz(-2.278573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7097968) q[0];
sx q[0];
rz(-1.7683832) q[0];
sx q[0];
rz(-0.61169949) q[0];
rz(-1.1071831) q[1];
sx q[1];
rz(-0.8586084) q[1];
sx q[1];
rz(0.5501737) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4470091) q[0];
sx q[0];
rz(-2.3059079) q[0];
sx q[0];
rz(1.6741333) q[0];
rz(-pi) q[1];
x q[1];
rz(0.30087154) q[2];
sx q[2];
rz(-2.8340325) q[2];
sx q[2];
rz(3.0645264) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.78754567) q[1];
sx q[1];
rz(-2.5552632) q[1];
sx q[1];
rz(-1.7802618) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6226193) q[3];
sx q[3];
rz(-1.341815) q[3];
sx q[3];
rz(2.1932972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.6198373) q[2];
sx q[2];
rz(-2.6553314) q[2];
sx q[2];
rz(2.8660529) q[2];
rz(3.0299305) q[3];
sx q[3];
rz(-1.9409981) q[3];
sx q[3];
rz(-2.6707941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4438641) q[0];
sx q[0];
rz(-2.0621018) q[0];
sx q[0];
rz(2.2648947) q[0];
rz(2.450401) q[1];
sx q[1];
rz(-2.2677939) q[1];
sx q[1];
rz(-0.91526389) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4782151) q[0];
sx q[0];
rz(-1.5236679) q[0];
sx q[0];
rz(0.94232725) q[0];
rz(-2.6164242) q[2];
sx q[2];
rz(-1.5319954) q[2];
sx q[2];
rz(2.2280488) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.76449672) q[1];
sx q[1];
rz(-0.15863523) q[1];
sx q[1];
rz(-1.3678958) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3239273) q[3];
sx q[3];
rz(-2.2125803) q[3];
sx q[3];
rz(1.4828009) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.9203732) q[2];
sx q[2];
rz(-2.129014) q[2];
sx q[2];
rz(-0.4635703) q[2];
rz(0.56435895) q[3];
sx q[3];
rz(-0.99269358) q[3];
sx q[3];
rz(0.92818964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1148465) q[0];
sx q[0];
rz(-0.62018728) q[0];
sx q[0];
rz(-2.0625431) q[0];
rz(-0.59533978) q[1];
sx q[1];
rz(-0.7535615) q[1];
sx q[1];
rz(-2.7450096) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88543372) q[0];
sx q[0];
rz(-2.0246756) q[0];
sx q[0];
rz(0.28215779) q[0];
x q[1];
rz(-2.5846892) q[2];
sx q[2];
rz(-2.5957426) q[2];
sx q[2];
rz(-0.57005537) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.086416883) q[1];
sx q[1];
rz(-2.2284818) q[1];
sx q[1];
rz(3.0228826) q[1];
rz(-0.58432213) q[3];
sx q[3];
rz(-1.7752247) q[3];
sx q[3];
rz(2.1401329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.98809272) q[2];
sx q[2];
rz(-0.92553878) q[2];
sx q[2];
rz(-1.5863824) q[2];
rz(-1.4554626) q[3];
sx q[3];
rz(-2.5374135) q[3];
sx q[3];
rz(0.82715183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
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
rz(-0.39847386) q[0];
sx q[0];
rz(-1.9819336) q[0];
sx q[0];
rz(2.356785) q[0];
rz(1.2706884) q[1];
sx q[1];
rz(-1.7653468) q[1];
sx q[1];
rz(1.126359) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9046017) q[0];
sx q[0];
rz(-0.31123529) q[0];
sx q[0];
rz(-2.8471332) q[0];
x q[1];
rz(1.8225841) q[2];
sx q[2];
rz(-2.6195824) q[2];
sx q[2];
rz(2.8587647) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.2396444) q[1];
sx q[1];
rz(-2.6522954) q[1];
sx q[1];
rz(2.7326665) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2231636) q[3];
sx q[3];
rz(-1.4813444) q[3];
sx q[3];
rz(1.4060494) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.9324947) q[2];
sx q[2];
rz(-2.035049) q[2];
sx q[2];
rz(-0.15360019) q[2];
rz(-2.8364654) q[3];
sx q[3];
rz(-2.025445) q[3];
sx q[3];
rz(-1.7677527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7711733) q[0];
sx q[0];
rz(-2.6035247) q[0];
sx q[0];
rz(3.0287108) q[0];
rz(-1.0007292) q[1];
sx q[1];
rz(-2.395144) q[1];
sx q[1];
rz(3.1088366) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8772802) q[0];
sx q[0];
rz(-0.41752975) q[0];
sx q[0];
rz(0.89950048) q[0];
rz(0.92808) q[2];
sx q[2];
rz(-0.66814458) q[2];
sx q[2];
rz(-1.966114) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.9084839) q[1];
sx q[1];
rz(-1.1540124) q[1];
sx q[1];
rz(3.0770739) q[1];
rz(-pi) q[2];
rz(-1.1912187) q[3];
sx q[3];
rz(-1.6536342) q[3];
sx q[3];
rz(-0.4656725) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.1774566) q[2];
sx q[2];
rz(-0.63825858) q[2];
sx q[2];
rz(-0.62210554) q[2];
rz(1.1671676) q[3];
sx q[3];
rz(-1.111235) q[3];
sx q[3];
rz(-2.7511403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.099982925) q[0];
sx q[0];
rz(-0.5031302) q[0];
sx q[0];
rz(1.5266248) q[0];
rz(2.408662) q[1];
sx q[1];
rz(-2.4885978) q[1];
sx q[1];
rz(2.656235) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55886666) q[0];
sx q[0];
rz(-2.9968046) q[0];
sx q[0];
rz(0.38082122) q[0];
rz(0.59818563) q[2];
sx q[2];
rz(-1.4387812) q[2];
sx q[2];
rz(2.8796632) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.24385142) q[1];
sx q[1];
rz(-0.85695367) q[1];
sx q[1];
rz(2.7276911) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0210578) q[3];
sx q[3];
rz(-2.1053208) q[3];
sx q[3];
rz(0.083732097) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.0316524) q[2];
sx q[2];
rz(-1.8376708) q[2];
sx q[2];
rz(-0.1594485) q[2];
rz(1.738328) q[3];
sx q[3];
rz(-0.38968971) q[3];
sx q[3];
rz(3.1260417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52206802) q[0];
sx q[0];
rz(-0.63729006) q[0];
sx q[0];
rz(-1.4341226) q[0];
rz(1.8822949) q[1];
sx q[1];
rz(-0.95016304) q[1];
sx q[1];
rz(2.5433345) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78865096) q[0];
sx q[0];
rz(-1.7876248) q[0];
sx q[0];
rz(1.3593332) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.60600772) q[2];
sx q[2];
rz(-2.8851644) q[2];
sx q[2];
rz(-0.83895776) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.2559402) q[1];
sx q[1];
rz(-1.9347895) q[1];
sx q[1];
rz(3.0401858) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.80501276) q[3];
sx q[3];
rz(-2.9152438) q[3];
sx q[3];
rz(1.602802) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0593807) q[2];
sx q[2];
rz(-1.8965315) q[2];
sx q[2];
rz(0.65199488) q[2];
rz(-0.59371289) q[3];
sx q[3];
rz(-1.1810602) q[3];
sx q[3];
rz(-1.1317071) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.839529) q[0];
sx q[0];
rz(-2.4840214) q[0];
sx q[0];
rz(-2.1485463) q[0];
rz(-1.9051753) q[1];
sx q[1];
rz(-2.1724783) q[1];
sx q[1];
rz(1.9289) q[1];
rz(0.61959735) q[2];
sx q[2];
rz(-1.6710812) q[2];
sx q[2];
rz(-2.8166213) q[2];
rz(2.2451154) q[3];
sx q[3];
rz(-1.258068) q[3];
sx q[3];
rz(0.51545959) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
