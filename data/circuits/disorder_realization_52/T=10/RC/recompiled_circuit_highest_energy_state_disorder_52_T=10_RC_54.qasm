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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5587529) q[0];
sx q[0];
rz(-2.2338998) q[0];
sx q[0];
rz(-0.21132529) q[0];
x q[1];
rz(-0.65324983) q[2];
sx q[2];
rz(-1.1949277) q[2];
sx q[2];
rz(-2.1819161) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.1728178) q[1];
sx q[1];
rz(-2.7451395) q[1];
sx q[1];
rz(2.5345567) q[1];
rz(2.7755819) q[3];
sx q[3];
rz(-1.5029534) q[3];
sx q[3];
rz(-2.3930166) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.6583307) q[2];
sx q[2];
rz(-0.61507812) q[2];
sx q[2];
rz(0.94493803) q[2];
rz(0.0065217892) q[3];
sx q[3];
rz(-0.76144731) q[3];
sx q[3];
rz(2.884088) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1061123) q[0];
sx q[0];
rz(-2.1529614) q[0];
sx q[0];
rz(-0.60580564) q[0];
rz(2.2094191) q[1];
sx q[1];
rz(-1.7180387) q[1];
sx q[1];
rz(3.0359643) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2759342) q[0];
sx q[0];
rz(-1.4167042) q[0];
sx q[0];
rz(-2.4529414) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0829686) q[2];
sx q[2];
rz(-1.8886856) q[2];
sx q[2];
rz(2.6018104) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.2534598) q[1];
sx q[1];
rz(-1.8646629) q[1];
sx q[1];
rz(0.47305686) q[1];
x q[2];
rz(-0.37704269) q[3];
sx q[3];
rz(-1.1867513) q[3];
sx q[3];
rz(2.3689601) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.85670829) q[2];
sx q[2];
rz(-1.7785037) q[2];
sx q[2];
rz(1.6178097) q[2];
rz(-0.81426042) q[3];
sx q[3];
rz(-1.3596478) q[3];
sx q[3];
rz(2.2566569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.043269) q[0];
sx q[0];
rz(-0.79541484) q[0];
sx q[0];
rz(1.5765618) q[0];
rz(0.99524975) q[1];
sx q[1];
rz(-2.1627656) q[1];
sx q[1];
rz(-2.3547122) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7835307) q[0];
sx q[0];
rz(-1.770505) q[0];
sx q[0];
rz(0.13808967) q[0];
rz(-pi) q[1];
rz(2.9126105) q[2];
sx q[2];
rz(-2.7858284) q[2];
sx q[2];
rz(0.34504978) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.7461024) q[1];
sx q[1];
rz(-1.2416035) q[1];
sx q[1];
rz(1.9026865) q[1];
x q[2];
rz(-0.059642369) q[3];
sx q[3];
rz(-1.630405) q[3];
sx q[3];
rz(0.27935057) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3028822) q[2];
sx q[2];
rz(-1.6532712) q[2];
sx q[2];
rz(-2.5035456) q[2];
rz(0.4979411) q[3];
sx q[3];
rz(-2.0697856) q[3];
sx q[3];
rz(-2.0518484) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(1.2447253) q[0];
sx q[0];
rz(-1.7843972) q[0];
sx q[0];
rz(-0.063752739) q[0];
rz(-1.4490734) q[1];
sx q[1];
rz(-1.8359102) q[1];
sx q[1];
rz(-1.8316899) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44353911) q[0];
sx q[0];
rz(-1.5334645) q[0];
sx q[0];
rz(-0.079825618) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9430423) q[2];
sx q[2];
rz(-1.7524002) q[2];
sx q[2];
rz(1.9924194) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.0367239) q[1];
sx q[1];
rz(-0.83170676) q[1];
sx q[1];
rz(-1.7387407) q[1];
rz(0.3569266) q[3];
sx q[3];
rz(-0.91536575) q[3];
sx q[3];
rz(-2.7259158) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-3.0706851) q[2];
sx q[2];
rz(-1.1068152) q[2];
sx q[2];
rz(-3.0359388) q[2];
rz(1.9581155) q[3];
sx q[3];
rz(-0.89619291) q[3];
sx q[3];
rz(0.67160523) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3463109) q[0];
sx q[0];
rz(-2.2310937) q[0];
sx q[0];
rz(1.7972535) q[0];
rz(-0.67289871) q[1];
sx q[1];
rz(-1.8310603) q[1];
sx q[1];
rz(-0.013462822) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5199842) q[0];
sx q[0];
rz(-0.88112393) q[0];
sx q[0];
rz(-0.033999835) q[0];
rz(-pi) q[1];
rz(3.1273753) q[2];
sx q[2];
rz(-1.1566312) q[2];
sx q[2];
rz(-0.80314512) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.8720189) q[1];
sx q[1];
rz(-2.7260352) q[1];
sx q[1];
rz(3.1088105) q[1];
rz(-pi) q[2];
rz(1.1290523) q[3];
sx q[3];
rz(-2.3886282) q[3];
sx q[3];
rz(2.4747965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.5992735) q[2];
sx q[2];
rz(-0.66212526) q[2];
sx q[2];
rz(1.2467965) q[2];
rz(3.0689734) q[3];
sx q[3];
rz(-2.9473372) q[3];
sx q[3];
rz(-2.8323925) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4319864) q[0];
sx q[0];
rz(-1.1259587) q[0];
sx q[0];
rz(0.11269888) q[0];
rz(2.9365149) q[1];
sx q[1];
rz(-1.8094742) q[1];
sx q[1];
rz(0.90788666) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68112198) q[0];
sx q[0];
rz(-2.2325071) q[0];
sx q[0];
rz(0.76028334) q[0];
rz(-1.142318) q[2];
sx q[2];
rz(-2.1399763) q[2];
sx q[2];
rz(-0.5165529) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.19872936) q[1];
sx q[1];
rz(-2.2009549) q[1];
sx q[1];
rz(-1.4575926) q[1];
rz(-2.870918) q[3];
sx q[3];
rz(-0.68289103) q[3];
sx q[3];
rz(1.0574993) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.9997361) q[2];
sx q[2];
rz(-0.33679589) q[2];
sx q[2];
rz(0.067961819) q[2];
rz(1.9939907) q[3];
sx q[3];
rz(-1.2679029) q[3];
sx q[3];
rz(-1.2609153) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9481908) q[0];
sx q[0];
rz(-1.8302487) q[0];
sx q[0];
rz(2.6780658) q[0];
rz(0.14687982) q[1];
sx q[1];
rz(-1.905922) q[1];
sx q[1];
rz(0.99259496) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50336058) q[0];
sx q[0];
rz(-1.322317) q[0];
sx q[0];
rz(0.72297024) q[0];
rz(3.0743344) q[2];
sx q[2];
rz(-1.2853299) q[2];
sx q[2];
rz(-0.87905264) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0907572) q[1];
sx q[1];
rz(-0.91716424) q[1];
sx q[1];
rz(0.98934503) q[1];
x q[2];
rz(-1.1168043) q[3];
sx q[3];
rz(-0.89056784) q[3];
sx q[3];
rz(-0.46453634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2831882) q[2];
sx q[2];
rz(-1.4294727) q[2];
sx q[2];
rz(2.6962213) q[2];
rz(-1.8177659) q[3];
sx q[3];
rz(-2.0711074) q[3];
sx q[3];
rz(2.0557192) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0367947) q[0];
sx q[0];
rz(-3.1322271) q[0];
sx q[0];
rz(2.985756) q[0];
rz(1.2771295) q[1];
sx q[1];
rz(-1.7946449) q[1];
sx q[1];
rz(3.1256622) q[1];
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
rz(1.2066293) q[2];
sx q[2];
rz(-0.81239919) q[2];
sx q[2];
rz(1.3398085) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.6813035) q[1];
sx q[1];
rz(-0.88424129) q[1];
sx q[1];
rz(0.89214561) q[1];
rz(-pi) q[2];
rz(0.63834493) q[3];
sx q[3];
rz(-1.9192291) q[3];
sx q[3];
rz(-1.750647) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.5126123) q[2];
sx q[2];
rz(-0.11233687) q[2];
sx q[2];
rz(-0.12574276) q[2];
rz(-0.9564774) q[3];
sx q[3];
rz(-1.4230909) q[3];
sx q[3];
rz(-2.8023348) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
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
rz(-2.1379037) q[1];
sx q[1];
rz(2.454954) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4545802) q[0];
sx q[0];
rz(-1.079315) q[0];
sx q[0];
rz(-0.67680741) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3977658) q[2];
sx q[2];
rz(-1.4464207) q[2];
sx q[2];
rz(-1.5501407) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.2575779) q[1];
sx q[1];
rz(-0.19022372) q[1];
sx q[1];
rz(-0.68489267) q[1];
rz(-pi) q[2];
rz(-1.500528) q[3];
sx q[3];
rz(-1.3965551) q[3];
sx q[3];
rz(-1.4605923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.3829284) q[2];
sx q[2];
rz(-2.0186581) q[2];
sx q[2];
rz(-2.0056966) q[2];
rz(-2.1618333) q[3];
sx q[3];
rz(-1.3330678) q[3];
sx q[3];
rz(2.4433344) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
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
rz(-1.2084992) q[1];
sx q[1];
rz(-0.6932238) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1826349) q[0];
sx q[0];
rz(-1.4719677) q[0];
sx q[0];
rz(-1.1846428) q[0];
rz(-pi) q[1];
rz(-2.3568527) q[2];
sx q[2];
rz(-2.5973136) q[2];
sx q[2];
rz(1.2960824) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.7267513) q[1];
sx q[1];
rz(-1.5637239) q[1];
sx q[1];
rz(-1.4727791) q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[1];
rz(-2.2263055) q[2];
sx q[2];
rz(-1.8733571) q[2];
sx q[2];
rz(-2.477296) q[2];
rz(-1.9281467) q[3];
sx q[3];
rz(-2.4197141) q[3];
sx q[3];
rz(0.87219316) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.450347) q[0];
sx q[0];
rz(-2.3007614) q[0];
sx q[0];
rz(1.7578516) q[0];
rz(-1.6830403) q[1];
sx q[1];
rz(-0.28266193) q[1];
sx q[1];
rz(-2.2255486) q[1];
rz(-0.37049313) q[2];
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
