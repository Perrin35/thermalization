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
rz(1.3820833) q[0];
sx q[0];
rz(-0.74892646) q[0];
sx q[0];
rz(1.7480667) q[0];
rz(0.8028318) q[1];
sx q[1];
rz(-0.79610151) q[1];
sx q[1];
rz(-2.610745) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7323468) q[0];
sx q[0];
rz(-2.0589925) q[0];
sx q[0];
rz(-0.40803473) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5535001) q[2];
sx q[2];
rz(-1.7375542) q[2];
sx q[2];
rz(1.1823428) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.2474607) q[1];
sx q[1];
rz(-1.0399721) q[1];
sx q[1];
rz(3.0546419) q[1];
rz(1.7479806) q[3];
sx q[3];
rz(-0.52588829) q[3];
sx q[3];
rz(2.2729659) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.73534766) q[2];
sx q[2];
rz(-2.2082059) q[2];
sx q[2];
rz(-0.99785844) q[2];
rz(-0.88456279) q[3];
sx q[3];
rz(-1.1793143) q[3];
sx q[3];
rz(-0.71200371) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74580055) q[0];
sx q[0];
rz(-2.6048248) q[0];
sx q[0];
rz(-2.0641548) q[0];
rz(-0.61915818) q[1];
sx q[1];
rz(-1.8749323) q[1];
sx q[1];
rz(-1.9177297) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9911988) q[0];
sx q[0];
rz(-1.7454892) q[0];
sx q[0];
rz(1.2758486) q[0];
rz(-pi) q[1];
rz(3.139758) q[2];
sx q[2];
rz(-1.6009838) q[2];
sx q[2];
rz(-1.0890397) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.89586414) q[1];
sx q[1];
rz(-0.99555496) q[1];
sx q[1];
rz(-1.0771345) q[1];
rz(0.72325752) q[3];
sx q[3];
rz(-0.92953909) q[3];
sx q[3];
rz(-2.0632405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.18664843) q[2];
sx q[2];
rz(-2.2961605) q[2];
sx q[2];
rz(2.1050982) q[2];
rz(-1.9519818) q[3];
sx q[3];
rz(-1.9018973) q[3];
sx q[3];
rz(-3.0724683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67967296) q[0];
sx q[0];
rz(-2.9000977) q[0];
sx q[0];
rz(1.4659721) q[0];
rz(-1.8849323) q[1];
sx q[1];
rz(-0.64777056) q[1];
sx q[1];
rz(-0.57674903) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5546122) q[0];
sx q[0];
rz(-2.1471427) q[0];
sx q[0];
rz(-2.2423241) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5336214) q[2];
sx q[2];
rz(-1.1122744) q[2];
sx q[2];
rz(-0.22500817) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-3.0054606) q[1];
sx q[1];
rz(-1.1578655) q[1];
sx q[1];
rz(-0.056483566) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4068094) q[3];
sx q[3];
rz(-1.6583558) q[3];
sx q[3];
rz(1.6370186) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.74167788) q[2];
sx q[2];
rz(-2.8027813) q[2];
sx q[2];
rz(-0.76876387) q[2];
rz(-2.9703043) q[3];
sx q[3];
rz(-1.7007217) q[3];
sx q[3];
rz(0.35458529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4549558) q[0];
sx q[0];
rz(-0.88382116) q[0];
sx q[0];
rz(-0.53632847) q[0];
rz(0.80279154) q[1];
sx q[1];
rz(-1.8575467) q[1];
sx q[1];
rz(0.098043052) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2778846) q[0];
sx q[0];
rz(-1.4914762) q[0];
sx q[0];
rz(-1.2427928) q[0];
rz(-pi) q[1];
rz(-0.78196193) q[2];
sx q[2];
rz(-1.1267917) q[2];
sx q[2];
rz(0.16829106) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.50419549) q[1];
sx q[1];
rz(-1.3679844) q[1];
sx q[1];
rz(-1.0493953) q[1];
x q[2];
rz(0.40880568) q[3];
sx q[3];
rz(-1.173436) q[3];
sx q[3];
rz(-1.8026343) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.1596277) q[2];
sx q[2];
rz(-2.1106796) q[2];
sx q[2];
rz(-0.31943303) q[2];
rz(2.5751513) q[3];
sx q[3];
rz(-0.74685493) q[3];
sx q[3];
rz(1.1613065) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7403858) q[0];
sx q[0];
rz(-1.8219319) q[0];
sx q[0];
rz(2.5886986) q[0];
rz(-2.3623908) q[1];
sx q[1];
rz(-2.3190277) q[1];
sx q[1];
rz(2.353277) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32892431) q[0];
sx q[0];
rz(-1.5018437) q[0];
sx q[0];
rz(2.8806995) q[0];
rz(-pi) q[1];
rz(-1.5469903) q[2];
sx q[2];
rz(-1.362065) q[2];
sx q[2];
rz(0.1114276) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.31438875) q[1];
sx q[1];
rz(-2.6174859) q[1];
sx q[1];
rz(0.57669611) q[1];
x q[2];
rz(-1.3230927) q[3];
sx q[3];
rz(-1.8616397) q[3];
sx q[3];
rz(1.820562) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.7439482) q[2];
sx q[2];
rz(-2.0136191) q[2];
sx q[2];
rz(2.9359342) q[2];
rz(-1.7913943) q[3];
sx q[3];
rz(-0.74068991) q[3];
sx q[3];
rz(-3.1309483) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3288997) q[0];
sx q[0];
rz(-2.6943272) q[0];
sx q[0];
rz(1.6269667) q[0];
rz(2.336592) q[1];
sx q[1];
rz(-1.4470419) q[1];
sx q[1];
rz(-2.8256493) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7132063) q[0];
sx q[0];
rz(-1.9613736) q[0];
sx q[0];
rz(1.2001363) q[0];
rz(-2.5300643) q[2];
sx q[2];
rz(-1.8793162) q[2];
sx q[2];
rz(-2.6055544) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.8711612) q[1];
sx q[1];
rz(-1.3140895) q[1];
sx q[1];
rz(1.8553599) q[1];
rz(-1.7234572) q[3];
sx q[3];
rz(-1.8777285) q[3];
sx q[3];
rz(-1.538601) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.67779764) q[2];
sx q[2];
rz(-2.4424876) q[2];
sx q[2];
rz(-0.15764906) q[2];
rz(-0.53358233) q[3];
sx q[3];
rz(-1.491051) q[3];
sx q[3];
rz(-1.6238448) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
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
rz(-0.59876281) q[0];
sx q[0];
rz(-2.6451126) q[0];
sx q[0];
rz(-0.98094034) q[0];
rz(-0.44278231) q[1];
sx q[1];
rz(-1.9552224) q[1];
sx q[1];
rz(-0.47169366) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1991292) q[0];
sx q[0];
rz(-0.53265306) q[0];
sx q[0];
rz(-0.96766242) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6006655) q[2];
sx q[2];
rz(-2.8188883) q[2];
sx q[2];
rz(2.2216715) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.065475796) q[1];
sx q[1];
rz(-0.13607506) q[1];
sx q[1];
rz(0.10175609) q[1];
rz(2.8011462) q[3];
sx q[3];
rz(-1.7448815) q[3];
sx q[3];
rz(0.3917087) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.9337351) q[2];
sx q[2];
rz(-2.2722878) q[2];
sx q[2];
rz(-0.11824879) q[2];
rz(-0.56784981) q[3];
sx q[3];
rz(-1.3563145) q[3];
sx q[3];
rz(0.80287272) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0826913) q[0];
sx q[0];
rz(-1.4014129) q[0];
sx q[0];
rz(-2.4543104) q[0];
rz(1.2673238) q[1];
sx q[1];
rz(-1.9272389) q[1];
sx q[1];
rz(0.6257239) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95803496) q[0];
sx q[0];
rz(-1.5661217) q[0];
sx q[0];
rz(1.5497394) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2801493) q[2];
sx q[2];
rz(-2.3683135) q[2];
sx q[2];
rz(2.803162) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.76624745) q[1];
sx q[1];
rz(-1.5583923) q[1];
sx q[1];
rz(0.9961025) q[1];
rz(-2.1504907) q[3];
sx q[3];
rz(-1.3091581) q[3];
sx q[3];
rz(-1.2708479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.22831336) q[2];
sx q[2];
rz(-1.8931171) q[2];
sx q[2];
rz(2.0466364) q[2];
rz(-2.6878808) q[3];
sx q[3];
rz(-2.3504421) q[3];
sx q[3];
rz(-2.9554844) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5472645) q[0];
sx q[0];
rz(-0.493395) q[0];
sx q[0];
rz(-0.067597978) q[0];
rz(-1.0293845) q[1];
sx q[1];
rz(-1.8600347) q[1];
sx q[1];
rz(0.13551113) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.23897509) q[0];
sx q[0];
rz(-1.0229361) q[0];
sx q[0];
rz(2.2051714) q[0];
x q[1];
rz(-2.4004647) q[2];
sx q[2];
rz(-1.4504045) q[2];
sx q[2];
rz(-0.84726221) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.77733946) q[1];
sx q[1];
rz(-2.4852264) q[1];
sx q[1];
rz(-0.47231625) q[1];
x q[2];
rz(-1.816808) q[3];
sx q[3];
rz(-1.8455077) q[3];
sx q[3];
rz(-1.6149855) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.4542666) q[2];
sx q[2];
rz(-2.8886075) q[2];
sx q[2];
rz(-0.48736408) q[2];
rz(0.38960114) q[3];
sx q[3];
rz(-0.95249683) q[3];
sx q[3];
rz(0.065936955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49366632) q[0];
sx q[0];
rz(-1.0856029) q[0];
sx q[0];
rz(-1.7568463) q[0];
rz(2.0866277) q[1];
sx q[1];
rz(-2.2672548) q[1];
sx q[1];
rz(-0.25996444) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6523165) q[0];
sx q[0];
rz(-2.2201009) q[0];
sx q[0];
rz(0.32353521) q[0];
rz(0.10037701) q[2];
sx q[2];
rz(-1.4401199) q[2];
sx q[2];
rz(1.584077) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.205906) q[1];
sx q[1];
rz(-2.5126713) q[1];
sx q[1];
rz(1.511093) q[1];
rz(-pi) q[2];
rz(2.9197844) q[3];
sx q[3];
rz(-1.7647499) q[3];
sx q[3];
rz(1.2916331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.42183033) q[2];
sx q[2];
rz(-2.7037342) q[2];
sx q[2];
rz(-1.3502632) q[2];
rz(-2.0566025) q[3];
sx q[3];
rz(-1.9271873) q[3];
sx q[3];
rz(2.0124281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0709025) q[0];
sx q[0];
rz(-2.4546843) q[0];
sx q[0];
rz(-0.51921459) q[0];
rz(1.4821953) q[1];
sx q[1];
rz(-1.2980325) q[1];
sx q[1];
rz(1.4549805) q[1];
rz(-0.71375511) q[2];
sx q[2];
rz(-0.71520765) q[2];
sx q[2];
rz(2.2463828) q[2];
rz(2.8419513) q[3];
sx q[3];
rz(-1.089923) q[3];
sx q[3];
rz(-3.0724473) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
