OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.20004162) q[0];
sx q[0];
rz(-2.9574432) q[0];
sx q[0];
rz(1.3567691) q[0];
rz(-1.6954724) q[1];
sx q[1];
rz(-1.5741916) q[1];
sx q[1];
rz(2.8653436) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41208273) q[0];
sx q[0];
rz(-0.19664054) q[0];
sx q[0];
rz(1.6565408) q[0];
rz(3.0577781) q[2];
sx q[2];
rz(-3.028125) q[2];
sx q[2];
rz(-0.49926114) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(3.0608643) q[1];
sx q[1];
rz(-1.5866188) q[1];
sx q[1];
rz(3.0982261) q[1];
rz(-pi) q[2];
rz(2.9570691) q[3];
sx q[3];
rz(-0.34466668) q[3];
sx q[3];
rz(2.7027948) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.19286284) q[2];
sx q[2];
rz(-0.36878815) q[2];
sx q[2];
rz(-0.31664872) q[2];
rz(-2.6266895) q[3];
sx q[3];
rz(-0.0080990214) q[3];
sx q[3];
rz(0.81075877) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0136451) q[0];
sx q[0];
rz(-0.34334308) q[0];
sx q[0];
rz(0.031033255) q[0];
rz(-1.4399928) q[1];
sx q[1];
rz(-0.80039918) q[1];
sx q[1];
rz(-1.4061141) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3701551) q[0];
sx q[0];
rz(-1.9906128) q[0];
sx q[0];
rz(0.19924723) q[0];
rz(-3.1028264) q[2];
sx q[2];
rz(-1.3666461) q[2];
sx q[2];
rz(0.37224712) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.5760562) q[1];
sx q[1];
rz(-1.4829552) q[1];
sx q[1];
rz(2.9599543) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4214899) q[3];
sx q[3];
rz(-2.341695) q[3];
sx q[3];
rz(-0.61198583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.42355737) q[2];
sx q[2];
rz(-0.2090629) q[2];
sx q[2];
rz(-0.62530708) q[2];
rz(-1.5710255) q[3];
sx q[3];
rz(-2.8795241) q[3];
sx q[3];
rz(3.1077793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5686947) q[0];
sx q[0];
rz(-2.608572) q[0];
sx q[0];
rz(2.0030588) q[0];
rz(-0.89433995) q[1];
sx q[1];
rz(-3.1411451) q[1];
sx q[1];
rz(-2.0785887) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.759623) q[0];
sx q[0];
rz(-2.2142525) q[0];
sx q[0];
rz(2.8354885) q[0];
rz(-pi) q[1];
rz(1.3787621) q[2];
sx q[2];
rz(-2.3096032) q[2];
sx q[2];
rz(2.2453824) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.56617006) q[1];
sx q[1];
rz(-0.1488007) q[1];
sx q[1];
rz(2.0729561) q[1];
rz(2.9195905) q[3];
sx q[3];
rz(-1.3679805) q[3];
sx q[3];
rz(-2.0012282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.87463093) q[2];
sx q[2];
rz(-2.3571372) q[2];
sx q[2];
rz(2.0060914) q[2];
rz(-1.1442432) q[3];
sx q[3];
rz(-3.0502697) q[3];
sx q[3];
rz(1.7448366) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8046232) q[0];
sx q[0];
rz(-2.7380044) q[0];
sx q[0];
rz(-1.1176156) q[0];
rz(3.0506253) q[1];
sx q[1];
rz(-0.03799835) q[1];
sx q[1];
rz(1.631558) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9695797) q[0];
sx q[0];
rz(-1.0019127) q[0];
sx q[0];
rz(-0.24636951) q[0];
x q[1];
rz(0.20524673) q[2];
sx q[2];
rz(-1.6778429) q[2];
sx q[2];
rz(2.7595124) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.7854033) q[1];
sx q[1];
rz(-1.8757687) q[1];
sx q[1];
rz(3.0126291) q[1];
rz(-2.4337188) q[3];
sx q[3];
rz(-2.3972569) q[3];
sx q[3];
rz(-0.78619781) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.9522004) q[2];
sx q[2];
rz(-0.11710937) q[2];
sx q[2];
rz(1.2516775) q[2];
rz(-0.39024726) q[3];
sx q[3];
rz(-0.016642112) q[3];
sx q[3];
rz(0.19002953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8613043) q[0];
sx q[0];
rz(-2.5087293) q[0];
sx q[0];
rz(-1.8508152) q[0];
rz(3.1351008) q[1];
sx q[1];
rz(-0.26632729) q[1];
sx q[1];
rz(1.806462) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4075942) q[0];
sx q[0];
rz(-1.7069792) q[0];
sx q[0];
rz(-2.7969317) q[0];
rz(-0.13258055) q[2];
sx q[2];
rz(-2.6040316) q[2];
sx q[2];
rz(0.81576306) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.3975828) q[1];
sx q[1];
rz(-2.8967794) q[1];
sx q[1];
rz(-0.53298612) q[1];
rz(-2.9560821) q[3];
sx q[3];
rz(-1.9855567) q[3];
sx q[3];
rz(-1.1832373) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.409965) q[2];
sx q[2];
rz(-0.16606398) q[2];
sx q[2];
rz(0.49327332) q[2];
rz(-1.0924245) q[3];
sx q[3];
rz(-3.137393) q[3];
sx q[3];
rz(2.5815729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6888371) q[0];
sx q[0];
rz(-1.9552564) q[0];
sx q[0];
rz(1.8) q[0];
rz(-1.5635368) q[1];
sx q[1];
rz(-2.9069803) q[1];
sx q[1];
rz(0.87306195) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.045220203) q[0];
sx q[0];
rz(-2.0866418) q[0];
sx q[0];
rz(1.7302897) q[0];
rz(3.0186924) q[2];
sx q[2];
rz(-2.6989938) q[2];
sx q[2];
rz(1.0198319) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.4640518) q[1];
sx q[1];
rz(-1.4960491) q[1];
sx q[1];
rz(-1.6932373) q[1];
rz(0.88718702) q[3];
sx q[3];
rz(-0.33970133) q[3];
sx q[3];
rz(-1.8357133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.0688429) q[2];
sx q[2];
rz(-1.9708956) q[2];
sx q[2];
rz(-1.076131) q[2];
rz(0.38232803) q[3];
sx q[3];
rz(-0.040150661) q[3];
sx q[3];
rz(0.7676777) q[3];
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
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.024993984) q[0];
sx q[0];
rz(-3.0285663) q[0];
sx q[0];
rz(2.9207927) q[0];
rz(-2.2258672) q[1];
sx q[1];
rz(-2.974739) q[1];
sx q[1];
rz(-1.6545666) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.090335719) q[0];
sx q[0];
rz(-1.1658586) q[0];
sx q[0];
rz(-1.1722408) q[0];
rz(-pi) q[1];
rz(-0.19796298) q[2];
sx q[2];
rz(-1.930365) q[2];
sx q[2];
rz(-0.61164226) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.4055192) q[1];
sx q[1];
rz(-0.87784518) q[1];
sx q[1];
rz(0.939721) q[1];
rz(-pi) q[2];
rz(-1.7127895) q[3];
sx q[3];
rz(-1.6902704) q[3];
sx q[3];
rz(2.9657397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.9080092) q[2];
sx q[2];
rz(-3.0685232) q[2];
sx q[2];
rz(1.118534) q[2];
rz(1.3633049) q[3];
sx q[3];
rz(-0.076863591) q[3];
sx q[3];
rz(-1.3479056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(1.3240647) q[0];
sx q[0];
rz(-0.20453608) q[0];
sx q[0];
rz(-1.8472141) q[0];
rz(1.4140465) q[1];
sx q[1];
rz(-0.050948016) q[1];
sx q[1];
rz(2.9185413) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4572499) q[0];
sx q[0];
rz(-0.97829223) q[0];
sx q[0];
rz(0.98513435) q[0];
rz(-pi) q[1];
rz(-2.9518218) q[2];
sx q[2];
rz(-1.4354354) q[2];
sx q[2];
rz(3.1369097) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.9893045) q[1];
sx q[1];
rz(-0.76123255) q[1];
sx q[1];
rz(-2.1851995) q[1];
x q[2];
rz(-3.024142) q[3];
sx q[3];
rz(-2.2625105) q[3];
sx q[3];
rz(2.1827616) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.53905067) q[2];
sx q[2];
rz(-1.9571303) q[2];
sx q[2];
rz(2.7089673) q[2];
rz(2.1116347) q[3];
sx q[3];
rz(-0.10051388) q[3];
sx q[3];
rz(2.4242134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3887727) q[0];
sx q[0];
rz(-2.6322375) q[0];
sx q[0];
rz(-0.21554047) q[0];
rz(1.132025) q[1];
sx q[1];
rz(-0.80765453) q[1];
sx q[1];
rz(1.4707123) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4819538) q[0];
sx q[0];
rz(-0.99876548) q[0];
sx q[0];
rz(-2.7670949) q[0];
rz(-pi) q[1];
rz(-2.7637761) q[2];
sx q[2];
rz(-1.9799798) q[2];
sx q[2];
rz(0.62407392) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(3.1241144) q[1];
sx q[1];
rz(-0.062659293) q[1];
sx q[1];
rz(-2.6397328) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6588259) q[3];
sx q[3];
rz(-1.996187) q[3];
sx q[3];
rz(1.9874043) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.97951621) q[2];
sx q[2];
rz(-3.1023878) q[2];
sx q[2];
rz(2.6230679) q[2];
rz(-2.3833852) q[3];
sx q[3];
rz(-2.3617305) q[3];
sx q[3];
rz(1.6583748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46490797) q[0];
sx q[0];
rz(-1.9187036) q[0];
sx q[0];
rz(1.884888) q[0];
rz(2.8769809) q[1];
sx q[1];
rz(-0.0022609641) q[1];
sx q[1];
rz(1.1983926) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0007858) q[0];
sx q[0];
rz(-0.0065751271) q[0];
sx q[0];
rz(3.0964609) q[0];
rz(1.9047059) q[2];
sx q[2];
rz(-1.5057202) q[2];
sx q[2];
rz(-1.2953616) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.8492412) q[1];
sx q[1];
rz(-1.8051984) q[1];
sx q[1];
rz(1.5567907) q[1];
rz(-1.7890187) q[3];
sx q[3];
rz(-0.17328158) q[3];
sx q[3];
rz(0.97001433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.4635224) q[2];
sx q[2];
rz(-0.0040201298) q[2];
sx q[2];
rz(-0.013997495) q[2];
rz(-2.8409581) q[3];
sx q[3];
rz(-0.0036573452) q[3];
sx q[3];
rz(1.5716871) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32776253) q[0];
sx q[0];
rz(-1.4284842) q[0];
sx q[0];
rz(1.6416657) q[0];
rz(0.035540237) q[1];
sx q[1];
rz(-0.4353558) q[1];
sx q[1];
rz(0.33886649) q[1];
rz(0.28066228) q[2];
sx q[2];
rz(-2.4246695) q[2];
sx q[2];
rz(-0.05447745) q[2];
rz(-2.9406459) q[3];
sx q[3];
rz(-2.3280716) q[3];
sx q[3];
rz(0.18542326) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
