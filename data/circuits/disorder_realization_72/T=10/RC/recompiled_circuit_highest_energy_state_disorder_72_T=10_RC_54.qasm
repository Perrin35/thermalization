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
rz(1.3067955) q[0];
sx q[0];
rz(5.4352221) q[0];
sx q[0];
rz(9.3873831) q[0];
rz(-1.0132064) q[1];
sx q[1];
rz(-0.4816882) q[1];
sx q[1];
rz(1.4451292) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2045528) q[0];
sx q[0];
rz(-1.3890966) q[0];
sx q[0];
rz(-2.7054525) q[0];
rz(-pi) q[1];
rz(2.1712473) q[2];
sx q[2];
rz(-3.0171347) q[2];
sx q[2];
rz(2.4336634) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.5132222) q[1];
sx q[1];
rz(-2.0525825) q[1];
sx q[1];
rz(-2.1867153) q[1];
rz(-pi) q[2];
rz(-1.3919034) q[3];
sx q[3];
rz(-0.57235347) q[3];
sx q[3];
rz(1.5327041) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.41363132) q[2];
sx q[2];
rz(-3.0484338) q[2];
sx q[2];
rz(1.5067345) q[2];
rz(2.9480751) q[3];
sx q[3];
rz(-0.7842803) q[3];
sx q[3];
rz(0.94582742) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.042260878) q[0];
sx q[0];
rz(-1.7758545) q[0];
sx q[0];
rz(2.9625764) q[0];
rz(-1.6049113) q[1];
sx q[1];
rz(-0.34114006) q[1];
sx q[1];
rz(1.590033) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59373271) q[0];
sx q[0];
rz(-1.4734125) q[0];
sx q[0];
rz(1.2800526) q[0];
rz(1.8297029) q[2];
sx q[2];
rz(-2.2635411) q[2];
sx q[2];
rz(1.9996114) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.68374709) q[1];
sx q[1];
rz(-0.48626712) q[1];
sx q[1];
rz(0.97944983) q[1];
rz(-2.670774) q[3];
sx q[3];
rz(-1.075346) q[3];
sx q[3];
rz(-1.9268056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4701074) q[2];
sx q[2];
rz(-1.476373) q[2];
sx q[2];
rz(0.022484953) q[2];
rz(-0.89618987) q[3];
sx q[3];
rz(-2.360207) q[3];
sx q[3];
rz(-1.5145068) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(0.80980587) q[0];
sx q[0];
rz(-0.2066732) q[0];
sx q[0];
rz(-1.5326387) q[0];
rz(1.1398075) q[1];
sx q[1];
rz(-2.8228357) q[1];
sx q[1];
rz(0.87361139) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4777268) q[0];
sx q[0];
rz(-0.88318825) q[0];
sx q[0];
rz(2.4640883) q[0];
rz(-1.0192835) q[2];
sx q[2];
rz(-0.9992558) q[2];
sx q[2];
rz(2.4311266) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.4406179) q[1];
sx q[1];
rz(-1.7243885) q[1];
sx q[1];
rz(-2.6540592) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8687137) q[3];
sx q[3];
rz(-1.498184) q[3];
sx q[3];
rz(2.2842294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.5736299) q[2];
sx q[2];
rz(-2.6831388) q[2];
sx q[2];
rz(-2.6652794) q[2];
rz(2.1161946) q[3];
sx q[3];
rz(-1.3566596) q[3];
sx q[3];
rz(2.8477113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3697701) q[0];
sx q[0];
rz(-2.2875468) q[0];
sx q[0];
rz(-2.5452132) q[0];
rz(2.3884933) q[1];
sx q[1];
rz(-1.3558847) q[1];
sx q[1];
rz(-2.7900043) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6839978) q[0];
sx q[0];
rz(-0.54349923) q[0];
sx q[0];
rz(0.25590055) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.48248191) q[2];
sx q[2];
rz(-0.94581476) q[2];
sx q[2];
rz(-3.0925054) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.56573136) q[1];
sx q[1];
rz(-1.3412151) q[1];
sx q[1];
rz(1.1665535) q[1];
x q[2];
rz(-2.5078689) q[3];
sx q[3];
rz(-2.124064) q[3];
sx q[3];
rz(-2.3330757) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.4190462) q[2];
sx q[2];
rz(-1.498035) q[2];
sx q[2];
rz(2.2124115) q[2];
rz(1.2292713) q[3];
sx q[3];
rz(-3.1049187) q[3];
sx q[3];
rz(-1.2693955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(0.96255985) q[0];
sx q[0];
rz(-2.5293009) q[0];
sx q[0];
rz(-2.6137733) q[0];
rz(0.57012308) q[1];
sx q[1];
rz(-0.56990439) q[1];
sx q[1];
rz(0.39400563) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4657426) q[0];
sx q[0];
rz(-0.52211863) q[0];
sx q[0];
rz(0.39710775) q[0];
rz(-pi) q[1];
x q[1];
rz(0.5929596) q[2];
sx q[2];
rz(-1.8793794) q[2];
sx q[2];
rz(-1.9401996) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.2886823) q[1];
sx q[1];
rz(-2.8894419) q[1];
sx q[1];
rz(2.0618477) q[1];
rz(-pi) q[2];
x q[2];
rz(0.72368716) q[3];
sx q[3];
rz(-2.1877648) q[3];
sx q[3];
rz(2.8599515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.7492619) q[2];
sx q[2];
rz(-1.8529842) q[2];
sx q[2];
rz(1.1345081) q[2];
rz(3.115861) q[3];
sx q[3];
rz(-2.0870356) q[3];
sx q[3];
rz(2.499685) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57539097) q[0];
sx q[0];
rz(-1.3302777) q[0];
sx q[0];
rz(0.45912418) q[0];
rz(-0.8210012) q[1];
sx q[1];
rz(-2.2586925) q[1];
sx q[1];
rz(0.51148907) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3200681) q[0];
sx q[0];
rz(-1.2787158) q[0];
sx q[0];
rz(-2.065737) q[0];
rz(-2.8452427) q[2];
sx q[2];
rz(-2.4561693) q[2];
sx q[2];
rz(2.1330331) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.59907179) q[1];
sx q[1];
rz(-0.59172809) q[1];
sx q[1];
rz(-0.27246957) q[1];
x q[2];
rz(-3.0036075) q[3];
sx q[3];
rz(-0.89179865) q[3];
sx q[3];
rz(0.0019055923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.56200999) q[2];
sx q[2];
rz(-1.1082114) q[2];
sx q[2];
rz(2.3902334) q[2];
rz(0.56685081) q[3];
sx q[3];
rz(-1.5375429) q[3];
sx q[3];
rz(-1.8142627) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71996561) q[0];
sx q[0];
rz(-0.1816853) q[0];
sx q[0];
rz(-0.0075465329) q[0];
rz(2.201572) q[1];
sx q[1];
rz(-0.66497856) q[1];
sx q[1];
rz(-0.13253458) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7791393) q[0];
sx q[0];
rz(-0.7374239) q[0];
sx q[0];
rz(2.2064184) q[0];
rz(2.7785702) q[2];
sx q[2];
rz(-1.6953812) q[2];
sx q[2];
rz(1.7302135) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.346437) q[1];
sx q[1];
rz(-2.7043685) q[1];
sx q[1];
rz(-0.45262419) q[1];
rz(-2.9671415) q[3];
sx q[3];
rz(-1.9088591) q[3];
sx q[3];
rz(1.6473339) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.14564766) q[2];
sx q[2];
rz(-1.4952679) q[2];
sx q[2];
rz(0.10124595) q[2];
rz(2.4920987) q[3];
sx q[3];
rz(-2.5443304) q[3];
sx q[3];
rz(-0.35972843) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7262909) q[0];
sx q[0];
rz(-1.8446209) q[0];
sx q[0];
rz(-1.2671965) q[0];
rz(-1.2665117) q[1];
sx q[1];
rz(-0.15788618) q[1];
sx q[1];
rz(-2.3983541) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0948461) q[0];
sx q[0];
rz(-1.0047997) q[0];
sx q[0];
rz(-1.8106242) q[0];
x q[1];
rz(-2.8617992) q[2];
sx q[2];
rz(-0.81088561) q[2];
sx q[2];
rz(1.5992129) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.4488277) q[1];
sx q[1];
rz(-1.8069441) q[1];
sx q[1];
rz(2.9797878) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.69265794) q[3];
sx q[3];
rz(-2.1941278) q[3];
sx q[3];
rz(-2.3828363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.0172952) q[2];
sx q[2];
rz(-2.9399019) q[2];
sx q[2];
rz(-1.7151493) q[2];
rz(2.1258449) q[3];
sx q[3];
rz(-1.4371212) q[3];
sx q[3];
rz(0.82971853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5468686) q[0];
sx q[0];
rz(-2.4877553) q[0];
sx q[0];
rz(-0.015901707) q[0];
rz(-1.5025899) q[1];
sx q[1];
rz(-0.67633164) q[1];
sx q[1];
rz(0.99448386) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5976006) q[0];
sx q[0];
rz(-2.3975537) q[0];
sx q[0];
rz(-0.97980325) q[0];
rz(2.7632668) q[2];
sx q[2];
rz(-1.6294453) q[2];
sx q[2];
rz(-2.891922) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.35895015) q[1];
sx q[1];
rz(-0.94767682) q[1];
sx q[1];
rz(2.5324178) q[1];
rz(-pi) q[2];
rz(-2.0845549) q[3];
sx q[3];
rz(-1.7026859) q[3];
sx q[3];
rz(-1.3472324) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.9859621) q[2];
sx q[2];
rz(-1.4302284) q[2];
sx q[2];
rz(2.963781) q[2];
rz(-1.6832247) q[3];
sx q[3];
rz(-2.7198313) q[3];
sx q[3];
rz(1.873707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47852248) q[0];
sx q[0];
rz(-1.2489742) q[0];
sx q[0];
rz(-1.2836237) q[0];
rz(2.3578857) q[1];
sx q[1];
rz(-2.3678534) q[1];
sx q[1];
rz(-1.8342038) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6790153) q[0];
sx q[0];
rz(-2.0082012) q[0];
sx q[0];
rz(-1.4694197) q[0];
rz(-3.1112017) q[2];
sx q[2];
rz(-1.9585397) q[2];
sx q[2];
rz(0.0645685) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(3.1124005) q[1];
sx q[1];
rz(-1.5123197) q[1];
sx q[1];
rz(-1.7450208) q[1];
rz(-2.7268355) q[3];
sx q[3];
rz(-1.0598174) q[3];
sx q[3];
rz(-1.3172883) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.54138294) q[2];
sx q[2];
rz(-1.4098488) q[2];
sx q[2];
rz(0.48995885) q[2];
rz(-1.1146891) q[3];
sx q[3];
rz(-1.9652818) q[3];
sx q[3];
rz(-2.8118242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1212696) q[0];
sx q[0];
rz(-1.9135495) q[0];
sx q[0];
rz(1.3187153) q[0];
rz(-1.2976788) q[1];
sx q[1];
rz(-0.73328016) q[1];
sx q[1];
rz(2.7136623) q[1];
rz(-0.78133493) q[2];
sx q[2];
rz(-0.35751783) q[2];
sx q[2];
rz(1.9903571) q[2];
rz(-1.8036203) q[3];
sx q[3];
rz(-2.4216087) q[3];
sx q[3];
rz(-2.5523228) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
