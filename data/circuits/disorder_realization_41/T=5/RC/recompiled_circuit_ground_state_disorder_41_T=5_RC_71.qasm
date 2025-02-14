OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.5557033) q[0];
sx q[0];
rz(4.3309431) q[0];
sx q[0];
rz(10.588607) q[0];
rz(2.3948506) q[1];
sx q[1];
rz(-2.876694) q[1];
sx q[1];
rz(2.9704111) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25397444) q[0];
sx q[0];
rz(-1.3518392) q[0];
sx q[0];
rz(-2.9848214) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.77734485) q[2];
sx q[2];
rz(-1.0413578) q[2];
sx q[2];
rz(2.6482761) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.2210126) q[1];
sx q[1];
rz(-1.5713931) q[1];
sx q[1];
rz(-1.8662054) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.1247507) q[3];
sx q[3];
rz(-1.5762268) q[3];
sx q[3];
rz(-0.098244103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(3.050194) q[2];
sx q[2];
rz(-1.564743) q[2];
sx q[2];
rz(2.0719299) q[2];
rz(-1.4403249) q[3];
sx q[3];
rz(-1.0245208) q[3];
sx q[3];
rz(1.9694156) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1759724) q[0];
sx q[0];
rz(-1.6930641) q[0];
sx q[0];
rz(-2.5536221) q[0];
rz(-1.2929471) q[1];
sx q[1];
rz(-2.1473532) q[1];
sx q[1];
rz(2.9867244) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17093824) q[0];
sx q[0];
rz(-1.0736671) q[0];
sx q[0];
rz(1.4488104) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9005342) q[2];
sx q[2];
rz(-1.7130642) q[2];
sx q[2];
rz(0.011474284) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.57666393) q[1];
sx q[1];
rz(-2.0368848) q[1];
sx q[1];
rz(-2.7647777) q[1];
rz(-2.9294176) q[3];
sx q[3];
rz(-1.7042233) q[3];
sx q[3];
rz(1.9047036) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.64314848) q[2];
sx q[2];
rz(-2.4109106) q[2];
sx q[2];
rz(-3.0291962) q[2];
rz(2.3183838) q[3];
sx q[3];
rz(-1.9197437) q[3];
sx q[3];
rz(0.52197758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4560029) q[0];
sx q[0];
rz(-1.0887479) q[0];
sx q[0];
rz(0.96813694) q[0];
rz(-2.0227506) q[1];
sx q[1];
rz(-1.5348996) q[1];
sx q[1];
rz(-1.8082089) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4844167) q[0];
sx q[0];
rz(-1.4098995) q[0];
sx q[0];
rz(-0.96477525) q[0];
x q[1];
rz(2.5381375) q[2];
sx q[2];
rz(-1.9656745) q[2];
sx q[2];
rz(-1.4970571) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3231877) q[1];
sx q[1];
rz(-1.0099704) q[1];
sx q[1];
rz(-0.11769275) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9288238) q[3];
sx q[3];
rz(-1.4910526) q[3];
sx q[3];
rz(2.2984576) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1916888) q[2];
sx q[2];
rz(-2.1046941) q[2];
sx q[2];
rz(2.8766768) q[2];
rz(-2.8252937) q[3];
sx q[3];
rz(-0.33360544) q[3];
sx q[3];
rz(1.3914289) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.815627) q[0];
sx q[0];
rz(-2.9209324) q[0];
sx q[0];
rz(-1.4937481) q[0];
rz(1.4119459) q[1];
sx q[1];
rz(-1.0271881) q[1];
sx q[1];
rz(0.6573917) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18497224) q[0];
sx q[0];
rz(-1.4669384) q[0];
sx q[0];
rz(1.64479) q[0];
rz(-pi) q[1];
rz(0.78423402) q[2];
sx q[2];
rz(-2.0080697) q[2];
sx q[2];
rz(2.0511049) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.3616391) q[1];
sx q[1];
rz(-1.7657537) q[1];
sx q[1];
rz(-0.037752823) q[1];
rz(-pi) q[2];
rz(2.3393624) q[3];
sx q[3];
rz(-2.1539237) q[3];
sx q[3];
rz(1.0232353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.21939453) q[2];
sx q[2];
rz(-1.9191091) q[2];
sx q[2];
rz(-2.738319) q[2];
rz(-1.203631) q[3];
sx q[3];
rz(-1.5091242) q[3];
sx q[3];
rz(-0.50886124) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4598684) q[0];
sx q[0];
rz(-0.42400703) q[0];
sx q[0];
rz(-2.5176609) q[0];
rz(2.4195747) q[1];
sx q[1];
rz(-1.7941509) q[1];
sx q[1];
rz(3.0375979) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34891221) q[0];
sx q[0];
rz(-0.70250547) q[0];
sx q[0];
rz(0.2200495) q[0];
rz(-pi) q[1];
x q[1];
rz(0.35997593) q[2];
sx q[2];
rz(-1.896654) q[2];
sx q[2];
rz(0.47285541) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6045058) q[1];
sx q[1];
rz(-2.1963122) q[1];
sx q[1];
rz(-2.4948859) q[1];
x q[2];
rz(-2.4173043) q[3];
sx q[3];
rz(-0.87832993) q[3];
sx q[3];
rz(-2.9812733) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.35088745) q[2];
sx q[2];
rz(-2.4217589) q[2];
sx q[2];
rz(-0.78901115) q[2];
rz(-1.2645432) q[3];
sx q[3];
rz(-2.0980554) q[3];
sx q[3];
rz(-2.1544382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8378545) q[0];
sx q[0];
rz(-0.721295) q[0];
sx q[0];
rz(-0.20027941) q[0];
rz(-1.6750977) q[1];
sx q[1];
rz(-1.8148345) q[1];
sx q[1];
rz(1.5178348) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87720931) q[0];
sx q[0];
rz(-1.0671927) q[0];
sx q[0];
rz(-1.992454) q[0];
rz(-3.0110961) q[2];
sx q[2];
rz(-2.0708212) q[2];
sx q[2];
rz(-0.44483003) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.955217) q[1];
sx q[1];
rz(-2.5506335) q[1];
sx q[1];
rz(-0.0062828961) q[1];
rz(-pi) q[2];
rz(-1.8143043) q[3];
sx q[3];
rz(-0.78726913) q[3];
sx q[3];
rz(0.7574429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.5184861) q[2];
sx q[2];
rz(-0.97272626) q[2];
sx q[2];
rz(2.8002807) q[2];
rz(2.2845279) q[3];
sx q[3];
rz(-1.2271481) q[3];
sx q[3];
rz(2.6763693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0905867) q[0];
sx q[0];
rz(-1.0336646) q[0];
sx q[0];
rz(1.8940014) q[0];
rz(3.0233851) q[1];
sx q[1];
rz(-1.8659614) q[1];
sx q[1];
rz(0.035621312) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44862952) q[0];
sx q[0];
rz(-1.5482731) q[0];
sx q[0];
rz(1.6883786) q[0];
rz(-0.18892498) q[2];
sx q[2];
rz(-1.1729203) q[2];
sx q[2];
rz(-0.16824761) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.9034152) q[1];
sx q[1];
rz(-2.1768355) q[1];
sx q[1];
rz(1.2708962) q[1];
rz(-pi) q[2];
x q[2];
rz(0.12322085) q[3];
sx q[3];
rz(-2.356592) q[3];
sx q[3];
rz(0.30320019) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.9901765) q[2];
sx q[2];
rz(-1.4750865) q[2];
sx q[2];
rz(-2.0659633) q[2];
rz(2.6208124) q[3];
sx q[3];
rz(-2.5213089) q[3];
sx q[3];
rz(-1.5889997) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2369279) q[0];
sx q[0];
rz(-2.1187145) q[0];
sx q[0];
rz(2.2383595) q[0];
rz(1.5029933) q[1];
sx q[1];
rz(-1.6939949) q[1];
sx q[1];
rz(-1.2618056) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55106581) q[0];
sx q[0];
rz(-1.8000949) q[0];
sx q[0];
rz(0.14044827) q[0];
rz(-pi) q[1];
rz(-2.4440358) q[2];
sx q[2];
rz(-0.64854014) q[2];
sx q[2];
rz(-1.7096138) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.8846036) q[1];
sx q[1];
rz(-1.6620364) q[1];
sx q[1];
rz(-1.6216519) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4724305) q[3];
sx q[3];
rz(-0.72297308) q[3];
sx q[3];
rz(-2.1694136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5624076) q[2];
sx q[2];
rz(-0.6520485) q[2];
sx q[2];
rz(-0.11422608) q[2];
rz(1.3769897) q[3];
sx q[3];
rz(-1.1403133) q[3];
sx q[3];
rz(0.50447869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1219516) q[0];
sx q[0];
rz(-2.1421102) q[0];
sx q[0];
rz(1.9695388) q[0];
rz(1.7578341) q[1];
sx q[1];
rz(-1.9969321) q[1];
sx q[1];
rz(-3.0269472) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8269236) q[0];
sx q[0];
rz(-2.0631587) q[0];
sx q[0];
rz(-0.63757245) q[0];
rz(-pi) q[1];
rz(-0.17342148) q[2];
sx q[2];
rz(-2.6422814) q[2];
sx q[2];
rz(2.2713646) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.23073611) q[1];
sx q[1];
rz(-1.641526) q[1];
sx q[1];
rz(-1.3225698) q[1];
rz(-pi) q[2];
rz(3.036384) q[3];
sx q[3];
rz(-0.79920879) q[3];
sx q[3];
rz(1.5463643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.5961479) q[2];
sx q[2];
rz(-0.20716509) q[2];
sx q[2];
rz(2.2607415) q[2];
rz(2.7018069) q[3];
sx q[3];
rz(-1.1652911) q[3];
sx q[3];
rz(-2.0161276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83774829) q[0];
sx q[0];
rz(-1.4864018) q[0];
sx q[0];
rz(-3.0434171) q[0];
rz(-0.57442609) q[1];
sx q[1];
rz(-1.4593294) q[1];
sx q[1];
rz(2.548545) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0837729) q[0];
sx q[0];
rz(-1.2811986) q[0];
sx q[0];
rz(-2.3168091) q[0];
rz(-pi) q[1];
rz(2.2111528) q[2];
sx q[2];
rz(-2.499492) q[2];
sx q[2];
rz(-3.0885027) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.8349466) q[1];
sx q[1];
rz(-0.017408522) q[1];
sx q[1];
rz(-0.86916344) q[1];
rz(-pi) q[2];
rz(1.9022868) q[3];
sx q[3];
rz(-0.96486366) q[3];
sx q[3];
rz(0.87406033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.9888931) q[2];
sx q[2];
rz(-1.9579192) q[2];
sx q[2];
rz(-2.2031671) q[2];
rz(-1.7484131) q[3];
sx q[3];
rz(-1.8311071) q[3];
sx q[3];
rz(-2.9032629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4599081) q[0];
sx q[0];
rz(-0.60612283) q[0];
sx q[0];
rz(-1.7932307) q[0];
rz(1.0472736) q[1];
sx q[1];
rz(-0.92887639) q[1];
sx q[1];
rz(2.1705719) q[1];
rz(3.0964839) q[2];
sx q[2];
rz(-1.5061629) q[2];
sx q[2];
rz(-0.20284222) q[2];
rz(0.88225928) q[3];
sx q[3];
rz(-0.282448) q[3];
sx q[3];
rz(-0.99568102) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
