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
rz(-0.64033163) q[0];
sx q[0];
rz(5.3618199) q[0];
sx q[0];
rz(11.034135) q[0];
rz(0.20707239) q[1];
sx q[1];
rz(4.3011811) q[1];
sx q[1];
rz(11.094697) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6651513) q[0];
sx q[0];
rz(-1.8855321) q[0];
sx q[0];
rz(1.2405618) q[0];
rz(0.81469131) q[2];
sx q[2];
rz(-2.1213795) q[2];
sx q[2];
rz(1.6901291) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.0589367) q[1];
sx q[1];
rz(-1.3322222) q[1];
sx q[1];
rz(-3.039837) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.26443414) q[3];
sx q[3];
rz(-0.61226058) q[3];
sx q[3];
rz(-2.0443288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.6472935) q[2];
sx q[2];
rz(-2.1493201) q[2];
sx q[2];
rz(-2.5956019) q[2];
rz(-2.2030988) q[3];
sx q[3];
rz(-0.5135082) q[3];
sx q[3];
rz(-2.688496) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.12478011) q[0];
sx q[0];
rz(-2.0240968) q[0];
sx q[0];
rz(2.5002531) q[0];
rz(-1.9524139) q[1];
sx q[1];
rz(-0.42354217) q[1];
sx q[1];
rz(2.0372527) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76726145) q[0];
sx q[0];
rz(-2.4924954) q[0];
sx q[0];
rz(-0.24477203) q[0];
rz(-pi) q[1];
rz(2.8139427) q[2];
sx q[2];
rz(-1.0585241) q[2];
sx q[2];
rz(-1.22359) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.0094665388) q[1];
sx q[1];
rz(-2.3915572) q[1];
sx q[1];
rz(-0.86923952) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6112174) q[3];
sx q[3];
rz(-1.4835525) q[3];
sx q[3];
rz(-0.62007346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.39893338) q[2];
sx q[2];
rz(-2.0519966) q[2];
sx q[2];
rz(2.147414) q[2];
rz(1.582877) q[3];
sx q[3];
rz(-2.4623058) q[3];
sx q[3];
rz(-2.5905632) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1347374) q[0];
sx q[0];
rz(-2.0515433) q[0];
sx q[0];
rz(-2.5110974) q[0];
rz(-0.14671239) q[1];
sx q[1];
rz(-1.881733) q[1];
sx q[1];
rz(-0.86348081) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1711581) q[0];
sx q[0];
rz(-0.92464952) q[0];
sx q[0];
rz(-0.37982725) q[0];
rz(-0.93117945) q[2];
sx q[2];
rz(-0.34602994) q[2];
sx q[2];
rz(-2.211386) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.32541986) q[1];
sx q[1];
rz(-0.4494986) q[1];
sx q[1];
rz(-1.3260496) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9894838) q[3];
sx q[3];
rz(-2.8242691) q[3];
sx q[3];
rz(0.0025686669) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.33519393) q[2];
sx q[2];
rz(-1.6890182) q[2];
sx q[2];
rz(2.2550968) q[2];
rz(0.42996201) q[3];
sx q[3];
rz(-2.6468247) q[3];
sx q[3];
rz(0.89216843) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-0.48915136) q[0];
sx q[0];
rz(-0.2944856) q[0];
sx q[0];
rz(-0.44892204) q[0];
rz(2.4849675) q[1];
sx q[1];
rz(-1.7078327) q[1];
sx q[1];
rz(-2.8112559) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9245488) q[0];
sx q[0];
rz(-1.7946436) q[0];
sx q[0];
rz(1.8987109) q[0];
rz(1.0895715) q[2];
sx q[2];
rz(-0.46483609) q[2];
sx q[2];
rz(-0.46457738) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.8790414) q[1];
sx q[1];
rz(-1.8617668) q[1];
sx q[1];
rz(-2.0249184) q[1];
rz(-pi) q[2];
rz(-1.4497595) q[3];
sx q[3];
rz(-1.5797628) q[3];
sx q[3];
rz(-3.1154273) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.0535447) q[2];
sx q[2];
rz(-1.6036754) q[2];
sx q[2];
rz(-2.9223082) q[2];
rz(-0.80026475) q[3];
sx q[3];
rz(-2.181668) q[3];
sx q[3];
rz(0.76103359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4735755) q[0];
sx q[0];
rz(-1.1290843) q[0];
sx q[0];
rz(-1.6582723) q[0];
rz(-0.56770101) q[1];
sx q[1];
rz(-1.2115819) q[1];
sx q[1];
rz(1.4989046) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.606584) q[0];
sx q[0];
rz(-2.5902915) q[0];
sx q[0];
rz(-2.8583741) q[0];
x q[1];
rz(1.5818059) q[2];
sx q[2];
rz(-1.4729452) q[2];
sx q[2];
rz(-0.40153533) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.066047) q[1];
sx q[1];
rz(-0.64295628) q[1];
sx q[1];
rz(0.33400225) q[1];
rz(-pi) q[2];
rz(-2.1338283) q[3];
sx q[3];
rz(-0.22528409) q[3];
sx q[3];
rz(2.99175) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.5119001) q[2];
sx q[2];
rz(-2.3834159) q[2];
sx q[2];
rz(-2.8153815) q[2];
rz(-0.1263667) q[3];
sx q[3];
rz(-0.56505239) q[3];
sx q[3];
rz(0.27157426) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7997953) q[0];
sx q[0];
rz(-1.837715) q[0];
sx q[0];
rz(-2.6190992) q[0];
rz(2.3937468) q[1];
sx q[1];
rz(-1.5380842) q[1];
sx q[1];
rz(-1.1572256) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9686988) q[0];
sx q[0];
rz(-1.1215498) q[0];
sx q[0];
rz(-0.32562738) q[0];
x q[1];
rz(1.1076449) q[2];
sx q[2];
rz(-1.2105296) q[2];
sx q[2];
rz(-0.86371326) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.5392235) q[1];
sx q[1];
rz(-1.4280978) q[1];
sx q[1];
rz(3.0485832) q[1];
rz(-pi) q[2];
rz(0.22208235) q[3];
sx q[3];
rz(-0.33470585) q[3];
sx q[3];
rz(0.10192733) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(3.130827) q[2];
sx q[2];
rz(-0.5126493) q[2];
sx q[2];
rz(1.958468) q[2];
rz(3.0986687) q[3];
sx q[3];
rz(-1.5019417) q[3];
sx q[3];
rz(-2.8973798) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
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
rz(-0.44508988) q[0];
sx q[0];
rz(-0.88560167) q[0];
sx q[0];
rz(-2.4687299) q[0];
rz(-0.54975763) q[1];
sx q[1];
rz(-1.2943228) q[1];
sx q[1];
rz(-1.4901644) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25355089) q[0];
sx q[0];
rz(-1.3924989) q[0];
sx q[0];
rz(-2.5744497) q[0];
rz(-pi) q[1];
rz(-2.6308833) q[2];
sx q[2];
rz(-2.0119397) q[2];
sx q[2];
rz(3.0808251) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0615599) q[1];
sx q[1];
rz(-1.6668238) q[1];
sx q[1];
rz(-0.43635861) q[1];
rz(-pi) q[2];
x q[2];
rz(3.1079632) q[3];
sx q[3];
rz(-1.6190268) q[3];
sx q[3];
rz(-1.5329211) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0979536) q[2];
sx q[2];
rz(-2.592228) q[2];
sx q[2];
rz(1.9926386) q[2];
rz(0.39135459) q[3];
sx q[3];
rz(-1.9214182) q[3];
sx q[3];
rz(-2.2933188) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5452165) q[0];
sx q[0];
rz(-1.2059809) q[0];
sx q[0];
rz(2.4429831) q[0];
rz(-2.2907603) q[1];
sx q[1];
rz(-0.60092503) q[1];
sx q[1];
rz(-0.94815475) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8671001) q[0];
sx q[0];
rz(-0.72625181) q[0];
sx q[0];
rz(-3.0175643) q[0];
rz(-pi) q[1];
rz(1.9004563) q[2];
sx q[2];
rz(-1.621581) q[2];
sx q[2];
rz(0.18731681) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.52719927) q[1];
sx q[1];
rz(-2.4233257) q[1];
sx q[1];
rz(2.2092256) q[1];
rz(-0.83017577) q[3];
sx q[3];
rz(-0.98168938) q[3];
sx q[3];
rz(-2.0996706) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.4725388) q[2];
sx q[2];
rz(-1.6927745) q[2];
sx q[2];
rz(1.4581468) q[2];
rz(-1.091188) q[3];
sx q[3];
rz(-2.2444221) q[3];
sx q[3];
rz(0.18479656) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3882196) q[0];
sx q[0];
rz(-2.9599157) q[0];
sx q[0];
rz(-1.8411807) q[0];
rz(2.6822207) q[1];
sx q[1];
rz(-1.6875024) q[1];
sx q[1];
rz(-0.5836817) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3193286) q[0];
sx q[0];
rz(-2.2177025) q[0];
sx q[0];
rz(-0.95563332) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5351717) q[2];
sx q[2];
rz(-1.6093264) q[2];
sx q[2];
rz(-2.5690479) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.65185279) q[1];
sx q[1];
rz(-2.0803343) q[1];
sx q[1];
rz(-0.40470985) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0784057) q[3];
sx q[3];
rz(-2.655683) q[3];
sx q[3];
rz(2.8775281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.5251081) q[2];
sx q[2];
rz(-0.27518347) q[2];
sx q[2];
rz(-0.52213651) q[2];
rz(0.76256847) q[3];
sx q[3];
rz(-1.4985761) q[3];
sx q[3];
rz(-2.7975119) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6793215) q[0];
sx q[0];
rz(-1.9336047) q[0];
sx q[0];
rz(-0.59360498) q[0];
rz(-0.69315928) q[1];
sx q[1];
rz(-0.79265541) q[1];
sx q[1];
rz(2.3902182) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2776141) q[0];
sx q[0];
rz(-1.1728012) q[0];
sx q[0];
rz(-1.2897183) q[0];
rz(-0.88971207) q[2];
sx q[2];
rz(-1.0719704) q[2];
sx q[2];
rz(-0.46173438) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.46111108) q[1];
sx q[1];
rz(-2.2348677) q[1];
sx q[1];
rz(1.4898438) q[1];
x q[2];
rz(0.16887197) q[3];
sx q[3];
rz(-1.6467301) q[3];
sx q[3];
rz(0.31496615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.0005325) q[2];
sx q[2];
rz(-0.70912051) q[2];
sx q[2];
rz(0.023690311) q[2];
rz(-1.116811) q[3];
sx q[3];
rz(-1.4477372) q[3];
sx q[3];
rz(-2.007133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
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
rz(2.5524207) q[0];
sx q[0];
rz(-1.2180653) q[0];
sx q[0];
rz(1.1267452) q[0];
rz(-1.3202271) q[1];
sx q[1];
rz(-1.2757433) q[1];
sx q[1];
rz(-2.129) q[1];
rz(2.597328) q[2];
sx q[2];
rz(-1.3514634) q[2];
sx q[2];
rz(0.8018871) q[2];
rz(-1.3397401) q[3];
sx q[3];
rz(-1.7920296) q[3];
sx q[3];
rz(-1.5248004) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
