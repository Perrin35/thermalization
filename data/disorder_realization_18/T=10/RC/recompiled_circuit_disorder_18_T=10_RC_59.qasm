OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.82233959) q[0];
sx q[0];
rz(5.3263721) q[0];
sx q[0];
rz(10.858067) q[0];
rz(-0.23694555) q[1];
sx q[1];
rz(4.3720923) q[1];
sx q[1];
rz(11.421539) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5231397) q[0];
sx q[0];
rz(-1.8123527) q[0];
sx q[0];
rz(-0.068058204) q[0];
rz(1.990591) q[2];
sx q[2];
rz(-1.3090054) q[2];
sx q[2];
rz(-0.37444886) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.5073587) q[1];
sx q[1];
rz(-0.45295742) q[1];
sx q[1];
rz(1.5276315) q[1];
x q[2];
rz(1.2727229) q[3];
sx q[3];
rz(-1.0469336) q[3];
sx q[3];
rz(3.1071172) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-3.0248489) q[2];
sx q[2];
rz(-1.8086834) q[2];
sx q[2];
rz(-1.9072745) q[2];
rz(-2.0862789) q[3];
sx q[3];
rz(-1.0588667) q[3];
sx q[3];
rz(1.1528667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18594436) q[0];
sx q[0];
rz(-1.9444822) q[0];
sx q[0];
rz(0.82988513) q[0];
rz(2.8886967) q[1];
sx q[1];
rz(-0.77650944) q[1];
sx q[1];
rz(-2.2944962) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69764187) q[0];
sx q[0];
rz(-1.1753433) q[0];
sx q[0];
rz(2.6614499) q[0];
rz(-pi) q[1];
x q[1];
rz(0.29909889) q[2];
sx q[2];
rz(-2.3397589) q[2];
sx q[2];
rz(-2.9849844) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3165247) q[1];
sx q[1];
rz(-1.0267338) q[1];
sx q[1];
rz(-1.6506509) q[1];
rz(-pi) q[2];
rz(-2.794572) q[3];
sx q[3];
rz(-1.9850529) q[3];
sx q[3];
rz(-1.3115209) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.117924) q[2];
sx q[2];
rz(-0.37332049) q[2];
sx q[2];
rz(-0.71933293) q[2];
rz(-1.8524648) q[3];
sx q[3];
rz(-2.9057861) q[3];
sx q[3];
rz(-1.6524564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8711202) q[0];
sx q[0];
rz(-1.2501165) q[0];
sx q[0];
rz(2.5701994) q[0];
rz(-0.96673036) q[1];
sx q[1];
rz(-1.8833908) q[1];
sx q[1];
rz(2.6142696) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8427211) q[0];
sx q[0];
rz(-1.5639595) q[0];
sx q[0];
rz(-2.8323035) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4685523) q[2];
sx q[2];
rz(-1.8290142) q[2];
sx q[2];
rz(2.8845127) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.1151162) q[1];
sx q[1];
rz(-2.3079702) q[1];
sx q[1];
rz(1.2181746) q[1];
rz(-pi) q[2];
x q[2];
rz(0.35964386) q[3];
sx q[3];
rz(-2.375964) q[3];
sx q[3];
rz(1.5639203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.8335235) q[2];
sx q[2];
rz(-1.4493194) q[2];
sx q[2];
rz(-2.4148338) q[2];
rz(0.99749342) q[3];
sx q[3];
rz(-0.62524978) q[3];
sx q[3];
rz(-1.6231026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8961287) q[0];
sx q[0];
rz(-1.6587057) q[0];
sx q[0];
rz(-2.1380651) q[0];
rz(3.1009122) q[1];
sx q[1];
rz(-1.129312) q[1];
sx q[1];
rz(2.3153268) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47047939) q[0];
sx q[0];
rz(-2.2320036) q[0];
sx q[0];
rz(0.83755042) q[0];
rz(-pi) q[1];
x q[1];
rz(0.70010186) q[2];
sx q[2];
rz(-1.6119909) q[2];
sx q[2];
rz(0.61566478) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.367825) q[1];
sx q[1];
rz(-2.7959679) q[1];
sx q[1];
rz(-1.3454382) q[1];
x q[2];
rz(0.63185933) q[3];
sx q[3];
rz(-1.8903036) q[3];
sx q[3];
rz(-2.7416122) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.56759175) q[2];
sx q[2];
rz(-1.2446128) q[2];
sx q[2];
rz(0.91439247) q[2];
rz(-0.10522035) q[3];
sx q[3];
rz(-1.5172232) q[3];
sx q[3];
rz(0.58925327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88702622) q[0];
sx q[0];
rz(-1.5019324) q[0];
sx q[0];
rz(-2.0078833) q[0];
rz(3.1359613) q[1];
sx q[1];
rz(-1.4375604) q[1];
sx q[1];
rz(-2.5240135) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78248258) q[0];
sx q[0];
rz(-1.6164653) q[0];
sx q[0];
rz(-3.111905) q[0];
x q[1];
rz(-0.56474833) q[2];
sx q[2];
rz(-2.4167633) q[2];
sx q[2];
rz(-2.2861779) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.87318201) q[1];
sx q[1];
rz(-2.1408484) q[1];
sx q[1];
rz(0.33318712) q[1];
rz(-pi) q[2];
rz(1.5579617) q[3];
sx q[3];
rz(-1.395263) q[3];
sx q[3];
rz(-0.83735355) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.79409838) q[2];
sx q[2];
rz(-1.8811052) q[2];
sx q[2];
rz(-0.23362544) q[2];
rz(2.1485093) q[3];
sx q[3];
rz(-2.1916094) q[3];
sx q[3];
rz(0.63794199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1722906) q[0];
sx q[0];
rz(-0.06047051) q[0];
sx q[0];
rz(0.0897952) q[0];
rz(2.2604997) q[1];
sx q[1];
rz(-2.6125364) q[1];
sx q[1];
rz(0.18009137) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0463379) q[0];
sx q[0];
rz(-2.1349847) q[0];
sx q[0];
rz(0.5395704) q[0];
rz(-2.1744556) q[2];
sx q[2];
rz(-1.7269616) q[2];
sx q[2];
rz(-0.16317633) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(3.0568911) q[1];
sx q[1];
rz(-1.4086205) q[1];
sx q[1];
rz(-0.2918891) q[1];
rz(-pi) q[2];
rz(-1.5363541) q[3];
sx q[3];
rz(-2.648571) q[3];
sx q[3];
rz(2.8480414) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.2580516) q[2];
sx q[2];
rz(-1.5774612) q[2];
sx q[2];
rz(-1.1759261) q[2];
rz(0.073143395) q[3];
sx q[3];
rz(-2.4172343) q[3];
sx q[3];
rz(-2.6426962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(2.183855) q[0];
sx q[0];
rz(-2.4860005) q[0];
sx q[0];
rz(-1.4181597) q[0];
rz(-1.9765967) q[1];
sx q[1];
rz(-2.5823451) q[1];
sx q[1];
rz(3.1013536) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2049853) q[0];
sx q[0];
rz(-1.4892704) q[0];
sx q[0];
rz(1.8100912) q[0];
rz(-1.6778498) q[2];
sx q[2];
rz(-0.40823001) q[2];
sx q[2];
rz(-1.6183311) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.38899598) q[1];
sx q[1];
rz(-0.95655671) q[1];
sx q[1];
rz(0.1715626) q[1];
x q[2];
rz(-1.2651029) q[3];
sx q[3];
rz(-1.2672658) q[3];
sx q[3];
rz(0.85206735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-3.0030901) q[2];
sx q[2];
rz(-2.3562727) q[2];
sx q[2];
rz(0.17383943) q[2];
rz(-1.396817) q[3];
sx q[3];
rz(-1.4649748) q[3];
sx q[3];
rz(-2.2824536) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1180856) q[0];
sx q[0];
rz(-2.9243587) q[0];
sx q[0];
rz(-1.7370976) q[0];
rz(-1.2069758) q[1];
sx q[1];
rz(-1.6563621) q[1];
sx q[1];
rz(-1.6361902) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0818427) q[0];
sx q[0];
rz(-1.8847701) q[0];
sx q[0];
rz(2.7008675) q[0];
rz(2.5152399) q[2];
sx q[2];
rz(-2.2144631) q[2];
sx q[2];
rz(-2.6477637) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.9951524) q[1];
sx q[1];
rz(-0.69680981) q[1];
sx q[1];
rz(-0.050683024) q[1];
rz(-2.8137915) q[3];
sx q[3];
rz(-2.1468688) q[3];
sx q[3];
rz(0.93186659) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.3796842) q[2];
sx q[2];
rz(-2.6613993) q[2];
sx q[2];
rz(0.15052477) q[2];
rz(-1.6020417) q[3];
sx q[3];
rz(-1.1952885) q[3];
sx q[3];
rz(2.5185744) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4432916) q[0];
sx q[0];
rz(-2.0069831) q[0];
sx q[0];
rz(-0.64569965) q[0];
rz(2.6422016) q[1];
sx q[1];
rz(-1.7130518) q[1];
sx q[1];
rz(1.2649149) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9547573) q[0];
sx q[0];
rz(-1.1823913) q[0];
sx q[0];
rz(-2.3507621) q[0];
rz(2.5876797) q[2];
sx q[2];
rz(-1.950112) q[2];
sx q[2];
rz(-1.6162789) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.7281108) q[1];
sx q[1];
rz(-0.15004798) q[1];
sx q[1];
rz(-0.24197443) q[1];
rz(-1.7725138) q[3];
sx q[3];
rz(-1.9864051) q[3];
sx q[3];
rz(-2.6430074) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.634793) q[2];
sx q[2];
rz(-0.96968499) q[2];
sx q[2];
rz(2.6169422) q[2];
rz(-2.5850885) q[3];
sx q[3];
rz(-1.5126712) q[3];
sx q[3];
rz(0.35486737) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5234579) q[0];
sx q[0];
rz(-2.3045461) q[0];
sx q[0];
rz(-2.8961704) q[0];
rz(0.55631176) q[1];
sx q[1];
rz(-1.5309155) q[1];
sx q[1];
rz(-1.4642749) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9368162) q[0];
sx q[0];
rz(-1.6813155) q[0];
sx q[0];
rz(-2.9226581) q[0];
rz(-pi) q[1];
rz(-0.98430888) q[2];
sx q[2];
rz(-1.3830292) q[2];
sx q[2];
rz(0.078660065) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.163584) q[1];
sx q[1];
rz(-0.89466909) q[1];
sx q[1];
rz(0.34983695) q[1];
x q[2];
rz(1.1226095) q[3];
sx q[3];
rz(-2.9467694) q[3];
sx q[3];
rz(0.71988718) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.8993373) q[2];
sx q[2];
rz(-1.238404) q[2];
sx q[2];
rz(2.396092) q[2];
rz(1.7177104) q[3];
sx q[3];
rz(-2.8427377) q[3];
sx q[3];
rz(2.0132813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2765008) q[0];
sx q[0];
rz(-1.4958953) q[0];
sx q[0];
rz(-2.4687742) q[0];
rz(1.8854234) q[1];
sx q[1];
rz(-2.3354463) q[1];
sx q[1];
rz(-1.068402) q[1];
rz(2.5901569) q[2];
sx q[2];
rz(-0.79179344) q[2];
sx q[2];
rz(1.0168016) q[2];
rz(1.9477378) q[3];
sx q[3];
rz(-1.0300763) q[3];
sx q[3];
rz(-1.1618617) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];