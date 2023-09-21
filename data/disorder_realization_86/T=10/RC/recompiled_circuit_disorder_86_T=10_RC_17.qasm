OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.3671626) q[0];
sx q[0];
rz(4.055152) q[0];
sx q[0];
rz(11.154296) q[0];
rz(-2.9867759) q[1];
sx q[1];
rz(-0.59564367) q[1];
sx q[1];
rz(-1.6593978) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1887814) q[0];
sx q[0];
rz(-2.6872098) q[0];
sx q[0];
rz(2.7812468) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6719195) q[2];
sx q[2];
rz(-2.854752) q[2];
sx q[2];
rz(-1.6490205) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.60018051) q[1];
sx q[1];
rz(-1.0804847) q[1];
sx q[1];
rz(2.6580826) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0406038) q[3];
sx q[3];
rz(-2.1312993) q[3];
sx q[3];
rz(1.8141754) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.98510629) q[2];
sx q[2];
rz(-0.50922314) q[2];
sx q[2];
rz(0.86581725) q[2];
rz(0.95430294) q[3];
sx q[3];
rz(-1.6031957) q[3];
sx q[3];
rz(-1.8538063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1433379) q[0];
sx q[0];
rz(-1.7049494) q[0];
sx q[0];
rz(-3.1153733) q[0];
rz(1.5401309) q[1];
sx q[1];
rz(-1.5988348) q[1];
sx q[1];
rz(-0.96347934) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.022097691) q[0];
sx q[0];
rz(-2.1848328) q[0];
sx q[0];
rz(3.1387781) q[0];
rz(-pi) q[1];
rz(2.0630402) q[2];
sx q[2];
rz(-2.5139367) q[2];
sx q[2];
rz(0.17222675) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-3.1184247) q[1];
sx q[1];
rz(-1.8903362) q[1];
sx q[1];
rz(2.097514) q[1];
x q[2];
rz(-1.2198592) q[3];
sx q[3];
rz(-1.3388472) q[3];
sx q[3];
rz(-1.9333145) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.5144689) q[2];
sx q[2];
rz(-2.0141979) q[2];
sx q[2];
rz(3.0070686) q[2];
rz(-0.7450122) q[3];
sx q[3];
rz(-0.22694215) q[3];
sx q[3];
rz(2.1988595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9298252) q[0];
sx q[0];
rz(-2.7524502) q[0];
sx q[0];
rz(-2.3441558) q[0];
rz(1.047661) q[1];
sx q[1];
rz(-0.14973775) q[1];
sx q[1];
rz(0.55999666) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3968351) q[0];
sx q[0];
rz(-2.7042537) q[0];
sx q[0];
rz(-2.0646981) q[0];
rz(1.5921002) q[2];
sx q[2];
rz(-1.8739803) q[2];
sx q[2];
rz(1.6456749) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.8289889) q[1];
sx q[1];
rz(-1.7692411) q[1];
sx q[1];
rz(0.36638422) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0873763) q[3];
sx q[3];
rz(-1.9994352) q[3];
sx q[3];
rz(2.7408858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.3893163) q[2];
sx q[2];
rz(-1.1976778) q[2];
sx q[2];
rz(-2.9690572) q[2];
rz(2.1595188) q[3];
sx q[3];
rz(-1.3970102) q[3];
sx q[3];
rz(-1.0579695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.003222) q[0];
sx q[0];
rz(-1.0976185) q[0];
sx q[0];
rz(2.8570783) q[0];
rz(2.8248887) q[1];
sx q[1];
rz(-2.7088294) q[1];
sx q[1];
rz(1.2987312) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7553058) q[0];
sx q[0];
rz(-0.55314976) q[0];
sx q[0];
rz(-1.132071) q[0];
rz(-3.1403055) q[2];
sx q[2];
rz(-0.80438559) q[2];
sx q[2];
rz(-0.019891642) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.97949308) q[1];
sx q[1];
rz(-0.59826189) q[1];
sx q[1];
rz(1.23566) q[1];
x q[2];
rz(1.8907359) q[3];
sx q[3];
rz(-2.7971929) q[3];
sx q[3];
rz(0.26564769) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.6056885) q[2];
sx q[2];
rz(-0.19583344) q[2];
sx q[2];
rz(0.38468012) q[2];
rz(-2.3875333) q[3];
sx q[3];
rz(-2.0850756) q[3];
sx q[3];
rz(-1.6872905) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5383179) q[0];
sx q[0];
rz(-2.1544927) q[0];
sx q[0];
rz(-1.7549365) q[0];
rz(0.23100135) q[1];
sx q[1];
rz(-1.8004386) q[1];
sx q[1];
rz(2.8447661) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7291527) q[0];
sx q[0];
rz(-1.6099596) q[0];
sx q[0];
rz(-2.5705283) q[0];
rz(-0.74283959) q[2];
sx q[2];
rz(-1.683871) q[2];
sx q[2];
rz(-0.46087056) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0823114) q[1];
sx q[1];
rz(-1.1059522) q[1];
sx q[1];
rz(-2.6962198) q[1];
x q[2];
rz(-0.10260251) q[3];
sx q[3];
rz(-0.71628621) q[3];
sx q[3];
rz(-2.7308381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.7632873) q[2];
sx q[2];
rz(-1.3102691) q[2];
sx q[2];
rz(2.7491167) q[2];
rz(-1.9893507) q[3];
sx q[3];
rz(-0.71458721) q[3];
sx q[3];
rz(0.31744441) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0157938) q[0];
sx q[0];
rz(-1.5725461) q[0];
sx q[0];
rz(2.3902067) q[0];
rz(1.8136576) q[1];
sx q[1];
rz(-1.2633879) q[1];
sx q[1];
rz(2.5352535) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6316815) q[0];
sx q[0];
rz(-1.5541549) q[0];
sx q[0];
rz(0.04583866) q[0];
rz(-2.0727856) q[2];
sx q[2];
rz(-2.398218) q[2];
sx q[2];
rz(0.547264) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.6341083) q[1];
sx q[1];
rz(-2.1180696) q[1];
sx q[1];
rz(0.21168153) q[1];
x q[2];
rz(1.353225) q[3];
sx q[3];
rz(-1.4975582) q[3];
sx q[3];
rz(2.6069802) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.63885826) q[2];
sx q[2];
rz(-2.0998462) q[2];
sx q[2];
rz(-1.139337) q[2];
rz(-1.4849439) q[3];
sx q[3];
rz(-1.9610201) q[3];
sx q[3];
rz(-0.10425723) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6095603) q[0];
sx q[0];
rz(-0.74815265) q[0];
sx q[0];
rz(0.50810057) q[0];
rz(1.5787026) q[1];
sx q[1];
rz(-2.0527614) q[1];
sx q[1];
rz(2.3513444) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.761844) q[0];
sx q[0];
rz(-1.8659235) q[0];
sx q[0];
rz(2.5248812) q[0];
rz(3.087567) q[2];
sx q[2];
rz(-2.9276491) q[2];
sx q[2];
rz(-2.8097048) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1192757) q[1];
sx q[1];
rz(-1.1806618) q[1];
sx q[1];
rz(-1.2783865) q[1];
rz(0.4831794) q[3];
sx q[3];
rz(-2.4348767) q[3];
sx q[3];
rz(0.28373517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.053085176) q[2];
sx q[2];
rz(-2.6997456) q[2];
sx q[2];
rz(1.4132168) q[2];
rz(-1.4767856) q[3];
sx q[3];
rz(-1.0352742) q[3];
sx q[3];
rz(3.0800381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7168032) q[0];
sx q[0];
rz(-0.030310832) q[0];
sx q[0];
rz(1.0472263) q[0];
rz(-0.60910243) q[1];
sx q[1];
rz(-1.4139688) q[1];
sx q[1];
rz(-1.75288) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4980709) q[0];
sx q[0];
rz(-0.56034351) q[0];
sx q[0];
rz(-1.5146921) q[0];
rz(-pi) q[1];
x q[1];
rz(0.65865626) q[2];
sx q[2];
rz(-0.63630644) q[2];
sx q[2];
rz(3.0898526) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.74275201) q[1];
sx q[1];
rz(-1.5594149) q[1];
sx q[1];
rz(-3.0497754) q[1];
x q[2];
rz(0.96447585) q[3];
sx q[3];
rz(-1.3250321) q[3];
sx q[3];
rz(2.5442459) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.1887112) q[2];
sx q[2];
rz(-2.7313576) q[2];
sx q[2];
rz(2.2593373) q[2];
rz(1.4011718) q[3];
sx q[3];
rz(-1.1663576) q[3];
sx q[3];
rz(-1.2004948) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27154487) q[0];
sx q[0];
rz(-0.4168059) q[0];
sx q[0];
rz(-1.4260938) q[0];
rz(0.081461279) q[1];
sx q[1];
rz(-1.9790244) q[1];
sx q[1];
rz(2.5833599) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4687913) q[0];
sx q[0];
rz(-1.9580012) q[0];
sx q[0];
rz(-2.0331435) q[0];
x q[1];
rz(1.4633281) q[2];
sx q[2];
rz(-1.5170013) q[2];
sx q[2];
rz(1.4449643) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.7575175) q[1];
sx q[1];
rz(-1.5069403) q[1];
sx q[1];
rz(0.82660316) q[1];
rz(0.48645143) q[3];
sx q[3];
rz(-1.7351741) q[3];
sx q[3];
rz(-2.1742976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.8490863) q[2];
sx q[2];
rz(-1.8722653) q[2];
sx q[2];
rz(-2.9294087) q[2];
rz(2.9296181) q[3];
sx q[3];
rz(-2.4583355) q[3];
sx q[3];
rz(-1.9395444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6367209) q[0];
sx q[0];
rz(-2.3175406) q[0];
sx q[0];
rz(-1.6037534) q[0];
rz(2.3161855) q[1];
sx q[1];
rz(-2.4688265) q[1];
sx q[1];
rz(2.6182981) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6494203) q[0];
sx q[0];
rz(-1.5086552) q[0];
sx q[0];
rz(2.5323575) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5027572) q[2];
sx q[2];
rz(-2.5685852) q[2];
sx q[2];
rz(-0.21493658) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1112422) q[1];
sx q[1];
rz(-0.48644629) q[1];
sx q[1];
rz(-2.4187947) q[1];
rz(2.8966122) q[3];
sx q[3];
rz(-1.7953201) q[3];
sx q[3];
rz(2.6905439) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.837073) q[2];
sx q[2];
rz(-0.7154811) q[2];
sx q[2];
rz(0.26930299) q[2];
rz(-0.4942016) q[3];
sx q[3];
rz(-2.2952357) q[3];
sx q[3];
rz(2.0555029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3257278) q[0];
sx q[0];
rz(-1.6115191) q[0];
sx q[0];
rz(1.4900526) q[0];
rz(1.6745463) q[1];
sx q[1];
rz(-0.29232262) q[1];
sx q[1];
rz(1.2437337) q[1];
rz(1.833563) q[2];
sx q[2];
rz(-1.5599712) q[2];
sx q[2];
rz(0.47906265) q[2];
rz(1.8105103) q[3];
sx q[3];
rz(-2.3532383) q[3];
sx q[3];
rz(-0.91010703) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];