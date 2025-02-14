OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.9118948) q[0];
sx q[0];
rz(-1.2461646) q[0];
sx q[0];
rz(-1.5204313) q[0];
rz(0.35960943) q[1];
sx q[1];
rz(-2.8728027) q[1];
sx q[1];
rz(-0.51528817) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20977192) q[0];
sx q[0];
rz(-1.2270667) q[0];
sx q[0];
rz(-1.3552257) q[0];
rz(-pi) q[1];
rz(-2.5955321) q[2];
sx q[2];
rz(-1.310643) q[2];
sx q[2];
rz(-0.85064155) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.88564116) q[1];
sx q[1];
rz(-2.6868539) q[1];
sx q[1];
rz(0.498147) q[1];
rz(0.73909594) q[3];
sx q[3];
rz(-1.0217654) q[3];
sx q[3];
rz(2.9573553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.9109853) q[2];
sx q[2];
rz(-1.5585941) q[2];
sx q[2];
rz(-2.4674463) q[2];
rz(2.9437183) q[3];
sx q[3];
rz(-1.9015692) q[3];
sx q[3];
rz(1.5052634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3182217) q[0];
sx q[0];
rz(-2.8405393) q[0];
sx q[0];
rz(0.2163042) q[0];
rz(-3.0220616) q[1];
sx q[1];
rz(-2.0152338) q[1];
sx q[1];
rz(-1.4215887) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3960806) q[0];
sx q[0];
rz(-2.0876309) q[0];
sx q[0];
rz(-3.050258) q[0];
x q[1];
rz(1.7825215) q[2];
sx q[2];
rz(-1.9634517) q[2];
sx q[2];
rz(-2.3711039) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6232963) q[1];
sx q[1];
rz(-0.5591679) q[1];
sx q[1];
rz(-2.1914047) q[1];
rz(0.98540636) q[3];
sx q[3];
rz(-2.2088361) q[3];
sx q[3];
rz(0.93321484) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.30430946) q[2];
sx q[2];
rz(-1.9459566) q[2];
sx q[2];
rz(2.1514814) q[2];
rz(2.385251) q[3];
sx q[3];
rz(-2.4939311) q[3];
sx q[3];
rz(-2.9653449) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15683098) q[0];
sx q[0];
rz(-0.75060833) q[0];
sx q[0];
rz(1.1358776) q[0];
rz(-2.1319977) q[1];
sx q[1];
rz(-1.986809) q[1];
sx q[1];
rz(-0.44581595) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8404482) q[0];
sx q[0];
rz(-3.0160286) q[0];
sx q[0];
rz(-2.156394) q[0];
rz(3.1059044) q[2];
sx q[2];
rz(-2.4378187) q[2];
sx q[2];
rz(1.4427217) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.43128428) q[1];
sx q[1];
rz(-0.7529707) q[1];
sx q[1];
rz(0.94628872) q[1];
rz(3.0296311) q[3];
sx q[3];
rz(-2.6918133) q[3];
sx q[3];
rz(-1.2686319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.77057114) q[2];
sx q[2];
rz(-2.1393445) q[2];
sx q[2];
rz(-2.501343) q[2];
rz(-2.443743) q[3];
sx q[3];
rz(-1.6777439) q[3];
sx q[3];
rz(2.8016134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(2.8401538) q[0];
sx q[0];
rz(-2.2541663) q[0];
sx q[0];
rz(1.812717) q[0];
rz(1.2170732) q[1];
sx q[1];
rz(-0.71958676) q[1];
sx q[1];
rz(2.365239) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76809873) q[0];
sx q[0];
rz(-1.2006309) q[0];
sx q[0];
rz(-0.64565701) q[0];
x q[1];
rz(2.0587875) q[2];
sx q[2];
rz(-1.0363724) q[2];
sx q[2];
rz(-2.7492409) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4763219) q[1];
sx q[1];
rz(-1.5821536) q[1];
sx q[1];
rz(-1.2922056) q[1];
rz(1.5723002) q[3];
sx q[3];
rz(-0.88620159) q[3];
sx q[3];
rz(-0.94438121) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.85294574) q[2];
sx q[2];
rz(-2.8198346) q[2];
sx q[2];
rz(1.9256437) q[2];
rz(-1.1207885) q[3];
sx q[3];
rz(-1.8782764) q[3];
sx q[3];
rz(-2.2585675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
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
rz(-2.9057366) q[0];
sx q[0];
rz(-2.164542) q[0];
sx q[0];
rz(-0.46347722) q[0];
rz(2.0526759) q[1];
sx q[1];
rz(-1.8810279) q[1];
sx q[1];
rz(1.3190528) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4396297) q[0];
sx q[0];
rz(-1.3101398) q[0];
sx q[0];
rz(-1.319665) q[0];
rz(1.7402788) q[2];
sx q[2];
rz(-2.3613644) q[2];
sx q[2];
rz(1.7192769) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6336167) q[1];
sx q[1];
rz(-2.1514261) q[1];
sx q[1];
rz(0.10426048) q[1];
rz(-pi) q[2];
rz(-1.9604881) q[3];
sx q[3];
rz(-2.6929085) q[3];
sx q[3];
rz(-1.9657024) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.9702381) q[2];
sx q[2];
rz(-2.3562608) q[2];
sx q[2];
rz(-0.63924092) q[2];
rz(-0.81131896) q[3];
sx q[3];
rz(-2.6526484) q[3];
sx q[3];
rz(1.5343687) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2163579) q[0];
sx q[0];
rz(-0.41330591) q[0];
sx q[0];
rz(-0.21892029) q[0];
rz(1.9256598) q[1];
sx q[1];
rz(-2.6775807) q[1];
sx q[1];
rz(-1.6485515) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7812778) q[0];
sx q[0];
rz(-2.6584315) q[0];
sx q[0];
rz(0.29668087) q[0];
rz(-pi) q[1];
rz(-1.5677388) q[2];
sx q[2];
rz(-1.1532461) q[2];
sx q[2];
rz(1.9034346) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.1472125) q[1];
sx q[1];
rz(-1.5541422) q[1];
sx q[1];
rz(-2.4092595) q[1];
x q[2];
rz(3.1102255) q[3];
sx q[3];
rz(-0.099197242) q[3];
sx q[3];
rz(-1.8435696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.0367937) q[2];
sx q[2];
rz(-1.3497738) q[2];
sx q[2];
rz(1.8297423) q[2];
rz(-1.6364243) q[3];
sx q[3];
rz(-1.2994095) q[3];
sx q[3];
rz(2.6817491) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5055607) q[0];
sx q[0];
rz(-0.4902896) q[0];
sx q[0];
rz(-1.7946515) q[0];
rz(-1.3268283) q[1];
sx q[1];
rz(-0.83419269) q[1];
sx q[1];
rz(0.38536513) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6198719) q[0];
sx q[0];
rz(-1.6024029) q[0];
sx q[0];
rz(2.7824529) q[0];
rz(1.0406251) q[2];
sx q[2];
rz(-1.1405924) q[2];
sx q[2];
rz(0.29597798) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.24033156) q[1];
sx q[1];
rz(-1.5973395) q[1];
sx q[1];
rz(-0.44601299) q[1];
rz(-pi) q[2];
rz(0.47355424) q[3];
sx q[3];
rz(-1.9723168) q[3];
sx q[3];
rz(0.36916379) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.68606004) q[2];
sx q[2];
rz(-2.3251688) q[2];
sx q[2];
rz(-1.4875937) q[2];
rz(-2.9758596) q[3];
sx q[3];
rz(-2.3159852) q[3];
sx q[3];
rz(2.9674271) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1626749) q[0];
sx q[0];
rz(-0.62541494) q[0];
sx q[0];
rz(-2.9130574) q[0];
rz(-2.8288815) q[1];
sx q[1];
rz(-2.2622908) q[1];
sx q[1];
rz(-1.762134) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1918743) q[0];
sx q[0];
rz(-2.7316425) q[0];
sx q[0];
rz(-0.19522841) q[0];
x q[1];
rz(2.8174322) q[2];
sx q[2];
rz(-1.7334786) q[2];
sx q[2];
rz(1.206347) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.5818565) q[1];
sx q[1];
rz(-0.35071555) q[1];
sx q[1];
rz(1.4697671) q[1];
rz(-0.70693858) q[3];
sx q[3];
rz(-1.9732554) q[3];
sx q[3];
rz(-1.5618531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.10869965) q[2];
sx q[2];
rz(-1.0445107) q[2];
sx q[2];
rz(1.0408939) q[2];
rz(0.4852455) q[3];
sx q[3];
rz(-1.3121366) q[3];
sx q[3];
rz(0.41788873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(0.27337209) q[0];
sx q[0];
rz(-0.98216787) q[0];
sx q[0];
rz(-0.81106538) q[0];
rz(1.9048994) q[1];
sx q[1];
rz(-0.86634723) q[1];
sx q[1];
rz(2.9467357) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70640608) q[0];
sx q[0];
rz(-1.8098477) q[0];
sx q[0];
rz(1.3066221) q[0];
x q[1];
rz(-0.56320517) q[2];
sx q[2];
rz(-0.99894542) q[2];
sx q[2];
rz(2.2541313) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.5681618) q[1];
sx q[1];
rz(-2.5622612) q[1];
sx q[1];
rz(-0.90888494) q[1];
x q[2];
rz(2.2236597) q[3];
sx q[3];
rz(-2.0313162) q[3];
sx q[3];
rz(1.1499109) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.98086944) q[2];
sx q[2];
rz(-1.8941433) q[2];
sx q[2];
rz(0.54171872) q[2];
rz(-1.8187652) q[3];
sx q[3];
rz(-1.3397237) q[3];
sx q[3];
rz(2.8790348) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5498891) q[0];
sx q[0];
rz(-1.5174958) q[0];
sx q[0];
rz(-0.24895915) q[0];
rz(1.9242363) q[1];
sx q[1];
rz(-1.8739506) q[1];
sx q[1];
rz(0.41044661) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7792203) q[0];
sx q[0];
rz(-1.1817314) q[0];
sx q[0];
rz(-2.7011407) q[0];
rz(-pi) q[1];
rz(1.4538953) q[2];
sx q[2];
rz(-0.091594897) q[2];
sx q[2];
rz(-0.26008886) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6140266) q[1];
sx q[1];
rz(-2.0260794) q[1];
sx q[1];
rz(-0.27099314) q[1];
rz(-pi) q[2];
rz(-2.5312891) q[3];
sx q[3];
rz(-1.9962976) q[3];
sx q[3];
rz(0.60275383) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.88811389) q[2];
sx q[2];
rz(-1.9733182) q[2];
sx q[2];
rz(-1.1178364) q[2];
rz(-3.0228293) q[3];
sx q[3];
rz(-0.6898841) q[3];
sx q[3];
rz(1.9411055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4111907) q[0];
sx q[0];
rz(-1.7350736) q[0];
sx q[0];
rz(-0.34179678) q[0];
rz(0.047601184) q[1];
sx q[1];
rz(-1.2536512) q[1];
sx q[1];
rz(1.8258078) q[1];
rz(-2.6834221) q[2];
sx q[2];
rz(-2.2728163) q[2];
sx q[2];
rz(-1.386281) q[2];
rz(0.94275311) q[3];
sx q[3];
rz(-0.92798622) q[3];
sx q[3];
rz(-0.29792132) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
