OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.5717918) q[0];
sx q[0];
rz(-0.90919149) q[0];
sx q[0];
rz(-1.6348913) q[0];
rz(1.2070967) q[1];
sx q[1];
rz(-2.6488882) q[1];
sx q[1];
rz(2.5821813) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4643505) q[0];
sx q[0];
rz(-1.4046245) q[0];
sx q[0];
rz(-0.31810036) q[0];
rz(-1.962553) q[2];
sx q[2];
rz(-1.9843352) q[2];
sx q[2];
rz(2.3512083) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.2845498) q[1];
sx q[1];
rz(-2.8618097) q[1];
sx q[1];
rz(-1.5955052) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1640638) q[3];
sx q[3];
rz(-1.1479605) q[3];
sx q[3];
rz(-0.55227597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.5514465) q[2];
sx q[2];
rz(-2.2821653) q[2];
sx q[2];
rz(-2.0476511) q[2];
rz(-1.9048196) q[3];
sx q[3];
rz(-1.9779132) q[3];
sx q[3];
rz(2.8804603) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(1.3742438) q[0];
sx q[0];
rz(-1.2232895) q[0];
sx q[0];
rz(-0.71794024) q[0];
rz(-1.1942489) q[1];
sx q[1];
rz(-0.4078882) q[1];
sx q[1];
rz(1.492122) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6787918) q[0];
sx q[0];
rz(-0.18175069) q[0];
sx q[0];
rz(0.75580718) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9487593) q[2];
sx q[2];
rz(-1.9828258) q[2];
sx q[2];
rz(0.90726024) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.1528984) q[1];
sx q[1];
rz(-1.1225268) q[1];
sx q[1];
rz(2.6397735) q[1];
rz(-pi) q[2];
x q[2];
rz(0.82573311) q[3];
sx q[3];
rz(-1.7307708) q[3];
sx q[3];
rz(-2.3478594) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.0933928) q[2];
sx q[2];
rz(-0.77646774) q[2];
sx q[2];
rz(0.36273599) q[2];
rz(2.2705966) q[3];
sx q[3];
rz(-1.769915) q[3];
sx q[3];
rz(1.3236275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.094630346) q[0];
sx q[0];
rz(-1.8383263) q[0];
sx q[0];
rz(2.5808425) q[0];
rz(-2.1669855) q[1];
sx q[1];
rz(-0.86169306) q[1];
sx q[1];
rz(-0.21751705) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4982078) q[0];
sx q[0];
rz(-1.5899993) q[0];
sx q[0];
rz(2.5570844) q[0];
x q[1];
rz(-1.0571805) q[2];
sx q[2];
rz(-1.660778) q[2];
sx q[2];
rz(-1.0124026) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.8629093) q[1];
sx q[1];
rz(-2.8562198) q[1];
sx q[1];
rz(1.1527658) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5325696) q[3];
sx q[3];
rz(-1.5191188) q[3];
sx q[3];
rz(-0.9059815) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.64323032) q[2];
sx q[2];
rz(-2.3537894) q[2];
sx q[2];
rz(-2.1882449) q[2];
rz(-2.0638454) q[3];
sx q[3];
rz(-1.4606303) q[3];
sx q[3];
rz(-2.6673711) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0039907) q[0];
sx q[0];
rz(-0.1739665) q[0];
sx q[0];
rz(2.2250788) q[0];
rz(-0.42463955) q[1];
sx q[1];
rz(-1.759513) q[1];
sx q[1];
rz(2.0133846) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5461499) q[0];
sx q[0];
rz(-0.91888035) q[0];
sx q[0];
rz(-0.94080399) q[0];
rz(-pi) q[1];
rz(-1.7776971) q[2];
sx q[2];
rz(-0.82447663) q[2];
sx q[2];
rz(1.7306223) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0227851) q[1];
sx q[1];
rz(-1.8703543) q[1];
sx q[1];
rz(-0.54356281) q[1];
rz(-2.1682285) q[3];
sx q[3];
rz(-1.7586244) q[3];
sx q[3];
rz(2.0300161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.8757223) q[2];
sx q[2];
rz(-1.6758726) q[2];
sx q[2];
rz(-1.9423368) q[2];
rz(-1.1572329) q[3];
sx q[3];
rz(-2.3374989) q[3];
sx q[3];
rz(-0.51586241) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98738325) q[0];
sx q[0];
rz(-0.043954285) q[0];
sx q[0];
rz(0.92535812) q[0];
rz(2.5788653) q[1];
sx q[1];
rz(-2.1268763) q[1];
sx q[1];
rz(-0.37511197) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1031418) q[0];
sx q[0];
rz(-1.4618357) q[0];
sx q[0];
rz(-2.5953963) q[0];
rz(2.0061467) q[2];
sx q[2];
rz(-1.5223961) q[2];
sx q[2];
rz(-1.8627) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.8598392) q[1];
sx q[1];
rz(-1.5956164) q[1];
sx q[1];
rz(2.9725705) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.971832) q[3];
sx q[3];
rz(-1.9936863) q[3];
sx q[3];
rz(-2.324375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.23842266) q[2];
sx q[2];
rz(-2.2245202) q[2];
sx q[2];
rz(-0.49752107) q[2];
rz(2.7072952) q[3];
sx q[3];
rz(-1.4562166) q[3];
sx q[3];
rz(0.98880497) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76181805) q[0];
sx q[0];
rz(-2.2803545) q[0];
sx q[0];
rz(0.27989835) q[0];
rz(2.1474536) q[1];
sx q[1];
rz(-2.2810664) q[1];
sx q[1];
rz(2.5631189) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6743603) q[0];
sx q[0];
rz(-2.0772837) q[0];
sx q[0];
rz(2.0048398) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6874993) q[2];
sx q[2];
rz(-2.3857949) q[2];
sx q[2];
rz(2.2735655) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.0117409) q[1];
sx q[1];
rz(-1.7659582) q[1];
sx q[1];
rz(-0.73634781) q[1];
rz(-pi) q[2];
rz(-1.6719477) q[3];
sx q[3];
rz(-1.3037852) q[3];
sx q[3];
rz(-1.8570516) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.0824739) q[2];
sx q[2];
rz(-2.5422577) q[2];
sx q[2];
rz(-0.88258755) q[2];
rz(0.62697083) q[3];
sx q[3];
rz(-0.65492237) q[3];
sx q[3];
rz(0.89491189) q[3];
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
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6249228) q[0];
sx q[0];
rz(-2.1997917) q[0];
sx q[0];
rz(3.140977) q[0];
rz(-2.2387538) q[1];
sx q[1];
rz(-2.2689464) q[1];
sx q[1];
rz(0.045086233) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49360291) q[0];
sx q[0];
rz(-1.6479074) q[0];
sx q[0];
rz(0.48992975) q[0];
rz(-pi) q[1];
rz(-2.3995724) q[2];
sx q[2];
rz(-1.4145523) q[2];
sx q[2];
rz(-0.39575037) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.5645091) q[1];
sx q[1];
rz(-2.2840743) q[1];
sx q[1];
rz(1.1611931) q[1];
x q[2];
rz(-0.61130133) q[3];
sx q[3];
rz(-2.0620278) q[3];
sx q[3];
rz(0.69913759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.2948598) q[2];
sx q[2];
rz(-2.5856057) q[2];
sx q[2];
rz(2.2247458) q[2];
rz(0.89546853) q[3];
sx q[3];
rz(-1.6296891) q[3];
sx q[3];
rz(-1.9303314) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13670066) q[0];
sx q[0];
rz(-3.0558375) q[0];
sx q[0];
rz(-0.46491796) q[0];
rz(1.9526941) q[1];
sx q[1];
rz(-1.7949972) q[1];
sx q[1];
rz(-2.7704923) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44921267) q[0];
sx q[0];
rz(-1.977773) q[0];
sx q[0];
rz(-1.1212342) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6109054) q[2];
sx q[2];
rz(-2.3903987) q[2];
sx q[2];
rz(1.7263279) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6976561) q[1];
sx q[1];
rz(-0.84007004) q[1];
sx q[1];
rz(2.1751462) q[1];
x q[2];
rz(-1.8598167) q[3];
sx q[3];
rz(-1.9403696) q[3];
sx q[3];
rz(0.7909382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.4798639) q[2];
sx q[2];
rz(-0.78833818) q[2];
sx q[2];
rz(-2.7313477) q[2];
rz(1.6346301) q[3];
sx q[3];
rz(-2.7362636) q[3];
sx q[3];
rz(-3.098367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0931382) q[0];
sx q[0];
rz(-2.6148836) q[0];
sx q[0];
rz(0.17573892) q[0];
rz(-2.3161092) q[1];
sx q[1];
rz(-1.8800507) q[1];
sx q[1];
rz(-2.3186191) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9776445) q[0];
sx q[0];
rz(-0.86714449) q[0];
sx q[0];
rz(-2.4022341) q[0];
rz(-pi) q[1];
rz(-0.31655471) q[2];
sx q[2];
rz(-1.9783894) q[2];
sx q[2];
rz(-2.0616693) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.3210248) q[1];
sx q[1];
rz(-1.9181644) q[1];
sx q[1];
rz(0.72183164) q[1];
rz(-pi) q[2];
rz(-3.104745) q[3];
sx q[3];
rz(-1.1603796) q[3];
sx q[3];
rz(-2.0582123) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.2814111) q[2];
sx q[2];
rz(-1.9965636) q[2];
sx q[2];
rz(-0.54110503) q[2];
rz(0.42630729) q[3];
sx q[3];
rz(-2.6796894) q[3];
sx q[3];
rz(2.2947252) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68778872) q[0];
sx q[0];
rz(-2.332088) q[0];
sx q[0];
rz(-0.33568207) q[0];
rz(0.022739284) q[1];
sx q[1];
rz(-0.72240654) q[1];
sx q[1];
rz(0.67363277) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.861872) q[0];
sx q[0];
rz(-2.4422944) q[0];
sx q[0];
rz(-2.5847816) q[0];
rz(-pi) q[1];
rz(-0.19825528) q[2];
sx q[2];
rz(-2.6645244) q[2];
sx q[2];
rz(-1.4587634) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.9424098) q[1];
sx q[1];
rz(-1.7734151) q[1];
sx q[1];
rz(-2.1658648) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7106179) q[3];
sx q[3];
rz(-1.9019902) q[3];
sx q[3];
rz(0.14023031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6966887) q[2];
sx q[2];
rz(-2.1477369) q[2];
sx q[2];
rz(-1.0059086) q[2];
rz(2.3575947) q[3];
sx q[3];
rz(-2.8204212) q[3];
sx q[3];
rz(-1.706749) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1226596) q[0];
sx q[0];
rz(-1.4228595) q[0];
sx q[0];
rz(2.0056437) q[0];
rz(0.38416531) q[1];
sx q[1];
rz(-1.1687678) q[1];
sx q[1];
rz(1.6484177) q[1];
rz(-1.775209) q[2];
sx q[2];
rz(-1.778893) q[2];
sx q[2];
rz(-1.9472029) q[2];
rz(0.85687153) q[3];
sx q[3];
rz(-1.8663434) q[3];
sx q[3];
rz(0.28874884) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
