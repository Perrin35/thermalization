OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-3.0394734) q[0];
sx q[0];
rz(-1.4730299) q[0];
sx q[0];
rz(0.13134512) q[0];
rz(0.036845358) q[1];
sx q[1];
rz(-0.52739066) q[1];
sx q[1];
rz(-1.3091458) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0378758) q[0];
sx q[0];
rz(-1.193422) q[0];
sx q[0];
rz(-2.3233633) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0264583) q[2];
sx q[2];
rz(-1.1725468) q[2];
sx q[2];
rz(1.8894701) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.4507452) q[1];
sx q[1];
rz(-1.8826797) q[1];
sx q[1];
rz(1.8212832) q[1];
rz(1.741949) q[3];
sx q[3];
rz(-1.7516836) q[3];
sx q[3];
rz(0.76598841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.5347791) q[2];
sx q[2];
rz(-2.7853192) q[2];
sx q[2];
rz(-0.66205364) q[2];
rz(-2.3946297) q[3];
sx q[3];
rz(-2.2253939) q[3];
sx q[3];
rz(1.0158739) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5544283) q[0];
sx q[0];
rz(-0.14475188) q[0];
sx q[0];
rz(1.1821049) q[0];
rz(-0.74554044) q[1];
sx q[1];
rz(-2.5423971) q[1];
sx q[1];
rz(-0.10339698) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7000843) q[0];
sx q[0];
rz(-1.9943976) q[0];
sx q[0];
rz(1.8258498) q[0];
rz(-pi) q[1];
x q[1];
rz(0.38455268) q[2];
sx q[2];
rz(-2.1332624) q[2];
sx q[2];
rz(-0.99316521) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.7999801) q[1];
sx q[1];
rz(-0.80474058) q[1];
sx q[1];
rz(1.1031723) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3820135) q[3];
sx q[3];
rz(-0.75136614) q[3];
sx q[3];
rz(2.7392859) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.55249247) q[2];
sx q[2];
rz(-1.8742259) q[2];
sx q[2];
rz(-2.6972771) q[2];
rz(1.0537423) q[3];
sx q[3];
rz(-0.070662347) q[3];
sx q[3];
rz(1.1963371) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(1.0788197) q[0];
sx q[0];
rz(-1.8906931) q[0];
sx q[0];
rz(-0.66251063) q[0];
rz(1.484681) q[1];
sx q[1];
rz(-1.0228446) q[1];
sx q[1];
rz(0.2167162) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9154926) q[0];
sx q[0];
rz(-2.2558172) q[0];
sx q[0];
rz(0.95695509) q[0];
rz(-2.6609801) q[2];
sx q[2];
rz(-2.3490902) q[2];
sx q[2];
rz(2.8414244) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.413447) q[1];
sx q[1];
rz(-1.6108509) q[1];
sx q[1];
rz(-0.89786454) q[1];
x q[2];
rz(1.7314579) q[3];
sx q[3];
rz(-1.0815797) q[3];
sx q[3];
rz(2.1833724) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.4776939) q[2];
sx q[2];
rz(-0.52637664) q[2];
sx q[2];
rz(-2.9065175) q[2];
rz(2.7663686) q[3];
sx q[3];
rz(-2.2347361) q[3];
sx q[3];
rz(2.4231329) q[3];
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
rz(pi/2) q[3];
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
rz(1.5821871) q[0];
sx q[0];
rz(-0.087787293) q[0];
sx q[0];
rz(1.0002332) q[0];
rz(1.9384711) q[1];
sx q[1];
rz(-1.5088046) q[1];
sx q[1];
rz(-0.54471725) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4017553) q[0];
sx q[0];
rz(-0.71439636) q[0];
sx q[0];
rz(-0.06846662) q[0];
rz(2.4248872) q[2];
sx q[2];
rz(-0.96895987) q[2];
sx q[2];
rz(-0.81317893) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.1724068) q[1];
sx q[1];
rz(-1.452407) q[1];
sx q[1];
rz(-1.7612996) q[1];
rz(-pi) q[2];
rz(-0.29981837) q[3];
sx q[3];
rz(-2.162419) q[3];
sx q[3];
rz(-3.0772131) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.41877052) q[2];
sx q[2];
rz(-0.79221574) q[2];
sx q[2];
rz(-2.344632) q[2];
rz(-0.289251) q[3];
sx q[3];
rz(-0.97024337) q[3];
sx q[3];
rz(2.7204035) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61409426) q[0];
sx q[0];
rz(-0.088936381) q[0];
sx q[0];
rz(1.8726789) q[0];
rz(0.96619636) q[1];
sx q[1];
rz(-1.393001) q[1];
sx q[1];
rz(0.46636811) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0206576) q[0];
sx q[0];
rz(-2.0657206) q[0];
sx q[0];
rz(1.0878956) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.63299243) q[2];
sx q[2];
rz(-0.55321732) q[2];
sx q[2];
rz(-2.4621682) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6154229) q[1];
sx q[1];
rz(-0.96994441) q[1];
sx q[1];
rz(2.7533965) q[1];
x q[2];
rz(2.2339646) q[3];
sx q[3];
rz(-1.8844127) q[3];
sx q[3];
rz(-2.7682891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.7352778) q[2];
sx q[2];
rz(-2.3390528) q[2];
sx q[2];
rz(2.3373513) q[2];
rz(-2.6546226) q[3];
sx q[3];
rz(-2.099497) q[3];
sx q[3];
rz(-3.13412) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93953472) q[0];
sx q[0];
rz(-0.29629961) q[0];
sx q[0];
rz(0.95788389) q[0];
rz(1.6288039) q[1];
sx q[1];
rz(-1.0572546) q[1];
sx q[1];
rz(-1.5527976) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68715019) q[0];
sx q[0];
rz(-1.5936562) q[0];
sx q[0];
rz(1.4909362) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3423296) q[2];
sx q[2];
rz(-2.7529319) q[2];
sx q[2];
rz(-2.6262019) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.1576288) q[1];
sx q[1];
rz(-2.2215054) q[1];
sx q[1];
rz(-1.7939066) q[1];
x q[2];
rz(-0.62726043) q[3];
sx q[3];
rz(-1.903844) q[3];
sx q[3];
rz(1.2675389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.2019041) q[2];
sx q[2];
rz(-2.0857911) q[2];
sx q[2];
rz(-2.4064257) q[2];
rz(0.9489263) q[3];
sx q[3];
rz(-2.2831423) q[3];
sx q[3];
rz(-2.2775876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3801124) q[0];
sx q[0];
rz(-2.1367456) q[0];
sx q[0];
rz(-0.6066221) q[0];
rz(2.3573549) q[1];
sx q[1];
rz(-1.8083068) q[1];
sx q[1];
rz(-1.7971136) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79216751) q[0];
sx q[0];
rz(-2.0443522) q[0];
sx q[0];
rz(2.4132632) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.49656252) q[2];
sx q[2];
rz(-1.7933868) q[2];
sx q[2];
rz(-2.6816161) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.5658907) q[1];
sx q[1];
rz(-1.1739578) q[1];
sx q[1];
rz(1.8116519) q[1];
rz(1.1056455) q[3];
sx q[3];
rz(-1.4896258) q[3];
sx q[3];
rz(1.2042883) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.532423) q[2];
sx q[2];
rz(-0.50749856) q[2];
sx q[2];
rz(1.3896821) q[2];
rz(-2.9621647) q[3];
sx q[3];
rz(-0.98589412) q[3];
sx q[3];
rz(2.7348886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0824025) q[0];
sx q[0];
rz(-2.817509) q[0];
sx q[0];
rz(2.3865336) q[0];
rz(0.45626196) q[1];
sx q[1];
rz(-1.4249964) q[1];
sx q[1];
rz(1.0770575) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9108601) q[0];
sx q[0];
rz(-0.60711475) q[0];
sx q[0];
rz(-0.40672983) q[0];
rz(-pi) q[1];
rz(-2.3577214) q[2];
sx q[2];
rz(-1.4848564) q[2];
sx q[2];
rz(0.29168561) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.3216599) q[1];
sx q[1];
rz(-2.0256307) q[1];
sx q[1];
rz(-1.8881892) q[1];
rz(-pi) q[2];
rz(0.15865008) q[3];
sx q[3];
rz(-1.0358184) q[3];
sx q[3];
rz(-1.8381929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.72851744) q[2];
sx q[2];
rz(-1.1001526) q[2];
sx q[2];
rz(-0.73823482) q[2];
rz(1.5089367) q[3];
sx q[3];
rz(-1.3239219) q[3];
sx q[3];
rz(-1.9807321) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64410925) q[0];
sx q[0];
rz(-1.4583541) q[0];
sx q[0];
rz(-0.64895502) q[0];
rz(0.68967462) q[1];
sx q[1];
rz(-0.70976218) q[1];
sx q[1];
rz(-2.7241657) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2608503) q[0];
sx q[0];
rz(-0.98310231) q[0];
sx q[0];
rz(-1.2686612) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8990977) q[2];
sx q[2];
rz(-2.8333377) q[2];
sx q[2];
rz(-0.85390845) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.46681225) q[1];
sx q[1];
rz(-1.9793452) q[1];
sx q[1];
rz(1.684434) q[1];
x q[2];
rz(-0.40006684) q[3];
sx q[3];
rz(-1.9377982) q[3];
sx q[3];
rz(0.36959592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.6680341) q[2];
sx q[2];
rz(-1.4213976) q[2];
sx q[2];
rz(-0.045698015) q[2];
rz(-0.33347305) q[3];
sx q[3];
rz(-2.5662751) q[3];
sx q[3];
rz(-3.0157109) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-2.4561653) q[0];
sx q[0];
rz(-2.6021155) q[0];
sx q[0];
rz(-1.9239377) q[0];
rz(-0.24670163) q[1];
sx q[1];
rz(-1.8496937) q[1];
sx q[1];
rz(-0.47952476) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2829153) q[0];
sx q[0];
rz(-2.0020131) q[0];
sx q[0];
rz(2.2515576) q[0];
rz(-0.87534753) q[2];
sx q[2];
rz(-1.9672868) q[2];
sx q[2];
rz(3.0975395) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.785276) q[1];
sx q[1];
rz(-2.0802167) q[1];
sx q[1];
rz(1.2356349) q[1];
rz(-pi) q[2];
rz(-2.6761618) q[3];
sx q[3];
rz(-1.3579391) q[3];
sx q[3];
rz(-1.2133017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.55629998) q[2];
sx q[2];
rz(-1.2703398) q[2];
sx q[2];
rz(-1.7362107) q[2];
rz(-2.4893238) q[3];
sx q[3];
rz(-0.26572078) q[3];
sx q[3];
rz(2.4238267) q[3];
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
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2300867) q[0];
sx q[0];
rz(-1.7510887) q[0];
sx q[0];
rz(-0.44816309) q[0];
rz(-2.3975092) q[1];
sx q[1];
rz(-2.4241445) q[1];
sx q[1];
rz(-1.4345899) q[1];
rz(-0.23575959) q[2];
sx q[2];
rz(-0.59638646) q[2];
sx q[2];
rz(-1.4183945) q[2];
rz(0.35318315) q[3];
sx q[3];
rz(-1.0241057) q[3];
sx q[3];
rz(3.1288341) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
