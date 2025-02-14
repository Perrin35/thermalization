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
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1037169) q[0];
sx q[0];
rz(-1.9481707) q[0];
sx q[0];
rz(-0.81822936) q[0];
rz(-pi) q[1];
rz(-0.88850682) q[2];
sx q[2];
rz(-2.4791988) q[2];
sx q[2];
rz(-0.88844013) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.6908474) q[1];
sx q[1];
rz(-1.2589129) q[1];
sx q[1];
rz(1.3203095) q[1];
x q[2];
rz(2.9580826) q[3];
sx q[3];
rz(-1.402463) q[3];
sx q[3];
rz(0.77372293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.6068136) q[2];
sx q[2];
rz(-0.35627347) q[2];
sx q[2];
rz(2.479539) q[2];
rz(-0.74696294) q[3];
sx q[3];
rz(-2.2253939) q[3];
sx q[3];
rz(-1.0158739) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58716431) q[0];
sx q[0];
rz(-0.14475188) q[0];
sx q[0];
rz(-1.9594877) q[0];
rz(-2.3960522) q[1];
sx q[1];
rz(-0.59919557) q[1];
sx q[1];
rz(-0.10339698) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7000843) q[0];
sx q[0];
rz(-1.147195) q[0];
sx q[0];
rz(-1.8258498) q[0];
rz(2.1680122) q[2];
sx q[2];
rz(-1.8937308) q[2];
sx q[2];
rz(-0.79018776) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.1704639) q[1];
sx q[1];
rz(-0.87201768) q[1];
sx q[1];
rz(0.43817313) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.96805) q[3];
sx q[3];
rz(-2.3056917) q[3];
sx q[3];
rz(-0.65803448) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.55249247) q[2];
sx q[2];
rz(-1.2673667) q[2];
sx q[2];
rz(-2.6972771) q[2];
rz(1.0537423) q[3];
sx q[3];
rz(-3.0709303) q[3];
sx q[3];
rz(1.9452555) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.062773) q[0];
sx q[0];
rz(-1.2508996) q[0];
sx q[0];
rz(-0.66251063) q[0];
rz(1.6569116) q[1];
sx q[1];
rz(-1.0228446) q[1];
sx q[1];
rz(2.9248765) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22610006) q[0];
sx q[0];
rz(-0.88577548) q[0];
sx q[0];
rz(-2.1846376) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6609801) q[2];
sx q[2];
rz(-0.79250249) q[2];
sx q[2];
rz(0.30016826) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.81074508) q[1];
sx q[1];
rz(-0.89850366) q[1];
sx q[1];
rz(3.0903893) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6469798) q[3];
sx q[3];
rz(-1.429116) q[3];
sx q[3];
rz(2.6050267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.66389877) q[2];
sx q[2];
rz(-2.615216) q[2];
sx q[2];
rz(2.9065175) q[2];
rz(-2.7663686) q[3];
sx q[3];
rz(-0.90685654) q[3];
sx q[3];
rz(-0.71845976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5594056) q[0];
sx q[0];
rz(-0.087787293) q[0];
sx q[0];
rz(2.1413595) q[0];
rz(-1.9384711) q[1];
sx q[1];
rz(-1.5088046) q[1];
sx q[1];
rz(0.54471725) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7398374) q[0];
sx q[0];
rz(-2.4271963) q[0];
sx q[0];
rz(0.06846662) q[0];
rz(0.71670549) q[2];
sx q[2];
rz(-0.96895987) q[2];
sx q[2];
rz(0.81317893) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.9904321) q[1];
sx q[1];
rz(-2.9176788) q[1];
sx q[1];
rz(-1.0099221) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.29981837) q[3];
sx q[3];
rz(-0.97917367) q[3];
sx q[3];
rz(3.0772131) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.41877052) q[2];
sx q[2];
rz(-2.3493769) q[2];
sx q[2];
rz(-2.344632) q[2];
rz(0.289251) q[3];
sx q[3];
rz(-2.1713493) q[3];
sx q[3];
rz(2.7204035) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5274984) q[0];
sx q[0];
rz(-3.0526563) q[0];
sx q[0];
rz(-1.2689137) q[0];
rz(-0.96619636) q[1];
sx q[1];
rz(-1.393001) q[1];
sx q[1];
rz(-0.46636811) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0206576) q[0];
sx q[0];
rz(-1.0758721) q[0];
sx q[0];
rz(1.0878956) q[0];
rz(-2.5086002) q[2];
sx q[2];
rz(-0.55321732) q[2];
sx q[2];
rz(-0.67942441) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.18257817) q[1];
sx q[1];
rz(-1.8883288) q[1];
sx q[1];
rz(2.2081801) q[1];
rz(-pi) q[2];
rz(0.39042274) q[3];
sx q[3];
rz(-2.1964034) q[3];
sx q[3];
rz(-2.1805891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.40631488) q[2];
sx q[2];
rz(-0.80253989) q[2];
sx q[2];
rz(2.3373513) q[2];
rz(-0.4869701) q[3];
sx q[3];
rz(-1.0420957) q[3];
sx q[3];
rz(0.0074726661) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93953472) q[0];
sx q[0];
rz(-0.29629961) q[0];
sx q[0];
rz(-0.95788389) q[0];
rz(-1.6288039) q[1];
sx q[1];
rz(-1.0572546) q[1];
sx q[1];
rz(1.5527976) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68715019) q[0];
sx q[0];
rz(-1.5479364) q[0];
sx q[0];
rz(-1.4909362) q[0];
x q[1];
rz(-0.79926305) q[2];
sx q[2];
rz(-0.38866079) q[2];
sx q[2];
rz(-2.6262019) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.342345) q[1];
sx q[1];
rz(-2.4589897) q[1];
sx q[1];
rz(2.8587539) q[1];
x q[2];
rz(2.6090066) q[3];
sx q[3];
rz(-2.4420693) q[3];
sx q[3];
rz(0.12040779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.9396886) q[2];
sx q[2];
rz(-1.0558015) q[2];
sx q[2];
rz(2.4064257) q[2];
rz(0.9489263) q[3];
sx q[3];
rz(-0.85845033) q[3];
sx q[3];
rz(2.2775876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-0.76148024) q[0];
sx q[0];
rz(-2.1367456) q[0];
sx q[0];
rz(2.5349706) q[0];
rz(2.3573549) q[1];
sx q[1];
rz(-1.8083068) q[1];
sx q[1];
rz(-1.7971136) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30596581) q[0];
sx q[0];
rz(-2.2971662) q[0];
sx q[0];
rz(2.4854922) q[0];
rz(-pi) q[1];
rz(0.49656252) q[2];
sx q[2];
rz(-1.3482058) q[2];
sx q[2];
rz(0.45997657) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.231338) q[1];
sx q[1];
rz(-1.3489854) q[1];
sx q[1];
rz(-2.7342058) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1056455) q[3];
sx q[3];
rz(-1.4896258) q[3];
sx q[3];
rz(-1.9373044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.6091696) q[2];
sx q[2];
rz(-0.50749856) q[2];
sx q[2];
rz(1.7519105) q[2];
rz(0.17942795) q[3];
sx q[3];
rz(-0.98589412) q[3];
sx q[3];
rz(-0.40670407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0824025) q[0];
sx q[0];
rz(-0.32408369) q[0];
sx q[0];
rz(2.3865336) q[0];
rz(-2.6853307) q[1];
sx q[1];
rz(-1.7165963) q[1];
sx q[1];
rz(-1.0770575) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74769831) q[0];
sx q[0];
rz(-1.0193045) q[0];
sx q[0];
rz(1.3026139) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3577214) q[2];
sx q[2];
rz(-1.6567363) q[2];
sx q[2];
rz(0.29168561) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.7491319) q[1];
sx q[1];
rz(-1.8549671) q[1];
sx q[1];
rz(0.47537132) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3102688) q[3];
sx q[3];
rz(-2.5857877) q[3];
sx q[3];
rz(2.1422841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.4130752) q[2];
sx q[2];
rz(-2.04144) q[2];
sx q[2];
rz(0.73823482) q[2];
rz(-1.5089367) q[3];
sx q[3];
rz(-1.3239219) q[3];
sx q[3];
rz(-1.1608605) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4974834) q[0];
sx q[0];
rz(-1.4583541) q[0];
sx q[0];
rz(0.64895502) q[0];
rz(-2.451918) q[1];
sx q[1];
rz(-0.70976218) q[1];
sx q[1];
rz(0.41742691) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74873105) q[0];
sx q[0];
rz(-2.4890206) q[0];
sx q[0];
rz(2.7215385) q[0];
rz(-pi) q[1];
rz(-1.6471049) q[2];
sx q[2];
rz(-1.869749) q[2];
sx q[2];
rz(1.1079009) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.0829186) q[1];
sx q[1];
rz(-1.6750458) q[1];
sx q[1];
rz(-0.41091316) q[1];
rz(-1.1754009) q[3];
sx q[3];
rz(-1.1987276) q[3];
sx q[3];
rz(-1.7898066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.6680341) q[2];
sx q[2];
rz(-1.7201951) q[2];
sx q[2];
rz(3.0958946) q[2];
rz(-2.8081196) q[3];
sx q[3];
rz(-0.57531753) q[3];
sx q[3];
rz(0.12588178) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68542737) q[0];
sx q[0];
rz(-2.6021155) q[0];
sx q[0];
rz(1.9239377) q[0];
rz(-0.24670163) q[1];
sx q[1];
rz(-1.291899) q[1];
sx q[1];
rz(-2.6620679) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1801301) q[0];
sx q[0];
rz(-0.96213522) q[0];
sx q[0];
rz(0.5345688) q[0];
rz(2.6423655) q[2];
sx q[2];
rz(-2.2031234) q[2];
sx q[2];
rz(1.8385173) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.61726924) q[1];
sx q[1];
rz(-1.8620544) q[1];
sx q[1];
rz(2.6074383) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3334951) q[3];
sx q[3];
rz(-1.1166683) q[3];
sx q[3];
rz(2.8898006) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.55629998) q[2];
sx q[2];
rz(-1.8712529) q[2];
sx q[2];
rz(-1.4053819) q[2];
rz(-0.65226883) q[3];
sx q[3];
rz(-0.26572078) q[3];
sx q[3];
rz(0.71776596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(2.9058331) q[2];
sx q[2];
rz(-0.59638646) q[2];
sx q[2];
rz(-1.4183945) q[2];
rz(2.1461829) q[3];
sx q[3];
rz(-1.8707471) q[3];
sx q[3];
rz(-1.7729014) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
