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
rz(-3.0102475) q[0];
rz(-3.1047473) q[1];
sx q[1];
rz(-2.614202) q[1];
sx q[1];
rz(-1.8324469) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15802424) q[0];
sx q[0];
rz(-2.3166172) q[0];
sx q[0];
rz(-2.0963066) q[0];
rz(-pi) q[1];
rz(0.88850682) q[2];
sx q[2];
rz(-2.4791988) q[2];
sx q[2];
rz(0.88844013) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.19840163) q[1];
sx q[1];
rz(-1.8089589) q[1];
sx q[1];
rz(-0.32126255) q[1];
rz(-pi) q[2];
rz(0.18351002) q[3];
sx q[3];
rz(-1.402463) q[3];
sx q[3];
rz(2.3678697) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.5347791) q[2];
sx q[2];
rz(-0.35627347) q[2];
sx q[2];
rz(0.66205364) q[2];
rz(-2.3946297) q[3];
sx q[3];
rz(-2.2253939) q[3];
sx q[3];
rz(1.0158739) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58716431) q[0];
sx q[0];
rz(-0.14475188) q[0];
sx q[0];
rz(1.9594877) q[0];
rz(-0.74554044) q[1];
sx q[1];
rz(-2.5423971) q[1];
sx q[1];
rz(3.0381957) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.11907) q[0];
sx q[0];
rz(-1.8028717) q[0];
sx q[0];
rz(0.43605767) q[0];
rz(0.97358046) q[2];
sx q[2];
rz(-1.8937308) q[2];
sx q[2];
rz(-2.3514049) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.34161257) q[1];
sx q[1];
rz(-0.80474058) q[1];
sx q[1];
rz(2.0384203) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3820135) q[3];
sx q[3];
rz(-2.3902265) q[3];
sx q[3];
rz(-2.7392859) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.55249247) q[2];
sx q[2];
rz(-1.8742259) q[2];
sx q[2];
rz(2.6972771) q[2];
rz(-1.0537423) q[3];
sx q[3];
rz(-0.070662347) q[3];
sx q[3];
rz(1.9452555) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.062773) q[0];
sx q[0];
rz(-1.2508996) q[0];
sx q[0];
rz(0.66251063) q[0];
rz(1.484681) q[1];
sx q[1];
rz(-1.0228446) q[1];
sx q[1];
rz(-2.9248765) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2162735) q[0];
sx q[0];
rz(-2.0331622) q[0];
sx q[0];
rz(-0.78512773) q[0];
rz(-2.4091085) q[2];
sx q[2];
rz(-1.9062796) q[2];
sx q[2];
rz(-1.5200638) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3308476) q[1];
sx q[1];
rz(-2.243089) q[1];
sx q[1];
rz(-0.051203392) q[1];
rz(1.4101348) q[3];
sx q[3];
rz(-1.0815797) q[3];
sx q[3];
rz(-2.1833724) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.4776939) q[2];
sx q[2];
rz(-2.615216) q[2];
sx q[2];
rz(-0.23507512) q[2];
rz(-2.7663686) q[3];
sx q[3];
rz(-0.90685654) q[3];
sx q[3];
rz(2.4231329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5594056) q[0];
sx q[0];
rz(-3.0538054) q[0];
sx q[0];
rz(-2.1413595) q[0];
rz(1.2031215) q[1];
sx q[1];
rz(-1.5088046) q[1];
sx q[1];
rz(0.54471725) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8303568) q[0];
sx q[0];
rz(-2.2831627) q[0];
sx q[0];
rz(1.6300549) q[0];
rz(-pi) q[1];
rz(-2.3339231) q[2];
sx q[2];
rz(-2.2414506) q[2];
sx q[2];
rz(0.18137056) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.9904321) q[1];
sx q[1];
rz(-2.9176788) q[1];
sx q[1];
rz(1.0099221) q[1];
x q[2];
rz(2.8417743) q[3];
sx q[3];
rz(-0.97917367) q[3];
sx q[3];
rz(-0.064379582) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.7228221) q[2];
sx q[2];
rz(-0.79221574) q[2];
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
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61409426) q[0];
sx q[0];
rz(-3.0526563) q[0];
sx q[0];
rz(1.8726789) q[0];
rz(-2.1753963) q[1];
sx q[1];
rz(-1.7485917) q[1];
sx q[1];
rz(-0.46636811) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.185925) q[0];
sx q[0];
rz(-0.67712443) q[0];
sx q[0];
rz(0.71046513) q[0];
rz(-pi) q[1];
rz(0.46196533) q[2];
sx q[2];
rz(-1.8868539) q[2];
sx q[2];
rz(-0.3332999) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.1524017) q[1];
sx q[1];
rz(-2.4394803) q[1];
sx q[1];
rz(-2.0753808) q[1];
rz(1.0859231) q[3];
sx q[3];
rz(-2.4182662) q[3];
sx q[3];
rz(-1.5679899) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7352778) q[2];
sx q[2];
rz(-2.3390528) q[2];
sx q[2];
rz(-0.80424133) q[2];
rz(2.6546226) q[3];
sx q[3];
rz(-2.099497) q[3];
sx q[3];
rz(-0.0074726661) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2020579) q[0];
sx q[0];
rz(-2.845293) q[0];
sx q[0];
rz(0.95788389) q[0];
rz(-1.5127888) q[1];
sx q[1];
rz(-1.0572546) q[1];
sx q[1];
rz(1.588795) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4544425) q[0];
sx q[0];
rz(-1.5936562) q[0];
sx q[0];
rz(-1.6506565) q[0];
rz(-pi) q[1];
x q[1];
rz(1.856316) q[2];
sx q[2];
rz(-1.3034045) q[2];
sx q[2];
rz(-1.7882787) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.342345) q[1];
sx q[1];
rz(-2.4589897) q[1];
sx q[1];
rz(0.28283878) q[1];
rz(-pi) q[2];
rz(0.53258606) q[3];
sx q[3];
rz(-0.69952337) q[3];
sx q[3];
rz(-3.0211849) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.2019041) q[2];
sx q[2];
rz(-2.0857911) q[2];
sx q[2];
rz(2.4064257) q[2];
rz(0.9489263) q[3];
sx q[3];
rz(-2.2831423) q[3];
sx q[3];
rz(0.86400509) q[3];
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
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76148024) q[0];
sx q[0];
rz(-1.004847) q[0];
sx q[0];
rz(-2.5349706) q[0];
rz(-2.3573549) q[1];
sx q[1];
rz(-1.3332858) q[1];
sx q[1];
rz(-1.7971136) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1649496) q[0];
sx q[0];
rz(-2.2047979) q[0];
sx q[0];
rz(-2.1725146) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.44353087) q[2];
sx q[2];
rz(-2.6012528) q[2];
sx q[2];
rz(-0.72393723) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.1320051) q[1];
sx q[1];
rz(-2.6807206) q[1];
sx q[1];
rz(0.51746093) q[1];
rz(-pi) q[2];
rz(0.090769692) q[3];
sx q[3];
rz(-2.0342954) q[3];
sx q[3];
rz(0.32583729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.6091696) q[2];
sx q[2];
rz(-2.6340941) q[2];
sx q[2];
rz(-1.7519105) q[2];
rz(-0.17942795) q[3];
sx q[3];
rz(-0.98589412) q[3];
sx q[3];
rz(0.40670407) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0824025) q[0];
sx q[0];
rz(-0.32408369) q[0];
sx q[0];
rz(0.75505906) q[0];
rz(-0.45626196) q[1];
sx q[1];
rz(-1.7165963) q[1];
sx q[1];
rz(1.0770575) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
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
rz(-2.3577214) q[2];
sx q[2];
rz(-1.4848564) q[2];
sx q[2];
rz(0.29168561) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.39246074) q[1];
sx q[1];
rz(-1.2866255) q[1];
sx q[1];
rz(0.47537132) q[1];
rz(-1.3102688) q[3];
sx q[3];
rz(-0.55580492) q[3];
sx q[3];
rz(2.1422841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.72851744) q[2];
sx q[2];
rz(-1.1001526) q[2];
sx q[2];
rz(0.73823482) q[2];
rz(-1.5089367) q[3];
sx q[3];
rz(-1.3239219) q[3];
sx q[3];
rz(1.9807321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64410925) q[0];
sx q[0];
rz(-1.4583541) q[0];
sx q[0];
rz(0.64895502) q[0];
rz(-0.68967462) q[1];
sx q[1];
rz(-0.70976218) q[1];
sx q[1];
rz(2.7241657) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48106347) q[0];
sx q[0];
rz(-1.821035) q[0];
sx q[0];
rz(0.60926837) q[0];
x q[1];
rz(-1.4944878) q[2];
sx q[2];
rz(-1.2718437) q[2];
sx q[2];
rz(1.1079009) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.74655731) q[1];
sx q[1];
rz(-2.7183924) q[1];
sx q[1];
rz(2.8854135) q[1];
rz(-pi) q[2];
rz(-1.9661918) q[3];
sx q[3];
rz(-1.9428651) q[3];
sx q[3];
rz(-1.7898066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.6680341) q[2];
sx q[2];
rz(-1.4213976) q[2];
sx q[2];
rz(0.045698015) q[2];
rz(-2.8081196) q[3];
sx q[3];
rz(-0.57531753) q[3];
sx q[3];
rz(0.12588178) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68542737) q[0];
sx q[0];
rz(-2.6021155) q[0];
sx q[0];
rz(-1.9239377) q[0];
rz(0.24670163) q[1];
sx q[1];
rz(-1.291899) q[1];
sx q[1];
rz(2.6620679) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8586774) q[0];
sx q[0];
rz(-1.1395795) q[0];
sx q[0];
rz(0.89003508) q[0];
rz(-0.87534753) q[2];
sx q[2];
rz(-1.9672868) q[2];
sx q[2];
rz(-0.04405313) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.4054786) q[1];
sx q[1];
rz(-0.60156721) q[1];
sx q[1];
rz(2.6094237) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.4654309) q[3];
sx q[3];
rz(-1.3579391) q[3];
sx q[3];
rz(-1.928291) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.55629998) q[2];
sx q[2];
rz(-1.8712529) q[2];
sx q[2];
rz(-1.4053819) q[2];
rz(2.4893238) q[3];
sx q[3];
rz(-0.26572078) q[3];
sx q[3];
rz(-2.4238267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2300867) q[0];
sx q[0];
rz(-1.390504) q[0];
sx q[0];
rz(2.6934296) q[0];
rz(-2.3975092) q[1];
sx q[1];
rz(-2.4241445) q[1];
sx q[1];
rz(-1.4345899) q[1];
rz(-0.58341917) q[2];
sx q[2];
rz(-1.4392244) q[2];
sx q[2];
rz(0.34860162) q[2];
rz(-1.0539609) q[3];
sx q[3];
rz(-2.500633) q[3];
sx q[3];
rz(2.5120203) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
