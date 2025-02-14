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
rz(0.024700392) q[0];
sx q[0];
rz(-1.2588809) q[0];
sx q[0];
rz(-1.8688686) q[0];
rz(1.3920353) q[1];
sx q[1];
rz(-1.4104383) q[1];
sx q[1];
rz(-0.8937723) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0451242) q[0];
sx q[0];
rz(-0.91334263) q[0];
sx q[0];
rz(-0.7492926) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8270764) q[2];
sx q[2];
rz(-1.8602784) q[2];
sx q[2];
rz(0.50965259) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.3113678) q[1];
sx q[1];
rz(-1.1347949) q[1];
sx q[1];
rz(0.59200759) q[1];
x q[2];
rz(-1.5529384) q[3];
sx q[3];
rz(-1.2027338) q[3];
sx q[3];
rz(-0.56572711) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.7749403) q[2];
sx q[2];
rz(-1.4619991) q[2];
sx q[2];
rz(-2.6095663) q[2];
rz(-2.4074647) q[3];
sx q[3];
rz(-2.8668154) q[3];
sx q[3];
rz(-1.1067357) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4183913) q[0];
sx q[0];
rz(-0.98576236) q[0];
sx q[0];
rz(-3.0611839) q[0];
rz(0.53781992) q[1];
sx q[1];
rz(-2.0729013) q[1];
sx q[1];
rz(-2.417876) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88392576) q[0];
sx q[0];
rz(-2.2307255) q[0];
sx q[0];
rz(-1.0038478) q[0];
rz(-0.89214708) q[2];
sx q[2];
rz(-0.79024678) q[2];
sx q[2];
rz(0.21871601) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.5840877) q[1];
sx q[1];
rz(-2.1592916) q[1];
sx q[1];
rz(0.001838847) q[1];
x q[2];
rz(-2.3238417) q[3];
sx q[3];
rz(-1.9649385) q[3];
sx q[3];
rz(-0.5241636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2878652) q[2];
sx q[2];
rz(-2.008805) q[2];
sx q[2];
rz(-2.2354324) q[2];
rz(-0.3624889) q[3];
sx q[3];
rz(-1.5972842) q[3];
sx q[3];
rz(3.0947065) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.872181) q[0];
sx q[0];
rz(-1.8159287) q[0];
sx q[0];
rz(2.0500702) q[0];
rz(0.12256924) q[1];
sx q[1];
rz(-1.4361607) q[1];
sx q[1];
rz(1.7950119) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99341644) q[0];
sx q[0];
rz(-1.5412586) q[0];
sx q[0];
rz(1.2146773) q[0];
rz(-pi) q[1];
rz(0.41287072) q[2];
sx q[2];
rz(-2.4390568) q[2];
sx q[2];
rz(-1.2324049) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.0071738) q[1];
sx q[1];
rz(-1.6225623) q[1];
sx q[1];
rz(1.365713) q[1];
rz(1.6835378) q[3];
sx q[3];
rz(-2.097297) q[3];
sx q[3];
rz(-0.91673382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.71318212) q[2];
sx q[2];
rz(-1.1298263) q[2];
sx q[2];
rz(2.909519) q[2];
rz(1.170916) q[3];
sx q[3];
rz(-0.62058312) q[3];
sx q[3];
rz(2.8569729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.93037) q[0];
sx q[0];
rz(-0.79293293) q[0];
sx q[0];
rz(-1.6388182) q[0];
rz(0.59208313) q[1];
sx q[1];
rz(-1.1555669) q[1];
sx q[1];
rz(-0.88964644) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9735721) q[0];
sx q[0];
rz(-2.4102927) q[0];
sx q[0];
rz(2.2941053) q[0];
rz(-pi) q[1];
rz(-2.0060894) q[2];
sx q[2];
rz(-0.5331299) q[2];
sx q[2];
rz(-1.0077623) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.0977323) q[1];
sx q[1];
rz(-1.3907038) q[1];
sx q[1];
rz(-1.8921683) q[1];
rz(1.9806421) q[3];
sx q[3];
rz(-1.92627) q[3];
sx q[3];
rz(0.4996757) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4952116) q[2];
sx q[2];
rz(-1.4633598) q[2];
sx q[2];
rz(2.4656673) q[2];
rz(-1.7050788) q[3];
sx q[3];
rz(-0.33994514) q[3];
sx q[3];
rz(-2.9743312) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2109569) q[0];
sx q[0];
rz(-2.8093331) q[0];
sx q[0];
rz(-2.9462872) q[0];
rz(-3.0209814) q[1];
sx q[1];
rz(-0.21574012) q[1];
sx q[1];
rz(-2.9511071) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57225655) q[0];
sx q[0];
rz(-1.4829206) q[0];
sx q[0];
rz(1.5059307) q[0];
x q[1];
rz(-0.10730524) q[2];
sx q[2];
rz(-0.679418) q[2];
sx q[2];
rz(-0.36919644) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.021564158) q[1];
sx q[1];
rz(-0.11332527) q[1];
sx q[1];
rz(0.47685195) q[1];
rz(-pi) q[2];
rz(1.9917914) q[3];
sx q[3];
rz(-2.1432869) q[3];
sx q[3];
rz(1.5508625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.5087937) q[2];
sx q[2];
rz(-1.4845279) q[2];
sx q[2];
rz(-1.18139) q[2];
rz(1.3440291) q[3];
sx q[3];
rz(-0.7414147) q[3];
sx q[3];
rz(2.369829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33264273) q[0];
sx q[0];
rz(-1.3430261) q[0];
sx q[0];
rz(-1.1395662) q[0];
rz(-0.66967669) q[1];
sx q[1];
rz(-2.2256336) q[1];
sx q[1];
rz(-2.0119827) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2655339) q[0];
sx q[0];
rz(-1.8647412) q[0];
sx q[0];
rz(-2.8493899) q[0];
rz(2.8390719) q[2];
sx q[2];
rz(-0.74069689) q[2];
sx q[2];
rz(-2.4940235) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.672513) q[1];
sx q[1];
rz(-1.7653078) q[1];
sx q[1];
rz(0.37176337) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7288759) q[3];
sx q[3];
rz(-1.3922833) q[3];
sx q[3];
rz(-1.3016537) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.7149522) q[2];
sx q[2];
rz(-1.4098097) q[2];
sx q[2];
rz(2.7362774) q[2];
rz(0.82695588) q[3];
sx q[3];
rz(-0.6260286) q[3];
sx q[3];
rz(-0.85957447) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7760794) q[0];
sx q[0];
rz(-3.1190393) q[0];
sx q[0];
rz(-0.93589163) q[0];
rz(0.86839688) q[1];
sx q[1];
rz(-0.52803841) q[1];
sx q[1];
rz(-1.4194999) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7129684) q[0];
sx q[0];
rz(-2.0789062) q[0];
sx q[0];
rz(1.1626121) q[0];
rz(1.1677891) q[2];
sx q[2];
rz(-2.5984077) q[2];
sx q[2];
rz(-0.78389558) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8384435) q[1];
sx q[1];
rz(-2.0334963) q[1];
sx q[1];
rz(-2.3118603) q[1];
x q[2];
rz(0.0017599864) q[3];
sx q[3];
rz(-1.5704182) q[3];
sx q[3];
rz(3.1217679) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.3203656) q[2];
sx q[2];
rz(-1.7419523) q[2];
sx q[2];
rz(1.0199176) q[2];
rz(0.50926456) q[3];
sx q[3];
rz(-1.7002707) q[3];
sx q[3];
rz(3.063108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7996063) q[0];
sx q[0];
rz(-1.6260363) q[0];
sx q[0];
rz(1.1588143) q[0];
rz(0.47053567) q[1];
sx q[1];
rz(-1.4152941) q[1];
sx q[1];
rz(-1.6758957) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62773317) q[0];
sx q[0];
rz(-3.0854221) q[0];
sx q[0];
rz(-2.7059731) q[0];
rz(-pi) q[1];
rz(-2.7297425) q[2];
sx q[2];
rz(-0.73397103) q[2];
sx q[2];
rz(-0.07871544) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.1701057) q[1];
sx q[1];
rz(-2.4772812) q[1];
sx q[1];
rz(-0.12477915) q[1];
rz(-0.37181446) q[3];
sx q[3];
rz(-0.68228693) q[3];
sx q[3];
rz(2.2390847) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.9965685) q[2];
sx q[2];
rz(-1.327876) q[2];
sx q[2];
rz(-0.81929755) q[2];
rz(-2.3451292) q[3];
sx q[3];
rz(-0.4070681) q[3];
sx q[3];
rz(1.6527294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9374989) q[0];
sx q[0];
rz(-0.096385328) q[0];
sx q[0];
rz(2.3199484) q[0];
rz(2.4687528) q[1];
sx q[1];
rz(-0.57603374) q[1];
sx q[1];
rz(-1.0427262) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.722114) q[0];
sx q[0];
rz(-1.6432149) q[0];
sx q[0];
rz(0.2054604) q[0];
x q[1];
rz(-0.64260245) q[2];
sx q[2];
rz(-2.2589141) q[2];
sx q[2];
rz(0.43660276) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0303749) q[1];
sx q[1];
rz(-1.0008604) q[1];
sx q[1];
rz(-2.6968677) q[1];
x q[2];
rz(-2.6058873) q[3];
sx q[3];
rz(-0.2920449) q[3];
sx q[3];
rz(0.10242187) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.3599856) q[2];
sx q[2];
rz(-1.2722509) q[2];
sx q[2];
rz(-2.2992415) q[2];
rz(0.53884566) q[3];
sx q[3];
rz(-1.4875965) q[3];
sx q[3];
rz(0.52411383) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46776029) q[0];
sx q[0];
rz(-2.9957275) q[0];
sx q[0];
rz(-1.9061331) q[0];
rz(2.1927059) q[1];
sx q[1];
rz(-0.83273879) q[1];
sx q[1];
rz(1.4804776) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7296071) q[0];
sx q[0];
rz(-0.74568891) q[0];
sx q[0];
rz(1.9870158) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3034978) q[2];
sx q[2];
rz(-2.3857255) q[2];
sx q[2];
rz(-2.2264293) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.21232852) q[1];
sx q[1];
rz(-2.9605299) q[1];
sx q[1];
rz(-0.4275112) q[1];
x q[2];
rz(-0.62678316) q[3];
sx q[3];
rz(-2.7434182) q[3];
sx q[3];
rz(-0.45792031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.71198717) q[2];
sx q[2];
rz(-1.7179855) q[2];
sx q[2];
rz(2.831366) q[2];
rz(-2.6817536) q[3];
sx q[3];
rz(-2.583677) q[3];
sx q[3];
rz(-1.3742113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
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
rz(-1.9643758) q[0];
sx q[0];
rz(-1.1342659) q[0];
sx q[0];
rz(-2.076617) q[0];
rz(0.89827697) q[1];
sx q[1];
rz(-2.5986462) q[1];
sx q[1];
rz(2.6959261) q[1];
rz(-2.73478) q[2];
sx q[2];
rz(-0.94712232) q[2];
sx q[2];
rz(0.85489544) q[2];
rz(-1.0887922) q[3];
sx q[3];
rz(-0.72127753) q[3];
sx q[3];
rz(0.75281561) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
