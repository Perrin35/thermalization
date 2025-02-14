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
rz(-2.1190116) q[0];
sx q[0];
rz(-2.1562161) q[0];
sx q[0];
rz(-1.1126385) q[0];
rz(-2.7467709) q[1];
sx q[1];
rz(-2.0070751) q[1];
sx q[1];
rz(-1.9505824) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9562839) q[0];
sx q[0];
rz(-0.6260159) q[0];
sx q[0];
rz(2.6176388) q[0];
rz(-pi) q[1];
rz(1.8687227) q[2];
sx q[2];
rz(-0.42537826) q[2];
sx q[2];
rz(0.59051248) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.0886695) q[1];
sx q[1];
rz(-3.0471339) q[1];
sx q[1];
rz(-0.7868305) q[1];
rz(-pi) q[2];
rz(2.7026579) q[3];
sx q[3];
rz(-1.1263444) q[3];
sx q[3];
rz(1.4107454) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.4224008) q[2];
sx q[2];
rz(-0.16133186) q[2];
sx q[2];
rz(-0.42616978) q[2];
rz(2.4424477) q[3];
sx q[3];
rz(-1.1933425) q[3];
sx q[3];
rz(-0.41081158) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8461269) q[0];
sx q[0];
rz(-2.182425) q[0];
sx q[0];
rz(0.83569431) q[0];
rz(0.524638) q[1];
sx q[1];
rz(-1.6181889) q[1];
sx q[1];
rz(1.5118648) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.152911) q[0];
sx q[0];
rz(-2.464474) q[0];
sx q[0];
rz(-2.0858913) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2316547) q[2];
sx q[2];
rz(-0.99915394) q[2];
sx q[2];
rz(2.2425368) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.32538336) q[1];
sx q[1];
rz(-2.0289501) q[1];
sx q[1];
rz(1.7857331) q[1];
rz(-pi) q[2];
rz(0.55142656) q[3];
sx q[3];
rz(-1.2987483) q[3];
sx q[3];
rz(-0.37032933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.52874804) q[2];
sx q[2];
rz(-1.560644) q[2];
sx q[2];
rz(0.35231248) q[2];
rz(-2.6490372) q[3];
sx q[3];
rz(-2.0217321) q[3];
sx q[3];
rz(-0.39920863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6414129) q[0];
sx q[0];
rz(-0.68920511) q[0];
sx q[0];
rz(2.8885762) q[0];
rz(-2.6192656) q[1];
sx q[1];
rz(-0.42799196) q[1];
sx q[1];
rz(-2.7740251) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51742348) q[0];
sx q[0];
rz(-1.6332024) q[0];
sx q[0];
rz(-3.1145658) q[0];
x q[1];
rz(2.1913826) q[2];
sx q[2];
rz(-1.2428811) q[2];
sx q[2];
rz(2.7958946) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.8525703) q[1];
sx q[1];
rz(-1.7779568) q[1];
sx q[1];
rz(3.0604355) q[1];
x q[2];
rz(3.0250164) q[3];
sx q[3];
rz(-1.7128782) q[3];
sx q[3];
rz(0.97654479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.14105178) q[2];
sx q[2];
rz(-0.29033574) q[2];
sx q[2];
rz(-0.42624897) q[2];
rz(-2.352377) q[3];
sx q[3];
rz(-1.1169249) q[3];
sx q[3];
rz(1.6023908) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57192794) q[0];
sx q[0];
rz(-0.33721384) q[0];
sx q[0];
rz(2.9544882) q[0];
rz(1.3693753) q[1];
sx q[1];
rz(-1.2792842) q[1];
sx q[1];
rz(-0.61417907) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3423796) q[0];
sx q[0];
rz(-0.49509096) q[0];
sx q[0];
rz(1.449388) q[0];
rz(-2.1151727) q[2];
sx q[2];
rz(-1.8649057) q[2];
sx q[2];
rz(-1.1687129) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.5803197) q[1];
sx q[1];
rz(-0.6918219) q[1];
sx q[1];
rz(-1.1279593) q[1];
rz(-pi) q[2];
rz(-1.6574347) q[3];
sx q[3];
rz(-1.9443439) q[3];
sx q[3];
rz(-3.0141351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.46089178) q[2];
sx q[2];
rz(-2.2368175) q[2];
sx q[2];
rz(-1.8756078) q[2];
rz(-0.2868109) q[3];
sx q[3];
rz(-2.3293142) q[3];
sx q[3];
rz(-3.0205145) q[3];
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
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1665961) q[0];
sx q[0];
rz(-2.2280405) q[0];
sx q[0];
rz(1.9972557) q[0];
rz(2.3727349) q[1];
sx q[1];
rz(-2.2298593) q[1];
sx q[1];
rz(-1.2164046) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74690565) q[0];
sx q[0];
rz(-1.1228859) q[0];
sx q[0];
rz(-2.7120717) q[0];
rz(-3.0450174) q[2];
sx q[2];
rz(-1.7032832) q[2];
sx q[2];
rz(-0.16667067) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.7256355) q[1];
sx q[1];
rz(-2.8548988) q[1];
sx q[1];
rz(1.2691203) q[1];
rz(-pi) q[2];
x q[2];
rz(0.22138314) q[3];
sx q[3];
rz(-0.37015823) q[3];
sx q[3];
rz(-2.5551318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.52648175) q[2];
sx q[2];
rz(-2.5418126) q[2];
sx q[2];
rz(-0.46670023) q[2];
rz(-2.1488721) q[3];
sx q[3];
rz(-1.4027184) q[3];
sx q[3];
rz(-1.7464975) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4863131) q[0];
sx q[0];
rz(-1.7236973) q[0];
sx q[0];
rz(-0.41743761) q[0];
rz(0.4256658) q[1];
sx q[1];
rz(-0.8368496) q[1];
sx q[1];
rz(-0.72798896) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8288475) q[0];
sx q[0];
rz(-2.0916601) q[0];
sx q[0];
rz(-2.62045) q[0];
rz(0.73982088) q[2];
sx q[2];
rz(-1.707381) q[2];
sx q[2];
rz(-0.53892577) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.14236808) q[1];
sx q[1];
rz(-1.8538239) q[1];
sx q[1];
rz(-1.9623161) q[1];
x q[2];
rz(-0.40188222) q[3];
sx q[3];
rz(-1.0392351) q[3];
sx q[3];
rz(1.2653923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.41842469) q[2];
sx q[2];
rz(-1.4801414) q[2];
sx q[2];
rz(-1.9815014) q[2];
rz(-0.55231071) q[3];
sx q[3];
rz(-2.4513125) q[3];
sx q[3];
rz(2.6753329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0277249) q[0];
sx q[0];
rz(-0.92806569) q[0];
sx q[0];
rz(-3.0357251) q[0];
rz(-1.8766807) q[1];
sx q[1];
rz(-2.5906339) q[1];
sx q[1];
rz(1.6884621) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47521771) q[0];
sx q[0];
rz(-1.2547818) q[0];
sx q[0];
rz(-2.3316335) q[0];
rz(1.3927473) q[2];
sx q[2];
rz(-1.6916923) q[2];
sx q[2];
rz(2.4089782) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.65583751) q[1];
sx q[1];
rz(-1.7110516) q[1];
sx q[1];
rz(1.8717688) q[1];
rz(-2.9960521) q[3];
sx q[3];
rz(-0.91314935) q[3];
sx q[3];
rz(1.6931134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.10446163) q[2];
sx q[2];
rz(-2.0115325) q[2];
sx q[2];
rz(-1.5360443) q[2];
rz(-1.3691085) q[3];
sx q[3];
rz(-0.57309279) q[3];
sx q[3];
rz(-1.8039186) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29984125) q[0];
sx q[0];
rz(-1.2513237) q[0];
sx q[0];
rz(-1.1773671) q[0];
rz(-1.8384701) q[1];
sx q[1];
rz(-1.6603371) q[1];
sx q[1];
rz(0.67828137) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4561397) q[0];
sx q[0];
rz(-1.8243921) q[0];
sx q[0];
rz(2.0965791) q[0];
rz(-1.6839333) q[2];
sx q[2];
rz(-1.8461421) q[2];
sx q[2];
rz(-2.6983698) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.57215124) q[1];
sx q[1];
rz(-1.8742233) q[1];
sx q[1];
rz(-1.7904758) q[1];
rz(-0.60538624) q[3];
sx q[3];
rz(-1.5380713) q[3];
sx q[3];
rz(2.721019) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9699041) q[2];
sx q[2];
rz(-0.50953141) q[2];
sx q[2];
rz(0.92863885) q[2];
rz(1.6396133) q[3];
sx q[3];
rz(-2.00311) q[3];
sx q[3];
rz(-1.5665215) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-2.0314727) q[0];
sx q[0];
rz(-2.0125772) q[0];
sx q[0];
rz(1.336115) q[0];
rz(-2.8307092) q[1];
sx q[1];
rz(-1.5300405) q[1];
sx q[1];
rz(1.8119887) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7696849) q[0];
sx q[0];
rz(-1.9084255) q[0];
sx q[0];
rz(-0.67358526) q[0];
rz(-pi) q[1];
x q[1];
rz(0.86756687) q[2];
sx q[2];
rz(-2.0793781) q[2];
sx q[2];
rz(0.17365467) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.80881608) q[1];
sx q[1];
rz(-3.0114698) q[1];
sx q[1];
rz(2.6511741) q[1];
x q[2];
rz(2.2317367) q[3];
sx q[3];
rz(-1.5136711) q[3];
sx q[3];
rz(1.2965681) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.10799321) q[2];
sx q[2];
rz(-0.42069837) q[2];
sx q[2];
rz(-1.1619953) q[2];
rz(-1.9977995) q[3];
sx q[3];
rz(-1.1868718) q[3];
sx q[3];
rz(-2.7403045) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4556295) q[0];
sx q[0];
rz(-2.3271053) q[0];
sx q[0];
rz(2.234835) q[0];
rz(-0.39457679) q[1];
sx q[1];
rz(-2.5386609) q[1];
sx q[1];
rz(-1.1366064) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6429657) q[0];
sx q[0];
rz(-1.8436448) q[0];
sx q[0];
rz(-2.8754995) q[0];
rz(-pi) q[1];
rz(-0.074176057) q[2];
sx q[2];
rz(-1.8674506) q[2];
sx q[2];
rz(-2.3919472) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.26173466) q[1];
sx q[1];
rz(-2.2764808) q[1];
sx q[1];
rz(2.9565587) q[1];
rz(1.2796938) q[3];
sx q[3];
rz(-1.8661235) q[3];
sx q[3];
rz(0.81084033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-3.0236728) q[2];
sx q[2];
rz(-1.4503786) q[2];
sx q[2];
rz(1.0517906) q[2];
rz(-3.1045095) q[3];
sx q[3];
rz(-0.6784234) q[3];
sx q[3];
rz(-1.7557433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7294075) q[0];
sx q[0];
rz(-0.89898983) q[0];
sx q[0];
rz(0.40869024) q[0];
rz(-2.9858934) q[1];
sx q[1];
rz(-1.0380048) q[1];
sx q[1];
rz(0.17539594) q[1];
rz(2.2536106) q[2];
sx q[2];
rz(-1.8941634) q[2];
sx q[2];
rz(-1.8695199) q[2];
rz(-0.59378271) q[3];
sx q[3];
rz(-1.8206222) q[3];
sx q[3];
rz(0.016236246) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
