OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-3.5716612) q[0];
sx q[0];
rz(-0.067458955) q[0];
sx q[0];
rz(10.098739) q[0];
rz(-0.31710467) q[1];
sx q[1];
rz(-1.6333406) q[1];
sx q[1];
rz(-0.79467264) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7793517) q[0];
sx q[0];
rz(-1.2138214) q[0];
sx q[0];
rz(-0.39874052) q[0];
rz(0.27172471) q[2];
sx q[2];
rz(-0.50479111) q[2];
sx q[2];
rz(-1.7287849) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.0104048) q[1];
sx q[1];
rz(-2.511575) q[1];
sx q[1];
rz(-0.6063993) q[1];
rz(-3.0029293) q[3];
sx q[3];
rz(-2.0141467) q[3];
sx q[3];
rz(1.5594304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.6661466) q[2];
sx q[2];
rz(-1.1300766) q[2];
sx q[2];
rz(1.3249935) q[2];
rz(-0.31630668) q[3];
sx q[3];
rz(-0.91747147) q[3];
sx q[3];
rz(-2.2345208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9544202) q[0];
sx q[0];
rz(-2.3311054) q[0];
sx q[0];
rz(2.1818838) q[0];
rz(-1.6540487) q[1];
sx q[1];
rz(-2.4580749) q[1];
sx q[1];
rz(-2.5141292) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7213631) q[0];
sx q[0];
rz(-3.1313167) q[0];
sx q[0];
rz(0.72364877) q[0];
rz(-pi) q[1];
rz(-0.21164125) q[2];
sx q[2];
rz(-1.5923119) q[2];
sx q[2];
rz(0.053552901) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1712449) q[1];
sx q[1];
rz(-1.8106318) q[1];
sx q[1];
rz(0.14648267) q[1];
rz(-pi) q[2];
rz(-1.596132) q[3];
sx q[3];
rz(-1.5832033) q[3];
sx q[3];
rz(-0.62904639) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.54005694) q[2];
sx q[2];
rz(-1.7404218) q[2];
sx q[2];
rz(-2.1616948) q[2];
rz(-1.6905789) q[3];
sx q[3];
rz(-0.7106978) q[3];
sx q[3];
rz(1.6987945) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5805805) q[0];
sx q[0];
rz(-3.09364) q[0];
sx q[0];
rz(-1.7618435) q[0];
rz(-0.17669949) q[1];
sx q[1];
rz(-1.2172164) q[1];
sx q[1];
rz(1.7938991) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1130226) q[0];
sx q[0];
rz(-2.0355814) q[0];
sx q[0];
rz(1.3010498) q[0];
x q[1];
rz(1.5930487) q[2];
sx q[2];
rz(-1.77102) q[2];
sx q[2];
rz(-2.3929838) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.1469517) q[1];
sx q[1];
rz(-2.9839582) q[1];
sx q[1];
rz(-2.4528802) q[1];
rz(0.090431902) q[3];
sx q[3];
rz(-1.3867497) q[3];
sx q[3];
rz(-0.66891608) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.0195062) q[2];
sx q[2];
rz(-0.85947376) q[2];
sx q[2];
rz(1.9879831) q[2];
rz(-2.7518318) q[3];
sx q[3];
rz(-0.4699769) q[3];
sx q[3];
rz(2.3526758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47238123) q[0];
sx q[0];
rz(-0.62494576) q[0];
sx q[0];
rz(-2.5770082) q[0];
rz(-0.083104221) q[1];
sx q[1];
rz(-1.9080947) q[1];
sx q[1];
rz(-0.21534236) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1911083) q[0];
sx q[0];
rz(-1.9035625) q[0];
sx q[0];
rz(2.3118408) q[0];
x q[1];
rz(-1.3232857) q[2];
sx q[2];
rz(-0.60192054) q[2];
sx q[2];
rz(-1.9844696) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.5986961) q[1];
sx q[1];
rz(-0.22194949) q[1];
sx q[1];
rz(1.0735682) q[1];
rz(-pi) q[2];
rz(-0.54939778) q[3];
sx q[3];
rz(-1.1664207) q[3];
sx q[3];
rz(-2.3466563) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.6354436) q[2];
sx q[2];
rz(-2.052867) q[2];
sx q[2];
rz(-2.5396458) q[2];
rz(2.4510032) q[3];
sx q[3];
rz(-0.75406051) q[3];
sx q[3];
rz(-0.62600342) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(0.85201207) q[0];
sx q[0];
rz(-2.1471725) q[0];
sx q[0];
rz(2.1767298) q[0];
rz(0.016618641) q[1];
sx q[1];
rz(-2.7162572) q[1];
sx q[1];
rz(0.66666493) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0913038) q[0];
sx q[0];
rz(-1.6717523) q[0];
sx q[0];
rz(-1.1703277) q[0];
rz(-pi) q[1];
rz(-2.7890116) q[2];
sx q[2];
rz(-1.0513628) q[2];
sx q[2];
rz(-0.40371603) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.40068377) q[1];
sx q[1];
rz(-0.92036696) q[1];
sx q[1];
rz(1.9294192) q[1];
rz(-0.22294754) q[3];
sx q[3];
rz(-1.2811536) q[3];
sx q[3];
rz(1.851561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.2255286) q[2];
sx q[2];
rz(-1.3801882) q[2];
sx q[2];
rz(1.2529681) q[2];
rz(-1.3145674) q[3];
sx q[3];
rz(-1.9606749) q[3];
sx q[3];
rz(2.0107646) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1154293) q[0];
sx q[0];
rz(-2.207343) q[0];
sx q[0];
rz(1.3622267) q[0];
rz(-1.7419787) q[1];
sx q[1];
rz(-1.5618961) q[1];
sx q[1];
rz(-1.5302352) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0479193) q[0];
sx q[0];
rz(-1.1512655) q[0];
sx q[0];
rz(-1.490926) q[0];
x q[1];
rz(1.2457232) q[2];
sx q[2];
rz(-1.9893861) q[2];
sx q[2];
rz(-2.3649154) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.1215026) q[1];
sx q[1];
rz(-1.1768747) q[1];
sx q[1];
rz(-2.9083088) q[1];
rz(3.0909782) q[3];
sx q[3];
rz(-1.317987) q[3];
sx q[3];
rz(-1.8265754) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.35649148) q[2];
sx q[2];
rz(-2.1337528) q[2];
sx q[2];
rz(0.029416857) q[2];
rz(-1.4908837) q[3];
sx q[3];
rz(-2.1865032) q[3];
sx q[3];
rz(1.9194549) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4239663) q[0];
sx q[0];
rz(-0.76537776) q[0];
sx q[0];
rz(2.9587865) q[0];
rz(2.706066) q[1];
sx q[1];
rz(-2.2874449) q[1];
sx q[1];
rz(0.51087728) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8144007) q[0];
sx q[0];
rz(-1.4609903) q[0];
sx q[0];
rz(-0.26392428) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3350231) q[2];
sx q[2];
rz(-0.91420805) q[2];
sx q[2];
rz(-0.29733959) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.6888914) q[1];
sx q[1];
rz(-0.93033067) q[1];
sx q[1];
rz(-2.074261) q[1];
x q[2];
rz(1.0501409) q[3];
sx q[3];
rz(-2.5556892) q[3];
sx q[3];
rz(-1.0019765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.2304948) q[2];
sx q[2];
rz(-2.348933) q[2];
sx q[2];
rz(-2.5210209) q[2];
rz(0.42256045) q[3];
sx q[3];
rz(-1.3603323) q[3];
sx q[3];
rz(-1.871199) q[3];
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
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7997416) q[0];
sx q[0];
rz(-0.36062476) q[0];
sx q[0];
rz(-2.0525232) q[0];
rz(2.0607121) q[1];
sx q[1];
rz(-1.2842478) q[1];
sx q[1];
rz(-2.506315) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1875293) q[0];
sx q[0];
rz(-2.7494193) q[0];
sx q[0];
rz(-1.9847045) q[0];
x q[1];
rz(0.17417553) q[2];
sx q[2];
rz(-1.0798228) q[2];
sx q[2];
rz(2.8665286) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.095070187) q[1];
sx q[1];
rz(-2.1696643) q[1];
sx q[1];
rz(1.8938766) q[1];
rz(-1.2653207) q[3];
sx q[3];
rz(-2.5442122) q[3];
sx q[3];
rz(-1.7978096) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.13715956) q[2];
sx q[2];
rz(-1.6798423) q[2];
sx q[2];
rz(-3.1414202) q[2];
rz(0.6862644) q[3];
sx q[3];
rz(-2.4133447) q[3];
sx q[3];
rz(2.2837158) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2962608) q[0];
sx q[0];
rz(-2.2289472) q[0];
sx q[0];
rz(2.4849179) q[0];
rz(2.7776921) q[1];
sx q[1];
rz(-0.44635043) q[1];
sx q[1];
rz(0.90973967) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37659392) q[0];
sx q[0];
rz(-1.24876) q[0];
sx q[0];
rz(1.8565208) q[0];
rz(-pi) q[1];
rz(-1.3526731) q[2];
sx q[2];
rz(-1.2303196) q[2];
sx q[2];
rz(0.097629658) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.4917131) q[1];
sx q[1];
rz(-2.8976106) q[1];
sx q[1];
rz(0.53423832) q[1];
rz(1.2551366) q[3];
sx q[3];
rz(-1.8161895) q[3];
sx q[3];
rz(1.0823702) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.3108814) q[2];
sx q[2];
rz(-0.41822663) q[2];
sx q[2];
rz(0.48661423) q[2];
rz(-3.0330372) q[3];
sx q[3];
rz(-0.71043772) q[3];
sx q[3];
rz(1.9177115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9412823) q[0];
sx q[0];
rz(-1.0815942) q[0];
sx q[0];
rz(0.51668984) q[0];
rz(-1.5453045) q[1];
sx q[1];
rz(-1.0524806) q[1];
sx q[1];
rz(2.2470078) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54297011) q[0];
sx q[0];
rz(-0.93790903) q[0];
sx q[0];
rz(-3.0892239) q[0];
rz(-pi) q[1];
rz(2.9355029) q[2];
sx q[2];
rz(-2.0317674) q[2];
sx q[2];
rz(3.0989625) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.6539508) q[1];
sx q[1];
rz(-2.7751121) q[1];
sx q[1];
rz(1.7483064) q[1];
x q[2];
rz(-1.7747202) q[3];
sx q[3];
rz(-0.66019928) q[3];
sx q[3];
rz(3.077293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.9937925) q[2];
sx q[2];
rz(-2.1788969) q[2];
sx q[2];
rz(1.0220035) q[2];
rz(-1.3379124) q[3];
sx q[3];
rz(-1.6513848) q[3];
sx q[3];
rz(-0.62906229) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5932896) q[0];
sx q[0];
rz(-1.5240482) q[0];
sx q[0];
rz(-0.89292009) q[0];
rz(-2.100636) q[1];
sx q[1];
rz(-0.092408471) q[1];
sx q[1];
rz(1.668781) q[1];
rz(-2.2330877) q[2];
sx q[2];
rz(-1.1237332) q[2];
sx q[2];
rz(-0.087469812) q[2];
rz(1.7265504) q[3];
sx q[3];
rz(-2.8372436) q[3];
sx q[3];
rz(-2.5929034) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
