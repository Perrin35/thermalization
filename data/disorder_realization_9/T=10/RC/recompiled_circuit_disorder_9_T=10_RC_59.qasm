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
rz(2.824488) q[1];
sx q[1];
rz(-1.5082521) q[1];
sx q[1];
rz(-2.34692) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90097839) q[0];
sx q[0];
rz(-0.52871791) q[0];
sx q[0];
rz(0.76529495) q[0];
x q[1];
rz(2.8698679) q[2];
sx q[2];
rz(-2.6368015) q[2];
sx q[2];
rz(1.4128078) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0104048) q[1];
sx q[1];
rz(-0.63001761) q[1];
sx q[1];
rz(0.6063993) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2875597) q[3];
sx q[3];
rz(-0.4631511) q[3];
sx q[3];
rz(-1.2676261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.47544605) q[2];
sx q[2];
rz(-1.1300766) q[2];
sx q[2];
rz(1.8165992) q[2];
rz(-0.31630668) q[3];
sx q[3];
rz(-2.2241212) q[3];
sx q[3];
rz(2.2345208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1871724) q[0];
sx q[0];
rz(-2.3311054) q[0];
sx q[0];
rz(-0.95970884) q[0];
rz(1.487544) q[1];
sx q[1];
rz(-2.4580749) q[1];
sx q[1];
rz(-2.5141292) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7146485) q[0];
sx q[0];
rz(-1.5776002) q[0];
sx q[0];
rz(0.0077008458) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.10208315) q[2];
sx q[2];
rz(-2.9288769) q[2];
sx q[2];
rz(1.724147) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.5761766) q[1];
sx q[1];
rz(-1.4285354) q[1];
sx q[1];
rz(1.3284645) q[1];
rz(-pi) q[2];
rz(-2.0262358) q[3];
sx q[3];
rz(-3.1133828) q[3];
sx q[3];
rz(-1.7445604) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.54005694) q[2];
sx q[2];
rz(-1.4011708) q[2];
sx q[2];
rz(2.1616948) q[2];
rz(-1.6905789) q[3];
sx q[3];
rz(-0.7106978) q[3];
sx q[3];
rz(1.6987945) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5805805) q[0];
sx q[0];
rz(-3.09364) q[0];
sx q[0];
rz(1.7618435) q[0];
rz(-2.9648932) q[1];
sx q[1];
rz(-1.9243762) q[1];
sx q[1];
rz(1.7938991) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4189258) q[0];
sx q[0];
rz(-1.8113266) q[0];
sx q[0];
rz(-0.47970432) q[0];
rz(-pi) q[1];
x q[1];
rz(1.548544) q[2];
sx q[2];
rz(-1.3705727) q[2];
sx q[2];
rz(-2.3929838) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.8417931) q[1];
sx q[1];
rz(-1.6922957) q[1];
sx q[1];
rz(-1.6714765) q[1];
rz(-pi) q[2];
rz(3.0511608) q[3];
sx q[3];
rz(-1.3867497) q[3];
sx q[3];
rz(0.66891608) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.1220864) q[2];
sx q[2];
rz(-2.2821189) q[2];
sx q[2];
rz(1.9879831) q[2];
rz(-0.38976088) q[3];
sx q[3];
rz(-0.4699769) q[3];
sx q[3];
rz(-2.3526758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47238123) q[0];
sx q[0];
rz(-0.62494576) q[0];
sx q[0];
rz(-2.5770082) q[0];
rz(3.0584884) q[1];
sx q[1];
rz(-1.233498) q[1];
sx q[1];
rz(-2.9262503) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72258654) q[0];
sx q[0];
rz(-0.79918396) q[0];
sx q[0];
rz(-1.0976085) q[0];
x q[1];
rz(-1.3232857) q[2];
sx q[2];
rz(-0.60192054) q[2];
sx q[2];
rz(1.1571231) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.5147869) q[1];
sx q[1];
rz(-1.6759911) q[1];
sx q[1];
rz(1.3749966) q[1];
x q[2];
rz(0.54939778) q[3];
sx q[3];
rz(-1.975172) q[3];
sx q[3];
rz(-2.3466563) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.6354436) q[2];
sx q[2];
rz(-2.052867) q[2];
sx q[2];
rz(2.5396458) q[2];
rz(-0.69058949) q[3];
sx q[3];
rz(-2.3875321) q[3];
sx q[3];
rz(-2.5155892) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2895806) q[0];
sx q[0];
rz(-2.1471725) q[0];
sx q[0];
rz(-2.1767298) q[0];
rz(3.124974) q[1];
sx q[1];
rz(-2.7162572) q[1];
sx q[1];
rz(2.4749277) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.704741) q[0];
sx q[0];
rz(-1.1724823) q[0];
sx q[0];
rz(-3.0320291) q[0];
x q[1];
rz(-0.35258106) q[2];
sx q[2];
rz(-1.0513628) q[2];
sx q[2];
rz(0.40371603) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.95496817) q[1];
sx q[1];
rz(-0.72997) q[1];
sx q[1];
rz(0.43220046) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9186451) q[3];
sx q[3];
rz(-1.2811536) q[3];
sx q[3];
rz(-1.2900316) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.2255286) q[2];
sx q[2];
rz(-1.3801882) q[2];
sx q[2];
rz(-1.8886245) q[2];
rz(-1.3145674) q[3];
sx q[3];
rz(-1.1809177) q[3];
sx q[3];
rz(-2.0107646) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0261633) q[0];
sx q[0];
rz(-2.207343) q[0];
sx q[0];
rz(-1.3622267) q[0];
rz(1.399614) q[1];
sx q[1];
rz(-1.5618961) q[1];
sx q[1];
rz(-1.5302352) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2419469) q[0];
sx q[0];
rz(-2.7149704) q[0];
sx q[0];
rz(-0.17701478) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.62263454) q[2];
sx q[2];
rz(-2.6175675) q[2];
sx q[2];
rz(3.0571836) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.64165243) q[1];
sx q[1];
rz(-1.3556726) q[1];
sx q[1];
rz(1.9745449) q[1];
rz(-0.050614428) q[3];
sx q[3];
rz(-1.317987) q[3];
sx q[3];
rz(1.3150172) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.35649148) q[2];
sx q[2];
rz(-1.0078398) q[2];
sx q[2];
rz(0.029416857) q[2];
rz(-1.4908837) q[3];
sx q[3];
rz(-0.95508948) q[3];
sx q[3];
rz(1.2221378) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71762639) q[0];
sx q[0];
rz(-2.3762149) q[0];
sx q[0];
rz(-0.18280612) q[0];
rz(-0.43552661) q[1];
sx q[1];
rz(-0.85414779) q[1];
sx q[1];
rz(2.6307154) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8683838) q[0];
sx q[0];
rz(-1.8330935) q[0];
sx q[0];
rz(-1.457085) q[0];
rz(-pi) q[1];
rz(-0.81803825) q[2];
sx q[2];
rz(-2.15089) q[2];
sx q[2];
rz(1.8028508) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.4335945) q[1];
sx q[1];
rz(-0.79213789) q[1];
sx q[1];
rz(0.57452332) q[1];
rz(-pi) q[2];
rz(-2.0931582) q[3];
sx q[3];
rz(-1.8494542) q[3];
sx q[3];
rz(0.1230965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.2304948) q[2];
sx q[2];
rz(-2.348933) q[2];
sx q[2];
rz(2.5210209) q[2];
rz(0.42256045) q[3];
sx q[3];
rz(-1.3603323) q[3];
sx q[3];
rz(1.2703936) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34185103) q[0];
sx q[0];
rz(-2.7809679) q[0];
sx q[0];
rz(2.0525232) q[0];
rz(2.0607121) q[1];
sx q[1];
rz(-1.2842478) q[1];
sx q[1];
rz(0.63527766) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95406336) q[0];
sx q[0];
rz(-2.7494193) q[0];
sx q[0];
rz(-1.1568882) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.073457) q[2];
sx q[2];
rz(-1.4173696) q[2];
sx q[2];
rz(-1.2129601) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.852408) q[1];
sx q[1];
rz(-1.305456) q[1];
sx q[1];
rz(2.5177588) q[1];
rz(2.9397804) q[3];
sx q[3];
rz(-1.0045856) q[3];
sx q[3];
rz(1.4334397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0044331) q[2];
sx q[2];
rz(-1.6798423) q[2];
sx q[2];
rz(0.00017246406) q[2];
rz(2.4553283) q[3];
sx q[3];
rz(-2.4133447) q[3];
sx q[3];
rz(-2.2837158) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2962608) q[0];
sx q[0];
rz(-0.91264549) q[0];
sx q[0];
rz(0.6566748) q[0];
rz(-0.36390057) q[1];
sx q[1];
rz(-0.44635043) q[1];
sx q[1];
rz(0.90973967) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37659392) q[0];
sx q[0];
rz(-1.24876) q[0];
sx q[0];
rz(-1.8565208) q[0];
x q[1];
rz(1.3526731) q[2];
sx q[2];
rz(-1.911273) q[2];
sx q[2];
rz(0.097629658) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.197387) q[1];
sx q[1];
rz(-1.3613609) q[1];
sx q[1];
rz(-1.6968813) q[1];
x q[2];
rz(2.2496201) q[3];
sx q[3];
rz(-2.744305) q[3];
sx q[3];
rz(2.9904423) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.83071128) q[2];
sx q[2];
rz(-0.41822663) q[2];
sx q[2];
rz(0.48661423) q[2];
rz(-0.1085554) q[3];
sx q[3];
rz(-0.71043772) q[3];
sx q[3];
rz(-1.9177115) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20031032) q[0];
sx q[0];
rz(-2.0599984) q[0];
sx q[0];
rz(0.51668984) q[0];
rz(1.5453045) q[1];
sx q[1];
rz(-2.089112) q[1];
sx q[1];
rz(2.2470078) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45458083) q[0];
sx q[0];
rz(-0.63475383) q[0];
sx q[0];
rz(1.4995585) q[0];
x q[1];
rz(2.9355029) q[2];
sx q[2];
rz(-2.0317674) q[2];
sx q[2];
rz(-0.042630171) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.67748653) q[1];
sx q[1];
rz(-1.9312526) q[1];
sx q[1];
rz(0.067671138) q[1];
rz(-2.9856332) q[3];
sx q[3];
rz(-0.92658639) q[3];
sx q[3];
rz(-2.9498266) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.14780012) q[2];
sx q[2];
rz(-2.1788969) q[2];
sx q[2];
rz(2.1195892) q[2];
rz(1.8036802) q[3];
sx q[3];
rz(-1.6513848) q[3];
sx q[3];
rz(2.5125304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.548303) q[0];
sx q[0];
rz(-1.6175445) q[0];
sx q[0];
rz(2.2486726) q[0];
rz(-2.100636) q[1];
sx q[1];
rz(-0.092408471) q[1];
sx q[1];
rz(1.668781) q[1];
rz(2.2330877) q[2];
sx q[2];
rz(-2.0178595) q[2];
sx q[2];
rz(3.0541228) q[2];
rz(-1.7265504) q[3];
sx q[3];
rz(-0.30434904) q[3];
sx q[3];
rz(0.54868922) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
