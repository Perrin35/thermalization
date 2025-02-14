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
rz(0.69831508) q[0];
sx q[0];
rz(-0.3787711) q[0];
sx q[0];
rz(1.1573428) q[0];
rz(10.209822) q[1];
sx q[1];
rz(0.68612376) q[1];
sx q[1];
rz(2.6148028) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8136381) q[0];
sx q[0];
rz(-1.0109954) q[0];
sx q[0];
rz(0.35388057) q[0];
rz(-0.73424642) q[2];
sx q[2];
rz(-0.99319211) q[2];
sx q[2];
rz(-0.89370382) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.4545645) q[1];
sx q[1];
rz(-2.6978081) q[1];
sx q[1];
rz(3.0566178) q[1];
rz(-pi) q[2];
rz(3.111638) q[3];
sx q[3];
rz(-2.1308793) q[3];
sx q[3];
rz(-2.4966979) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.36246768) q[2];
sx q[2];
rz(-0.96131009) q[2];
sx q[2];
rz(-0.15164068) q[2];
rz(-0.14196299) q[3];
sx q[3];
rz(-1.4700593) q[3];
sx q[3];
rz(-1.1375455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66459429) q[0];
sx q[0];
rz(-2.4095896) q[0];
sx q[0];
rz(2.3240996) q[0];
rz(0.11257653) q[1];
sx q[1];
rz(-1.45603) q[1];
sx q[1];
rz(-0.89961019) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0257638) q[0];
sx q[0];
rz(-1.2983822) q[0];
sx q[0];
rz(-1.3221571) q[0];
x q[1];
rz(-2.1915273) q[2];
sx q[2];
rz(-2.143444) q[2];
sx q[2];
rz(3.0814296) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.039520415) q[1];
sx q[1];
rz(-1.6270492) q[1];
sx q[1];
rz(0.7489167) q[1];
rz(-pi) q[2];
rz(-2.2475776) q[3];
sx q[3];
rz(-1.1543659) q[3];
sx q[3];
rz(-2.3125966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.67109913) q[2];
sx q[2];
rz(-1.495139) q[2];
sx q[2];
rz(-1.0279083) q[2];
rz(0.73728621) q[3];
sx q[3];
rz(-1.6643915) q[3];
sx q[3];
rz(-0.024287311) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0055493) q[0];
sx q[0];
rz(-2.5856954) q[0];
sx q[0];
rz(-1.9480202) q[0];
rz(-1.2981124) q[1];
sx q[1];
rz(-1.4793414) q[1];
sx q[1];
rz(1.4847635) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1910716) q[0];
sx q[0];
rz(-1.6170039) q[0];
sx q[0];
rz(0.025006983) q[0];
x q[1];
rz(1.3593) q[2];
sx q[2];
rz(-1.4471608) q[2];
sx q[2];
rz(-0.57963003) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.72276593) q[1];
sx q[1];
rz(-1.2294985) q[1];
sx q[1];
rz(-1.4272293) q[1];
rz(1.9336352) q[3];
sx q[3];
rz(-1.4936183) q[3];
sx q[3];
rz(3.0953323) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.9803311) q[2];
sx q[2];
rz(-1.9183466) q[2];
sx q[2];
rz(-0.70183357) q[2];
rz(-3.0724604) q[3];
sx q[3];
rz(-2.2528503) q[3];
sx q[3];
rz(2.4338636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0858916) q[0];
sx q[0];
rz(-2.0022855) q[0];
sx q[0];
rz(2.5162146) q[0];
rz(2.0471795) q[1];
sx q[1];
rz(-1.5233327) q[1];
sx q[1];
rz(-1.3166924) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1392446) q[0];
sx q[0];
rz(-1.7138094) q[0];
sx q[0];
rz(2.1190179) q[0];
rz(-1.473963) q[2];
sx q[2];
rz(-1.8928573) q[2];
sx q[2];
rz(-1.9063973) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.55344501) q[1];
sx q[1];
rz(-0.7281174) q[1];
sx q[1];
rz(-0.27256984) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5685124) q[3];
sx q[3];
rz(-0.88643688) q[3];
sx q[3];
rz(2.7671368) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.40654287) q[2];
sx q[2];
rz(-0.64261618) q[2];
sx q[2];
rz(3.1114846) q[2];
rz(-0.41664577) q[3];
sx q[3];
rz(-1.7130339) q[3];
sx q[3];
rz(-0.055195181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77867126) q[0];
sx q[0];
rz(-1.8632977) q[0];
sx q[0];
rz(-1.2177421) q[0];
rz(0.22449224) q[1];
sx q[1];
rz(-1.4935363) q[1];
sx q[1];
rz(-0.40103689) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7571791) q[0];
sx q[0];
rz(-1.9862439) q[0];
sx q[0];
rz(-2.7336304) q[0];
rz(-pi) q[1];
x q[1];
rz(0.65481477) q[2];
sx q[2];
rz(-1.8698955) q[2];
sx q[2];
rz(1.6950184) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.68919888) q[1];
sx q[1];
rz(-1.3340063) q[1];
sx q[1];
rz(1.522032) q[1];
rz(1.6717853) q[3];
sx q[3];
rz(-1.4091531) q[3];
sx q[3];
rz(0.18835959) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.2991221) q[2];
sx q[2];
rz(-1.9186019) q[2];
sx q[2];
rz(0.50745884) q[2];
rz(0.92042813) q[3];
sx q[3];
rz(-2.7046552) q[3];
sx q[3];
rz(-0.18032716) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
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
rz(-2.4949263) q[0];
sx q[0];
rz(-3.086402) q[0];
sx q[0];
rz(-1.0194417) q[0];
rz(1.1722209) q[1];
sx q[1];
rz(-1.9688789) q[1];
sx q[1];
rz(-1.9196462) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2766314) q[0];
sx q[0];
rz(-0.91827938) q[0];
sx q[0];
rz(-0.44256532) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6331633) q[2];
sx q[2];
rz(-0.51566873) q[2];
sx q[2];
rz(-0.0002278681) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.5355365) q[1];
sx q[1];
rz(-2.5185985) q[1];
sx q[1];
rz(0.48320233) q[1];
rz(-2.0133408) q[3];
sx q[3];
rz(-2.5199515) q[3];
sx q[3];
rz(0.78237247) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.4814066) q[2];
sx q[2];
rz(-0.6654827) q[2];
sx q[2];
rz(-1.0931724) q[2];
rz(2.6045351) q[3];
sx q[3];
rz(-1.4635181) q[3];
sx q[3];
rz(-0.66948906) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90401232) q[0];
sx q[0];
rz(-2.8222988) q[0];
sx q[0];
rz(-1.681666) q[0];
rz(-1.6319252) q[1];
sx q[1];
rz(-0.62848148) q[1];
sx q[1];
rz(-3.0622838) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58188841) q[0];
sx q[0];
rz(-1.4900582) q[0];
sx q[0];
rz(1.8084007) q[0];
rz(-2.9379876) q[2];
sx q[2];
rz(-1.7053268) q[2];
sx q[2];
rz(2.8367538) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.7497158) q[1];
sx q[1];
rz(-1.8444711) q[1];
sx q[1];
rz(2.7065212) q[1];
x q[2];
rz(-2.4238911) q[3];
sx q[3];
rz(-1.4677748) q[3];
sx q[3];
rz(2.3915402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.0038393) q[2];
sx q[2];
rz(-1.9766108) q[2];
sx q[2];
rz(-0.4129146) q[2];
rz(1.2952992) q[3];
sx q[3];
rz(-2.2173939) q[3];
sx q[3];
rz(-2.7220461) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15692784) q[0];
sx q[0];
rz(-0.65021896) q[0];
sx q[0];
rz(0.28537634) q[0];
rz(-0.05263075) q[1];
sx q[1];
rz(-1.7335408) q[1];
sx q[1];
rz(0.1836798) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.015332) q[0];
sx q[0];
rz(-1.53526) q[0];
sx q[0];
rz(1.4592341) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6551774) q[2];
sx q[2];
rz(-0.32054311) q[2];
sx q[2];
rz(2.965791) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.3425828) q[1];
sx q[1];
rz(-1.6470121) q[1];
sx q[1];
rz(0.59554312) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7643572) q[3];
sx q[3];
rz(-2.4580015) q[3];
sx q[3];
rz(2.205276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.24667428) q[2];
sx q[2];
rz(-1.4072714) q[2];
sx q[2];
rz(2.0764009) q[2];
rz(-1.6861606) q[3];
sx q[3];
rz(-1.8543517) q[3];
sx q[3];
rz(0.58568946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0378549) q[0];
sx q[0];
rz(-2.3268564) q[0];
sx q[0];
rz(-1.5392342) q[0];
rz(-0.94929758) q[1];
sx q[1];
rz(-2.0930591) q[1];
sx q[1];
rz(1.4303713) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2968144) q[0];
sx q[0];
rz(-2.0955293) q[0];
sx q[0];
rz(2.120976) q[0];
x q[1];
rz(1.4865506) q[2];
sx q[2];
rz(-2.4110594) q[2];
sx q[2];
rz(0.73405311) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.6257244) q[1];
sx q[1];
rz(-1.6004381) q[1];
sx q[1];
rz(0.29582204) q[1];
rz(-pi) q[2];
x q[2];
rz(0.28258459) q[3];
sx q[3];
rz(-1.1573829) q[3];
sx q[3];
rz(1.9081209) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.1324233) q[2];
sx q[2];
rz(-1.3045661) q[2];
sx q[2];
rz(-0.12651786) q[2];
rz(-2.8070519) q[3];
sx q[3];
rz(-0.25944513) q[3];
sx q[3];
rz(-0.87219605) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8406521) q[0];
sx q[0];
rz(-0.88467389) q[0];
sx q[0];
rz(0.51710039) q[0];
rz(-3.0866947) q[1];
sx q[1];
rz(-1.6322735) q[1];
sx q[1];
rz(0.089769207) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9909279) q[0];
sx q[0];
rz(-2.1164843) q[0];
sx q[0];
rz(2.664458) q[0];
rz(-pi) q[1];
rz(-1.6705065) q[2];
sx q[2];
rz(-2.5784628) q[2];
sx q[2];
rz(2.827284) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.0193737) q[1];
sx q[1];
rz(-1.3189338) q[1];
sx q[1];
rz(-1.4068961) q[1];
x q[2];
rz(0.2281727) q[3];
sx q[3];
rz(-1.2163278) q[3];
sx q[3];
rz(-1.4972403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.53139293) q[2];
sx q[2];
rz(-1.1343845) q[2];
sx q[2];
rz(2.4143977) q[2];
rz(2.3860892) q[3];
sx q[3];
rz(-0.38755363) q[3];
sx q[3];
rz(0.21272794) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.94201921) q[0];
sx q[0];
rz(-1.6006391) q[0];
sx q[0];
rz(-2.0508456) q[0];
rz(1.749281) q[1];
sx q[1];
rz(-1.5131469) q[1];
sx q[1];
rz(-1.6319235) q[1];
rz(-2.2752599) q[2];
sx q[2];
rz(-1.6523747) q[2];
sx q[2];
rz(-2.0518641) q[2];
rz(-0.94115067) q[3];
sx q[3];
rz(-2.0654021) q[3];
sx q[3];
rz(0.13765814) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
