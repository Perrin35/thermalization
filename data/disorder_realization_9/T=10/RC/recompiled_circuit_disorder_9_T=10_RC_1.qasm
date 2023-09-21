OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.7115241) q[0];
sx q[0];
rz(-0.067458955) q[0];
sx q[0];
rz(0.67396069) q[0];
rz(2.824488) q[1];
sx q[1];
rz(-1.5082521) q[1];
sx q[1];
rz(-2.34692) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0792102) q[0];
sx q[0];
rz(-1.9431207) q[0];
sx q[0];
rz(1.9553493) q[0];
x q[1];
rz(-2.8698679) q[2];
sx q[2];
rz(-2.6368015) q[2];
sx q[2];
rz(1.7287849) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.0104048) q[1];
sx q[1];
rz(-0.63001761) q[1];
sx q[1];
rz(-0.6063993) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0029293) q[3];
sx q[3];
rz(-2.0141467) q[3];
sx q[3];
rz(-1.5821622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.47544605) q[2];
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
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9544202) q[0];
sx q[0];
rz(-2.3311054) q[0];
sx q[0];
rz(2.1818838) q[0];
rz(1.487544) q[1];
sx q[1];
rz(-0.68351775) q[1];
sx q[1];
rz(-0.62746343) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9976881) q[0];
sx q[0];
rz(-1.578497) q[0];
sx q[0];
rz(-1.5639923) q[0];
rz(-pi) q[1];
rz(1.5928028) q[2];
sx q[2];
rz(-1.3592048) q[2];
sx q[2];
rz(-1.6197268) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.5761766) q[1];
sx q[1];
rz(-1.7130573) q[1];
sx q[1];
rz(-1.8131282) q[1];
rz(-pi) q[2];
x q[2];
rz(1.596132) q[3];
sx q[3];
rz(-1.5583894) q[3];
sx q[3];
rz(2.5125463) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.6015357) q[2];
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
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56101218) q[0];
sx q[0];
rz(-0.0479527) q[0];
sx q[0];
rz(-1.3797492) q[0];
rz(0.17669949) q[1];
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
rz(1.7226669) q[0];
sx q[0];
rz(-1.330266) q[0];
sx q[0];
rz(0.47970432) q[0];
x q[1];
rz(1.548544) q[2];
sx q[2];
rz(-1.3705727) q[2];
sx q[2];
rz(-2.3929838) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.2587535) q[1];
sx q[1];
rz(-1.4708609) q[1];
sx q[1];
rz(-0.12211166) q[1];
x q[2];
rz(1.3860116) q[3];
sx q[3];
rz(-1.4818958) q[3];
sx q[3];
rz(0.88528663) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.0195062) q[2];
sx q[2];
rz(-0.85947376) q[2];
sx q[2];
rz(-1.9879831) q[2];
rz(-2.7518318) q[3];
sx q[3];
rz(-0.4699769) q[3];
sx q[3];
rz(2.3526758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47238123) q[0];
sx q[0];
rz(-2.5166469) q[0];
sx q[0];
rz(0.56458449) q[0];
rz(-0.083104221) q[1];
sx q[1];
rz(-1.9080947) q[1];
sx q[1];
rz(-0.21534236) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1911083) q[0];
sx q[0];
rz(-1.9035625) q[0];
sx q[0];
rz(2.3118408) q[0];
rz(-2.9748561) q[2];
sx q[2];
rz(-0.98966375) q[2];
sx q[2];
rz(-0.8596479) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.035186471) q[1];
sx q[1];
rz(-1.3760929) q[1];
sx q[1];
rz(-3.0343642) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.54939778) q[3];
sx q[3];
rz(-1.975172) q[3];
sx q[3];
rz(-0.79493633) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.6354436) q[2];
sx q[2];
rz(-2.052867) q[2];
sx q[2];
rz(2.5396458) q[2];
rz(0.69058949) q[3];
sx q[3];
rz(-0.75406051) q[3];
sx q[3];
rz(0.62600342) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85201207) q[0];
sx q[0];
rz(-2.1471725) q[0];
sx q[0];
rz(2.1767298) q[0];
rz(-3.124974) q[1];
sx q[1];
rz(-2.7162572) q[1];
sx q[1];
rz(-2.4749277) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.704741) q[0];
sx q[0];
rz(-1.9691103) q[0];
sx q[0];
rz(3.0320291) q[0];
x q[1];
rz(-2.7890116) q[2];
sx q[2];
rz(-1.0513628) q[2];
sx q[2];
rz(-0.40371603) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.40068377) q[1];
sx q[1];
rz(-2.2212257) q[1];
sx q[1];
rz(-1.2121735) q[1];
x q[2];
rz(-0.22294754) q[3];
sx q[3];
rz(-1.2811536) q[3];
sx q[3];
rz(-1.2900316) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.2255286) q[2];
sx q[2];
rz(-1.7614044) q[2];
sx q[2];
rz(1.2529681) q[2];
rz(-1.3145674) q[3];
sx q[3];
rz(-1.1809177) q[3];
sx q[3];
rz(-2.0107646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1154293) q[0];
sx q[0];
rz(-0.9342497) q[0];
sx q[0];
rz(-1.3622267) q[0];
rz(1.399614) q[1];
sx q[1];
rz(-1.5618961) q[1];
sx q[1];
rz(1.6113575) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0936733) q[0];
sx q[0];
rz(-1.1512655) q[0];
sx q[0];
rz(-1.6506667) q[0];
x q[1];
rz(-2.7026664) q[2];
sx q[2];
rz(-1.274684) q[2];
sx q[2];
rz(0.93026464) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.46576071) q[1];
sx q[1];
rz(-0.45468802) q[1];
sx q[1];
rz(-1.0632221) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8239162) q[3];
sx q[3];
rz(-1.6198006) q[3];
sx q[3];
rz(-0.24310902) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.35649148) q[2];
sx q[2];
rz(-1.0078398) q[2];
sx q[2];
rz(-0.029416857) q[2];
rz(-1.4908837) q[3];
sx q[3];
rz(-2.1865032) q[3];
sx q[3];
rz(1.9194549) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71762639) q[0];
sx q[0];
rz(-0.76537776) q[0];
sx q[0];
rz(-2.9587865) q[0];
rz(0.43552661) q[1];
sx q[1];
rz(-2.2874449) q[1];
sx q[1];
rz(2.6307154) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8144007) q[0];
sx q[0];
rz(-1.6806024) q[0];
sx q[0];
rz(-2.8776684) q[0];
rz(-pi) q[1];
rz(0.73166087) q[2];
sx q[2];
rz(-0.96207843) q[2];
sx q[2];
rz(-2.4351956) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.43607298) q[1];
sx q[1];
rz(-1.1735859) q[1];
sx q[1];
rz(-2.436609) q[1];
rz(-1.0484344) q[3];
sx q[3];
rz(-1.8494542) q[3];
sx q[3];
rz(-0.1230965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.2304948) q[2];
sx q[2];
rz(-0.7926597) q[2];
sx q[2];
rz(2.5210209) q[2];
rz(2.7190322) q[3];
sx q[3];
rz(-1.7812604) q[3];
sx q[3];
rz(-1.871199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34185103) q[0];
sx q[0];
rz(-0.36062476) q[0];
sx q[0];
rz(-1.0890695) q[0];
rz(1.0808806) q[1];
sx q[1];
rz(-1.2842478) q[1];
sx q[1];
rz(2.506315) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7437744) q[0];
sx q[0];
rz(-1.213307) q[0];
sx q[0];
rz(2.9767569) q[0];
rz(-2.0681357) q[2];
sx q[2];
rz(-1.724223) q[2];
sx q[2];
rz(1.9286326) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(3.0465225) q[1];
sx q[1];
rz(-0.97192837) q[1];
sx q[1];
rz(1.8938766) q[1];
x q[2];
rz(-1.876272) q[3];
sx q[3];
rz(-2.5442122) q[3];
sx q[3];
rz(-1.3437831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.13715956) q[2];
sx q[2];
rz(-1.4617504) q[2];
sx q[2];
rz(0.00017246406) q[2];
rz(2.4553283) q[3];
sx q[3];
rz(-2.4133447) q[3];
sx q[3];
rz(0.85787684) q[3];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2962608) q[0];
sx q[0];
rz(-0.91264549) q[0];
sx q[0];
rz(-2.4849179) q[0];
rz(0.36390057) q[1];
sx q[1];
rz(-0.44635043) q[1];
sx q[1];
rz(2.231853) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7649987) q[0];
sx q[0];
rz(-1.8928327) q[0];
sx q[0];
rz(-1.8565208) q[0];
x q[1];
rz(-0.34809525) q[2];
sx q[2];
rz(-1.7762134) q[2];
sx q[2];
rz(1.5470488) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.197387) q[1];
sx q[1];
rz(-1.7802317) q[1];
sx q[1];
rz(1.6968813) q[1];
x q[2];
rz(-1.886456) q[3];
sx q[3];
rz(-1.8161895) q[3];
sx q[3];
rz(1.0823702) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.3108814) q[2];
sx q[2];
rz(-2.723366) q[2];
sx q[2];
rz(2.6549784) q[2];
rz(0.1085554) q[3];
sx q[3];
rz(-0.71043772) q[3];
sx q[3];
rz(-1.2238812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9412823) q[0];
sx q[0];
rz(-2.0599984) q[0];
sx q[0];
rz(2.6249028) q[0];
rz(-1.5962881) q[1];
sx q[1];
rz(-1.0524806) q[1];
sx q[1];
rz(-2.2470078) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6870118) q[0];
sx q[0];
rz(-2.5068388) q[0];
sx q[0];
rz(-1.4995585) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9355029) q[2];
sx q[2];
rz(-2.0317674) q[2];
sx q[2];
rz(3.0989625) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.67748653) q[1];
sx q[1];
rz(-1.9312526) q[1];
sx q[1];
rz(3.0739215) q[1];
x q[2];
rz(-1.3668725) q[3];
sx q[3];
rz(-2.4813934) q[3];
sx q[3];
rz(-0.064299671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
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
rz(-2.5953318) q[2];
sx q[2];
rz(-0.98304521) q[2];
sx q[2];
rz(1.8084768) q[2];
rz(3.0929052) q[3];
sx q[3];
rz(-1.8713453) q[3];
sx q[3];
rz(0.71181675) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];