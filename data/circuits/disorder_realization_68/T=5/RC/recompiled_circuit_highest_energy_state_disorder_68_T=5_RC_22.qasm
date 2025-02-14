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
rz(0.8402549) q[0];
sx q[0];
rz(4.1794887) q[0];
sx q[0];
rz(10.269796) q[0];
rz(0.71212274) q[1];
sx q[1];
rz(4.1432015) q[1];
sx q[1];
rz(10.92034) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5069731) q[0];
sx q[0];
rz(-1.15043) q[0];
sx q[0];
rz(0.93376055) q[0];
rz(-pi) q[1];
rz(-0.38365429) q[2];
sx q[2];
rz(-1.4315413) q[2];
sx q[2];
rz(-0.68472199) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.6604267) q[1];
sx q[1];
rz(-1.1200953) q[1];
sx q[1];
rz(0.53498241) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.59652416) q[3];
sx q[3];
rz(-2.6748195) q[3];
sx q[3];
rz(-1.0214361) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.738395) q[2];
sx q[2];
rz(-1.9591816) q[2];
sx q[2];
rz(0.5564059) q[2];
rz(2.6465936) q[3];
sx q[3];
rz(-2.7904816) q[3];
sx q[3];
rz(2.2569412) q[3];
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
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86968386) q[0];
sx q[0];
rz(-0.10232919) q[0];
sx q[0];
rz(-2.5129357) q[0];
rz(0.87720811) q[1];
sx q[1];
rz(-2.6429206) q[1];
sx q[1];
rz(-2.8884851) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0639187) q[0];
sx q[0];
rz(-1.5787933) q[0];
sx q[0];
rz(-1.5551989) q[0];
rz(2.5649389) q[2];
sx q[2];
rz(-2.6284784) q[2];
sx q[2];
rz(2.6000044) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.69324609) q[1];
sx q[1];
rz(-2.9102059) q[1];
sx q[1];
rz(-0.53129249) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3023002) q[3];
sx q[3];
rz(-1.7613693) q[3];
sx q[3];
rz(0.96123248) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.46193281) q[2];
sx q[2];
rz(-1.2075281) q[2];
sx q[2];
rz(-1.0665464) q[2];
rz(-1.8343532) q[3];
sx q[3];
rz(-2.6756838) q[3];
sx q[3];
rz(-2.2333142) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9268554) q[0];
sx q[0];
rz(-2.1934788) q[0];
sx q[0];
rz(0.98534775) q[0];
rz(-0.39380479) q[1];
sx q[1];
rz(-0.99291283) q[1];
sx q[1];
rz(-2.6148112) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7486289) q[0];
sx q[0];
rz(-2.606578) q[0];
sx q[0];
rz(1.7345384) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4220139) q[2];
sx q[2];
rz(-2.4256896) q[2];
sx q[2];
rz(0.86052513) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.5722428) q[1];
sx q[1];
rz(-1.4016289) q[1];
sx q[1];
rz(1.1430686) q[1];
rz(-0.97350307) q[3];
sx q[3];
rz(-1.1949348) q[3];
sx q[3];
rz(1.5099389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3452611) q[2];
sx q[2];
rz(-0.76097208) q[2];
sx q[2];
rz(0.14656466) q[2];
rz(2.2740299) q[3];
sx q[3];
rz(-2.2680794) q[3];
sx q[3];
rz(-0.7578907) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72916156) q[0];
sx q[0];
rz(-0.49598345) q[0];
sx q[0];
rz(-2.5878986) q[0];
rz(1.6665392) q[1];
sx q[1];
rz(-0.29860425) q[1];
sx q[1];
rz(2.8972304) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9825443) q[0];
sx q[0];
rz(-1.2384909) q[0];
sx q[0];
rz(-2.766702) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.30714037) q[2];
sx q[2];
rz(-2.2835586) q[2];
sx q[2];
rz(-2.3736726) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1235413) q[1];
sx q[1];
rz(-1.9885716) q[1];
sx q[1];
rz(-1.17735) q[1];
rz(-3.1126106) q[3];
sx q[3];
rz(-2.4292415) q[3];
sx q[3];
rz(-2.1952352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.0814521) q[2];
sx q[2];
rz(-2.0144561) q[2];
sx q[2];
rz(0.78090182) q[2];
rz(0.39685708) q[3];
sx q[3];
rz(-0.01440993) q[3];
sx q[3];
rz(2.0675596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2395372) q[0];
sx q[0];
rz(-2.9750415) q[0];
sx q[0];
rz(-2.8848414) q[0];
rz(2.9638839) q[1];
sx q[1];
rz(-0.54568988) q[1];
sx q[1];
rz(0.94388747) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0487918) q[0];
sx q[0];
rz(-1.9328797) q[0];
sx q[0];
rz(-2.0440153) q[0];
rz(-pi) q[1];
rz(2.947784) q[2];
sx q[2];
rz(-1.8558981) q[2];
sx q[2];
rz(-1.5799892) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.9967794) q[1];
sx q[1];
rz(-2.6018892) q[1];
sx q[1];
rz(-0.0085047095) q[1];
x q[2];
rz(-2.2511399) q[3];
sx q[3];
rz(-2.4665006) q[3];
sx q[3];
rz(-2.4321041) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.2029767) q[2];
sx q[2];
rz(-2.1111033) q[2];
sx q[2];
rz(-0.20982783) q[2];
rz(0.046791568) q[3];
sx q[3];
rz(-1.4593461) q[3];
sx q[3];
rz(-0.34745026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.059134722) q[0];
sx q[0];
rz(-2.3248398) q[0];
sx q[0];
rz(1.9385852) q[0];
rz(-1.6251534) q[1];
sx q[1];
rz(-2.7639183) q[1];
sx q[1];
rz(-3.1086521) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0432949) q[0];
sx q[0];
rz(-1.2970719) q[0];
sx q[0];
rz(2.0541463) q[0];
rz(2.158098) q[2];
sx q[2];
rz(-1.1193395) q[2];
sx q[2];
rz(1.2304103) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.1748993) q[1];
sx q[1];
rz(-1.0224578) q[1];
sx q[1];
rz(-1.95005) q[1];
rz(-pi) q[2];
rz(2.2516656) q[3];
sx q[3];
rz(-1.0491199) q[3];
sx q[3];
rz(-1.3104749) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.5861627) q[2];
sx q[2];
rz(-2.1383492) q[2];
sx q[2];
rz(0.45247751) q[2];
rz(-2.9943976) q[3];
sx q[3];
rz(-0.23284027) q[3];
sx q[3];
rz(1.8424621) q[3];
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
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.066605695) q[0];
sx q[0];
rz(-0.43656483) q[0];
sx q[0];
rz(2.7845352) q[0];
rz(-1.0128516) q[1];
sx q[1];
rz(-1.169299) q[1];
sx q[1];
rz(0.036651932) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79082876) q[0];
sx q[0];
rz(-1.1393271) q[0];
sx q[0];
rz(-2.55326) q[0];
x q[1];
rz(-0.22256644) q[2];
sx q[2];
rz(-2.7024726) q[2];
sx q[2];
rz(-0.47638461) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.086163047) q[1];
sx q[1];
rz(-2.2698672) q[1];
sx q[1];
rz(3.0208339) q[1];
rz(-pi) q[2];
rz(-1.5753079) q[3];
sx q[3];
rz(-1.5399714) q[3];
sx q[3];
rz(-0.84336057) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.4770294) q[2];
sx q[2];
rz(-2.1881115) q[2];
sx q[2];
rz(-0.089740962) q[2];
rz(-1.2612032) q[3];
sx q[3];
rz(-1.3124876) q[3];
sx q[3];
rz(-1.7173654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
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
rz(2.4795714) q[0];
sx q[0];
rz(-2.5630072) q[0];
sx q[0];
rz(-1.4805967) q[0];
rz(1.4501694) q[1];
sx q[1];
rz(-2.4822576) q[1];
sx q[1];
rz(2.3816542) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1131033) q[0];
sx q[0];
rz(-3.0249502) q[0];
sx q[0];
rz(1.4820497) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0133661) q[2];
sx q[2];
rz(-0.63214801) q[2];
sx q[2];
rz(1.8116443) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.41553283) q[1];
sx q[1];
rz(-0.59846717) q[1];
sx q[1];
rz(-2.3319671) q[1];
rz(-pi) q[2];
rz(-2.6344243) q[3];
sx q[3];
rz(-1.7970603) q[3];
sx q[3];
rz(1.7698225) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.8370168) q[2];
sx q[2];
rz(-1.872007) q[2];
sx q[2];
rz(3.0820091) q[2];
rz(-0.42851055) q[3];
sx q[3];
rz(-0.84463745) q[3];
sx q[3];
rz(-2.5394411) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8932327) q[0];
sx q[0];
rz(-0.28286523) q[0];
sx q[0];
rz(0.34348139) q[0];
rz(0.19613014) q[1];
sx q[1];
rz(-2.2562512) q[1];
sx q[1];
rz(2.3416065) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4206652) q[0];
sx q[0];
rz(-1.1625097) q[0];
sx q[0];
rz(2.3022006) q[0];
rz(-pi) q[1];
rz(1.3197199) q[2];
sx q[2];
rz(-2.1188687) q[2];
sx q[2];
rz(-2.1679116) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.0240805) q[1];
sx q[1];
rz(-1.5264866) q[1];
sx q[1];
rz(-2.5818392) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0723582) q[3];
sx q[3];
rz(-0.88919176) q[3];
sx q[3];
rz(-0.72661663) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.34667748) q[2];
sx q[2];
rz(-0.76378834) q[2];
sx q[2];
rz(-2.9128892) q[2];
rz(-1.7250569) q[3];
sx q[3];
rz(-1.4226457) q[3];
sx q[3];
rz(-0.67030877) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5423841) q[0];
sx q[0];
rz(-0.51775652) q[0];
sx q[0];
rz(-1.1826578) q[0];
rz(-0.54284894) q[1];
sx q[1];
rz(-0.12195568) q[1];
sx q[1];
rz(-0.45546946) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96129721) q[0];
sx q[0];
rz(-1.1396798) q[0];
sx q[0];
rz(-2.6102553) q[0];
x q[1];
rz(2.9331839) q[2];
sx q[2];
rz(-1.128607) q[2];
sx q[2];
rz(0.65763523) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.9046136) q[1];
sx q[1];
rz(-1.1061258) q[1];
sx q[1];
rz(-1.418271) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1703759) q[3];
sx q[3];
rz(-1.9846623) q[3];
sx q[3];
rz(-0.032776959) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.21696422) q[2];
sx q[2];
rz(-2.2483726) q[2];
sx q[2];
rz(-2.7157937) q[2];
rz(-2.9567772) q[3];
sx q[3];
rz(-0.63822377) q[3];
sx q[3];
rz(3.1137915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7040779) q[0];
sx q[0];
rz(-1.5298433) q[0];
sx q[0];
rz(-0.035506305) q[0];
rz(-0.47954814) q[1];
sx q[1];
rz(-1.5621114) q[1];
sx q[1];
rz(1.5893804) q[1];
rz(-2.5465847) q[2];
sx q[2];
rz(-1.7562661) q[2];
sx q[2];
rz(-1.4944639) q[2];
rz(2.1816523) q[3];
sx q[3];
rz(-1.6889986) q[3];
sx q[3];
rz(-2.1903174) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
