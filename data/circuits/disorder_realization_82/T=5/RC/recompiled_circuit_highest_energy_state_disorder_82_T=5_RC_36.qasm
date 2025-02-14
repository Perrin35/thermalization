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
rz(1.0225811) q[0];
sx q[0];
rz(5.2978088) q[0];
sx q[0];
rz(10.537416) q[0];
rz(-2.7467709) q[1];
sx q[1];
rz(-2.0070751) q[1];
sx q[1];
rz(1.1910103) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1853088) q[0];
sx q[0];
rz(-2.5155768) q[0];
sx q[0];
rz(-0.52395384) q[0];
rz(0.13220867) q[2];
sx q[2];
rz(-1.9763051) q[2];
sx q[2];
rz(-2.8762238) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.2996051) q[1];
sx q[1];
rz(-1.6374432) q[1];
sx q[1];
rz(1.6377836) q[1];
rz(-pi) q[2];
rz(-0.43893473) q[3];
sx q[3];
rz(-1.1263444) q[3];
sx q[3];
rz(-1.7308472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.4224008) q[2];
sx q[2];
rz(-0.16133186) q[2];
sx q[2];
rz(2.7154229) q[2];
rz(-0.69914493) q[3];
sx q[3];
rz(-1.1933425) q[3];
sx q[3];
rz(-0.41081158) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
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
rz(0.29546577) q[0];
sx q[0];
rz(-2.182425) q[0];
sx q[0];
rz(2.3058983) q[0];
rz(0.524638) q[1];
sx q[1];
rz(-1.6181889) q[1];
sx q[1];
rz(1.5118648) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.152911) q[0];
sx q[0];
rz(-0.67711867) q[0];
sx q[0];
rz(-2.0858913) q[0];
rz(0.68372441) q[2];
sx q[2];
rz(-1.0283768) q[2];
sx q[2];
rz(2.8680141) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.8162093) q[1];
sx q[1];
rz(-2.0289501) q[1];
sx q[1];
rz(1.7857331) q[1];
rz(-2.6523051) q[3];
sx q[3];
rz(-0.60859546) q[3];
sx q[3];
rz(1.6123475) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.52874804) q[2];
sx q[2];
rz(-1.5809487) q[2];
sx q[2];
rz(0.35231248) q[2];
rz(0.49255541) q[3];
sx q[3];
rz(-1.1198606) q[3];
sx q[3];
rz(0.39920863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6414129) q[0];
sx q[0];
rz(-2.4523875) q[0];
sx q[0];
rz(-2.8885762) q[0];
rz(0.5223271) q[1];
sx q[1];
rz(-2.7136007) q[1];
sx q[1];
rz(2.7740251) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92644405) q[0];
sx q[0];
rz(-0.068000168) q[0];
sx q[0];
rz(-1.1626194) q[0];
x q[1];
rz(2.7455212) q[2];
sx q[2];
rz(-2.1537915) q[2];
sx q[2];
rz(-0.99882674) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.8525703) q[1];
sx q[1];
rz(-1.7779568) q[1];
sx q[1];
rz(3.0604355) q[1];
rz(-3.0250164) q[3];
sx q[3];
rz(-1.4287144) q[3];
sx q[3];
rz(0.97654479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(3.0005409) q[2];
sx q[2];
rz(-2.8512569) q[2];
sx q[2];
rz(-2.7153437) q[2];
rz(0.78921562) q[3];
sx q[3];
rz(-1.1169249) q[3];
sx q[3];
rz(-1.5392019) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57192794) q[0];
sx q[0];
rz(-2.8043788) q[0];
sx q[0];
rz(-0.1871044) q[0];
rz(-1.7722173) q[1];
sx q[1];
rz(-1.2792842) q[1];
sx q[1];
rz(2.5274136) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2046004) q[0];
sx q[0];
rz(-2.061917) q[0];
sx q[0];
rz(3.0762927) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.34030045) q[2];
sx q[2];
rz(-1.0522166) q[2];
sx q[2];
rz(0.22835635) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.56127292) q[1];
sx q[1];
rz(-0.6918219) q[1];
sx q[1];
rz(-2.0136334) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.21728269) q[3];
sx q[3];
rz(-0.38300316) q[3];
sx q[3];
rz(-0.10620761) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.46089178) q[2];
sx q[2];
rz(-0.9047752) q[2];
sx q[2];
rz(1.2659849) q[2];
rz(2.8547817) q[3];
sx q[3];
rz(-2.3293142) q[3];
sx q[3];
rz(0.12107818) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1665961) q[0];
sx q[0];
rz(-0.91355211) q[0];
sx q[0];
rz(-1.1443369) q[0];
rz(0.76885778) q[1];
sx q[1];
rz(-2.2298593) q[1];
sx q[1];
rz(1.2164046) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74690565) q[0];
sx q[0];
rz(-2.0187067) q[0];
sx q[0];
rz(0.42952092) q[0];
rz(-1.7038962) q[2];
sx q[2];
rz(-1.6665227) q[2];
sx q[2];
rz(-1.7246703) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.2768086) q[1];
sx q[1];
rz(-1.4866765) q[1];
sx q[1];
rz(-1.8451971) q[1];
rz(2.9202095) q[3];
sx q[3];
rz(-2.7714344) q[3];
sx q[3];
rz(-2.5551318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.6151109) q[2];
sx q[2];
rz(-0.59978008) q[2];
sx q[2];
rz(-0.46670023) q[2];
rz(-2.1488721) q[3];
sx q[3];
rz(-1.4027184) q[3];
sx q[3];
rz(1.3950951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65527958) q[0];
sx q[0];
rz(-1.7236973) q[0];
sx q[0];
rz(0.41743761) q[0];
rz(-2.7159269) q[1];
sx q[1];
rz(-0.8368496) q[1];
sx q[1];
rz(2.4136037) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54386747) q[0];
sx q[0];
rz(-2.4222582) q[0];
sx q[0];
rz(-0.85605232) q[0];
rz(-pi) q[1];
x q[1];
rz(0.73982088) q[2];
sx q[2];
rz(-1.4342116) q[2];
sx q[2];
rz(0.53892577) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.14236808) q[1];
sx q[1];
rz(-1.2877688) q[1];
sx q[1];
rz(1.9623161) q[1];
rz(-pi) q[2];
rz(-2.7397104) q[3];
sx q[3];
rz(-1.0392351) q[3];
sx q[3];
rz(1.8762003) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.723168) q[2];
sx q[2];
rz(-1.4801414) q[2];
sx q[2];
rz(1.9815014) q[2];
rz(0.55231071) q[3];
sx q[3];
rz(-0.6902802) q[3];
sx q[3];
rz(2.6753329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1138678) q[0];
sx q[0];
rz(-0.92806569) q[0];
sx q[0];
rz(3.0357251) q[0];
rz(1.8766807) q[1];
sx q[1];
rz(-2.5906339) q[1];
sx q[1];
rz(1.4531306) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47521771) q[0];
sx q[0];
rz(-1.8868108) q[0];
sx q[0];
rz(2.3316335) q[0];
rz(-pi) q[1];
rz(-1.7488453) q[2];
sx q[2];
rz(-1.4499003) q[2];
sx q[2];
rz(0.73261443) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.2699996) q[1];
sx q[1];
rz(-1.8687222) q[1];
sx q[1];
rz(2.9948283) q[1];
rz(-pi) q[2];
rz(0.14554056) q[3];
sx q[3];
rz(-0.91314935) q[3];
sx q[3];
rz(-1.4484792) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.10446163) q[2];
sx q[2];
rz(-1.1300602) q[2];
sx q[2];
rz(-1.5360443) q[2];
rz(1.7724841) q[3];
sx q[3];
rz(-0.57309279) q[3];
sx q[3];
rz(-1.8039186) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8417514) q[0];
sx q[0];
rz(-1.2513237) q[0];
sx q[0];
rz(1.9642255) q[0];
rz(-1.8384701) q[1];
sx q[1];
rz(-1.4812555) q[1];
sx q[1];
rz(2.4633113) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2934351) q[0];
sx q[0];
rz(-0.57852902) q[0];
sx q[0];
rz(-1.0941157) q[0];
rz(-pi) q[1];
rz(0.2770284) q[2];
sx q[2];
rz(-1.4619383) q[2];
sx q[2];
rz(-1.9831374) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.5694414) q[1];
sx q[1];
rz(-1.8742233) q[1];
sx q[1];
rz(1.3511168) q[1];
rz(-pi) q[2];
rz(-0.60538624) q[3];
sx q[3];
rz(-1.5380713) q[3];
sx q[3];
rz(2.721019) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.9699041) q[2];
sx q[2];
rz(-2.6320612) q[2];
sx q[2];
rz(0.92863885) q[2];
rz(1.6396133) q[3];
sx q[3];
rz(-2.00311) q[3];
sx q[3];
rz(1.5750711) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0314727) q[0];
sx q[0];
rz(-1.1290154) q[0];
sx q[0];
rz(-1.336115) q[0];
rz(-0.31088343) q[1];
sx q[1];
rz(-1.6115522) q[1];
sx q[1];
rz(1.8119887) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80547896) q[0];
sx q[0];
rz(-0.74148899) q[0];
sx q[0];
rz(-0.51261897) q[0];
x q[1];
rz(-2.5104292) q[2];
sx q[2];
rz(-0.97056015) q[2];
sx q[2];
rz(1.7886666) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.2488798) q[1];
sx q[1];
rz(-1.5096438) q[1];
sx q[1];
rz(-0.11492954) q[1];
rz(-0.90985591) q[3];
sx q[3];
rz(-1.5136711) q[3];
sx q[3];
rz(-1.8450246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.10799321) q[2];
sx q[2];
rz(-0.42069837) q[2];
sx q[2];
rz(-1.1619953) q[2];
rz(1.1437931) q[3];
sx q[3];
rz(-1.1868718) q[3];
sx q[3];
rz(-2.7403045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6859632) q[0];
sx q[0];
rz(-0.81448737) q[0];
sx q[0];
rz(2.234835) q[0];
rz(2.7470159) q[1];
sx q[1];
rz(-2.5386609) q[1];
sx q[1];
rz(2.0049863) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4986269) q[0];
sx q[0];
rz(-1.8436448) q[0];
sx q[0];
rz(0.26609315) q[0];
x q[1];
rz(-1.8682212) q[2];
sx q[2];
rz(-1.6417268) q[2];
sx q[2];
rz(-2.298722) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9533331) q[1];
sx q[1];
rz(-1.430295) q[1];
sx q[1];
rz(-2.2849915) q[1];
rz(-pi) q[2];
rz(1.2796938) q[3];
sx q[3];
rz(-1.2754692) q[3];
sx q[3];
rz(-0.81084033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-3.0236728) q[2];
sx q[2];
rz(-1.6912141) q[2];
sx q[2];
rz(-1.0517906) q[2];
rz(3.1045095) q[3];
sx q[3];
rz(-0.6784234) q[3];
sx q[3];
rz(1.7557433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
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
rz(-0.15569923) q[1];
sx q[1];
rz(-2.1035879) q[1];
sx q[1];
rz(-2.9661967) q[1];
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
