OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.1190342) q[0];
sx q[0];
rz(-2.0329539) q[0];
sx q[0];
rz(-0.57060567) q[0];
rz(1.624149) q[1];
sx q[1];
rz(-0.94944209) q[1];
sx q[1];
rz(-1.5631262) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9304912) q[0];
sx q[0];
rz(-1.6864713) q[0];
sx q[0];
rz(-2.0086847) q[0];
rz(0.3508829) q[2];
sx q[2];
rz(-2.299438) q[2];
sx q[2];
rz(-0.0082727783) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0277869) q[1];
sx q[1];
rz(-1.1996197) q[1];
sx q[1];
rz(1.3813301) q[1];
rz(-pi) q[2];
rz(-2.6622979) q[3];
sx q[3];
rz(-2.473978) q[3];
sx q[3];
rz(1.7874174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.7382875) q[2];
sx q[2];
rz(-2.3052576) q[2];
sx q[2];
rz(2.9264012) q[2];
rz(0.0067986851) q[3];
sx q[3];
rz(-2.2577622) q[3];
sx q[3];
rz(2.1782037) q[3];
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
rz(-2.6007518) q[0];
sx q[0];
rz(-2.231926) q[0];
sx q[0];
rz(1.2700861) q[0];
rz(-2.0974244) q[1];
sx q[1];
rz(-0.41193286) q[1];
sx q[1];
rz(2.5792436) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.017496271) q[0];
sx q[0];
rz(-1.3096333) q[0];
sx q[0];
rz(-0.11696741) q[0];
rz(-2.5452209) q[2];
sx q[2];
rz(-1.2105701) q[2];
sx q[2];
rz(0.63186592) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.54004242) q[1];
sx q[1];
rz(-1.6311495) q[1];
sx q[1];
rz(0.73436952) q[1];
rz(-1.2581732) q[3];
sx q[3];
rz(-1.2762831) q[3];
sx q[3];
rz(-2.9972863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(3.1336141) q[2];
sx q[2];
rz(-2.9972711) q[2];
sx q[2];
rz(-3.0001384) q[2];
rz(1.0072297) q[3];
sx q[3];
rz(-1.6601345) q[3];
sx q[3];
rz(-3.0032515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7815642) q[0];
sx q[0];
rz(-1.6369632) q[0];
sx q[0];
rz(0.052852782) q[0];
rz(-1.3797034) q[1];
sx q[1];
rz(-0.77769891) q[1];
sx q[1];
rz(-2.8676829) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8192379) q[0];
sx q[0];
rz(-2.9606426) q[0];
sx q[0];
rz(1.7767679) q[0];
x q[1];
rz(-3.048298) q[2];
sx q[2];
rz(-0.81893605) q[2];
sx q[2];
rz(-1.3064177) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.1812836) q[1];
sx q[1];
rz(-2.3707254) q[1];
sx q[1];
rz(0.038357448) q[1];
x q[2];
rz(1.307358) q[3];
sx q[3];
rz(-1.4589849) q[3];
sx q[3];
rz(-1.5105351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2770047) q[2];
sx q[2];
rz(-1.6463582) q[2];
sx q[2];
rz(-2.541339) q[2];
rz(0.65435919) q[3];
sx q[3];
rz(-2.8221966) q[3];
sx q[3];
rz(-1.4734242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7905423) q[0];
sx q[0];
rz(-1.496614) q[0];
sx q[0];
rz(0.45781621) q[0];
rz(0.59902016) q[1];
sx q[1];
rz(-2.6287703) q[1];
sx q[1];
rz(0.27214989) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3783589) q[0];
sx q[0];
rz(-2.1714179) q[0];
sx q[0];
rz(-1.8778223) q[0];
rz(-1.0402914) q[2];
sx q[2];
rz(-1.1018424) q[2];
sx q[2];
rz(-0.95382849) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.052201) q[1];
sx q[1];
rz(-1.7551345) q[1];
sx q[1];
rz(-1.6910529) q[1];
rz(-pi) q[2];
rz(0.64025819) q[3];
sx q[3];
rz(-1.4319515) q[3];
sx q[3];
rz(1.1306896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.6279471) q[2];
sx q[2];
rz(-2.7529035) q[2];
sx q[2];
rz(1.2316616) q[2];
rz(-1.4218467) q[3];
sx q[3];
rz(-1.3323517) q[3];
sx q[3];
rz(2.1739912) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.642639) q[0];
sx q[0];
rz(-3.0727144) q[0];
sx q[0];
rz(-2.4054085) q[0];
rz(0.63756293) q[1];
sx q[1];
rz(-1.2022737) q[1];
sx q[1];
rz(-2.7307302) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81911892) q[0];
sx q[0];
rz(-1.1000191) q[0];
sx q[0];
rz(-2.0841925) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.246894) q[2];
sx q[2];
rz(-2.5482087) q[2];
sx q[2];
rz(1.3932799) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.84437925) q[1];
sx q[1];
rz(-2.2003268) q[1];
sx q[1];
rz(1.5624863) q[1];
rz(-pi) q[2];
rz(-2.210064) q[3];
sx q[3];
rz(-1.0088822) q[3];
sx q[3];
rz(-0.84980295) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.8132482) q[2];
sx q[2];
rz(-1.5466377) q[2];
sx q[2];
rz(-0.52097875) q[2];
rz(-3.0327435) q[3];
sx q[3];
rz(-2.1284926) q[3];
sx q[3];
rz(-0.61645761) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9209038) q[0];
sx q[0];
rz(-1.2257129) q[0];
sx q[0];
rz(-1.2589681) q[0];
rz(1.0733696) q[1];
sx q[1];
rz(-1.4732889) q[1];
sx q[1];
rz(0.64356709) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63976906) q[0];
sx q[0];
rz(-1.5657525) q[0];
sx q[0];
rz(-2.6905294) q[0];
rz(-pi) q[1];
x q[1];
rz(0.50711716) q[2];
sx q[2];
rz(-2.7165453) q[2];
sx q[2];
rz(0.054946446) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.3583329) q[1];
sx q[1];
rz(-1.734501) q[1];
sx q[1];
rz(-0.81061426) q[1];
rz(0.57585133) q[3];
sx q[3];
rz(-2.0178049) q[3];
sx q[3];
rz(0.67547638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.7966938) q[2];
sx q[2];
rz(-1.0641791) q[2];
sx q[2];
rz(-1.8033484) q[2];
rz(-1.6898588) q[3];
sx q[3];
rz(-1.9386049) q[3];
sx q[3];
rz(-0.10045997) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.34375) q[0];
sx q[0];
rz(-1.8666973) q[0];
sx q[0];
rz(-2.1882958) q[0];
rz(0.15996179) q[1];
sx q[1];
rz(-2.5123031) q[1];
sx q[1];
rz(2.2728641) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1686483) q[0];
sx q[0];
rz(-0.7059083) q[0];
sx q[0];
rz(-1.9078518) q[0];
rz(-0.80256497) q[2];
sx q[2];
rz(-2.2266529) q[2];
sx q[2];
rz(-0.64284113) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.2273457) q[1];
sx q[1];
rz(-0.29623756) q[1];
sx q[1];
rz(0.072320894) q[1];
rz(2.8817435) q[3];
sx q[3];
rz(-0.9161549) q[3];
sx q[3];
rz(0.25370592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.4538883) q[2];
sx q[2];
rz(-1.1331465) q[2];
sx q[2];
rz(0.49393168) q[2];
rz(1.6535053) q[3];
sx q[3];
rz(-0.86797124) q[3];
sx q[3];
rz(-3.119829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22773753) q[0];
sx q[0];
rz(-1.7182925) q[0];
sx q[0];
rz(2.8505351) q[0];
rz(-0.13045467) q[1];
sx q[1];
rz(-2.7619402) q[1];
sx q[1];
rz(1.5694654) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1324738) q[0];
sx q[0];
rz(-1.5838916) q[0];
sx q[0];
rz(1.6214236) q[0];
rz(-pi) q[1];
rz(-1.7849017) q[2];
sx q[2];
rz(-2.4567025) q[2];
sx q[2];
rz(0.30944007) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2184184) q[1];
sx q[1];
rz(-1.4552792) q[1];
sx q[1];
rz(-0.61489132) q[1];
x q[2];
rz(-2.8329141) q[3];
sx q[3];
rz(-1.4490845) q[3];
sx q[3];
rz(-2.6409923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.7743885) q[2];
sx q[2];
rz(-2.9640894) q[2];
sx q[2];
rz(3.0299178) q[2];
rz(2.3292134) q[3];
sx q[3];
rz(-1.8492536) q[3];
sx q[3];
rz(1.5264548) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3647301) q[0];
sx q[0];
rz(-0.22219816) q[0];
sx q[0];
rz(-0.36866933) q[0];
rz(-0.61019507) q[1];
sx q[1];
rz(-1.5851603) q[1];
sx q[1];
rz(-1.0367702) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7198349) q[0];
sx q[0];
rz(-1.0204691) q[0];
sx q[0];
rz(-1.3070413) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5620932) q[2];
sx q[2];
rz(-2.6285951) q[2];
sx q[2];
rz(2.7287116) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.887275) q[1];
sx q[1];
rz(-0.84158449) q[1];
sx q[1];
rz(-2.9542854) q[1];
rz(-pi) q[2];
rz(0.16513326) q[3];
sx q[3];
rz(-2.0516178) q[3];
sx q[3];
rz(-0.1827249) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.6250352) q[2];
sx q[2];
rz(-0.80356193) q[2];
sx q[2];
rz(0.45300031) q[2];
rz(2.1895444) q[3];
sx q[3];
rz(-1.1166078) q[3];
sx q[3];
rz(3.0208352) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20671885) q[0];
sx q[0];
rz(-3.0937338) q[0];
sx q[0];
rz(0.55024838) q[0];
rz(1.0847367) q[1];
sx q[1];
rz(-2.5526498) q[1];
sx q[1];
rz(-1.3909856) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9828331) q[0];
sx q[0];
rz(-1.6102432) q[0];
sx q[0];
rz(2.0723913) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1065305) q[2];
sx q[2];
rz(-1.4260056) q[2];
sx q[2];
rz(2.5404251) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.15074888) q[1];
sx q[1];
rz(-2.1865784) q[1];
sx q[1];
rz(0.94373871) q[1];
rz(1.5143366) q[3];
sx q[3];
rz(-1.5096231) q[3];
sx q[3];
rz(0.45555761) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7636106) q[2];
sx q[2];
rz(-0.72145975) q[2];
sx q[2];
rz(1.1627496) q[2];
rz(0.55443305) q[3];
sx q[3];
rz(-0.88397637) q[3];
sx q[3];
rz(-0.46816167) q[3];
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
x q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8155895) q[0];
sx q[0];
rz(-1.8392039) q[0];
sx q[0];
rz(-1.4890672) q[0];
rz(-0.21492699) q[1];
sx q[1];
rz(-0.88354127) q[1];
sx q[1];
rz(3.0149928) q[1];
rz(2.8561572) q[2];
sx q[2];
rz(-2.4531721) q[2];
sx q[2];
rz(2.5317818) q[2];
rz(-2.4195803) q[3];
sx q[3];
rz(-1.5215254) q[3];
sx q[3];
rz(1.4302614) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
