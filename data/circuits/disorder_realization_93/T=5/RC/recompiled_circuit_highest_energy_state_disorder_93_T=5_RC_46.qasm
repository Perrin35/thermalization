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
rz(-0.45577249) q[0];
sx q[0];
rz(5.1626212) q[0];
sx q[0];
rz(10.838312) q[0];
rz(-2.7780374) q[1];
sx q[1];
rz(-1.8283365) q[1];
sx q[1];
rz(0.29120905) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.04714) q[0];
sx q[0];
rz(-1.5908851) q[0];
sx q[0];
rz(-1.2506586) q[0];
x q[1];
rz(0.12367308) q[2];
sx q[2];
rz(-1.9078662) q[2];
sx q[2];
rz(2.6840984) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.5346966) q[1];
sx q[1];
rz(-1.7217858) q[1];
sx q[1];
rz(-2.5642337) q[1];
rz(0.81756567) q[3];
sx q[3];
rz(-1.0888087) q[3];
sx q[3];
rz(2.7369902) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(3.1353961) q[2];
sx q[2];
rz(-1.5634894) q[2];
sx q[2];
rz(-3.1175933) q[2];
rz(0.35605797) q[3];
sx q[3];
rz(-0.67432109) q[3];
sx q[3];
rz(0.02478987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0559219) q[0];
sx q[0];
rz(-0.69127685) q[0];
sx q[0];
rz(-0.63496494) q[0];
rz(0.7629281) q[1];
sx q[1];
rz(-1.4632016) q[1];
sx q[1];
rz(-0.91711226) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.791009) q[0];
sx q[0];
rz(-1.6054594) q[0];
sx q[0];
rz(3.123218) q[0];
x q[1];
rz(2.1371827) q[2];
sx q[2];
rz(-2.3766563) q[2];
sx q[2];
rz(1.1980008) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.32815427) q[1];
sx q[1];
rz(-2.0106689) q[1];
sx q[1];
rz(3.0564223) q[1];
rz(-pi) q[2];
rz(0.95083745) q[3];
sx q[3];
rz(-1.7449494) q[3];
sx q[3];
rz(-3.1319428) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.94747535) q[2];
sx q[2];
rz(-2.6185991) q[2];
sx q[2];
rz(2.4083162) q[2];
rz(-2.0745847) q[3];
sx q[3];
rz(-1.7206444) q[3];
sx q[3];
rz(-0.1799306) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5718403) q[0];
sx q[0];
rz(-1.952992) q[0];
sx q[0];
rz(-0.52823129) q[0];
rz(-2.275548) q[1];
sx q[1];
rz(-1.3287611) q[1];
sx q[1];
rz(2.6285062) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7932803) q[0];
sx q[0];
rz(-1.6133119) q[0];
sx q[0];
rz(2.0094064) q[0];
rz(-pi) q[1];
rz(0.92136356) q[2];
sx q[2];
rz(-1.8214392) q[2];
sx q[2];
rz(-1.373137) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.60050234) q[1];
sx q[1];
rz(-1.2230495) q[1];
sx q[1];
rz(3.0060124) q[1];
rz(3.0894075) q[3];
sx q[3];
rz(-0.5002678) q[3];
sx q[3];
rz(1.0104928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-3.1055866) q[2];
sx q[2];
rz(-2.0522108) q[2];
sx q[2];
rz(-0.038724381) q[2];
rz(0.45768467) q[3];
sx q[3];
rz(-0.80515146) q[3];
sx q[3];
rz(1.3409486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7078581) q[0];
sx q[0];
rz(-2.7271294) q[0];
sx q[0];
rz(-2.0204954) q[0];
rz(-0.7577678) q[1];
sx q[1];
rz(-1.6242177) q[1];
sx q[1];
rz(0.68414348) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.093903784) q[0];
sx q[0];
rz(-2.2950116) q[0];
sx q[0];
rz(-2.8709477) q[0];
x q[1];
rz(1.5779547) q[2];
sx q[2];
rz(-1.560692) q[2];
sx q[2];
rz(-1.2105389) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.0348548) q[1];
sx q[1];
rz(-1.8632006) q[1];
sx q[1];
rz(0.3180983) q[1];
rz(-2.2627801) q[3];
sx q[3];
rz(-2.4185816) q[3];
sx q[3];
rz(-2.3634152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.6571558) q[2];
sx q[2];
rz(-1.5269205) q[2];
sx q[2];
rz(0.26051513) q[2];
rz(-0.83812964) q[3];
sx q[3];
rz(-1.1476436) q[3];
sx q[3];
rz(-0.014613541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92358661) q[0];
sx q[0];
rz(-1.7638693) q[0];
sx q[0];
rz(2.6309784) q[0];
rz(-1.249373) q[1];
sx q[1];
rz(-1.9793972) q[1];
sx q[1];
rz(1.6872663) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73390612) q[0];
sx q[0];
rz(-0.21557325) q[0];
sx q[0];
rz(3.0093497) q[0];
x q[1];
rz(-2.8122105) q[2];
sx q[2];
rz(-2.9940146) q[2];
sx q[2];
rz(1.7640698) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.92125398) q[1];
sx q[1];
rz(-0.61437139) q[1];
sx q[1];
rz(-2.6831616) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.733414) q[3];
sx q[3];
rz(-1.4282886) q[3];
sx q[3];
rz(2.5385025) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.55383596) q[2];
sx q[2];
rz(-2.7601056) q[2];
sx q[2];
rz(-0.39056632) q[2];
rz(-0.25911123) q[3];
sx q[3];
rz(-1.6389537) q[3];
sx q[3];
rz(2.794877) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3542341) q[0];
sx q[0];
rz(-0.43208313) q[0];
sx q[0];
rz(2.4969192) q[0];
rz(1.8395754) q[1];
sx q[1];
rz(-1.527486) q[1];
sx q[1];
rz(-0.34861809) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.020418305) q[0];
sx q[0];
rz(-2.0266337) q[0];
sx q[0];
rz(-2.0860079) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.604282) q[2];
sx q[2];
rz(-1.9331353) q[2];
sx q[2];
rz(-2.5700931) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.7346361) q[1];
sx q[1];
rz(-2.6183081) q[1];
sx q[1];
rz(0.95013817) q[1];
x q[2];
rz(-2.8556267) q[3];
sx q[3];
rz(-3.0197201) q[3];
sx q[3];
rz(-1.7274477) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.0384486) q[2];
sx q[2];
rz(-2.1861173) q[2];
sx q[2];
rz(3.1032584) q[2];
rz(2.1740055) q[3];
sx q[3];
rz(-1.7318232) q[3];
sx q[3];
rz(-2.7819395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2751145) q[0];
sx q[0];
rz(-2.622128) q[0];
sx q[0];
rz(3.0679829) q[0];
rz(1.7346802) q[1];
sx q[1];
rz(-0.55714566) q[1];
sx q[1];
rz(0.73307347) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0778962) q[0];
sx q[0];
rz(-2.1025582) q[0];
sx q[0];
rz(-2.1666906) q[0];
rz(-pi) q[1];
rz(-0.91821155) q[2];
sx q[2];
rz(-2.2951153) q[2];
sx q[2];
rz(1.0670964) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.0508557) q[1];
sx q[1];
rz(-0.49492087) q[1];
sx q[1];
rz(-0.75116091) q[1];
rz(0.81871512) q[3];
sx q[3];
rz(-2.8041841) q[3];
sx q[3];
rz(2.329987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.0343895) q[2];
sx q[2];
rz(-1.4185536) q[2];
sx q[2];
rz(-0.011073152) q[2];
rz(-2.8360046) q[3];
sx q[3];
rz(-2.0162851) q[3];
sx q[3];
rz(1.8886458) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8082064) q[0];
sx q[0];
rz(-0.60788637) q[0];
sx q[0];
rz(-0.73390865) q[0];
rz(-0.58087307) q[1];
sx q[1];
rz(-1.0092694) q[1];
sx q[1];
rz(3.0427921) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2487468) q[0];
sx q[0];
rz(-0.98631421) q[0];
sx q[0];
rz(-0.92846034) q[0];
x q[1];
rz(0.67811857) q[2];
sx q[2];
rz(-1.0762012) q[2];
sx q[2];
rz(-1.2505609) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.1360201) q[1];
sx q[1];
rz(-0.51678777) q[1];
sx q[1];
rz(2.6173475) q[1];
rz(-pi) q[2];
rz(2.5902254) q[3];
sx q[3];
rz(-0.15227642) q[3];
sx q[3];
rz(2.4874529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.29844478) q[2];
sx q[2];
rz(-1.2848022) q[2];
sx q[2];
rz(1.0996381) q[2];
rz(3.1067276) q[3];
sx q[3];
rz(-0.74368447) q[3];
sx q[3];
rz(0.59665027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.027503969) q[0];
sx q[0];
rz(-1.1570258) q[0];
sx q[0];
rz(1.9644568) q[0];
rz(1.6674532) q[1];
sx q[1];
rz(-2.2087704) q[1];
sx q[1];
rz(1.4449545) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5367573) q[0];
sx q[0];
rz(-2.6492181) q[0];
sx q[0];
rz(2.5144666) q[0];
rz(-pi) q[1];
rz(-2.3722052) q[2];
sx q[2];
rz(-1.1379998) q[2];
sx q[2];
rz(2.3355049) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.34945378) q[1];
sx q[1];
rz(-1.0568403) q[1];
sx q[1];
rz(3.0626014) q[1];
rz(-0.8952462) q[3];
sx q[3];
rz(-1.352868) q[3];
sx q[3];
rz(2.982614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.16704796) q[2];
sx q[2];
rz(-1.5844774) q[2];
sx q[2];
rz(-0.17189279) q[2];
rz(-0.78957549) q[3];
sx q[3];
rz(-2.004345) q[3];
sx q[3];
rz(-1.3353698) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0852785) q[0];
sx q[0];
rz(-2.9264937) q[0];
sx q[0];
rz(2.6000182) q[0];
rz(-1.9821292) q[1];
sx q[1];
rz(-1.8340725) q[1];
sx q[1];
rz(2.3084739) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3726927) q[0];
sx q[0];
rz(-1.8463377) q[0];
sx q[0];
rz(2.5252456) q[0];
x q[1];
rz(0.54260268) q[2];
sx q[2];
rz(-1.2719947) q[2];
sx q[2];
rz(-1.8712559) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.45634746) q[1];
sx q[1];
rz(-0.73328555) q[1];
sx q[1];
rz(-1.228778) q[1];
rz(-pi) q[2];
x q[2];
rz(0.0019440193) q[3];
sx q[3];
rz(-2.2984347) q[3];
sx q[3];
rz(1.9950642) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.1200166) q[2];
sx q[2];
rz(-0.77777255) q[2];
sx q[2];
rz(-0.94338256) q[2];
rz(1.0776445) q[3];
sx q[3];
rz(-1.5543944) q[3];
sx q[3];
rz(-0.32850346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42644603) q[0];
sx q[0];
rz(-2.1052512) q[0];
sx q[0];
rz(-0.25682009) q[0];
rz(2.9019451) q[1];
sx q[1];
rz(-0.9050723) q[1];
sx q[1];
rz(-1.1368652) q[1];
rz(2.4478365) q[2];
sx q[2];
rz(-1.0508595) q[2];
sx q[2];
rz(-1.6169242) q[2];
rz(0.31044337) q[3];
sx q[3];
rz(-1.8127999) q[3];
sx q[3];
rz(-1.8565069) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
