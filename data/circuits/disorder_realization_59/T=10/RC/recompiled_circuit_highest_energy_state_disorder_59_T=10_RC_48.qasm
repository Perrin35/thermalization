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
rz(1.264313) q[0];
sx q[0];
rz(-2.8180583) q[0];
sx q[0];
rz(-2.6927595) q[0];
rz(1.9857061) q[1];
sx q[1];
rz(-1.356025) q[1];
sx q[1];
rz(1.1745656) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6960775) q[0];
sx q[0];
rz(-2.4598274) q[0];
sx q[0];
rz(2.5755139) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5023793) q[2];
sx q[2];
rz(-1.8937773) q[2];
sx q[2];
rz(0.078106192) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.2043031) q[1];
sx q[1];
rz(-0.95033344) q[1];
sx q[1];
rz(0.043067769) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6714736) q[3];
sx q[3];
rz(-1.0591456) q[3];
sx q[3];
rz(-0.6938405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.590098) q[2];
sx q[2];
rz(-1.8415201) q[2];
sx q[2];
rz(-0.91280118) q[2];
rz(1.270795) q[3];
sx q[3];
rz(-1.0400892) q[3];
sx q[3];
rz(0.84996581) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1368644) q[0];
sx q[0];
rz(-2.5623463) q[0];
sx q[0];
rz(2.7194523) q[0];
rz(-2.7719851) q[1];
sx q[1];
rz(-1.0969176) q[1];
sx q[1];
rz(2.2014528) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9418798) q[0];
sx q[0];
rz(-1.5732978) q[0];
sx q[0];
rz(-1.9002761) q[0];
rz(-1.0375848) q[2];
sx q[2];
rz(-1.4681446) q[2];
sx q[2];
rz(0.29666049) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.861189) q[1];
sx q[1];
rz(-2.197816) q[1];
sx q[1];
rz(-2.1908733) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2359911) q[3];
sx q[3];
rz(-0.70476092) q[3];
sx q[3];
rz(-2.4089264) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.1993572) q[2];
sx q[2];
rz(-2.4713559) q[2];
sx q[2];
rz(-0.90636903) q[2];
rz(-0.97936112) q[3];
sx q[3];
rz(-0.40458471) q[3];
sx q[3];
rz(-1.2649068) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
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
rz(-1.1449428) q[0];
sx q[0];
rz(-3.0746089) q[0];
sx q[0];
rz(-1.2445194) q[0];
rz(2.7527346) q[1];
sx q[1];
rz(-1.2797979) q[1];
sx q[1];
rz(-2.8188474) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8334482) q[0];
sx q[0];
rz(-2.2287773) q[0];
sx q[0];
rz(-1.7127772) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1305442) q[2];
sx q[2];
rz(-1.8343535) q[2];
sx q[2];
rz(1.6853058) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.1908993) q[1];
sx q[1];
rz(-0.79792332) q[1];
sx q[1];
rz(-1.2446333) q[1];
rz(0.1348008) q[3];
sx q[3];
rz(-0.35405891) q[3];
sx q[3];
rz(-1.1969138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.9129703) q[2];
sx q[2];
rz(-1.3801489) q[2];
sx q[2];
rz(1.3529533) q[2];
rz(0.27255034) q[3];
sx q[3];
rz(-1.291178) q[3];
sx q[3];
rz(-1.4636309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6005741) q[0];
sx q[0];
rz(-2.8193642) q[0];
sx q[0];
rz(-1.898265) q[0];
rz(2.2728032) q[1];
sx q[1];
rz(-2.181535) q[1];
sx q[1];
rz(1.0702081) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5015429) q[0];
sx q[0];
rz(-0.86210712) q[0];
sx q[0];
rz(-2.2556644) q[0];
x q[1];
rz(1.1203499) q[2];
sx q[2];
rz(-2.3131158) q[2];
sx q[2];
rz(2.7591005) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.8319885) q[1];
sx q[1];
rz(-1.2491396) q[1];
sx q[1];
rz(3.020806) q[1];
rz(-2.7989705) q[3];
sx q[3];
rz(-1.3903769) q[3];
sx q[3];
rz(-0.29553771) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.59226817) q[2];
sx q[2];
rz(-0.78846875) q[2];
sx q[2];
rz(2.2198524) q[2];
rz(-2.8978469) q[3];
sx q[3];
rz(-1.4344401) q[3];
sx q[3];
rz(-0.29269472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.155596) q[0];
sx q[0];
rz(-1.3256185) q[0];
sx q[0];
rz(-2.6407114) q[0];
rz(-2.0241375) q[1];
sx q[1];
rz(-1.3331579) q[1];
sx q[1];
rz(2.8270328) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9524744) q[0];
sx q[0];
rz(-2.6388612) q[0];
sx q[0];
rz(-2.377654) q[0];
rz(-pi) q[1];
x q[1];
rz(3.1371181) q[2];
sx q[2];
rz(-1.5273849) q[2];
sx q[2];
rz(-2.2501066) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.3479473) q[1];
sx q[1];
rz(-1.0255019) q[1];
sx q[1];
rz(2.6752363) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0231832) q[3];
sx q[3];
rz(-1.1573912) q[3];
sx q[3];
rz(-1.6550696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.8353117) q[2];
sx q[2];
rz(-1.1113144) q[2];
sx q[2];
rz(-1.296898) q[2];
rz(-0.48526192) q[3];
sx q[3];
rz(-1.5819712) q[3];
sx q[3];
rz(2.0280793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49419633) q[0];
sx q[0];
rz(-1.8504471) q[0];
sx q[0];
rz(3.0418292) q[0];
rz(1.8708723) q[1];
sx q[1];
rz(-2.2427509) q[1];
sx q[1];
rz(-1.3894003) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64369624) q[0];
sx q[0];
rz(-2.8262016) q[0];
sx q[0];
rz(2.5853378) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6011141) q[2];
sx q[2];
rz(-0.20217824) q[2];
sx q[2];
rz(1.0982738) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.278881) q[1];
sx q[1];
rz(-2.2219255) q[1];
sx q[1];
rz(-2.34823) q[1];
rz(-pi) q[2];
rz(-0.12677315) q[3];
sx q[3];
rz(-2.354458) q[3];
sx q[3];
rz(-1.8925557) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6483267) q[2];
sx q[2];
rz(-0.47097012) q[2];
sx q[2];
rz(-0.89034447) q[2];
rz(-1.1912311) q[3];
sx q[3];
rz(-1.8343364) q[3];
sx q[3];
rz(0.90203917) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4297727) q[0];
sx q[0];
rz(-3.0968102) q[0];
sx q[0];
rz(3.1088767) q[0];
rz(0.50321594) q[1];
sx q[1];
rz(-1.0992173) q[1];
sx q[1];
rz(0.12282664) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4361849) q[0];
sx q[0];
rz(-2.2741973) q[0];
sx q[0];
rz(-1.1421655) q[0];
rz(-pi) q[1];
rz(1.9101479) q[2];
sx q[2];
rz(-2.0464346) q[2];
sx q[2];
rz(0.95107691) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.5090774) q[1];
sx q[1];
rz(-1.6700524) q[1];
sx q[1];
rz(-2.0429913) q[1];
rz(-pi) q[2];
x q[2];
rz(0.35251725) q[3];
sx q[3];
rz(-1.4488701) q[3];
sx q[3];
rz(0.33559775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.34933019) q[2];
sx q[2];
rz(-2.1749039) q[2];
sx q[2];
rz(1.1518504) q[2];
rz(-2.2564015) q[3];
sx q[3];
rz(-1.3693634) q[3];
sx q[3];
rz(1.7159897) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.952878) q[0];
sx q[0];
rz(-2.2080244) q[0];
sx q[0];
rz(0.92920148) q[0];
rz(0.62905351) q[1];
sx q[1];
rz(-1.355143) q[1];
sx q[1];
rz(-0.74310511) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4152195) q[0];
sx q[0];
rz(-2.6872244) q[0];
sx q[0];
rz(-1.3903862) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.17878187) q[2];
sx q[2];
rz(-0.4950854) q[2];
sx q[2];
rz(0.21735969) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.656132) q[1];
sx q[1];
rz(-0.4194057) q[1];
sx q[1];
rz(-2.9784883) q[1];
rz(-2.6499641) q[3];
sx q[3];
rz(-2.1475882) q[3];
sx q[3];
rz(-0.73161179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.5653845) q[2];
sx q[2];
rz(-0.78018633) q[2];
sx q[2];
rz(-2.4156477) q[2];
rz(-2.7706326) q[3];
sx q[3];
rz(-1.5434664) q[3];
sx q[3];
rz(-0.36271873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
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
rz(-0.87042701) q[0];
sx q[0];
rz(-2.3112264) q[0];
sx q[0];
rz(2.5776432) q[0];
rz(-1.9922527) q[1];
sx q[1];
rz(-2.1795858) q[1];
sx q[1];
rz(2.3990778) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19631736) q[0];
sx q[0];
rz(-0.68500297) q[0];
sx q[0];
rz(-0.79010583) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7269823) q[2];
sx q[2];
rz(-1.3002204) q[2];
sx q[2];
rz(-1.0502349) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.2143605) q[1];
sx q[1];
rz(-0.4402658) q[1];
sx q[1];
rz(-0.10378598) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.94007317) q[3];
sx q[3];
rz(-0.82118644) q[3];
sx q[3];
rz(-1.8668777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5138862) q[2];
sx q[2];
rz(-2.0676925) q[2];
sx q[2];
rz(-2.7970496) q[2];
rz(2.6978317) q[3];
sx q[3];
rz(-2.0659476) q[3];
sx q[3];
rz(-1.8189323) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
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
rz(-0.76170707) q[0];
sx q[0];
rz(-0.2825309) q[0];
sx q[0];
rz(-2.479582) q[0];
rz(-2.8061891) q[1];
sx q[1];
rz(-1.4619724) q[1];
sx q[1];
rz(-0.9368771) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10231384) q[0];
sx q[0];
rz(-2.1697308) q[0];
sx q[0];
rz(0.44012578) q[0];
rz(-pi) q[1];
rz(0.023597864) q[2];
sx q[2];
rz(-1.840971) q[2];
sx q[2];
rz(-2.9058356) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.7171927) q[1];
sx q[1];
rz(-2.6584627) q[1];
sx q[1];
rz(-0.65349726) q[1];
x q[2];
rz(0.85031894) q[3];
sx q[3];
rz(-2.8788044) q[3];
sx q[3];
rz(0.14498728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.718049) q[2];
sx q[2];
rz(-0.23423883) q[2];
sx q[2];
rz(2.6424109) q[2];
rz(1.6792123) q[3];
sx q[3];
rz(-1.1464109) q[3];
sx q[3];
rz(1.6818989) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.725631) q[0];
sx q[0];
rz(-1.220663) q[0];
sx q[0];
rz(2.2882373) q[0];
rz(0.60494963) q[1];
sx q[1];
rz(-1.9974983) q[1];
sx q[1];
rz(-2.6430184) q[1];
rz(0.22357438) q[2];
sx q[2];
rz(-0.58453544) q[2];
sx q[2];
rz(0.75605308) q[2];
rz(2.6400186) q[3];
sx q[3];
rz(-1.4091103) q[3];
sx q[3];
rz(1.4179358) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
