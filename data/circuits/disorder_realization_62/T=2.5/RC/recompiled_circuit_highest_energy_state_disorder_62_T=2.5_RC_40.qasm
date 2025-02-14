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
rz(-0.38464889) q[0];
sx q[0];
rz(0.83847133) q[0];
sx q[0];
rz(5.6738927) q[0];
rz(0.72120178) q[1];
sx q[1];
rz(5.3546049) q[1];
sx q[1];
rz(7.3794853) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.012714931) q[0];
sx q[0];
rz(-2.4372074) q[0];
sx q[0];
rz(-1.7449858) q[0];
rz(-pi) q[1];
rz(-2.6651971) q[2];
sx q[2];
rz(-0.98586776) q[2];
sx q[2];
rz(0.75854075) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.1692449) q[1];
sx q[1];
rz(-1.2903288) q[1];
sx q[1];
rz(2.601956) q[1];
rz(2.9864086) q[3];
sx q[3];
rz(-1.461457) q[3];
sx q[3];
rz(0.21462277) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.91813749) q[2];
sx q[2];
rz(-1.3104985) q[2];
sx q[2];
rz(-1.1892595) q[2];
rz(2.7506645) q[3];
sx q[3];
rz(-1.9783741) q[3];
sx q[3];
rz(2.0093567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.696058) q[0];
sx q[0];
rz(-1.7867418) q[0];
sx q[0];
rz(-1.6433486) q[0];
rz(2.4636726) q[1];
sx q[1];
rz(-1.3005715) q[1];
sx q[1];
rz(-0.71151412) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2145793) q[0];
sx q[0];
rz(-1.2138146) q[0];
sx q[0];
rz(-1.0735608) q[0];
x q[1];
rz(-0.46320148) q[2];
sx q[2];
rz(-1.3579662) q[2];
sx q[2];
rz(-2.6821399) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.1700279) q[1];
sx q[1];
rz(-1.4907279) q[1];
sx q[1];
rz(-2.7951334) q[1];
rz(-2.1944216) q[3];
sx q[3];
rz(-1.5452562) q[3];
sx q[3];
rz(-0.21835777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.70637643) q[2];
sx q[2];
rz(-1.7855568) q[2];
sx q[2];
rz(-2.0785275) q[2];
rz(-1.6478018) q[3];
sx q[3];
rz(-2.6726275) q[3];
sx q[3];
rz(2.4013605) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4110334) q[0];
sx q[0];
rz(-2.7506802) q[0];
sx q[0];
rz(2.539769) q[0];
rz(-0.14566323) q[1];
sx q[1];
rz(-1.0258976) q[1];
sx q[1];
rz(-2.513733) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2780298) q[0];
sx q[0];
rz(-0.13736996) q[0];
sx q[0];
rz(-2.6111215) q[0];
rz(-1.0192515) q[2];
sx q[2];
rz(-1.2639234) q[2];
sx q[2];
rz(0.17099525) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.2035297) q[1];
sx q[1];
rz(-2.7157686) q[1];
sx q[1];
rz(-0.99217207) q[1];
rz(0.46247356) q[3];
sx q[3];
rz(-2.1331027) q[3];
sx q[3];
rz(3.1218046) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.7566028) q[2];
sx q[2];
rz(-2.0774697) q[2];
sx q[2];
rz(-0.21558726) q[2];
rz(-2.0032517) q[3];
sx q[3];
rz(-1.3201069) q[3];
sx q[3];
rz(2.5707572) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0067921) q[0];
sx q[0];
rz(-1.726806) q[0];
sx q[0];
rz(0.11882812) q[0];
rz(1.4745332) q[1];
sx q[1];
rz(-0.88913616) q[1];
sx q[1];
rz(-2.3337591) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0519179) q[0];
sx q[0];
rz(-2.1362944) q[0];
sx q[0];
rz(-2.2321108) q[0];
x q[1];
rz(-0.045305552) q[2];
sx q[2];
rz(-1.8938365) q[2];
sx q[2];
rz(1.7957866) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4588759) q[1];
sx q[1];
rz(-1.4035051) q[1];
sx q[1];
rz(0.95750995) q[1];
rz(-pi) q[2];
rz(-0.019185493) q[3];
sx q[3];
rz(-0.67836715) q[3];
sx q[3];
rz(-2.8205591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.64041758) q[2];
sx q[2];
rz(-1.7265604) q[2];
sx q[2];
rz(0.51587063) q[2];
rz(3.0473895) q[3];
sx q[3];
rz(-2.3661864) q[3];
sx q[3];
rz(1.5883821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77189389) q[0];
sx q[0];
rz(-1.1006681) q[0];
sx q[0];
rz(-1.1871673) q[0];
rz(2.5502603) q[1];
sx q[1];
rz(-1.8449123) q[1];
sx q[1];
rz(2.3147413) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2255838) q[0];
sx q[0];
rz(-1.1904612) q[0];
sx q[0];
rz(1.4404141) q[0];
x q[1];
rz(-2.128878) q[2];
sx q[2];
rz(-1.6396171) q[2];
sx q[2];
rz(1.1281013) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.34870711) q[1];
sx q[1];
rz(-1.4408083) q[1];
sx q[1];
rz(-2.2568364) q[1];
rz(-pi) q[2];
x q[2];
rz(2.639797) q[3];
sx q[3];
rz(-2.083598) q[3];
sx q[3];
rz(0.38450188) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.7190711) q[2];
sx q[2];
rz(-2.8870388) q[2];
sx q[2];
rz(-1.0151601) q[2];
rz(0.90967956) q[3];
sx q[3];
rz(-1.9010952) q[3];
sx q[3];
rz(2.4801586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61843094) q[0];
sx q[0];
rz(-1.2105415) q[0];
sx q[0];
rz(-2.8316408) q[0];
rz(-0.018772086) q[1];
sx q[1];
rz(-1.680178) q[1];
sx q[1];
rz(-0.087336691) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76742889) q[0];
sx q[0];
rz(-1.6032251) q[0];
sx q[0];
rz(2.5686766) q[0];
rz(2.2703253) q[2];
sx q[2];
rz(-2.0972898) q[2];
sx q[2];
rz(0.73824182) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(3.0951335) q[1];
sx q[1];
rz(-1.1756056) q[1];
sx q[1];
rz(2.1161152) q[1];
rz(-pi) q[2];
rz(-1.0424006) q[3];
sx q[3];
rz(-1.434552) q[3];
sx q[3];
rz(2.6353177) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.35149082) q[2];
sx q[2];
rz(-2.6678706) q[2];
sx q[2];
rz(2.6215485) q[2];
rz(0.13488787) q[3];
sx q[3];
rz(-1.3210195) q[3];
sx q[3];
rz(1.409387) q[3];
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
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0670369) q[0];
sx q[0];
rz(-1.073607) q[0];
sx q[0];
rz(0.38597646) q[0];
rz(-1.0999673) q[1];
sx q[1];
rz(-2.4620582) q[1];
sx q[1];
rz(2.3540672) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58819492) q[0];
sx q[0];
rz(-0.46632354) q[0];
sx q[0];
rz(-2.2539) q[0];
rz(-pi) q[1];
rz(1.9454657) q[2];
sx q[2];
rz(-0.9866937) q[2];
sx q[2];
rz(-1.0874334) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.7907775) q[1];
sx q[1];
rz(-2.3190089) q[1];
sx q[1];
rz(-2.050553) q[1];
x q[2];
rz(0.94905628) q[3];
sx q[3];
rz(-0.28655401) q[3];
sx q[3];
rz(-1.7915981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.0218574) q[2];
sx q[2];
rz(-1.8334374) q[2];
sx q[2];
rz(1.5137399) q[2];
rz(1.2210023) q[3];
sx q[3];
rz(-1.1648213) q[3];
sx q[3];
rz(-0.47202078) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30962238) q[0];
sx q[0];
rz(-1.4291052) q[0];
sx q[0];
rz(0.23304644) q[0];
rz(-1.4876935) q[1];
sx q[1];
rz(-1.033604) q[1];
sx q[1];
rz(1.0071365) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8965573) q[0];
sx q[0];
rz(-2.5063337) q[0];
sx q[0];
rz(0.59478514) q[0];
rz(1.8677668) q[2];
sx q[2];
rz(-0.96189068) q[2];
sx q[2];
rz(1.1394274) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9527275) q[1];
sx q[1];
rz(-0.42546526) q[1];
sx q[1];
rz(-1.9317606) q[1];
rz(-2.8975186) q[3];
sx q[3];
rz(-1.5904038) q[3];
sx q[3];
rz(1.5527035) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2795589) q[2];
sx q[2];
rz(-2.0290012) q[2];
sx q[2];
rz(2.2588008) q[2];
rz(-0.27859303) q[3];
sx q[3];
rz(-1.168058) q[3];
sx q[3];
rz(0.45894233) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9423264) q[0];
sx q[0];
rz(-0.47573221) q[0];
sx q[0];
rz(-1.5698154) q[0];
rz(0.78872952) q[1];
sx q[1];
rz(-0.9659583) q[1];
sx q[1];
rz(0.68005651) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59916436) q[0];
sx q[0];
rz(-0.76272805) q[0];
sx q[0];
rz(0.086407089) q[0];
rz(-pi) q[1];
rz(-0.11876688) q[2];
sx q[2];
rz(-1.6965908) q[2];
sx q[2];
rz(-2.072352) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.9788821) q[1];
sx q[1];
rz(-1.0274402) q[1];
sx q[1];
rz(-1.1239978) q[1];
x q[2];
rz(-2.5933411) q[3];
sx q[3];
rz(-2.2202949) q[3];
sx q[3];
rz(-2.4626061) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.1532229) q[2];
sx q[2];
rz(-2.220927) q[2];
sx q[2];
rz(1.1663158) q[2];
rz(-1.9370646) q[3];
sx q[3];
rz(-1.322999) q[3];
sx q[3];
rz(-0.62674633) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6176497) q[0];
sx q[0];
rz(-2.2594422) q[0];
sx q[0];
rz(0.43706617) q[0];
rz(-2.3433459) q[1];
sx q[1];
rz(-0.50004807) q[1];
sx q[1];
rz(2.9214568) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1967314) q[0];
sx q[0];
rz(-1.5240977) q[0];
sx q[0];
rz(0.80094211) q[0];
rz(-2.9371168) q[2];
sx q[2];
rz(-1.8796953) q[2];
sx q[2];
rz(-2.1155865) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.70500532) q[1];
sx q[1];
rz(-1.844075) q[1];
sx q[1];
rz(0.80555861) q[1];
rz(-pi) q[2];
rz(0.80159558) q[3];
sx q[3];
rz(-1.6022046) q[3];
sx q[3];
rz(3.1289986) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2685214) q[2];
sx q[2];
rz(-1.2348509) q[2];
sx q[2];
rz(0.29042563) q[2];
rz(0.63646603) q[3];
sx q[3];
rz(-1.3822184) q[3];
sx q[3];
rz(1.9071473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-0.81454043) q[0];
sx q[0];
rz(-2.4335813) q[0];
sx q[0];
rz(-2.0306564) q[0];
rz(-0.064432714) q[1];
sx q[1];
rz(-1.671052) q[1];
sx q[1];
rz(-0.64238092) q[1];
rz(1.6173132) q[2];
sx q[2];
rz(-2.5462493) q[2];
sx q[2];
rz(0.52296849) q[2];
rz(-1.3051582) q[3];
sx q[3];
rz(-2.0943741) q[3];
sx q[3];
rz(1.5423273) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
