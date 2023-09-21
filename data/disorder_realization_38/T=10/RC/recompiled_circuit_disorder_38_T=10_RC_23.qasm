OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.71198553) q[0];
sx q[0];
rz(-2.7349732) q[0];
sx q[0];
rz(-0.24917319) q[0];
rz(3.0781526) q[1];
sx q[1];
rz(-0.97172207) q[1];
sx q[1];
rz(2.5914153) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7101947) q[0];
sx q[0];
rz(-0.22127998) q[0];
sx q[0];
rz(-1.6943323) q[0];
x q[1];
rz(1.0371738) q[2];
sx q[2];
rz(-1.954477) q[2];
sx q[2];
rz(1.8431611) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.6018511) q[1];
sx q[1];
rz(-1.1907693) q[1];
sx q[1];
rz(2.3439581) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1816741) q[3];
sx q[3];
rz(-1.5023408) q[3];
sx q[3];
rz(-1.5271036) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.41574079) q[2];
sx q[2];
rz(-2.6927413) q[2];
sx q[2];
rz(0.63981167) q[2];
rz(2.2885627) q[3];
sx q[3];
rz(-2.5363688) q[3];
sx q[3];
rz(-2.7602592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8752276) q[0];
sx q[0];
rz(-1.9748283) q[0];
sx q[0];
rz(0.27045989) q[0];
rz(-2.4282783) q[1];
sx q[1];
rz(-1.0353054) q[1];
sx q[1];
rz(-1.5126022) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0892031) q[0];
sx q[0];
rz(-2.4501778) q[0];
sx q[0];
rz(-1.0538488) q[0];
rz(-pi) q[1];
rz(-0.011629148) q[2];
sx q[2];
rz(-0.28380576) q[2];
sx q[2];
rz(-1.1139882) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.70817972) q[1];
sx q[1];
rz(-0.69060329) q[1];
sx q[1];
rz(-1.3604926) q[1];
rz(-pi) q[2];
rz(2.9019722) q[3];
sx q[3];
rz(-1.8488956) q[3];
sx q[3];
rz(-2.3688189) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.9397395) q[2];
sx q[2];
rz(-2.8994603) q[2];
sx q[2];
rz(2.3382323) q[2];
rz(-1.057829) q[3];
sx q[3];
rz(-1.4927031) q[3];
sx q[3];
rz(0.025618205) q[3];
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
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2597044) q[0];
sx q[0];
rz(-1.0972247) q[0];
sx q[0];
rz(2.846068) q[0];
rz(-2.9064536) q[1];
sx q[1];
rz(-1.7087015) q[1];
sx q[1];
rz(-2.3957516) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.196516) q[0];
sx q[0];
rz(-1.7917504) q[0];
sx q[0];
rz(3.0483079) q[0];
rz(1.7350115) q[2];
sx q[2];
rz(-2.0914372) q[2];
sx q[2];
rz(2.8525713) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.97400372) q[1];
sx q[1];
rz(-2.0473192) q[1];
sx q[1];
rz(1.2567026) q[1];
rz(-0.45332076) q[3];
sx q[3];
rz(-1.7954149) q[3];
sx q[3];
rz(2.2090705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.088034257) q[2];
sx q[2];
rz(-0.025667889) q[2];
sx q[2];
rz(0.68874613) q[2];
rz(-0.050343242) q[3];
sx q[3];
rz(-2.2294932) q[3];
sx q[3];
rz(1.5215727) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-0.21928366) q[0];
sx q[0];
rz(-2.121121) q[0];
sx q[0];
rz(3.0396089) q[0];
rz(0.12022262) q[1];
sx q[1];
rz(-0.52821237) q[1];
sx q[1];
rz(-0.27329683) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1647987) q[0];
sx q[0];
rz(-1.8935888) q[0];
sx q[0];
rz(0.093683634) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0068514) q[2];
sx q[2];
rz(-1.0014921) q[2];
sx q[2];
rz(-0.068892613) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.30444333) q[1];
sx q[1];
rz(-1.7177561) q[1];
sx q[1];
rz(-1.8069581) q[1];
rz(1.7498383) q[3];
sx q[3];
rz(-0.92902196) q[3];
sx q[3];
rz(-2.2620576) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3499202) q[2];
sx q[2];
rz(-2.4035954) q[2];
sx q[2];
rz(0.50393528) q[2];
rz(3.062011) q[3];
sx q[3];
rz(-1.9888398) q[3];
sx q[3];
rz(2.7664405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2824317) q[0];
sx q[0];
rz(-2.0895884) q[0];
sx q[0];
rz(3.058847) q[0];
rz(0.67963183) q[1];
sx q[1];
rz(-1.6487164) q[1];
sx q[1];
rz(0.98714978) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7894831) q[0];
sx q[0];
rz(-1.1687359) q[0];
sx q[0];
rz(0.91870086) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.64381386) q[2];
sx q[2];
rz(-1.5783196) q[2];
sx q[2];
rz(1.3155589) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.0070237006) q[1];
sx q[1];
rz(-2.1937074) q[1];
sx q[1];
rz(-1.9593777) q[1];
x q[2];
rz(1.4694674) q[3];
sx q[3];
rz(-0.39468995) q[3];
sx q[3];
rz(-2.8916388) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.8905028) q[2];
sx q[2];
rz(-0.7901929) q[2];
sx q[2];
rz(0.21128543) q[2];
rz(-2.7206897) q[3];
sx q[3];
rz(-0.55287164) q[3];
sx q[3];
rz(3.0781854) q[3];
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
rz(-1.485065) q[0];
sx q[0];
rz(-1.8873029) q[0];
sx q[0];
rz(0.791839) q[0];
rz(0.99545288) q[1];
sx q[1];
rz(-2.1738926) q[1];
sx q[1];
rz(-1.2794367) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6364481) q[0];
sx q[0];
rz(-2.8594058) q[0];
sx q[0];
rz(1.6237153) q[0];
rz(-pi) q[1];
rz(-1.860414) q[2];
sx q[2];
rz(-2.1589303) q[2];
sx q[2];
rz(-0.16298018) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.4931603) q[1];
sx q[1];
rz(-1.5409768) q[1];
sx q[1];
rz(2.0926507) q[1];
x q[2];
rz(-1.7518696) q[3];
sx q[3];
rz(-1.4317792) q[3];
sx q[3];
rz(0.74336038) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.8906158) q[2];
sx q[2];
rz(-2.1112517) q[2];
sx q[2];
rz(-2.3708169) q[2];
rz(-1.4701014) q[3];
sx q[3];
rz(-0.41430587) q[3];
sx q[3];
rz(-2.9582086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32271785) q[0];
sx q[0];
rz(-2.4920431) q[0];
sx q[0];
rz(-3.0859257) q[0];
rz(0.21559134) q[1];
sx q[1];
rz(-2.3781653) q[1];
sx q[1];
rz(-3.1380222) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5033747) q[0];
sx q[0];
rz(-0.57045454) q[0];
sx q[0];
rz(-2.4128782) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1793169) q[2];
sx q[2];
rz(-0.85748312) q[2];
sx q[2];
rz(1.6187402) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.6047302) q[1];
sx q[1];
rz(-2.4559048) q[1];
sx q[1];
rz(2.0525949) q[1];
x q[2];
rz(-1.9332063) q[3];
sx q[3];
rz(-0.21167314) q[3];
sx q[3];
rz(-0.025312245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5443762) q[2];
sx q[2];
rz(-2.1901013) q[2];
sx q[2];
rz(-0.92010951) q[2];
rz(0.19872935) q[3];
sx q[3];
rz(-1.2343497) q[3];
sx q[3];
rz(-2.7238817) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51628095) q[0];
sx q[0];
rz(-1.5100864) q[0];
sx q[0];
rz(-1.7238808) q[0];
rz(0.40813804) q[1];
sx q[1];
rz(-2.1022271) q[1];
sx q[1];
rz(-0.67869854) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40020254) q[0];
sx q[0];
rz(-2.7149902) q[0];
sx q[0];
rz(0.68047561) q[0];
rz(-pi) q[1];
rz(-2.3085824) q[2];
sx q[2];
rz(-1.9265124) q[2];
sx q[2];
rz(-0.21258159) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.24819599) q[1];
sx q[1];
rz(-2.0940603) q[1];
sx q[1];
rz(3.0688973) q[1];
rz(1.4123165) q[3];
sx q[3];
rz(-1.9915808) q[3];
sx q[3];
rz(0.63384038) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.9644908) q[2];
sx q[2];
rz(-2.047838) q[2];
sx q[2];
rz(2.2237681) q[2];
rz(1.142189) q[3];
sx q[3];
rz(-0.85421383) q[3];
sx q[3];
rz(2.2043622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.98638242) q[0];
sx q[0];
rz(-0.6739524) q[0];
sx q[0];
rz(-2.881799) q[0];
rz(0.70867509) q[1];
sx q[1];
rz(-2.8627113) q[1];
sx q[1];
rz(-2.6146467) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6078867) q[0];
sx q[0];
rz(-0.41782802) q[0];
sx q[0];
rz(-2.582344) q[0];
x q[1];
rz(-1.2277463) q[2];
sx q[2];
rz(-0.68708778) q[2];
sx q[2];
rz(0.81817852) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.2745797) q[1];
sx q[1];
rz(-1.8362987) q[1];
sx q[1];
rz(-0.32151476) q[1];
rz(-pi) q[2];
rz(-0.44550495) q[3];
sx q[3];
rz(-1.8337436) q[3];
sx q[3];
rz(-1.477369) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.8156585) q[2];
sx q[2];
rz(-2.2951173) q[2];
sx q[2];
rz(2.3507067) q[2];
rz(0.38665006) q[3];
sx q[3];
rz(-1.0145885) q[3];
sx q[3];
rz(-0.33716831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40795046) q[0];
sx q[0];
rz(-2.9691417) q[0];
sx q[0];
rz(-2.1561484) q[0];
rz(0.5685637) q[1];
sx q[1];
rz(-1.111258) q[1];
sx q[1];
rz(-0.3607761) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0246968) q[0];
sx q[0];
rz(-1.403406) q[0];
sx q[0];
rz(-0.025339729) q[0];
x q[1];
rz(-1.0656409) q[2];
sx q[2];
rz(-0.2444707) q[2];
sx q[2];
rz(0.62015647) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.1211179) q[1];
sx q[1];
rz(-1.4282244) q[1];
sx q[1];
rz(2.3974182) q[1];
x q[2];
rz(-1.4233227) q[3];
sx q[3];
rz(-1.4031938) q[3];
sx q[3];
rz(-0.22351219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.0976022) q[2];
sx q[2];
rz(-2.4520935) q[2];
sx q[2];
rz(-0.94341755) q[2];
rz(-0.57389456) q[3];
sx q[3];
rz(-2.6608163) q[3];
sx q[3];
rz(-1.3963612) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5744793) q[0];
sx q[0];
rz(-1.6945101) q[0];
sx q[0];
rz(2.2816818) q[0];
rz(1.7815331) q[1];
sx q[1];
rz(-2.139745) q[1];
sx q[1];
rz(2.3812961) q[1];
rz(1.3654937) q[2];
sx q[2];
rz(-1.6663972) q[2];
sx q[2];
rz(-1.321928) q[2];
rz(-0.73137024) q[3];
sx q[3];
rz(-2.2435355) q[3];
sx q[3];
rz(0.13767195) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
