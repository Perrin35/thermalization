OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-3.0225585) q[0];
sx q[0];
rz(-1.1086388) q[0];
sx q[0];
rz(-2.570987) q[0];
rz(-1.5174436) q[1];
sx q[1];
rz(-2.1921506) q[1];
sx q[1];
rz(1.5631262) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7279116) q[0];
sx q[0];
rz(-1.1360347) q[0];
sx q[0];
rz(-0.12760166) q[0];
rz(-pi) q[1];
x q[1];
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
rz(-pi) q[0];
x q[0];
rz(1.6001498) q[1];
sx q[1];
rz(-0.41473103) q[1];
sx q[1];
rz(0.45066898) q[1];
rz(-pi) q[2];
rz(-2.6622979) q[3];
sx q[3];
rz(-2.473978) q[3];
sx q[3];
rz(-1.3541753) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.40330517) q[2];
sx q[2];
rz(-0.83633509) q[2];
sx q[2];
rz(-2.9264012) q[2];
rz(-3.134794) q[3];
sx q[3];
rz(-2.2577622) q[3];
sx q[3];
rz(2.1782037) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6007518) q[0];
sx q[0];
rz(-0.90966666) q[0];
sx q[0];
rz(-1.2700861) q[0];
rz(-1.0441682) q[1];
sx q[1];
rz(-2.7296598) q[1];
sx q[1];
rz(-0.56234908) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5836307) q[0];
sx q[0];
rz(-1.6837801) q[0];
sx q[0];
rz(1.8336748) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9980045) q[2];
sx q[2];
rz(-2.1242122) q[2];
sx q[2];
rz(-1.9678022) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.085155) q[1];
sx q[1];
rz(-0.83806935) q[1];
sx q[1];
rz(-1.4895658) q[1];
rz(-pi) q[2];
rz(-2.8329912) q[3];
sx q[3];
rz(-1.2720593) q[3];
sx q[3];
rz(1.5200391) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.1336141) q[2];
sx q[2];
rz(-0.14432159) q[2];
sx q[2];
rz(3.0001384) q[2];
rz(2.134363) q[3];
sx q[3];
rz(-1.4814582) q[3];
sx q[3];
rz(-3.0032515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7815642) q[0];
sx q[0];
rz(-1.5046295) q[0];
sx q[0];
rz(-3.0887399) q[0];
rz(-1.7618893) q[1];
sx q[1];
rz(-2.3638937) q[1];
sx q[1];
rz(-2.8676829) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6099458) q[0];
sx q[0];
rz(-1.393712) q[0];
sx q[0];
rz(-3.1041932) q[0];
x q[1];
rz(1.4714965) q[2];
sx q[2];
rz(-0.75649951) q[2];
sx q[2];
rz(-1.4425636) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.2347331) q[1];
sx q[1];
rz(-2.3409493) q[1];
sx q[1];
rz(-1.6080286) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9778887) q[3];
sx q[3];
rz(-2.8559167) q[3];
sx q[3];
rz(0.33200246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.86458796) q[2];
sx q[2];
rz(-1.4952345) q[2];
sx q[2];
rz(0.60025364) q[2];
rz(2.4872335) q[3];
sx q[3];
rz(-2.8221966) q[3];
sx q[3];
rz(-1.6681685) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7905423) q[0];
sx q[0];
rz(-1.496614) q[0];
sx q[0];
rz(-0.45781621) q[0];
rz(-2.5425725) q[1];
sx q[1];
rz(-0.51282239) q[1];
sx q[1];
rz(-0.27214989) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.015135678) q[0];
sx q[0];
rz(-1.8227856) q[0];
sx q[0];
rz(0.62313764) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1013012) q[2];
sx q[2];
rz(-2.0397503) q[2];
sx q[2];
rz(-2.1877642) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6451431) q[1];
sx q[1];
rz(-1.4525868) q[1];
sx q[1];
rz(0.18564863) q[1];
rz(-pi) q[2];
rz(1.3982716) q[3];
sx q[3];
rz(-0.93768812) q[3];
sx q[3];
rz(2.5987491) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.6279471) q[2];
sx q[2];
rz(-2.7529035) q[2];
sx q[2];
rz(1.909931) q[2];
rz(1.7197459) q[3];
sx q[3];
rz(-1.809241) q[3];
sx q[3];
rz(0.96760145) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49895367) q[0];
sx q[0];
rz(-3.0727144) q[0];
sx q[0];
rz(2.4054085) q[0];
rz(-2.5040297) q[1];
sx q[1];
rz(-1.9393189) q[1];
sx q[1];
rz(-0.41086248) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6402886) q[0];
sx q[0];
rz(-2.023847) q[0];
sx q[0];
rz(-0.52877063) q[0];
rz(-2.9301398) q[2];
sx q[2];
rz(-1.0120858) q[2];
sx q[2];
rz(-1.3635456) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.3113271) q[1];
sx q[1];
rz(-0.62957785) q[1];
sx q[1];
rz(-0.011407995) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.75848364) q[3];
sx q[3];
rz(-2.3173769) q[3];
sx q[3];
rz(0.099140204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.32834443) q[2];
sx q[2];
rz(-1.5466377) q[2];
sx q[2];
rz(0.52097875) q[2];
rz(0.10884918) q[3];
sx q[3];
rz(-1.0131001) q[3];
sx q[3];
rz(-2.525135) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2206889) q[0];
sx q[0];
rz(-1.9158798) q[0];
sx q[0];
rz(-1.8826245) q[0];
rz(2.068223) q[1];
sx q[1];
rz(-1.4732889) q[1];
sx q[1];
rz(2.4980256) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9285842) q[0];
sx q[0];
rz(-1.1197392) q[0];
sx q[0];
rz(-1.5764007) q[0];
rz(-0.50711716) q[2];
sx q[2];
rz(-2.7165453) q[2];
sx q[2];
rz(-0.054946446) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.7593708) q[1];
sx q[1];
rz(-2.367451) q[1];
sx q[1];
rz(-1.8060807) q[1];
x q[2];
rz(0.57585133) q[3];
sx q[3];
rz(-1.1237877) q[3];
sx q[3];
rz(-0.67547638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.7966938) q[2];
sx q[2];
rz(-1.0641791) q[2];
sx q[2];
rz(1.8033484) q[2];
rz(-1.4517339) q[3];
sx q[3];
rz(-1.2029878) q[3];
sx q[3];
rz(3.0411327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79784262) q[0];
sx q[0];
rz(-1.2748953) q[0];
sx q[0];
rz(-0.95329681) q[0];
rz(-2.9816309) q[1];
sx q[1];
rz(-2.5123031) q[1];
sx q[1];
rz(2.2728641) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1686483) q[0];
sx q[0];
rz(-0.7059083) q[0];
sx q[0];
rz(-1.2337408) q[0];
rz(-pi) q[1];
x q[1];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9142469) q[1];
sx q[1];
rz(-0.29623756) q[1];
sx q[1];
rz(-0.072320894) q[1];
rz(-2.2420011) q[3];
sx q[3];
rz(-1.7760522) q[3];
sx q[3];
rz(-1.156607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.4538883) q[2];
sx q[2];
rz(-1.1331465) q[2];
sx q[2];
rz(0.49393168) q[2];
rz(-1.4880873) q[3];
sx q[3];
rz(-2.2736214) q[3];
sx q[3];
rz(-0.021763703) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9138551) q[0];
sx q[0];
rz(-1.7182925) q[0];
sx q[0];
rz(2.8505351) q[0];
rz(-3.011138) q[1];
sx q[1];
rz(-0.37965241) q[1];
sx q[1];
rz(-1.5721273) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0091189) q[0];
sx q[0];
rz(-1.557701) q[0];
sx q[0];
rz(1.520169) q[0];
x q[1];
rz(-1.7849017) q[2];
sx q[2];
rz(-2.4567025) q[2];
sx q[2];
rz(0.30944007) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.56640139) q[1];
sx q[1];
rz(-0.960604) q[1];
sx q[1];
rz(-1.7119049) q[1];
x q[2];
rz(0.30867851) q[3];
sx q[3];
rz(-1.6925081) q[3];
sx q[3];
rz(2.6409923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.7743885) q[2];
sx q[2];
rz(-2.9640894) q[2];
sx q[2];
rz(-3.0299178) q[2];
rz(0.81237927) q[3];
sx q[3];
rz(-1.2923391) q[3];
sx q[3];
rz(-1.6151379) q[3];
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
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77686253) q[0];
sx q[0];
rz(-0.22219816) q[0];
sx q[0];
rz(0.36866933) q[0];
rz(2.5313976) q[1];
sx q[1];
rz(-1.5564324) q[1];
sx q[1];
rz(-2.1048224) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0087413) q[0];
sx q[0];
rz(-1.794882) q[0];
sx q[0];
rz(0.56613825) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5794994) q[2];
sx q[2];
rz(-2.6285951) q[2];
sx q[2];
rz(2.7287116) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.887275) q[1];
sx q[1];
rz(-2.3000082) q[1];
sx q[1];
rz(-2.9542854) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0572564) q[3];
sx q[3];
rz(-1.7170625) q[3];
sx q[3];
rz(-1.4649966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.5165575) q[2];
sx q[2];
rz(-2.3380307) q[2];
sx q[2];
rz(2.6885923) q[2];
rz(0.9520483) q[3];
sx q[3];
rz(-2.0249849) q[3];
sx q[3];
rz(3.0208352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9348738) q[0];
sx q[0];
rz(-0.047858866) q[0];
sx q[0];
rz(-2.5913443) q[0];
rz(2.0568559) q[1];
sx q[1];
rz(-2.5526498) q[1];
sx q[1];
rz(-1.750607) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7511786) q[0];
sx q[0];
rz(-2.0719647) q[0];
sx q[0];
rz(0.044981754) q[0];
x q[1];
rz(-2.1065305) q[2];
sx q[2];
rz(-1.7155871) q[2];
sx q[2];
rz(2.5404251) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.7470737) q[1];
sx q[1];
rz(-0.84861524) q[1];
sx q[1];
rz(-2.4492521) q[1];
rz(-pi) q[2];
rz(-0.74446328) q[3];
sx q[3];
rz(-0.083221952) q[3];
sx q[3];
rz(-1.2017488) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.7636106) q[2];
sx q[2];
rz(-0.72145975) q[2];
sx q[2];
rz(1.9788431) q[2];
rz(-0.55443305) q[3];
sx q[3];
rz(-0.88397637) q[3];
sx q[3];
rz(-2.673431) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32600317) q[0];
sx q[0];
rz(-1.3023888) q[0];
sx q[0];
rz(1.6525255) q[0];
rz(-0.21492699) q[1];
sx q[1];
rz(-0.88354127) q[1];
sx q[1];
rz(3.0149928) q[1];
rz(1.3431637) q[2];
sx q[2];
rz(-0.91522436) q[2];
sx q[2];
rz(2.1686423) q[2];
rz(-0.074474143) q[3];
sx q[3];
rz(-2.4182033) q[3];
sx q[3];
rz(-0.084666336) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
