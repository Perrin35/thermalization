OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.6346729) q[0];
sx q[0];
rz(1.4022175) q[0];
sx q[0];
rz(13.316857) q[0];
rz(0.063272417) q[1];
sx q[1];
rz(-2.3228391) q[1];
sx q[1];
rz(-2.1225345) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69438231) q[0];
sx q[0];
rz(-1.3827494) q[0];
sx q[0];
rz(-0.5699372) q[0];
x q[1];
rz(1.8797714) q[2];
sx q[2];
rz(-0.7337386) q[2];
sx q[2];
rz(-0.66622615) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4338389) q[1];
sx q[1];
rz(-2.1496668) q[1];
sx q[1];
rz(1.7339091) q[1];
x q[2];
rz(-1.9492703) q[3];
sx q[3];
rz(-1.0590388) q[3];
sx q[3];
rz(0.056223362) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.2941306) q[2];
sx q[2];
rz(-0.99404603) q[2];
sx q[2];
rz(-1.6729986) q[2];
rz(-0.20644203) q[3];
sx q[3];
rz(-1.3279746) q[3];
sx q[3];
rz(2.4895721) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1714529) q[0];
sx q[0];
rz(-1.7296706) q[0];
sx q[0];
rz(-0.10301244) q[0];
rz(-0.39762321) q[1];
sx q[1];
rz(-0.42639521) q[1];
sx q[1];
rz(-2.2893589) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21450248) q[0];
sx q[0];
rz(-2.6047241) q[0];
sx q[0];
rz(-1.6887299) q[0];
rz(-pi) q[1];
rz(-1.909341) q[2];
sx q[2];
rz(-1.3024835) q[2];
sx q[2];
rz(1.2635096) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5168165) q[1];
sx q[1];
rz(-1.0472391) q[1];
sx q[1];
rz(0.015990301) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9217124) q[3];
sx q[3];
rz(-2.9438955) q[3];
sx q[3];
rz(1.3330595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.50920606) q[2];
sx q[2];
rz(-1.0078603) q[2];
sx q[2];
rz(2.3243135) q[2];
rz(-2.9362074) q[3];
sx q[3];
rz(-0.21586625) q[3];
sx q[3];
rz(-0.5425905) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83194724) q[0];
sx q[0];
rz(-1.0878071) q[0];
sx q[0];
rz(0.52016869) q[0];
rz(0.7236411) q[1];
sx q[1];
rz(-1.8579204) q[1];
sx q[1];
rz(2.9626194) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0332645) q[0];
sx q[0];
rz(-1.1305362) q[0];
sx q[0];
rz(2.5785793) q[0];
rz(2.6232324) q[2];
sx q[2];
rz(-0.39769618) q[2];
sx q[2];
rz(-0.30425554) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.8047898) q[1];
sx q[1];
rz(-1.8095142) q[1];
sx q[1];
rz(-1.6126339) q[1];
rz(-pi) q[2];
rz(-1.2988899) q[3];
sx q[3];
rz(-2.2103709) q[3];
sx q[3];
rz(3.0323882) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.8861063) q[2];
sx q[2];
rz(-0.28033689) q[2];
sx q[2];
rz(-2.172016) q[2];
rz(0.3420091) q[3];
sx q[3];
rz(-2.1186782) q[3];
sx q[3];
rz(-1.8146993) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2372357) q[0];
sx q[0];
rz(-0.80191737) q[0];
sx q[0];
rz(0.075415762) q[0];
rz(-0.30741179) q[1];
sx q[1];
rz(-1.9785898) q[1];
sx q[1];
rz(-0.42617646) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3097174) q[0];
sx q[0];
rz(-2.2906036) q[0];
sx q[0];
rz(-0.72711522) q[0];
x q[1];
rz(2.6226042) q[2];
sx q[2];
rz(-2.6140169) q[2];
sx q[2];
rz(-1.0949539) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.4170917) q[1];
sx q[1];
rz(-1.0001273) q[1];
sx q[1];
rz(-2.9076734) q[1];
rz(-pi) q[2];
rz(-0.05699031) q[3];
sx q[3];
rz(-0.40613031) q[3];
sx q[3];
rz(0.95088357) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.9093466) q[2];
sx q[2];
rz(-1.1541977) q[2];
sx q[2];
rz(-0.92440355) q[2];
rz(-2.038548) q[3];
sx q[3];
rz(-1.1140991) q[3];
sx q[3];
rz(-1.2801142) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0019317146) q[0];
sx q[0];
rz(-2.6502471) q[0];
sx q[0];
rz(-0.77144462) q[0];
rz(1.3123243) q[1];
sx q[1];
rz(-1.7701365) q[1];
sx q[1];
rz(1.9171453) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79662844) q[0];
sx q[0];
rz(-1.7078236) q[0];
sx q[0];
rz(-3.1132878) q[0];
rz(-1.5727051) q[2];
sx q[2];
rz(-0.95237007) q[2];
sx q[2];
rz(1.2313953) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.33358303) q[1];
sx q[1];
rz(-1.6621272) q[1];
sx q[1];
rz(-1.1846945) q[1];
x q[2];
rz(2.0290852) q[3];
sx q[3];
rz(-2.3009514) q[3];
sx q[3];
rz(-0.86220137) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9839342) q[2];
sx q[2];
rz(-2.5658786) q[2];
sx q[2];
rz(-2.0470587) q[2];
rz(-3.0648699) q[3];
sx q[3];
rz(-0.50944296) q[3];
sx q[3];
rz(0.69242394) q[3];
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
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7102966) q[0];
sx q[0];
rz(-0.65839473) q[0];
sx q[0];
rz(-2.9265213) q[0];
rz(-2.1521425) q[1];
sx q[1];
rz(-2.1727401) q[1];
sx q[1];
rz(-1.6542124) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7152255) q[0];
sx q[0];
rz(-1.0578814) q[0];
sx q[0];
rz(1.7576393) q[0];
rz(-0.39535661) q[2];
sx q[2];
rz(-1.8430508) q[2];
sx q[2];
rz(-0.56955063) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.3645059) q[1];
sx q[1];
rz(-1.1544656) q[1];
sx q[1];
rz(-0.39752605) q[1];
rz(-pi) q[2];
rz(-1.0808252) q[3];
sx q[3];
rz(-2.1741011) q[3];
sx q[3];
rz(-1.5321466) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.4789077) q[2];
sx q[2];
rz(-1.4143133) q[2];
sx q[2];
rz(1.5642081) q[2];
rz(-2.1936737) q[3];
sx q[3];
rz(-1.0736059) q[3];
sx q[3];
rz(1.7387559) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74744144) q[0];
sx q[0];
rz(-1.489137) q[0];
sx q[0];
rz(-0.58962756) q[0];
rz(-1.6472752) q[1];
sx q[1];
rz(-1.3811771) q[1];
sx q[1];
rz(3.0628915) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37109659) q[0];
sx q[0];
rz(-1.1770304) q[0];
sx q[0];
rz(1.257819) q[0];
x q[1];
rz(-2.6888618) q[2];
sx q[2];
rz(-0.50846487) q[2];
sx q[2];
rz(-1.1217211) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.4884686) q[1];
sx q[1];
rz(-0.32652509) q[1];
sx q[1];
rz(2.4420183) q[1];
rz(-pi) q[2];
rz(-1.0534473) q[3];
sx q[3];
rz(-2.0387408) q[3];
sx q[3];
rz(3.057923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.9239203) q[2];
sx q[2];
rz(-2.2675026) q[2];
sx q[2];
rz(0.11202845) q[2];
rz(-0.40515408) q[3];
sx q[3];
rz(-1.3978037) q[3];
sx q[3];
rz(0.20598327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9236295) q[0];
sx q[0];
rz(-0.97342747) q[0];
sx q[0];
rz(1.7226891) q[0];
rz(-2.0619242) q[1];
sx q[1];
rz(-0.38366145) q[1];
sx q[1];
rz(-1.4071646) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2699958) q[0];
sx q[0];
rz(-1.9752335) q[0];
sx q[0];
rz(1.2258421) q[0];
rz(1.4905995) q[2];
sx q[2];
rz(-2.0010304) q[2];
sx q[2];
rz(2.4938817) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.0041173) q[1];
sx q[1];
rz(-1.1761929) q[1];
sx q[1];
rz(2.7156107) q[1];
rz(-pi) q[2];
rz(-1.2958762) q[3];
sx q[3];
rz(-1.3310199) q[3];
sx q[3];
rz(1.4942684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5088523) q[2];
sx q[2];
rz(-0.036157046) q[2];
sx q[2];
rz(-0.73842326) q[2];
rz(0.19668713) q[3];
sx q[3];
rz(-2.5236712) q[3];
sx q[3];
rz(0.64016199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7554373) q[0];
sx q[0];
rz(-2.8897987) q[0];
sx q[0];
rz(3.1016896) q[0];
rz(-2.9207322) q[1];
sx q[1];
rz(-1.5233636) q[1];
sx q[1];
rz(1.9853282) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6823718) q[0];
sx q[0];
rz(-0.30293884) q[0];
sx q[0];
rz(-1.0653063) q[0];
rz(1.7650911) q[2];
sx q[2];
rz(-1.7826323) q[2];
sx q[2];
rz(1.719033) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.5584692) q[1];
sx q[1];
rz(-1.236318) q[1];
sx q[1];
rz(2.1534647) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4651894) q[3];
sx q[3];
rz(-1.6855414) q[3];
sx q[3];
rz(0.70699837) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.3952737) q[2];
sx q[2];
rz(-2.7776182) q[2];
sx q[2];
rz(0.054595646) q[2];
rz(-2.3553081) q[3];
sx q[3];
rz(-1.1966642) q[3];
sx q[3];
rz(1.1597077) q[3];
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
rz(0.38128024) q[0];
sx q[0];
rz(-0.78802839) q[0];
sx q[0];
rz(2.330761) q[0];
rz(-0.52931085) q[1];
sx q[1];
rz(-1.7050754) q[1];
sx q[1];
rz(-2.0307821) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62880077) q[0];
sx q[0];
rz(-1.1847825) q[0];
sx q[0];
rz(2.5284472) q[0];
x q[1];
rz(-1.0093685) q[2];
sx q[2];
rz(-1.9931442) q[2];
sx q[2];
rz(2.7627711) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.1326434) q[1];
sx q[1];
rz(-0.81671444) q[1];
sx q[1];
rz(1.7945748) q[1];
x q[2];
rz(-0.30382321) q[3];
sx q[3];
rz(-2.15832) q[3];
sx q[3];
rz(1.2161906) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.79798737) q[2];
sx q[2];
rz(-2.5729749) q[2];
sx q[2];
rz(-2.9000751) q[2];
rz(2.8999515) q[3];
sx q[3];
rz(-1.6499949) q[3];
sx q[3];
rz(2.9901166) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3342313) q[0];
sx q[0];
rz(-1.3483028) q[0];
sx q[0];
rz(-2.469113) q[0];
rz(-2.6895651) q[1];
sx q[1];
rz(-2.1990135) q[1];
sx q[1];
rz(2.6262851) q[1];
rz(3.0816513) q[2];
sx q[2];
rz(-1.1287545) q[2];
sx q[2];
rz(-1.9311369) q[2];
rz(-0.84108055) q[3];
sx q[3];
rz(-0.65541808) q[3];
sx q[3];
rz(-2.2578166) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
