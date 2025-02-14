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
rz(1.709197) q[0];
sx q[0];
rz(3.4442918) q[0];
sx q[0];
rz(11.533307) q[0];
rz(1.9167702) q[1];
sx q[1];
rz(-1.6515825) q[1];
sx q[1];
rz(-2.9472247) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2107735) q[0];
sx q[0];
rz(-1.418485) q[0];
sx q[0];
rz(2.5126626) q[0];
rz(-0.35440234) q[2];
sx q[2];
rz(-2.4131219) q[2];
sx q[2];
rz(-3.1284863) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.3857419) q[1];
sx q[1];
rz(-0.057690851) q[1];
sx q[1];
rz(-2.2123446) q[1];
rz(-pi) q[2];
rz(-2.5159149) q[3];
sx q[3];
rz(-1.9301658) q[3];
sx q[3];
rz(0.43454447) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.73244652) q[2];
sx q[2];
rz(-1.1017841) q[2];
sx q[2];
rz(2.5771602) q[2];
rz(1.268528) q[3];
sx q[3];
rz(-0.0081491834) q[3];
sx q[3];
rz(-0.90659365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.059435189) q[0];
sx q[0];
rz(-3.1340288) q[0];
sx q[0];
rz(-0.91307688) q[0];
rz(0.68436855) q[1];
sx q[1];
rz(-0.00033907779) q[1];
sx q[1];
rz(2.2025462) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2850239) q[0];
sx q[0];
rz(-1.9459581) q[0];
sx q[0];
rz(0.60174353) q[0];
x q[1];
rz(-0.026264391) q[2];
sx q[2];
rz(-2.677478) q[2];
sx q[2];
rz(0.26192203) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.9833656) q[1];
sx q[1];
rz(-1.905113) q[1];
sx q[1];
rz(-1.1329705) q[1];
rz(-pi) q[2];
rz(-1.9775852) q[3];
sx q[3];
rz(-1.6322246) q[3];
sx q[3];
rz(-0.2964501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3188476) q[2];
sx q[2];
rz(-3.1343967) q[2];
sx q[2];
rz(-0.89246559) q[2];
rz(-1.5626296) q[3];
sx q[3];
rz(-0.024450863) q[3];
sx q[3];
rz(1.310937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
rz(2.7398294) q[0];
sx q[0];
rz(-0.50245291) q[0];
sx q[0];
rz(0.40973642) q[0];
rz(-0.01092625) q[1];
sx q[1];
rz(-2.9210563) q[1];
sx q[1];
rz(1.3440557) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15191244) q[0];
sx q[0];
rz(-1.5912959) q[0];
sx q[0];
rz(2.0818039) q[0];
x q[1];
rz(3.0777099) q[2];
sx q[2];
rz(-1.4555571) q[2];
sx q[2];
rz(-1.921953) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.2933767) q[1];
sx q[1];
rz(-0.266222) q[1];
sx q[1];
rz(1.8195527) q[1];
rz(1.6860754) q[3];
sx q[3];
rz(-1.56101) q[3];
sx q[3];
rz(0.20231314) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.4695796) q[2];
sx q[2];
rz(-1.4521705) q[2];
sx q[2];
rz(-0.043188485) q[2];
rz(1.1399266) q[3];
sx q[3];
rz(-2.98525) q[3];
sx q[3];
rz(1.9388916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7286872) q[0];
sx q[0];
rz(-2.7396956) q[0];
sx q[0];
rz(-1.3380949) q[0];
rz(0.69522229) q[1];
sx q[1];
rz(-0.1278563) q[1];
sx q[1];
rz(-2.7241838) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.886717) q[0];
sx q[0];
rz(-0.69807295) q[0];
sx q[0];
rz(1.3017734) q[0];
rz(-pi) q[1];
rz(1.9193255) q[2];
sx q[2];
rz(-0.67643316) q[2];
sx q[2];
rz(1.2471022) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.96046916) q[1];
sx q[1];
rz(-2.4669018) q[1];
sx q[1];
rz(-1.2441649) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.14162161) q[3];
sx q[3];
rz(-1.9998261) q[3];
sx q[3];
rz(-0.98875586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.43856105) q[2];
sx q[2];
rz(-0.87623864) q[2];
sx q[2];
rz(-1.7523127) q[2];
rz(2.321068) q[3];
sx q[3];
rz(-1.5503784) q[3];
sx q[3];
rz(2.4337721) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97838068) q[0];
sx q[0];
rz(-3.1040525) q[0];
sx q[0];
rz(2.1782844) q[0];
rz(-2.9523201) q[1];
sx q[1];
rz(-3.1258588) q[1];
sx q[1];
rz(0.16545573) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4625797) q[0];
sx q[0];
rz(-1.6143342) q[0];
sx q[0];
rz(1.5971558) q[0];
x q[1];
rz(2.8420429) q[2];
sx q[2];
rz(-2.2431157) q[2];
sx q[2];
rz(-1.7393665) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.8973863) q[1];
sx q[1];
rz(-1.4546977) q[1];
sx q[1];
rz(-0.80516025) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1375764) q[3];
sx q[3];
rz(-1.2804083) q[3];
sx q[3];
rz(0.079513915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.307622) q[2];
sx q[2];
rz(-0.51887363) q[2];
sx q[2];
rz(-2.0818254) q[2];
rz(-2.4506532) q[3];
sx q[3];
rz(-0.86425224) q[3];
sx q[3];
rz(-1.8701657) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5237913) q[0];
sx q[0];
rz(-0.093136223) q[0];
sx q[0];
rz(1.6049438) q[0];
rz(0.45283428) q[1];
sx q[1];
rz(-0.0080527877) q[1];
sx q[1];
rz(1.4401999) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8062167) q[0];
sx q[0];
rz(-1.4037689) q[0];
sx q[0];
rz(1.7352292) q[0];
rz(0.78078713) q[2];
sx q[2];
rz(-0.30715725) q[2];
sx q[2];
rz(0.13328341) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.837484) q[1];
sx q[1];
rz(-2.9867801) q[1];
sx q[1];
rz(-3.0938838) q[1];
rz(-pi) q[2];
rz(-1.7415288) q[3];
sx q[3];
rz(-1.472982) q[3];
sx q[3];
rz(-0.82966891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3692533) q[2];
sx q[2];
rz(-0.99268308) q[2];
sx q[2];
rz(-3.0625647) q[2];
rz(0.8655656) q[3];
sx q[3];
rz(-0.95439684) q[3];
sx q[3];
rz(-1.0237833) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9158151) q[0];
sx q[0];
rz(-3.1362035) q[0];
sx q[0];
rz(-0.22501568) q[0];
rz(-2.8354538) q[1];
sx q[1];
rz(-0.016751079) q[1];
sx q[1];
rz(-0.87047815) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7891312) q[0];
sx q[0];
rz(-1.5232183) q[0];
sx q[0];
rz(-3.070288) q[0];
rz(-pi) q[1];
rz(1.00433) q[2];
sx q[2];
rz(-2.2029467) q[2];
sx q[2];
rz(1.8514148) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.5600109) q[1];
sx q[1];
rz(-2.3696179) q[1];
sx q[1];
rz(-3.0549269) q[1];
x q[2];
rz(2.0567749) q[3];
sx q[3];
rz(-1.2013913) q[3];
sx q[3];
rz(1.891275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.6656782) q[2];
sx q[2];
rz(-3.0147538) q[2];
sx q[2];
rz(2.1062984) q[2];
rz(2.3546442) q[3];
sx q[3];
rz(-3.0988155) q[3];
sx q[3];
rz(0.27534819) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.03054522) q[0];
sx q[0];
rz(-0.011140911) q[0];
sx q[0];
rz(-0.025644843) q[0];
rz(1.2643087) q[1];
sx q[1];
rz(-3.1176716) q[1];
sx q[1];
rz(2.4425676) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0920588) q[0];
sx q[0];
rz(-1.3483241) q[0];
sx q[0];
rz(-0.27357863) q[0];
x q[1];
rz(0.46866663) q[2];
sx q[2];
rz(-1.4184409) q[2];
sx q[2];
rz(0.050616905) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(3.0009202) q[1];
sx q[1];
rz(-0.39127054) q[1];
sx q[1];
rz(-0.19473445) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2243629) q[3];
sx q[3];
rz(-1.1207657) q[3];
sx q[3];
rz(-0.91564028) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.25253025) q[2];
sx q[2];
rz(-2.3928596) q[2];
sx q[2];
rz(2.8177281) q[2];
rz(2.9845386) q[3];
sx q[3];
rz(-1.2374977) q[3];
sx q[3];
rz(-0.15373716) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57482982) q[0];
sx q[0];
rz(-3.1269508) q[0];
sx q[0];
rz(-2.5515442) q[0];
rz(-2.4216962) q[1];
sx q[1];
rz(-3.0820334) q[1];
sx q[1];
rz(-0.86729008) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.86751) q[0];
sx q[0];
rz(-1.6757586) q[0];
sx q[0];
rz(-0.63469736) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2008823) q[2];
sx q[2];
rz(-2.5061766) q[2];
sx q[2];
rz(0.70178343) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.1220951) q[1];
sx q[1];
rz(-1.4585134) q[1];
sx q[1];
rz(0.068536802) q[1];
rz(-pi) q[2];
rz(1.7884364) q[3];
sx q[3];
rz(-10*pi/13) q[3];
sx q[3];
rz(-2.9293037) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.1934293) q[2];
sx q[2];
rz(-0.89274222) q[2];
sx q[2];
rz(-0.11290045) q[2];
rz(1.0490949) q[3];
sx q[3];
rz(-1.5594522) q[3];
sx q[3];
rz(0.73672867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0650487) q[0];
sx q[0];
rz(-1.5600518) q[0];
sx q[0];
rz(1.6381868) q[0];
rz(0.63649559) q[1];
sx q[1];
rz(-0.46000767) q[1];
sx q[1];
rz(1.5696625) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.628807) q[0];
sx q[0];
rz(-0.11462002) q[0];
sx q[0];
rz(0.8091457) q[0];
rz(-1.5701446) q[2];
sx q[2];
rz(-1.5739292) q[2];
sx q[2];
rz(-0.76515406) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.0850141) q[1];
sx q[1];
rz(-0.002925891) q[1];
sx q[1];
rz(3.0614733) q[1];
x q[2];
rz(-0.51497634) q[3];
sx q[3];
rz(-1.4287717) q[3];
sx q[3];
rz(1.9076626) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.2703209) q[2];
sx q[2];
rz(-1.511829) q[2];
sx q[2];
rz(2.9809269) q[2];
rz(1.2284944) q[3];
sx q[3];
rz(-3.1077423) q[3];
sx q[3];
rz(2.9401949) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0495618) q[0];
sx q[0];
rz(-1.1783538) q[0];
sx q[0];
rz(0.023068064) q[0];
rz(1.6231712) q[1];
sx q[1];
rz(-0.36133125) q[1];
sx q[1];
rz(0.28892118) q[1];
rz(1.4949746) q[2];
sx q[2];
rz(-0.29360148) q[2];
sx q[2];
rz(0.32099024) q[2];
rz(1.2613983) q[3];
sx q[3];
rz(-2.6831476) q[3];
sx q[3];
rz(0.89098709) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
