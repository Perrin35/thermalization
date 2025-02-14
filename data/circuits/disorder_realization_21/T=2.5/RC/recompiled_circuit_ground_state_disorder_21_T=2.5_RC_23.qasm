OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.8706239) q[0];
sx q[0];
rz(-2.5854817) q[0];
sx q[0];
rz(-2.1882353) q[0];
rz(0.019729992) q[1];
sx q[1];
rz(-2.1058197) q[1];
sx q[1];
rz(0.9486202) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0355546) q[0];
sx q[0];
rz(-0.69879888) q[0];
sx q[0];
rz(-2.4731556) q[0];
x q[1];
rz(0.23748246) q[2];
sx q[2];
rz(-1.5735224) q[2];
sx q[2];
rz(3.1196371) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1745156) q[1];
sx q[1];
rz(-2.3070455) q[1];
sx q[1];
rz(-0.59474919) q[1];
rz(2.0998276) q[3];
sx q[3];
rz(-2.5210125) q[3];
sx q[3];
rz(2.4247501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.76401508) q[2];
sx q[2];
rz(-3.0674051) q[2];
sx q[2];
rz(0.78262502) q[2];
rz(0.11893663) q[3];
sx q[3];
rz(-2.1075893) q[3];
sx q[3];
rz(-2.9940166) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19621944) q[0];
sx q[0];
rz(-2.0202899) q[0];
sx q[0];
rz(0.25783208) q[0];
rz(-0.091015426) q[1];
sx q[1];
rz(-1.0992522) q[1];
sx q[1];
rz(1.4965422) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.634406) q[0];
sx q[0];
rz(-1.0189462) q[0];
sx q[0];
rz(-0.080291434) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7690954) q[2];
sx q[2];
rz(-1.2476693) q[2];
sx q[2];
rz(1.1748287) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.0219258) q[1];
sx q[1];
rz(-2.7283784) q[1];
sx q[1];
rz(1.6456804) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4739059) q[3];
sx q[3];
rz(-1.2419309) q[3];
sx q[3];
rz(-1.937695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.1220793) q[2];
sx q[2];
rz(-1.8668819) q[2];
sx q[2];
rz(-0.40461928) q[2];
rz(-0.98958611) q[3];
sx q[3];
rz(-1.3000969) q[3];
sx q[3];
rz(-0.62304455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0290381) q[0];
sx q[0];
rz(-2.9974388) q[0];
sx q[0];
rz(-2.3032904) q[0];
rz(2.2932032) q[1];
sx q[1];
rz(-1.219607) q[1];
sx q[1];
rz(-2.1489876) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58439764) q[0];
sx q[0];
rz(-2.2879172) q[0];
sx q[0];
rz(-1.3846057) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6279334) q[2];
sx q[2];
rz(-1.2369725) q[2];
sx q[2];
rz(-1.8360857) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.4731208) q[1];
sx q[1];
rz(-1.8278012) q[1];
sx q[1];
rz(0.42196749) q[1];
x q[2];
rz(1.6818524) q[3];
sx q[3];
rz(-1.2821615) q[3];
sx q[3];
rz(1.9185324) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2720211) q[2];
sx q[2];
rz(-1.973899) q[2];
sx q[2];
rz(-2.0167548) q[2];
rz(-2.520842) q[3];
sx q[3];
rz(-2.1706332) q[3];
sx q[3];
rz(0.11515215) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(0.89897412) q[0];
sx q[0];
rz(-1.1599351) q[0];
sx q[0];
rz(2.9677891) q[0];
rz(0.68570343) q[1];
sx q[1];
rz(-1.655429) q[1];
sx q[1];
rz(2.3033843) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7156775) q[0];
sx q[0];
rz(-1.1155778) q[0];
sx q[0];
rz(-1.0550189) q[0];
rz(-2.4160552) q[2];
sx q[2];
rz(-0.71070403) q[2];
sx q[2];
rz(0.85884554) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.3633109) q[1];
sx q[1];
rz(-1.7125336) q[1];
sx q[1];
rz(-2.4635386) q[1];
rz(-1.6082786) q[3];
sx q[3];
rz(-2.0333383) q[3];
sx q[3];
rz(1.9328062) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.59914261) q[2];
sx q[2];
rz(-1.6833545) q[2];
sx q[2];
rz(-0.35169265) q[2];
rz(1.1951949) q[3];
sx q[3];
rz(-1.9545133) q[3];
sx q[3];
rz(-3.1174507) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46566063) q[0];
sx q[0];
rz(-2.6213578) q[0];
sx q[0];
rz(-2.7886673) q[0];
rz(-2.832761) q[1];
sx q[1];
rz(-2.1052723) q[1];
sx q[1];
rz(-3.1130863) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7362979) q[0];
sx q[0];
rz(-0.60966821) q[0];
sx q[0];
rz(1.8788615) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1813004) q[2];
sx q[2];
rz(-1.2151698) q[2];
sx q[2];
rz(-0.86330044) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5529768) q[1];
sx q[1];
rz(-1.3541095) q[1];
sx q[1];
rz(0.79766794) q[1];
rz(-1.1358374) q[3];
sx q[3];
rz(-1.3490145) q[3];
sx q[3];
rz(-1.3835761) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.21294022) q[2];
sx q[2];
rz(-0.063455909) q[2];
sx q[2];
rz(-2.1707936) q[2];
rz(1.0468696) q[3];
sx q[3];
rz(-1.6887083) q[3];
sx q[3];
rz(2.651732) q[3];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8397119) q[0];
sx q[0];
rz(-1.7827001) q[0];
sx q[0];
rz(0.090959892) q[0];
rz(-1.2184527) q[1];
sx q[1];
rz(-2.8107042) q[1];
sx q[1];
rz(2.5023696) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42933336) q[0];
sx q[0];
rz(-0.3540701) q[0];
sx q[0];
rz(2.7430915) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2940862) q[2];
sx q[2];
rz(-0.18346645) q[2];
sx q[2];
rz(-3.1364721) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.6248577) q[1];
sx q[1];
rz(-1.0565041) q[1];
sx q[1];
rz(-2.3964797) q[1];
rz(-1.5766673) q[3];
sx q[3];
rz(-2.8150055) q[3];
sx q[3];
rz(2.011855) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-3.0794534) q[2];
sx q[2];
rz(-2.1663351) q[2];
sx q[2];
rz(2.6677168) q[2];
rz(1.8114932) q[3];
sx q[3];
rz(-2.2606943) q[3];
sx q[3];
rz(-0.37548319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57629958) q[0];
sx q[0];
rz(-2.0602891) q[0];
sx q[0];
rz(2.9190049) q[0];
rz(2.0780156) q[1];
sx q[1];
rz(-1.5779326) q[1];
sx q[1];
rz(1.1933914) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9759068) q[0];
sx q[0];
rz(-0.75729174) q[0];
sx q[0];
rz(-0.16974561) q[0];
rz(-pi) q[1];
rz(-3.0203825) q[2];
sx q[2];
rz(-1.3451936) q[2];
sx q[2];
rz(0.76162213) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.8377922) q[1];
sx q[1];
rz(-1.1322316) q[1];
sx q[1];
rz(1.4213417) q[1];
rz(-pi) q[2];
rz(2.1126595) q[3];
sx q[3];
rz(-0.30908424) q[3];
sx q[3];
rz(2.2268471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-3.0003164) q[2];
sx q[2];
rz(-2.2681984) q[2];
sx q[2];
rz(0.6401965) q[2];
rz(0.021942465) q[3];
sx q[3];
rz(-1.7317737) q[3];
sx q[3];
rz(2.791361) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.491965) q[0];
sx q[0];
rz(-1.7438629) q[0];
sx q[0];
rz(-2.6212027) q[0];
rz(-2.261816) q[1];
sx q[1];
rz(-1.2326515) q[1];
sx q[1];
rz(-0.89281503) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.797282) q[0];
sx q[0];
rz(-1.6160674) q[0];
sx q[0];
rz(-0.08749732) q[0];
rz(-pi) q[1];
rz(-2.3310066) q[2];
sx q[2];
rz(-0.91601585) q[2];
sx q[2];
rz(2.9952632) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0860111) q[1];
sx q[1];
rz(-2.308929) q[1];
sx q[1];
rz(-2.7559381) q[1];
rz(-0.51629169) q[3];
sx q[3];
rz(-1.9441351) q[3];
sx q[3];
rz(3.1331799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.2616547) q[2];
sx q[2];
rz(-2.0377906) q[2];
sx q[2];
rz(3.0547764) q[2];
rz(-0.9451198) q[3];
sx q[3];
rz(-0.34543959) q[3];
sx q[3];
rz(-1.1300348) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9574808) q[0];
sx q[0];
rz(-3.0511973) q[0];
sx q[0];
rz(-0.064067319) q[0];
rz(-2.0059026) q[1];
sx q[1];
rz(-1.4402729) q[1];
sx q[1];
rz(2.6293829) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6653888) q[0];
sx q[0];
rz(-1.246016) q[0];
sx q[0];
rz(0.43412125) q[0];
x q[1];
rz(-1.9718599) q[2];
sx q[2];
rz(-0.74359119) q[2];
sx q[2];
rz(-1.0917853) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.4254949) q[1];
sx q[1];
rz(-1.6563882) q[1];
sx q[1];
rz(-0.25984989) q[1];
rz(-pi) q[2];
rz(3.104761) q[3];
sx q[3];
rz(-0.40912155) q[3];
sx q[3];
rz(-0.44884071) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.8168489) q[2];
sx q[2];
rz(-0.73255676) q[2];
sx q[2];
rz(1.0653488) q[2];
rz(1.4893074) q[3];
sx q[3];
rz(-0.95732147) q[3];
sx q[3];
rz(-0.092441946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5866933) q[0];
sx q[0];
rz(-2.7296992) q[0];
sx q[0];
rz(0.33083415) q[0];
rz(1.9031485) q[1];
sx q[1];
rz(-0.60791433) q[1];
sx q[1];
rz(2.8668561) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58178066) q[0];
sx q[0];
rz(-1.1413478) q[0];
sx q[0];
rz(2.7665124) q[0];
x q[1];
rz(0.27336911) q[2];
sx q[2];
rz(-1.6022575) q[2];
sx q[2];
rz(-0.69845573) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.292784) q[1];
sx q[1];
rz(-2.074179) q[1];
sx q[1];
rz(-0.15979691) q[1];
rz(1.4624743) q[3];
sx q[3];
rz(-2.5511955) q[3];
sx q[3];
rz(0.27449671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.1642509) q[2];
sx q[2];
rz(-2.3873316) q[2];
sx q[2];
rz(-1.784262) q[2];
rz(-2.6409798) q[3];
sx q[3];
rz(-1.5784135) q[3];
sx q[3];
rz(-1.4969426) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5464583) q[0];
sx q[0];
rz(-1.5820137) q[0];
sx q[0];
rz(-0.59376846) q[0];
rz(-0.068269923) q[1];
sx q[1];
rz(-1.2208114) q[1];
sx q[1];
rz(0.023718871) q[1];
rz(2.7082607) q[2];
sx q[2];
rz(-2.7637252) q[2];
sx q[2];
rz(-2.6890474) q[2];
rz(-0.81238644) q[3];
sx q[3];
rz(-2.8858622) q[3];
sx q[3];
rz(1.1271911) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
