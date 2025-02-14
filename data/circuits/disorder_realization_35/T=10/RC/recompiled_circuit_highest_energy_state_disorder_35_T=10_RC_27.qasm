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
rz(-2.45911) q[0];
sx q[0];
rz(-0.64829666) q[0];
sx q[0];
rz(0.034870235) q[0];
rz(1.3999445) q[1];
sx q[1];
rz(6.016461) q[1];
sx q[1];
rz(7.4363554) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54915743) q[0];
sx q[0];
rz(-1.9218635) q[0];
sx q[0];
rz(0.017819509) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.2049293) q[2];
sx q[2];
rz(-1.7589169) q[2];
sx q[2];
rz(-2.8486696) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.68426484) q[1];
sx q[1];
rz(-0.65456355) q[1];
sx q[1];
rz(1.2689204) q[1];
rz(-pi) q[2];
rz(1.5409914) q[3];
sx q[3];
rz(-0.4510703) q[3];
sx q[3];
rz(0.20584895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(3.1109041) q[2];
sx q[2];
rz(-1.8762174) q[2];
sx q[2];
rz(3.0008924) q[2];
rz(1.582102) q[3];
sx q[3];
rz(-2.9086106) q[3];
sx q[3];
rz(1.337602) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2577308) q[0];
sx q[0];
rz(-2.1900539) q[0];
sx q[0];
rz(2.0146433) q[0];
rz(1.1832712) q[1];
sx q[1];
rz(-2.0014346) q[1];
sx q[1];
rz(0.3068876) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81079262) q[0];
sx q[0];
rz(-1.9527657) q[0];
sx q[0];
rz(-1.6195787) q[0];
rz(-pi) q[1];
rz(0.65019239) q[2];
sx q[2];
rz(-0.99417328) q[2];
sx q[2];
rz(1.2675831) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.3908071) q[1];
sx q[1];
rz(-1.7089426) q[1];
sx q[1];
rz(2.9259794) q[1];
rz(-pi) q[2];
rz(1.0829665) q[3];
sx q[3];
rz(-1.4734419) q[3];
sx q[3];
rz(0.1258345) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7221476) q[2];
sx q[2];
rz(-2.751613) q[2];
sx q[2];
rz(1.6411714) q[2];
rz(-1.8309343) q[3];
sx q[3];
rz(-0.981841) q[3];
sx q[3];
rz(-0.76154861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0797043) q[0];
sx q[0];
rz(-0.56489983) q[0];
sx q[0];
rz(1.0796219) q[0];
rz(2.9234746) q[1];
sx q[1];
rz(-1.4061385) q[1];
sx q[1];
rz(-1.9535779) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69258286) q[0];
sx q[0];
rz(-0.94845686) q[0];
sx q[0];
rz(1.5507109) q[0];
rz(-pi) q[1];
x q[1];
rz(0.25323694) q[2];
sx q[2];
rz(-0.5917158) q[2];
sx q[2];
rz(-2.5160421) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0741785) q[1];
sx q[1];
rz(-1.6700831) q[1];
sx q[1];
rz(1.5669785) q[1];
rz(-0.025026021) q[3];
sx q[3];
rz(-2.9200168) q[3];
sx q[3];
rz(3.0794249) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.5334566) q[2];
sx q[2];
rz(-0.44508219) q[2];
sx q[2];
rz(-2.377887) q[2];
rz(-2.8869827) q[3];
sx q[3];
rz(-2.2201241) q[3];
sx q[3];
rz(-2.7931131) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3449645) q[0];
sx q[0];
rz(-2.1944955) q[0];
sx q[0];
rz(0.011888591) q[0];
rz(-0.97870007) q[1];
sx q[1];
rz(-1.3589166) q[1];
sx q[1];
rz(-1.1241283) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5886712) q[0];
sx q[0];
rz(-1.3359937) q[0];
sx q[0];
rz(0.42490439) q[0];
rz(-pi) q[1];
x q[1];
rz(0.61907436) q[2];
sx q[2];
rz(-2.4644669) q[2];
sx q[2];
rz(-0.4437333) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.2024514) q[1];
sx q[1];
rz(-1.141046) q[1];
sx q[1];
rz(-0.38265574) q[1];
rz(0.72953485) q[3];
sx q[3];
rz(-1.9771061) q[3];
sx q[3];
rz(2.0168347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-3.0594788) q[2];
sx q[2];
rz(-0.79136005) q[2];
sx q[2];
rz(-0.7938844) q[2];
rz(1.7700899) q[3];
sx q[3];
rz(-0.81727782) q[3];
sx q[3];
rz(-0.98384682) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3802948) q[0];
sx q[0];
rz(-1.607432) q[0];
sx q[0];
rz(0.5436596) q[0];
rz(-1.9937953) q[1];
sx q[1];
rz(-2.5416083) q[1];
sx q[1];
rz(-2.0701087) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9990014) q[0];
sx q[0];
rz(-1.1650225) q[0];
sx q[0];
rz(-0.31082992) q[0];
rz(-pi) q[1];
rz(2.8473994) q[2];
sx q[2];
rz(-1.718866) q[2];
sx q[2];
rz(1.8235056) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.0746911) q[1];
sx q[1];
rz(-1.2558535) q[1];
sx q[1];
rz(-2.6661162) q[1];
x q[2];
rz(-2.33458) q[3];
sx q[3];
rz(-1.1928946) q[3];
sx q[3];
rz(-0.81951362) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.078309623) q[2];
sx q[2];
rz(-0.50476837) q[2];
sx q[2];
rz(2.2341991) q[2];
rz(0.0026957938) q[3];
sx q[3];
rz(-1.4009652) q[3];
sx q[3];
rz(1.9122972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21830782) q[0];
sx q[0];
rz(-1.8508428) q[0];
sx q[0];
rz(1.1894591) q[0];
rz(2.9749191) q[1];
sx q[1];
rz(-2.1107626) q[1];
sx q[1];
rz(-1.6698242) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20118344) q[0];
sx q[0];
rz(-1.1795949) q[0];
sx q[0];
rz(0.67232491) q[0];
x q[1];
rz(1.5277063) q[2];
sx q[2];
rz(-0.22858563) q[2];
sx q[2];
rz(-3.0800386) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.7547524) q[1];
sx q[1];
rz(-1.0472782) q[1];
sx q[1];
rz(1.0922111) q[1];
rz(-pi) q[2];
rz(-0.81780602) q[3];
sx q[3];
rz(-1.3997625) q[3];
sx q[3];
rz(0.59137646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.2344096) q[2];
sx q[2];
rz(-0.11203354) q[2];
sx q[2];
rz(0.33667931) q[2];
rz(-1.9628717) q[3];
sx q[3];
rz(-1.4320635) q[3];
sx q[3];
rz(-2.7746576) q[3];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1692899) q[0];
sx q[0];
rz(-0.65880913) q[0];
sx q[0];
rz(1.7377874) q[0];
rz(0.16109578) q[1];
sx q[1];
rz(-0.95268551) q[1];
sx q[1];
rz(2.5162627) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17574425) q[0];
sx q[0];
rz(-2.1702156) q[0];
sx q[0];
rz(0.61300399) q[0];
x q[1];
rz(-2.1834945) q[2];
sx q[2];
rz(-0.73054689) q[2];
sx q[2];
rz(2.036866) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.8323601) q[1];
sx q[1];
rz(-2.6200868) q[1];
sx q[1];
rz(1.7108365) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1745542) q[3];
sx q[3];
rz(-1.4449287) q[3];
sx q[3];
rz(0.87115067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.27955851) q[2];
sx q[2];
rz(-1.6509849) q[2];
sx q[2];
rz(0.030869182) q[2];
rz(-0.75012642) q[3];
sx q[3];
rz(-0.98837367) q[3];
sx q[3];
rz(-1.3039024) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52791643) q[0];
sx q[0];
rz(-0.91687098) q[0];
sx q[0];
rz(3.0314639) q[0];
rz(0.88422173) q[1];
sx q[1];
rz(-1.5343816) q[1];
sx q[1];
rz(-1.0482739) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3582402) q[0];
sx q[0];
rz(-1.8723346) q[0];
sx q[0];
rz(1.8389971) q[0];
x q[1];
rz(-1.3051304) q[2];
sx q[2];
rz(-1.507826) q[2];
sx q[2];
rz(2.8612325) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.5267262) q[1];
sx q[1];
rz(-1.5121361) q[1];
sx q[1];
rz(-0.23414302) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0804668) q[3];
sx q[3];
rz(-1.8009652) q[3];
sx q[3];
rz(1.597061) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.5683762) q[2];
sx q[2];
rz(-1.5244851) q[2];
sx q[2];
rz(3.1061843) q[2];
rz(3.0951485) q[3];
sx q[3];
rz(-1.8567663) q[3];
sx q[3];
rz(0.057028381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34920084) q[0];
sx q[0];
rz(-2.2430895) q[0];
sx q[0];
rz(3.0080556) q[0];
rz(1.4117433) q[1];
sx q[1];
rz(-1.1213419) q[1];
sx q[1];
rz(2.0917361) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8658757) q[0];
sx q[0];
rz(-1.9483074) q[0];
sx q[0];
rz(-1.4488245) q[0];
x q[1];
rz(2.1112415) q[2];
sx q[2];
rz(-1.5771416) q[2];
sx q[2];
rz(-0.66101569) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.404001) q[1];
sx q[1];
rz(-1.6114514) q[1];
sx q[1];
rz(-2.3126751) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7383835) q[3];
sx q[3];
rz(-1.3018738) q[3];
sx q[3];
rz(1.2300242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.92020804) q[2];
sx q[2];
rz(-2.460026) q[2];
sx q[2];
rz(-1.3725012) q[2];
rz(-1.8364604) q[3];
sx q[3];
rz(-1.2499481) q[3];
sx q[3];
rz(-0.60177747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9735499) q[0];
sx q[0];
rz(-1.9815227) q[0];
sx q[0];
rz(2.2450182) q[0];
rz(-0.98094455) q[1];
sx q[1];
rz(-0.70471057) q[1];
sx q[1];
rz(-3.0006345) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94146848) q[0];
sx q[0];
rz(-0.57766047) q[0];
sx q[0];
rz(-2.0117652) q[0];
rz(-pi) q[1];
x q[1];
rz(3.1257764) q[2];
sx q[2];
rz(-2.2770303) q[2];
sx q[2];
rz(-3.0045403) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.837622) q[1];
sx q[1];
rz(-2.6065278) q[1];
sx q[1];
rz(3.1206888) q[1];
x q[2];
rz(1.3178918) q[3];
sx q[3];
rz(-1.0974972) q[3];
sx q[3];
rz(-0.85609037) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.36833802) q[2];
sx q[2];
rz(-2.1835623) q[2];
sx q[2];
rz(0.21993302) q[2];
rz(2.1851165) q[3];
sx q[3];
rz(-0.82436162) q[3];
sx q[3];
rz(-2.1970356) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3750951) q[0];
sx q[0];
rz(-1.5506333) q[0];
sx q[0];
rz(-0.63313708) q[0];
rz(-1.3491058) q[1];
sx q[1];
rz(-0.92020412) q[1];
sx q[1];
rz(-0.79457582) q[1];
rz(2.6075543) q[2];
sx q[2];
rz(-1.8587458) q[2];
sx q[2];
rz(0.15482422) q[2];
rz(1.7344162) q[3];
sx q[3];
rz(-2.4118578) q[3];
sx q[3];
rz(2.2414997) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
