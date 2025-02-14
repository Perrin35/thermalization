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
rz(-2.3031213) q[0];
sx q[0];
rz(0.60929259) q[0];
rz(-2.4203909) q[1];
sx q[1];
rz(-2.2130122) q[1];
sx q[1];
rz(-1.0963) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.012714931) q[0];
sx q[0];
rz(-0.70438526) q[0];
sx q[0];
rz(-1.3966068) q[0];
rz(-0.9303665) q[2];
sx q[2];
rz(-1.1784679) q[2];
sx q[2];
rz(1.0898255) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.83447245) q[1];
sx q[1];
rz(-0.60170071) q[1];
sx q[1];
rz(-0.51096075) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.61756977) q[3];
sx q[3];
rz(-0.18958108) q[3];
sx q[3];
rz(-1.965234) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.2234552) q[2];
sx q[2];
rz(-1.8310941) q[2];
sx q[2];
rz(-1.9523331) q[2];
rz(-0.39092815) q[3];
sx q[3];
rz(-1.1632185) q[3];
sx q[3];
rz(1.1322359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44553462) q[0];
sx q[0];
rz(-1.3548509) q[0];
sx q[0];
rz(1.6433486) q[0];
rz(0.6779201) q[1];
sx q[1];
rz(-1.8410212) q[1];
sx q[1];
rz(2.4300785) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2133117) q[0];
sx q[0];
rz(-0.60325256) q[0];
sx q[0];
rz(-2.2343982) q[0];
rz(-pi) q[1];
rz(1.3337801) q[2];
sx q[2];
rz(-1.1188325) q[2];
sx q[2];
rz(-2.1353561) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.7603858) q[1];
sx q[1];
rz(-0.35522705) q[1];
sx q[1];
rz(-2.9095501) q[1];
rz(-pi) q[2];
rz(-1.5270803) q[3];
sx q[3];
rz(-0.62407848) q[3];
sx q[3];
rz(1.8246375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.70637643) q[2];
sx q[2];
rz(-1.7855568) q[2];
sx q[2];
rz(-2.0785275) q[2];
rz(-1.4937909) q[3];
sx q[3];
rz(-0.46896514) q[3];
sx q[3];
rz(2.4013605) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4110334) q[0];
sx q[0];
rz(-0.39091245) q[0];
sx q[0];
rz(-2.539769) q[0];
rz(-0.14566323) q[1];
sx q[1];
rz(-1.0258976) q[1];
sx q[1];
rz(0.62785968) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3981741) q[0];
sx q[0];
rz(-1.6891915) q[0];
sx q[0];
rz(-1.640624) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0192515) q[2];
sx q[2];
rz(-1.8776692) q[2];
sx q[2];
rz(0.17099525) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.3158495) q[1];
sx q[1];
rz(-1.2176732) q[1];
sx q[1];
rz(2.8984515) q[1];
x q[2];
rz(2.1842923) q[3];
sx q[3];
rz(-1.9578551) q[3];
sx q[3];
rz(1.330803) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.38498983) q[2];
sx q[2];
rz(-2.0774697) q[2];
sx q[2];
rz(2.9260054) q[2];
rz(2.0032517) q[3];
sx q[3];
rz(-1.8214858) q[3];
sx q[3];
rz(2.5707572) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0067921) q[0];
sx q[0];
rz(-1.4147867) q[0];
sx q[0];
rz(0.11882812) q[0];
rz(-1.6670594) q[1];
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
rz(-0.87619239) q[0];
sx q[0];
rz(-2.1159439) q[0];
sx q[0];
rz(0.67727526) q[0];
rz(-0.045305552) q[2];
sx q[2];
rz(-1.2477562) q[2];
sx q[2];
rz(-1.7957866) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.1203907) q[1];
sx q[1];
rz(-0.63285108) q[1];
sx q[1];
rz(-1.8561884) q[1];
rz(-pi) q[2];
x q[2];
rz(1.586257) q[3];
sx q[3];
rz(-0.8925775) q[3];
sx q[3];
rz(0.34567269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.64041758) q[2];
sx q[2];
rz(-1.4150323) q[2];
sx q[2];
rz(-0.51587063) q[2];
rz(3.0473895) q[3];
sx q[3];
rz(-0.7754063) q[3];
sx q[3];
rz(1.5532106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3696988) q[0];
sx q[0];
rz(-2.0409245) q[0];
sx q[0];
rz(1.9544253) q[0];
rz(0.59133235) q[1];
sx q[1];
rz(-1.2966803) q[1];
sx q[1];
rz(2.3147413) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2255838) q[0];
sx q[0];
rz(-1.1904612) q[0];
sx q[0];
rz(1.4404141) q[0];
rz(1.0127147) q[2];
sx q[2];
rz(-1.6396171) q[2];
sx q[2];
rz(-2.0134913) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0650658) q[1];
sx q[1];
rz(-2.4453116) q[1];
sx q[1];
rz(-1.7743006) q[1];
x q[2];
rz(0.50179568) q[3];
sx q[3];
rz(-1.0579946) q[3];
sx q[3];
rz(-2.7570908) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.7190711) q[2];
sx q[2];
rz(-0.25455385) q[2];
sx q[2];
rz(-2.1264326) q[2];
rz(-2.2319131) q[3];
sx q[3];
rz(-1.9010952) q[3];
sx q[3];
rz(2.4801586) q[3];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61843094) q[0];
sx q[0];
rz(-1.9310512) q[0];
sx q[0];
rz(-0.30995187) q[0];
rz(-3.1228206) q[1];
sx q[1];
rz(-1.4614146) q[1];
sx q[1];
rz(-0.087336691) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7531484) q[0];
sx q[0];
rz(-0.57373057) q[0];
sx q[0];
rz(3.0818207) q[0];
rz(-pi) q[1];
rz(2.3051065) q[2];
sx q[2];
rz(-0.84794551) q[2];
sx q[2];
rz(1.370887) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.8467056) q[1];
sx q[1];
rz(-2.069983) q[1];
sx q[1];
rz(-0.45392068) q[1];
rz(-pi) q[2];
rz(1.8363073) q[3];
sx q[3];
rz(-0.54406057) q[3];
sx q[3];
rz(1.2931223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7901018) q[2];
sx q[2];
rz(-0.47372207) q[2];
sx q[2];
rz(2.6215485) q[2];
rz(0.13488787) q[3];
sx q[3];
rz(-1.3210195) q[3];
sx q[3];
rz(-1.7322056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(1.0670369) q[0];
sx q[0];
rz(-1.073607) q[0];
sx q[0];
rz(0.38597646) q[0];
rz(-2.0416253) q[1];
sx q[1];
rz(-2.4620582) q[1];
sx q[1];
rz(0.78752548) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35403901) q[0];
sx q[0];
rz(-1.2830495) q[0];
sx q[0];
rz(-1.1985874) q[0];
rz(-pi) q[1];
rz(-0.50561302) q[2];
sx q[2];
rz(-2.4596301) q[2];
sx q[2];
rz(1.7069034) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.7907775) q[1];
sx q[1];
rz(-2.3190089) q[1];
sx q[1];
rz(-2.050553) q[1];
rz(-pi) q[2];
rz(2.9716216) q[3];
sx q[3];
rz(-1.8026226) q[3];
sx q[3];
rz(-1.1500203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.1197352) q[2];
sx q[2];
rz(-1.3081552) q[2];
sx q[2];
rz(-1.5137399) q[2];
rz(1.2210023) q[3];
sx q[3];
rz(-1.1648213) q[3];
sx q[3];
rz(2.6695719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30962238) q[0];
sx q[0];
rz(-1.7124875) q[0];
sx q[0];
rz(-2.9085462) q[0];
rz(1.6538992) q[1];
sx q[1];
rz(-2.1079886) q[1];
sx q[1];
rz(-1.0071365) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3171948) q[0];
sx q[0];
rz(-1.2318512) q[0];
sx q[0];
rz(0.54817731) q[0];
rz(2.7442544) q[2];
sx q[2];
rz(-2.4724744) q[2];
sx q[2];
rz(2.4934409) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5817695) q[1];
sx q[1];
rz(-1.1743465) q[1];
sx q[1];
rz(-2.9828986) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.080970244) q[3];
sx q[3];
rz(-0.24484466) q[3];
sx q[3];
rz(3.0449344) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.2795589) q[2];
sx q[2];
rz(-1.1125914) q[2];
sx q[2];
rz(-2.2588008) q[2];
rz(-0.27859303) q[3];
sx q[3];
rz(-1.9735347) q[3];
sx q[3];
rz(-0.45894233) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19926628) q[0];
sx q[0];
rz(-0.47573221) q[0];
sx q[0];
rz(1.5698154) q[0];
rz(2.3528631) q[1];
sx q[1];
rz(-2.1756344) q[1];
sx q[1];
rz(0.68005651) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90908937) q[0];
sx q[0];
rz(-1.6304558) q[0];
sx q[0];
rz(-2.3807314) q[0];
rz(-pi) q[1];
rz(0.11876688) q[2];
sx q[2];
rz(-1.4450018) q[2];
sx q[2];
rz(-2.072352) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.16526651) q[1];
sx q[1];
rz(-1.191947) q[1];
sx q[1];
rz(2.5514609) q[1];
rz(-pi) q[2];
rz(-0.54825154) q[3];
sx q[3];
rz(-2.2202949) q[3];
sx q[3];
rz(2.4626061) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1532229) q[2];
sx q[2];
rz(-2.220927) q[2];
sx q[2];
rz(-1.9752768) q[2];
rz(1.9370646) q[3];
sx q[3];
rz(-1.322999) q[3];
sx q[3];
rz(-2.5148463) q[3];
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
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.523943) q[0];
sx q[0];
rz(-2.2594422) q[0];
sx q[0];
rz(-2.7045265) q[0];
rz(2.3433459) q[1];
sx q[1];
rz(-2.6415446) q[1];
sx q[1];
rz(-0.22013586) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1967314) q[0];
sx q[0];
rz(-1.6174949) q[0];
sx q[0];
rz(2.3406505) q[0];
rz(-pi) q[1];
rz(1.2557477) q[2];
sx q[2];
rz(-1.3761259) q[2];
sx q[2];
rz(-2.6597629) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0221631) q[1];
sx q[1];
rz(-2.3009752) q[1];
sx q[1];
rz(-0.37067159) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5256568) q[3];
sx q[3];
rz(-2.3718826) q[3];
sx q[3];
rz(1.5509645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.2685214) q[2];
sx q[2];
rz(-1.2348509) q[2];
sx q[2];
rz(2.851167) q[2];
rz(2.5051266) q[3];
sx q[3];
rz(-1.3822184) q[3];
sx q[3];
rz(-1.9071473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(-2.3270522) q[0];
sx q[0];
rz(-0.70801133) q[0];
sx q[0];
rz(1.1109362) q[0];
rz(-3.0771599) q[1];
sx q[1];
rz(-1.4705407) q[1];
sx q[1];
rz(2.4992117) q[1];
rz(-3.1101075) q[2];
sx q[2];
rz(-0.97618547) q[2];
sx q[2];
rz(0.57913274) q[2];
rz(-0.53917428) q[3];
sx q[3];
rz(-1.8001582) q[3];
sx q[3];
rz(-3.0348626) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
