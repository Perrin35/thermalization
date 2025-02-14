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
rz(1.4577515) q[0];
sx q[0];
rz(3.2447002) q[0];
sx q[0];
rz(8.6839747) q[0];
rz(-0.15751547) q[1];
sx q[1];
rz(3.8771602) q[1];
sx q[1];
rz(9.1280042) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65771253) q[0];
sx q[0];
rz(-1.9498511) q[0];
sx q[0];
rz(-1.1227648) q[0];
rz(-2.8662017) q[2];
sx q[2];
rz(-2.4045334) q[2];
sx q[2];
rz(0.136497) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.709503) q[1];
sx q[1];
rz(-0.55050921) q[1];
sx q[1];
rz(0.014660346) q[1];
rz(-pi) q[2];
x q[2];
rz(0.17916174) q[3];
sx q[3];
rz(-1.1350313) q[3];
sx q[3];
rz(-3.1181538) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.890581) q[2];
sx q[2];
rz(-0.0958395) q[2];
sx q[2];
rz(1.4061692) q[2];
rz(-2.3928394) q[3];
sx q[3];
rz(-2.483832) q[3];
sx q[3];
rz(0.21193084) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9894079) q[0];
sx q[0];
rz(-0.87225544) q[0];
sx q[0];
rz(2.8023791) q[0];
rz(2.5665414) q[1];
sx q[1];
rz(-1.1851858) q[1];
sx q[1];
rz(2.7599879) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0038892) q[0];
sx q[0];
rz(-2.6582672) q[0];
sx q[0];
rz(1.0522196) q[0];
rz(-0.91825466) q[2];
sx q[2];
rz(-2.7236669) q[2];
sx q[2];
rz(-1.9422836) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.4424628) q[1];
sx q[1];
rz(-0.3990631) q[1];
sx q[1];
rz(-1.5459874) q[1];
rz(-pi) q[2];
x q[2];
rz(0.80581237) q[3];
sx q[3];
rz(-1.1639813) q[3];
sx q[3];
rz(2.665887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.40994689) q[2];
sx q[2];
rz(-0.68845981) q[2];
sx q[2];
rz(-2.6551969) q[2];
rz(-2.5977123) q[3];
sx q[3];
rz(-1.0433652) q[3];
sx q[3];
rz(-0.16415183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51089066) q[0];
sx q[0];
rz(-2.7257901) q[0];
sx q[0];
rz(0.99935943) q[0];
rz(-2.4103145) q[1];
sx q[1];
rz(-0.46842289) q[1];
sx q[1];
rz(-0.17459248) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0431049) q[0];
sx q[0];
rz(-0.44666651) q[0];
sx q[0];
rz(1.5763603) q[0];
x q[1];
rz(2.3056957) q[2];
sx q[2];
rz(-0.23698254) q[2];
sx q[2];
rz(-1.3224755) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.0831956) q[1];
sx q[1];
rz(-0.10279142) q[1];
sx q[1];
rz(-2.1071069) q[1];
rz(2.9231151) q[3];
sx q[3];
rz(-1.1214773) q[3];
sx q[3];
rz(1.1547417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.752855) q[2];
sx q[2];
rz(-2.2896705) q[2];
sx q[2];
rz(1.277415) q[2];
rz(-0.82469213) q[3];
sx q[3];
rz(-0.84452355) q[3];
sx q[3];
rz(0.1674913) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81553066) q[0];
sx q[0];
rz(-2.675246) q[0];
sx q[0];
rz(1.9933568) q[0];
rz(2.8094021) q[1];
sx q[1];
rz(-0.59993184) q[1];
sx q[1];
rz(1.0848328) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6713632) q[0];
sx q[0];
rz(-1.4309023) q[0];
sx q[0];
rz(-1.8663919) q[0];
rz(0.33948989) q[2];
sx q[2];
rz(-0.84622806) q[2];
sx q[2];
rz(-0.17865114) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0177637) q[1];
sx q[1];
rz(-2.4014086) q[1];
sx q[1];
rz(-3.0637118) q[1];
rz(1.1365436) q[3];
sx q[3];
rz(-1.1449185) q[3];
sx q[3];
rz(1.9445462) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.79160249) q[2];
sx q[2];
rz(-2.377066) q[2];
sx q[2];
rz(-3.1349658) q[2];
rz(-2.918112) q[3];
sx q[3];
rz(-1.0379182) q[3];
sx q[3];
rz(2.3124783) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70581907) q[0];
sx q[0];
rz(-0.75278246) q[0];
sx q[0];
rz(1.9594132) q[0];
rz(-1.301282) q[1];
sx q[1];
rz(-2.2186406) q[1];
sx q[1];
rz(0.15929793) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2625339) q[0];
sx q[0];
rz(-0.28994432) q[0];
sx q[0];
rz(-0.13061173) q[0];
rz(1.8865239) q[2];
sx q[2];
rz(-1.3906533) q[2];
sx q[2];
rz(-0.9100998) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.22896773) q[1];
sx q[1];
rz(-2.0361613) q[1];
sx q[1];
rz(-2.1123337) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.42493762) q[3];
sx q[3];
rz(-0.94858381) q[3];
sx q[3];
rz(-3.034301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.84676877) q[2];
sx q[2];
rz(-2.939665) q[2];
sx q[2];
rz(-1.3429886) q[2];
rz(-2.790847) q[3];
sx q[3];
rz(-2.0028159) q[3];
sx q[3];
rz(0.32575592) q[3];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2419234) q[0];
sx q[0];
rz(-2.3133008) q[0];
sx q[0];
rz(2.9747466) q[0];
rz(2.9150561) q[1];
sx q[1];
rz(-1.3275361) q[1];
sx q[1];
rz(-2.66364) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83376965) q[0];
sx q[0];
rz(-1.9109028) q[0];
sx q[0];
rz(-0.97619636) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7331393) q[2];
sx q[2];
rz(-0.7614927) q[2];
sx q[2];
rz(3.048427) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8357827) q[1];
sx q[1];
rz(-1.0210345) q[1];
sx q[1];
rz(-0.98819403) q[1];
rz(-2.7109954) q[3];
sx q[3];
rz(-1.9816508) q[3];
sx q[3];
rz(-1.2263067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.53092521) q[2];
sx q[2];
rz(-1.3767367) q[2];
sx q[2];
rz(1.8492071) q[2];
rz(2.2409706) q[3];
sx q[3];
rz(-0.77018046) q[3];
sx q[3];
rz(1.5463411) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
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
rz(2.6677299) q[0];
sx q[0];
rz(-0.97309363) q[0];
sx q[0];
rz(-1.5245755) q[0];
rz(-1.4432817) q[1];
sx q[1];
rz(-2.2459005) q[1];
sx q[1];
rz(-0.5009833) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1847785) q[0];
sx q[0];
rz(-1.2010472) q[0];
sx q[0];
rz(-2.2201204) q[0];
rz(-pi) q[1];
rz(0.47246859) q[2];
sx q[2];
rz(-2.5295265) q[2];
sx q[2];
rz(-1.9111567) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.1827636) q[1];
sx q[1];
rz(-2.0685968) q[1];
sx q[1];
rz(2.2341841) q[1];
rz(-pi) q[2];
rz(2.768211) q[3];
sx q[3];
rz(-0.33003673) q[3];
sx q[3];
rz(-2.9756143) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.3524126) q[2];
sx q[2];
rz(-3.0173306) q[2];
sx q[2];
rz(0.46335709) q[2];
rz(3.0739259) q[3];
sx q[3];
rz(-1.3031518) q[3];
sx q[3];
rz(0.1629924) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8891334) q[0];
sx q[0];
rz(-1.2754138) q[0];
sx q[0];
rz(-0.5109936) q[0];
rz(-0.64396089) q[1];
sx q[1];
rz(-1.124758) q[1];
sx q[1];
rz(-0.28800979) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0181684) q[0];
sx q[0];
rz(-2.9459758) q[0];
sx q[0];
rz(-2.9430812) q[0];
x q[1];
rz(2.2637469) q[2];
sx q[2];
rz(-1.0509509) q[2];
sx q[2];
rz(0.23990384) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.6161313) q[1];
sx q[1];
rz(-1.0122293) q[1];
sx q[1];
rz(1.8957047) q[1];
rz(-1.6828641) q[3];
sx q[3];
rz(-1.9470511) q[3];
sx q[3];
rz(2.3628836) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.1945232) q[2];
sx q[2];
rz(-1.1627407) q[2];
sx q[2];
rz(2.9266749) q[2];
rz(-1.3373218) q[3];
sx q[3];
rz(-2.603172) q[3];
sx q[3];
rz(2.9609093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27338481) q[0];
sx q[0];
rz(-0.93623638) q[0];
sx q[0];
rz(-2.0599763) q[0];
rz(-2.1514905) q[1];
sx q[1];
rz(-1.5568045) q[1];
sx q[1];
rz(0.33499151) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8222447) q[0];
sx q[0];
rz(-1.2480422) q[0];
sx q[0];
rz(0.035495338) q[0];
rz(-pi) q[1];
rz(-0.0063517687) q[2];
sx q[2];
rz(-1.1489023) q[2];
sx q[2];
rz(-0.61074257) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.44030467) q[1];
sx q[1];
rz(-0.13092834) q[1];
sx q[1];
rz(0.71016117) q[1];
rz(-3.1331691) q[3];
sx q[3];
rz(-1.971774) q[3];
sx q[3];
rz(3.0110735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.058502402) q[2];
sx q[2];
rz(-0.52512705) q[2];
sx q[2];
rz(2.7527909) q[2];
rz(0.36924103) q[3];
sx q[3];
rz(-2.8856314) q[3];
sx q[3];
rz(0.84454876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19287547) q[0];
sx q[0];
rz(-0.95272869) q[0];
sx q[0];
rz(-0.70190758) q[0];
rz(1.5319872) q[1];
sx q[1];
rz(-2.46789) q[1];
sx q[1];
rz(-2.7598377) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4629342) q[0];
sx q[0];
rz(-1.4890428) q[0];
sx q[0];
rz(-1.5075179) q[0];
rz(-pi) q[1];
rz(0.41093801) q[2];
sx q[2];
rz(-1.4559846) q[2];
sx q[2];
rz(-1.4474335) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4873969) q[1];
sx q[1];
rz(-1.8113448) q[1];
sx q[1];
rz(1.3122561) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.89729805) q[3];
sx q[3];
rz(-0.89400154) q[3];
sx q[3];
rz(-0.38161665) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.6751042) q[2];
sx q[2];
rz(-2.1110822) q[2];
sx q[2];
rz(-2.4934736) q[2];
rz(2.134792) q[3];
sx q[3];
rz(-1.9961793) q[3];
sx q[3];
rz(0.072331585) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5299912) q[0];
sx q[0];
rz(-1.4016822) q[0];
sx q[0];
rz(0.66521426) q[0];
rz(2.8499659) q[1];
sx q[1];
rz(-1.7886152) q[1];
sx q[1];
rz(2.0033966) q[1];
rz(-0.57403471) q[2];
sx q[2];
rz(-0.59805255) q[2];
sx q[2];
rz(-2.6274353) q[2];
rz(-2.9838647) q[3];
sx q[3];
rz(-1.1554416) q[3];
sx q[3];
rz(-0.7994061) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
