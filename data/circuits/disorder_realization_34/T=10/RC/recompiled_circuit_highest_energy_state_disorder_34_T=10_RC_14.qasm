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
rz(1.0518987) q[0];
sx q[0];
rz(2.9419152) q[0];
sx q[0];
rz(10.022104) q[0];
rz(-0.53961331) q[1];
sx q[1];
rz(-0.55066723) q[1];
sx q[1];
rz(0.82672969) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.054928314) q[0];
sx q[0];
rz(-2.5205527) q[0];
sx q[0];
rz(-2.4989456) q[0];
x q[1];
rz(2.2924001) q[2];
sx q[2];
rz(-1.7623644) q[2];
sx q[2];
rz(-1.9208197) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.2079426) q[1];
sx q[1];
rz(-2.4550555) q[1];
sx q[1];
rz(1.0892434) q[1];
rz(-1.7021322) q[3];
sx q[3];
rz(-1.0447096) q[3];
sx q[3];
rz(0.50695723) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.4248767) q[2];
sx q[2];
rz(-1.1307319) q[2];
sx q[2];
rz(2.7190599) q[2];
rz(-1.2381964) q[3];
sx q[3];
rz(-2.7178552) q[3];
sx q[3];
rz(-0.94366664) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1084522) q[0];
sx q[0];
rz(-2.5863681) q[0];
sx q[0];
rz(-2.2174477) q[0];
rz(-0.67584258) q[1];
sx q[1];
rz(-1.5166413) q[1];
sx q[1];
rz(-2.1088375) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.402814) q[0];
sx q[0];
rz(-2.1354851) q[0];
sx q[0];
rz(-2.0069471) q[0];
rz(-pi) q[1];
x q[1];
rz(0.17272213) q[2];
sx q[2];
rz(-1.5859436) q[2];
sx q[2];
rz(0.33345371) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.7016587) q[1];
sx q[1];
rz(-0.55393078) q[1];
sx q[1];
rz(2.1560378) q[1];
rz(0.7224222) q[3];
sx q[3];
rz(-2.0075433) q[3];
sx q[3];
rz(-1.8590761) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.57130259) q[2];
sx q[2];
rz(-1.4455659) q[2];
sx q[2];
rz(-2.9147713) q[2];
rz(-1.4443385) q[3];
sx q[3];
rz(-3.0310013) q[3];
sx q[3];
rz(-2.159923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0572166) q[0];
sx q[0];
rz(-0.20350525) q[0];
sx q[0];
rz(-2.2454026) q[0];
rz(-0.33992386) q[1];
sx q[1];
rz(-1.2596005) q[1];
sx q[1];
rz(-2.6673754) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6723876) q[0];
sx q[0];
rz(-0.46142277) q[0];
sx q[0];
rz(-2.049033) q[0];
rz(1.7540054) q[2];
sx q[2];
rz(-0.18826655) q[2];
sx q[2];
rz(2.6383924) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.82638559) q[1];
sx q[1];
rz(-0.97123324) q[1];
sx q[1];
rz(1.0365328) q[1];
rz(0.83714788) q[3];
sx q[3];
rz(-1.4185393) q[3];
sx q[3];
rz(-0.6619795) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.077896) q[2];
sx q[2];
rz(-1.5301957) q[2];
sx q[2];
rz(1.873675) q[2];
rz(2.7481713) q[3];
sx q[3];
rz(-0.370341) q[3];
sx q[3];
rz(-2.2441277) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1657408) q[0];
sx q[0];
rz(-0.535088) q[0];
sx q[0];
rz(-2.72056) q[0];
rz(-2.9718705) q[1];
sx q[1];
rz(-1.7531027) q[1];
sx q[1];
rz(-1.9974744) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73070684) q[0];
sx q[0];
rz(-3.0907486) q[0];
sx q[0];
rz(-0.4836785) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0277191) q[2];
sx q[2];
rz(-1.9673229) q[2];
sx q[2];
rz(-0.20908326) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.01466929) q[1];
sx q[1];
rz(-1.7982535) q[1];
sx q[1];
rz(2.3421351) q[1];
x q[2];
rz(1.1781663) q[3];
sx q[3];
rz(-1.9646137) q[3];
sx q[3];
rz(-1.7497241) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.0890395) q[2];
sx q[2];
rz(-1.7567987) q[2];
sx q[2];
rz(2.4844737) q[2];
rz(2.0402015) q[3];
sx q[3];
rz(-1.907932) q[3];
sx q[3];
rz(1.5737981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9927486) q[0];
sx q[0];
rz(-2.086047) q[0];
sx q[0];
rz(-1.8788991) q[0];
rz(2.8008723) q[1];
sx q[1];
rz(-1.4422528) q[1];
sx q[1];
rz(0.54720324) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4870905) q[0];
sx q[0];
rz(-1.3939314) q[0];
sx q[0];
rz(-3.0804033) q[0];
x q[1];
rz(2.7589655) q[2];
sx q[2];
rz(-3.0225261) q[2];
sx q[2];
rz(-2.5061945) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.9090404) q[1];
sx q[1];
rz(-2.4989933) q[1];
sx q[1];
rz(-2.0765096) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9447127) q[3];
sx q[3];
rz(-1.8555879) q[3];
sx q[3];
rz(2.2367331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.0146694) q[2];
sx q[2];
rz(-1.9540484) q[2];
sx q[2];
rz(-0.32153258) q[2];
rz(2.1899636) q[3];
sx q[3];
rz(-1.245433) q[3];
sx q[3];
rz(-1.2661701) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5132009) q[0];
sx q[0];
rz(-2.7096847) q[0];
sx q[0];
rz(0.54650724) q[0];
rz(-2.6642117) q[1];
sx q[1];
rz(-2.7196306) q[1];
sx q[1];
rz(-1.2786) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7200658) q[0];
sx q[0];
rz(-1.2328487) q[0];
sx q[0];
rz(-0.30915156) q[0];
rz(-pi) q[1];
rz(-0.11083229) q[2];
sx q[2];
rz(-2.5907907) q[2];
sx q[2];
rz(-3.0299195) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9467612) q[1];
sx q[1];
rz(-1.9748747) q[1];
sx q[1];
rz(1.1412918) q[1];
rz(-pi) q[2];
rz(-2.8987721) q[3];
sx q[3];
rz(-2.2698463) q[3];
sx q[3];
rz(-2.9065913) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.4930341) q[2];
sx q[2];
rz(-0.15668046) q[2];
sx q[2];
rz(-2.1873761) q[2];
rz(-0.75718015) q[3];
sx q[3];
rz(-1.4933519) q[3];
sx q[3];
rz(2.2455588) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2467982) q[0];
sx q[0];
rz(-3.09258) q[0];
sx q[0];
rz(-3.0931296) q[0];
rz(-2.5080644) q[1];
sx q[1];
rz(-0.79094473) q[1];
sx q[1];
rz(1.6652426) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7325981) q[0];
sx q[0];
rz(-2.2267003) q[0];
sx q[0];
rz(0.31314416) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7557275) q[2];
sx q[2];
rz(-1.382916) q[2];
sx q[2];
rz(-0.056294346) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.7411091) q[1];
sx q[1];
rz(-1.2874075) q[1];
sx q[1];
rz(2.7251935) q[1];
rz(-2.879056) q[3];
sx q[3];
rz(-1.3634014) q[3];
sx q[3];
rz(-1.9241662) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.1518636) q[2];
sx q[2];
rz(-2.180474) q[2];
sx q[2];
rz(0.1667008) q[2];
rz(1.1866331) q[3];
sx q[3];
rz(-1.9125166) q[3];
sx q[3];
rz(2.1738539) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1879021) q[0];
sx q[0];
rz(-2.9832276) q[0];
sx q[0];
rz(2.4358391) q[0];
rz(-1.7965652) q[1];
sx q[1];
rz(-1.2860362) q[1];
sx q[1];
rz(-2.8154624) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7862575) q[0];
sx q[0];
rz(-1.8584492) q[0];
sx q[0];
rz(0.83922235) q[0];
x q[1];
rz(-1.5920611) q[2];
sx q[2];
rz(-0.51409634) q[2];
sx q[2];
rz(1.7618881) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.62319726) q[1];
sx q[1];
rz(-2.9120486) q[1];
sx q[1];
rz(-2.3238682) q[1];
rz(2.9860849) q[3];
sx q[3];
rz(-1.9846791) q[3];
sx q[3];
rz(-2.3055616) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.2976133) q[2];
sx q[2];
rz(-1.16301) q[2];
sx q[2];
rz(0.45912099) q[2];
rz(-1.2911568) q[3];
sx q[3];
rz(-1.2192817) q[3];
sx q[3];
rz(3.0572157) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5463663) q[0];
sx q[0];
rz(-1.3331174) q[0];
sx q[0];
rz(-0.03431933) q[0];
rz(-1.307084) q[1];
sx q[1];
rz(-2.391075) q[1];
sx q[1];
rz(1.7274571) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.14489) q[0];
sx q[0];
rz(-1.7976947) q[0];
sx q[0];
rz(-0.73510304) q[0];
rz(-pi) q[1];
x q[1];
rz(0.66811647) q[2];
sx q[2];
rz(-1.2496867) q[2];
sx q[2];
rz(2.5510066) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.203277) q[1];
sx q[1];
rz(-1.1585981) q[1];
sx q[1];
rz(-2.5286872) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.121641) q[3];
sx q[3];
rz(-2.8932778) q[3];
sx q[3];
rz(0.90398247) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2109334) q[2];
sx q[2];
rz(-2.2948269) q[2];
sx q[2];
rz(2.6625114) q[2];
rz(-2.4129698) q[3];
sx q[3];
rz(-1.9727547) q[3];
sx q[3];
rz(-0.59520477) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(1.0625192) q[0];
sx q[0];
rz(-2.4010824) q[0];
sx q[0];
rz(-2.8950574) q[0];
rz(1.1277699) q[1];
sx q[1];
rz(-2.0083387) q[1];
sx q[1];
rz(2.9182428) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52106954) q[0];
sx q[0];
rz(-0.4979254) q[0];
sx q[0];
rz(-3.0715406) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.29495802) q[2];
sx q[2];
rz(-0.56023894) q[2];
sx q[2];
rz(2.4192724) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.9586693) q[1];
sx q[1];
rz(-0.30411094) q[1];
sx q[1];
rz(-1.3728549) q[1];
rz(-pi) q[2];
rz(0.58227964) q[3];
sx q[3];
rz(-2.6447372) q[3];
sx q[3];
rz(0.41801807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8862306) q[2];
sx q[2];
rz(-2.2575111) q[2];
sx q[2];
rz(-3.0112126) q[2];
rz(2.5731795) q[3];
sx q[3];
rz(-1.7805028) q[3];
sx q[3];
rz(-0.4140678) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2100621) q[0];
sx q[0];
rz(-0.93183403) q[0];
sx q[0];
rz(0.43781042) q[0];
rz(1.1424278) q[1];
sx q[1];
rz(-1.2429968) q[1];
sx q[1];
rz(-1.4243781) q[1];
rz(2.7140638) q[2];
sx q[2];
rz(-1.9631546) q[2];
sx q[2];
rz(-3.0468804) q[2];
rz(1.2503446) q[3];
sx q[3];
rz(-1.9449825) q[3];
sx q[3];
rz(-2.2016761) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
