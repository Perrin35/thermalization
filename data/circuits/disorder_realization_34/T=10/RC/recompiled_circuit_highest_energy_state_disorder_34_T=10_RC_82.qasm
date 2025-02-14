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
rz(-0.19967747) q[0];
sx q[0];
rz(-0.59732616) q[0];
rz(-0.53961331) q[1];
sx q[1];
rz(-0.55066723) q[1];
sx q[1];
rz(0.82672969) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0866643) q[0];
sx q[0];
rz(-2.5205527) q[0];
sx q[0];
rz(-0.6426471) q[0];
rz(-pi) q[1];
rz(1.8563675) q[2];
sx q[2];
rz(-0.74215496) q[2];
sx q[2];
rz(0.56312219) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.1203314) q[1];
sx q[1];
rz(-1.2728294) q[1];
sx q[1];
rz(2.1989766) q[1];
x q[2];
rz(-2.6117418) q[3];
sx q[3];
rz(-1.4573026) q[3];
sx q[3];
rz(-0.99760357) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.716716) q[2];
sx q[2];
rz(-2.0108607) q[2];
sx q[2];
rz(-2.7190599) q[2];
rz(-1.2381964) q[3];
sx q[3];
rz(-0.4237375) q[3];
sx q[3];
rz(-2.197926) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0331405) q[0];
sx q[0];
rz(-2.5863681) q[0];
sx q[0];
rz(-0.92414498) q[0];
rz(-0.67584258) q[1];
sx q[1];
rz(-1.6249514) q[1];
sx q[1];
rz(-1.0327551) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.076457523) q[0];
sx q[0];
rz(-1.2058812) q[0];
sx q[0];
rz(-0.61000843) q[0];
rz(-pi) q[1];
rz(3.0536781) q[2];
sx q[2];
rz(-2.9682142) q[2];
sx q[2];
rz(-1.8176469) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.38235462) q[1];
sx q[1];
rz(-1.2759616) q[1];
sx q[1];
rz(1.0947541) q[1];
rz(0.61474425) q[3];
sx q[3];
rz(-0.82332506) q[3];
sx q[3];
rz(0.73574443) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5702901) q[2];
sx q[2];
rz(-1.4455659) q[2];
sx q[2];
rz(-2.9147713) q[2];
rz(-1.6972542) q[3];
sx q[3];
rz(-3.0310013) q[3];
sx q[3];
rz(2.159923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.084376) q[0];
sx q[0];
rz(-0.20350525) q[0];
sx q[0];
rz(2.2454026) q[0];
rz(2.8016688) q[1];
sx q[1];
rz(-1.8819921) q[1];
sx q[1];
rz(-0.47421727) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33297577) q[0];
sx q[0];
rz(-1.7771557) q[0];
sx q[0];
rz(1.986507) q[0];
rz(-pi) q[1];
rz(-1.3856085) q[2];
sx q[2];
rz(-1.6049002) q[2];
sx q[2];
rz(2.2540384) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.7194083) q[1];
sx q[1];
rz(-2.0046742) q[1];
sx q[1];
rz(-0.67117454) q[1];
rz(0.83714788) q[3];
sx q[3];
rz(-1.7230534) q[3];
sx q[3];
rz(0.6619795) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.0636966) q[2];
sx q[2];
rz(-1.6113969) q[2];
sx q[2];
rz(-1.873675) q[2];
rz(0.39342132) q[3];
sx q[3];
rz(-0.370341) q[3];
sx q[3];
rz(-0.89746499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1657408) q[0];
sx q[0];
rz(-0.535088) q[0];
sx q[0];
rz(-2.72056) q[0];
rz(-0.16972217) q[1];
sx q[1];
rz(-1.7531027) q[1];
sx q[1];
rz(-1.1441182) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3232358) q[0];
sx q[0];
rz(-1.5471598) q[0];
sx q[0];
rz(-0.045020176) q[0];
x q[1];
rz(3.0277191) q[2];
sx q[2];
rz(-1.1742697) q[2];
sx q[2];
rz(2.9325094) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.01466929) q[1];
sx q[1];
rz(-1.7982535) q[1];
sx q[1];
rz(-2.3421351) q[1];
x q[2];
rz(1.1781663) q[3];
sx q[3];
rz(-1.9646137) q[3];
sx q[3];
rz(1.3918685) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-3.0890395) q[2];
sx q[2];
rz(-1.384794) q[2];
sx q[2];
rz(0.65711898) q[2];
rz(2.0402015) q[3];
sx q[3];
rz(-1.907932) q[3];
sx q[3];
rz(1.5737981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14884406) q[0];
sx q[0];
rz(-1.0555457) q[0];
sx q[0];
rz(-1.8788991) q[0];
rz(-0.34072033) q[1];
sx q[1];
rz(-1.4422528) q[1];
sx q[1];
rz(-2.5943894) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.094484821) q[0];
sx q[0];
rz(-1.63103) q[0];
sx q[0];
rz(-1.3936067) q[0];
x q[1];
rz(-0.38262719) q[2];
sx q[2];
rz(-0.11906653) q[2];
sx q[2];
rz(2.5061945) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.5141827) q[1];
sx q[1];
rz(-2.1226494) q[1];
sx q[1];
rz(-0.34789209) q[1];
rz(-pi) q[2];
rz(2.2464348) q[3];
sx q[3];
rz(-0.46592281) q[3];
sx q[3];
rz(-1.8542445) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.0146694) q[2];
sx q[2];
rz(-1.9540484) q[2];
sx q[2];
rz(2.8200601) q[2];
rz(2.1899636) q[3];
sx q[3];
rz(-1.8961597) q[3];
sx q[3];
rz(-1.8754225) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5132009) q[0];
sx q[0];
rz(-0.43190792) q[0];
sx q[0];
rz(-0.54650724) q[0];
rz(0.47738099) q[1];
sx q[1];
rz(-2.7196306) q[1];
sx q[1];
rz(-1.2786) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9533449) q[0];
sx q[0];
rz(-0.45397568) q[0];
sx q[0];
rz(-2.2843642) q[0];
x q[1];
rz(1.5029656) q[2];
sx q[2];
rz(-2.1178341) q[2];
sx q[2];
rz(0.018195823) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.9467612) q[1];
sx q[1];
rz(-1.166718) q[1];
sx q[1];
rz(2.0003009) q[1];
x q[2];
rz(2.2845479) q[3];
sx q[3];
rz(-1.3856944) q[3];
sx q[3];
rz(-1.1777267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4930341) q[2];
sx q[2];
rz(-2.9849122) q[2];
sx q[2];
rz(-0.9542166) q[2];
rz(2.3844125) q[3];
sx q[3];
rz(-1.6482407) q[3];
sx q[3];
rz(0.89603388) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8947944) q[0];
sx q[0];
rz(-0.049012683) q[0];
sx q[0];
rz(0.04846305) q[0];
rz(-2.5080644) q[1];
sx q[1];
rz(-2.3506479) q[1];
sx q[1];
rz(1.4763501) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4089946) q[0];
sx q[0];
rz(-2.2267003) q[0];
sx q[0];
rz(-0.31314416) q[0];
rz(2.3728875) q[2];
sx q[2];
rz(-2.8787347) q[2];
sx q[2];
rz(-2.4117087) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.29334904) q[1];
sx q[1];
rz(-1.9696225) q[1];
sx q[1];
rz(-1.8790808) q[1];
rz(-1.3562629) q[3];
sx q[3];
rz(-1.8275785) q[3];
sx q[3];
rz(2.7329426) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.9897291) q[2];
sx q[2];
rz(-2.180474) q[2];
sx q[2];
rz(-2.9748919) q[2];
rz(-1.9549595) q[3];
sx q[3];
rz(-1.229076) q[3];
sx q[3];
rz(0.96773875) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95369059) q[0];
sx q[0];
rz(-0.1583651) q[0];
sx q[0];
rz(-0.70575356) q[0];
rz(1.7965652) q[1];
sx q[1];
rz(-1.2860362) q[1];
sx q[1];
rz(2.8154624) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46485422) q[0];
sx q[0];
rz(-2.2660672) q[0];
sx q[0];
rz(-0.37842964) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5920611) q[2];
sx q[2];
rz(-0.51409634) q[2];
sx q[2];
rz(-1.7618881) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3895411) q[1];
sx q[1];
rz(-1.4040177) q[1];
sx q[1];
rz(-0.15845815) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9891677) q[3];
sx q[3];
rz(-1.4285123) q[3];
sx q[3];
rz(-2.3438615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2976133) q[2];
sx q[2];
rz(-1.9785827) q[2];
sx q[2];
rz(-2.6824717) q[2];
rz(-1.8504359) q[3];
sx q[3];
rz(-1.2192817) q[3];
sx q[3];
rz(0.084376924) q[3];
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
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5952263) q[0];
sx q[0];
rz(-1.8084753) q[0];
sx q[0];
rz(0.03431933) q[0];
rz(-1.8345087) q[1];
sx q[1];
rz(-2.391075) q[1];
sx q[1];
rz(1.4141356) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.14489) q[0];
sx q[0];
rz(-1.7976947) q[0];
sx q[0];
rz(2.4064896) q[0];
x q[1];
rz(-0.66811647) q[2];
sx q[2];
rz(-1.2496867) q[2];
sx q[2];
rz(0.5905861) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.092897753) q[1];
sx q[1];
rz(-1.0156529) q[1];
sx q[1];
rz(1.0798961) q[1];
rz(-pi) q[2];
rz(3.121641) q[3];
sx q[3];
rz(-2.8932778) q[3];
sx q[3];
rz(-0.90398247) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.93065921) q[2];
sx q[2];
rz(-0.84676576) q[2];
sx q[2];
rz(-0.47908121) q[2];
rz(-0.72862285) q[3];
sx q[3];
rz(-1.9727547) q[3];
sx q[3];
rz(-2.5463879) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0790734) q[0];
sx q[0];
rz(-0.74051028) q[0];
sx q[0];
rz(-0.24653521) q[0];
rz(1.1277699) q[1];
sx q[1];
rz(-1.133254) q[1];
sx q[1];
rz(-2.9182428) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44137433) q[0];
sx q[0];
rz(-1.0742037) q[0];
sx q[0];
rz(1.5327647) q[0];
rz(-pi) q[1];
rz(2.8466346) q[2];
sx q[2];
rz(-2.5813537) q[2];
sx q[2];
rz(0.72232027) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.1829233) q[1];
sx q[1];
rz(-2.8374817) q[1];
sx q[1];
rz(1.7687377) q[1];
rz(-2.559313) q[3];
sx q[3];
rz(-2.6447372) q[3];
sx q[3];
rz(-2.7235746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.8862306) q[2];
sx q[2];
rz(-2.2575111) q[2];
sx q[2];
rz(-0.13038005) q[2];
rz(2.5731795) q[3];
sx q[3];
rz(-1.7805028) q[3];
sx q[3];
rz(2.7275248) q[3];
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
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9315306) q[0];
sx q[0];
rz(-0.93183403) q[0];
sx q[0];
rz(0.43781042) q[0];
rz(1.9991649) q[1];
sx q[1];
rz(-1.8985959) q[1];
sx q[1];
rz(1.7172145) q[1];
rz(0.78442153) q[2];
sx q[2];
rz(-0.57195819) q[2];
sx q[2];
rz(2.3637003) q[2];
rz(1.891248) q[3];
sx q[3];
rz(-1.1966101) q[3];
sx q[3];
rz(0.93991652) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
