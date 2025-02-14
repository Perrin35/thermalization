OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.7665793) q[0];
sx q[0];
rz(-1.2393247) q[0];
sx q[0];
rz(1.0898606) q[0];
rz(1.0640979) q[1];
sx q[1];
rz(4.394726) q[1];
sx q[1];
rz(10.328007) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3147557) q[0];
sx q[0];
rz(-2.1776548) q[0];
sx q[0];
rz(-2.3640102) q[0];
rz(-pi) q[1];
rz(1.5527234) q[2];
sx q[2];
rz(-0.68328349) q[2];
sx q[2];
rz(-2.9073213) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-3.0923903) q[1];
sx q[1];
rz(-2.4396439) q[1];
sx q[1];
rz(-0.54767056) q[1];
rz(1.3525659) q[3];
sx q[3];
rz(-0.50980836) q[3];
sx q[3];
rz(0.17579432) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.4452867) q[2];
sx q[2];
rz(-3.004965) q[2];
sx q[2];
rz(-0.43329263) q[2];
rz(-2.8516234) q[3];
sx q[3];
rz(-2.2765997) q[3];
sx q[3];
rz(0.32052952) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1746154) q[0];
sx q[0];
rz(-2.2791635) q[0];
sx q[0];
rz(1.8509266) q[0];
rz(2.4621452) q[1];
sx q[1];
rz(-2.3652855) q[1];
sx q[1];
rz(-0.89250934) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9226869) q[0];
sx q[0];
rz(-1.7480047) q[0];
sx q[0];
rz(0.025434504) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.65373924) q[2];
sx q[2];
rz(-1.8934665) q[2];
sx q[2];
rz(1.8316837) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.8645382) q[1];
sx q[1];
rz(-1.8450341) q[1];
sx q[1];
rz(0.67233337) q[1];
rz(2.1492747) q[3];
sx q[3];
rz(-0.43799339) q[3];
sx q[3];
rz(-2.7051444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9925925) q[2];
sx q[2];
rz(-1.0789824) q[2];
sx q[2];
rz(-2.0089669) q[2];
rz(0.66631404) q[3];
sx q[3];
rz(-1.3601114) q[3];
sx q[3];
rz(1.2247491) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
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
rz(-0.76729524) q[0];
sx q[0];
rz(-0.39300028) q[0];
sx q[0];
rz(-2.1597916) q[0];
rz(2.1994195) q[1];
sx q[1];
rz(-1.3092382) q[1];
sx q[1];
rz(-2.2312677) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.019429723) q[0];
sx q[0];
rz(-1.0896026) q[0];
sx q[0];
rz(-1.0896171) q[0];
x q[1];
rz(1.4508171) q[2];
sx q[2];
rz(-1.2960805) q[2];
sx q[2];
rz(-3.0697696) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.83214789) q[1];
sx q[1];
rz(-2.6845686) q[1];
sx q[1];
rz(0.8330658) q[1];
rz(-pi) q[2];
rz(1.0873454) q[3];
sx q[3];
rz(-2.6539565) q[3];
sx q[3];
rz(0.82043649) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.2008449) q[2];
sx q[2];
rz(-0.36483279) q[2];
sx q[2];
rz(2.9546837) q[2];
rz(-0.16380353) q[3];
sx q[3];
rz(-0.87023321) q[3];
sx q[3];
rz(3.0189309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
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
rz(2.0201037) q[0];
sx q[0];
rz(-1.7262456) q[0];
sx q[0];
rz(0.69766587) q[0];
rz(2.5949219) q[1];
sx q[1];
rz(-2.8488939) q[1];
sx q[1];
rz(1.1950511) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9055525) q[0];
sx q[0];
rz(-2.4396606) q[0];
sx q[0];
rz(-2.1208956) q[0];
rz(0.61364737) q[2];
sx q[2];
rz(-0.23229182) q[2];
sx q[2];
rz(0.42343806) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4711485) q[1];
sx q[1];
rz(-0.19757825) q[1];
sx q[1];
rz(-1.0001282) q[1];
rz(-pi) q[2];
rz(-1.2863808) q[3];
sx q[3];
rz(-1.5002325) q[3];
sx q[3];
rz(-1.1190991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.4635072) q[2];
sx q[2];
rz(-1.7184075) q[2];
sx q[2];
rz(-2.9441693) q[2];
rz(3.0366963) q[3];
sx q[3];
rz(-1.1095122) q[3];
sx q[3];
rz(-3.0535898) q[3];
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
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5941493) q[0];
sx q[0];
rz(-2.6056885) q[0];
sx q[0];
rz(1.0092258) q[0];
rz(-3.1127473) q[1];
sx q[1];
rz(-1.4776769) q[1];
sx q[1];
rz(0.65863329) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26732609) q[0];
sx q[0];
rz(-1.1972885) q[0];
sx q[0];
rz(2.0167355) q[0];
rz(1.1971742) q[2];
sx q[2];
rz(-0.51643096) q[2];
sx q[2];
rz(2.2852798) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5168136) q[1];
sx q[1];
rz(-0.1923696) q[1];
sx q[1];
rz(-2.6386847) q[1];
rz(-pi) q[2];
rz(2.8675788) q[3];
sx q[3];
rz(-0.7139132) q[3];
sx q[3];
rz(-1.7108202) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.7625526) q[2];
sx q[2];
rz(-2.5950044) q[2];
sx q[2];
rz(2.9075882) q[2];
rz(-2.0512569) q[3];
sx q[3];
rz(-1.5666311) q[3];
sx q[3];
rz(2.0719297) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91319084) q[0];
sx q[0];
rz(-2.985432) q[0];
sx q[0];
rz(-1.8432023) q[0];
rz(-0.78650728) q[1];
sx q[1];
rz(-0.85580099) q[1];
sx q[1];
rz(-1.3589121) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8131065) q[0];
sx q[0];
rz(-2.4755602) q[0];
sx q[0];
rz(-0.12909992) q[0];
rz(-pi) q[1];
rz(-0.26769079) q[2];
sx q[2];
rz(-1.7867733) q[2];
sx q[2];
rz(-2.9962073) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.7715367) q[1];
sx q[1];
rz(-3.1395457) q[1];
sx q[1];
rz(2.5084346) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5547736) q[3];
sx q[3];
rz(-2.6027461) q[3];
sx q[3];
rz(-2.9485045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.7274373) q[2];
sx q[2];
rz(-1.4667908) q[2];
sx q[2];
rz(-0.1725014) q[2];
rz(0.48834673) q[3];
sx q[3];
rz(-1.9988873) q[3];
sx q[3];
rz(-2.02777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97335029) q[0];
sx q[0];
rz(-0.872648) q[0];
sx q[0];
rz(0.33997047) q[0];
rz(0.32644692) q[1];
sx q[1];
rz(-2.4531334) q[1];
sx q[1];
rz(1.2350157) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8064855) q[0];
sx q[0];
rz(-2.3705784) q[0];
sx q[0];
rz(-2.5946027) q[0];
rz(-pi) q[1];
rz(2.6831362) q[2];
sx q[2];
rz(-1.8962911) q[2];
sx q[2];
rz(-2.1874962) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.75720471) q[1];
sx q[1];
rz(-1.5761307) q[1];
sx q[1];
rz(-2.6621006) q[1];
rz(-pi) q[2];
rz(-0.49576393) q[3];
sx q[3];
rz(-0.44124441) q[3];
sx q[3];
rz(0.10291162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.1116703) q[2];
sx q[2];
rz(-0.80283529) q[2];
sx q[2];
rz(2.5862582) q[2];
rz(0.16767821) q[3];
sx q[3];
rz(-1.5372814) q[3];
sx q[3];
rz(2.663747) q[3];
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
rz(-0.70917201) q[0];
sx q[0];
rz(-0.92996159) q[0];
sx q[0];
rz(-2.0322556) q[0];
rz(-1.1322016) q[1];
sx q[1];
rz(-2.7448476) q[1];
sx q[1];
rz(-0.94165492) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1718194) q[0];
sx q[0];
rz(-2.1226774) q[0];
sx q[0];
rz(-2.9540747) q[0];
rz(-pi) q[1];
rz(1.8563111) q[2];
sx q[2];
rz(-1.7196349) q[2];
sx q[2];
rz(-0.51644737) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.4577427) q[1];
sx q[1];
rz(-2.2096429) q[1];
sx q[1];
rz(-2.9838377) q[1];
x q[2];
rz(-2.7558318) q[3];
sx q[3];
rz(-1.3732135) q[3];
sx q[3];
rz(-2.945874) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7534916) q[2];
sx q[2];
rz(-1.0270303) q[2];
sx q[2];
rz(1.6027742) q[2];
rz(-1.3711035) q[3];
sx q[3];
rz(-1.95581) q[3];
sx q[3];
rz(-0.051430844) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(0.082315363) q[0];
sx q[0];
rz(-2.5264854) q[0];
sx q[0];
rz(0.96784651) q[0];
rz(1.8702501) q[1];
sx q[1];
rz(-1.8495193) q[1];
sx q[1];
rz(3.079792) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74154749) q[0];
sx q[0];
rz(-1.1843268) q[0];
sx q[0];
rz(0.58048141) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.25904663) q[2];
sx q[2];
rz(-1.7648089) q[2];
sx q[2];
rz(-1.2736748) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.73229746) q[1];
sx q[1];
rz(-1.6662855) q[1];
sx q[1];
rz(0.68640253) q[1];
rz(-2.8093407) q[3];
sx q[3];
rz(-1.8470754) q[3];
sx q[3];
rz(-1.075923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.9515848) q[2];
sx q[2];
rz(-0.20415674) q[2];
sx q[2];
rz(3.0657213) q[2];
rz(1.9742981) q[3];
sx q[3];
rz(-1.1018402) q[3];
sx q[3];
rz(2.8372138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17700125) q[0];
sx q[0];
rz(-2.8622506) q[0];
sx q[0];
rz(-0.28468537) q[0];
rz(-3.0668861) q[1];
sx q[1];
rz(-1.9120522) q[1];
sx q[1];
rz(-1.4178735) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1847294) q[0];
sx q[0];
rz(-1.3780344) q[0];
sx q[0];
rz(2.3817029) q[0];
x q[1];
rz(-2.9921164) q[2];
sx q[2];
rz(-2.1227395) q[2];
sx q[2];
rz(1.9665444) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.7602882) q[1];
sx q[1];
rz(-0.78236303) q[1];
sx q[1];
rz(0.60790498) q[1];
rz(-pi) q[2];
rz(-1.6154434) q[3];
sx q[3];
rz(-0.64638019) q[3];
sx q[3];
rz(-2.165739) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8031926) q[2];
sx q[2];
rz(-2.000838) q[2];
sx q[2];
rz(-1.8624064) q[2];
rz(0.98215669) q[3];
sx q[3];
rz(-1.6833865) q[3];
sx q[3];
rz(0.88472432) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2411156) q[0];
sx q[0];
rz(-1.3237088) q[0];
sx q[0];
rz(1.5644912) q[0];
rz(1.0152394) q[1];
sx q[1];
rz(-1.992234) q[1];
sx q[1];
rz(2.0043859) q[1];
rz(-1.8335291) q[2];
sx q[2];
rz(-1.6752401) q[2];
sx q[2];
rz(-0.40950767) q[2];
rz(0.99526631) q[3];
sx q[3];
rz(-2.223816) q[3];
sx q[3];
rz(1.73903) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
