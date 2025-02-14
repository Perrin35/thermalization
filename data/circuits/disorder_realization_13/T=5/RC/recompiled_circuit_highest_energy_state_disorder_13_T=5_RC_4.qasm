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
rz(0.19652551) q[0];
sx q[0];
rz(-0.82634574) q[0];
sx q[0];
rz(-2.2235121) q[0];
rz(3.0300568) q[1];
sx q[1];
rz(-1.7938951) q[1];
sx q[1];
rz(1.5741875) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.471431) q[0];
sx q[0];
rz(-1.5802556) q[0];
sx q[0];
rz(-0.0034250101) q[0];
rz(-pi) q[1];
rz(-0.56212933) q[2];
sx q[2];
rz(-1.8888753) q[2];
sx q[2];
rz(-2.0462917) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5985377) q[1];
sx q[1];
rz(-0.32078136) q[1];
sx q[1];
rz(1.9108652) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.918675) q[3];
sx q[3];
rz(-0.50361982) q[3];
sx q[3];
rz(0.38583392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.427318) q[2];
sx q[2];
rz(-1.486342) q[2];
sx q[2];
rz(1.8753258) q[2];
rz(-1.0424987) q[3];
sx q[3];
rz(-1.6008585) q[3];
sx q[3];
rz(-1.7460167) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.105044) q[0];
sx q[0];
rz(-2.5140913) q[0];
sx q[0];
rz(0.096916048) q[0];
rz(-2.3333683) q[1];
sx q[1];
rz(-0.40882912) q[1];
sx q[1];
rz(1.854863) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7185583) q[0];
sx q[0];
rz(-1.1211532) q[0];
sx q[0];
rz(3.1021523) q[0];
rz(-pi) q[1];
x q[1];
rz(2.497235) q[2];
sx q[2];
rz(-2.3080809) q[2];
sx q[2];
rz(-0.65284158) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3073439) q[1];
sx q[1];
rz(-1.1888224) q[1];
sx q[1];
rz(-0.81484199) q[1];
x q[2];
rz(-2.3318491) q[3];
sx q[3];
rz(-0.47172037) q[3];
sx q[3];
rz(-2.7641486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.6404932) q[2];
sx q[2];
rz(-2.2576136) q[2];
sx q[2];
rz(0.38197771) q[2];
rz(2.6169418) q[3];
sx q[3];
rz(-1.4068406) q[3];
sx q[3];
rz(0.14373556) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76661888) q[0];
sx q[0];
rz(-1.5559649) q[0];
sx q[0];
rz(-0.68614706) q[0];
rz(-1.067591) q[1];
sx q[1];
rz(-2.144404) q[1];
sx q[1];
rz(3.0608665) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9322745) q[0];
sx q[0];
rz(-2.9146412) q[0];
sx q[0];
rz(-1.7448241) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.92442928) q[2];
sx q[2];
rz(-1.0253128) q[2];
sx q[2];
rz(1.1875064) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.9517652) q[1];
sx q[1];
rz(-1.8009342) q[1];
sx q[1];
rz(-1.5946424) q[1];
rz(0.46960652) q[3];
sx q[3];
rz(-0.84495832) q[3];
sx q[3];
rz(-0.71401789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.8778235) q[2];
sx q[2];
rz(-0.70085415) q[2];
sx q[2];
rz(-0.43886718) q[2];
rz(-1.8435439) q[3];
sx q[3];
rz(-0.38709199) q[3];
sx q[3];
rz(-2.7264061) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-0.91442672) q[0];
sx q[0];
rz(-0.54410797) q[0];
sx q[0];
rz(2.4650204) q[0];
rz(0.1380955) q[1];
sx q[1];
rz(-2.0510249) q[1];
sx q[1];
rz(0.20060435) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56708589) q[0];
sx q[0];
rz(-2.4100523) q[0];
sx q[0];
rz(-1.0163496) q[0];
rz(2.8851231) q[2];
sx q[2];
rz(-2.4962728) q[2];
sx q[2];
rz(-2.8706474) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.8183867) q[1];
sx q[1];
rz(-0.96508677) q[1];
sx q[1];
rz(1.731864) q[1];
rz(0.225876) q[3];
sx q[3];
rz(-2.568733) q[3];
sx q[3];
rz(1.7879888) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.6382008) q[2];
sx q[2];
rz(-1.2752504) q[2];
sx q[2];
rz(1.7207883) q[2];
rz(0.11451379) q[3];
sx q[3];
rz(-1.6679461) q[3];
sx q[3];
rz(0.99944559) q[3];
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
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3223709) q[0];
sx q[0];
rz(-0.79077661) q[0];
sx q[0];
rz(0.97212273) q[0];
rz(-0.32360336) q[1];
sx q[1];
rz(-1.4515896) q[1];
sx q[1];
rz(1.1824898) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33859461) q[0];
sx q[0];
rz(-1.3756244) q[0];
sx q[0];
rz(0.040935733) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.57020541) q[2];
sx q[2];
rz(-1.5657164) q[2];
sx q[2];
rz(3.0722116) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.2731635) q[1];
sx q[1];
rz(-1.7291964) q[1];
sx q[1];
rz(3.0636154) q[1];
rz(-2.5706815) q[3];
sx q[3];
rz(-2.0200673) q[3];
sx q[3];
rz(0.8415287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.16367308) q[2];
sx q[2];
rz(-1.369643) q[2];
sx q[2];
rz(2.6046806) q[2];
rz(0.72530693) q[3];
sx q[3];
rz(-0.0497497) q[3];
sx q[3];
rz(-2.3427826) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2742915) q[0];
sx q[0];
rz(-2.0766356) q[0];
sx q[0];
rz(1.6987479) q[0];
rz(0.2105712) q[1];
sx q[1];
rz(-1.6981533) q[1];
sx q[1];
rz(0.64705667) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7919901) q[0];
sx q[0];
rz(-2.6812883) q[0];
sx q[0];
rz(-0.60636284) q[0];
rz(1.4616632) q[2];
sx q[2];
rz(-2.1246548) q[2];
sx q[2];
rz(0.46968974) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.3270204) q[1];
sx q[1];
rz(-1.2271735) q[1];
sx q[1];
rz(2.6966641) q[1];
rz(2.4652867) q[3];
sx q[3];
rz(-2.7190373) q[3];
sx q[3];
rz(0.12862118) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.1232542) q[2];
sx q[2];
rz(-2.047057) q[2];
sx q[2];
rz(1.9617762) q[2];
rz(1.2849464) q[3];
sx q[3];
rz(-0.40950567) q[3];
sx q[3];
rz(0.063260945) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1061851) q[0];
sx q[0];
rz(-1.4466865) q[0];
sx q[0];
rz(-2.2440198) q[0];
rz(-1.3342185) q[1];
sx q[1];
rz(-1.7709657) q[1];
sx q[1];
rz(2.1468377) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9809231) q[0];
sx q[0];
rz(-0.86359016) q[0];
sx q[0];
rz(-0.41236931) q[0];
rz(2.1734235) q[2];
sx q[2];
rz(-1.1873086) q[2];
sx q[2];
rz(-0.17959514) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9092332) q[1];
sx q[1];
rz(-2.5933752) q[1];
sx q[1];
rz(-2.1562804) q[1];
rz(-1.8514824) q[3];
sx q[3];
rz(-1.5475905) q[3];
sx q[3];
rz(-0.9965903) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.0576386) q[2];
sx q[2];
rz(-1.3292162) q[2];
sx q[2];
rz(2.8482021) q[2];
rz(-3.0883664) q[3];
sx q[3];
rz(-0.86881995) q[3];
sx q[3];
rz(1.4330385) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5720125) q[0];
sx q[0];
rz(-1.7520289) q[0];
sx q[0];
rz(0.30297512) q[0];
rz(1.6527893) q[1];
sx q[1];
rz(-1.3158512) q[1];
sx q[1];
rz(0.66910076) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44620546) q[0];
sx q[0];
rz(-1.900937) q[0];
sx q[0];
rz(2.6690525) q[0];
rz(-1.9225408) q[2];
sx q[2];
rz(-1.3430077) q[2];
sx q[2];
rz(-2.5643189) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.90780202) q[1];
sx q[1];
rz(-1.0668584) q[1];
sx q[1];
rz(-2.6556117) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7176526) q[3];
sx q[3];
rz(-1.3518847) q[3];
sx q[3];
rz(3.0242137) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.1423433) q[2];
sx q[2];
rz(-1.3331022) q[2];
sx q[2];
rz(0.86714253) q[2];
rz(-2.0182746) q[3];
sx q[3];
rz(-2.2893548) q[3];
sx q[3];
rz(1.1423906) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.860723) q[0];
sx q[0];
rz(-0.46361247) q[0];
sx q[0];
rz(3.0238357) q[0];
rz(-2.2640758) q[1];
sx q[1];
rz(-2.120647) q[1];
sx q[1];
rz(2.1868475) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3316393) q[0];
sx q[0];
rz(-2.3970986) q[0];
sx q[0];
rz(-0.5793504) q[0];
rz(1.0756726) q[2];
sx q[2];
rz(-1.1469764) q[2];
sx q[2];
rz(-2.0002805) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2429208) q[1];
sx q[1];
rz(-1.3551215) q[1];
sx q[1];
rz(2.2663003) q[1];
rz(-2.4525663) q[3];
sx q[3];
rz(-0.9456767) q[3];
sx q[3];
rz(0.90009119) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.2473601) q[2];
sx q[2];
rz(-0.33984137) q[2];
sx q[2];
rz(1.4840688) q[2];
rz(-1.5225211) q[3];
sx q[3];
rz(-1.7584453) q[3];
sx q[3];
rz(1.9671666) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5840983) q[0];
sx q[0];
rz(-1.431594) q[0];
sx q[0];
rz(-2.3884921) q[0];
rz(1.9901216) q[1];
sx q[1];
rz(-2.617372) q[1];
sx q[1];
rz(-1.6023191) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95349089) q[0];
sx q[0];
rz(-1.7537781) q[0];
sx q[0];
rz(2.0273927) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.50984211) q[2];
sx q[2];
rz(-0.93889369) q[2];
sx q[2];
rz(2.2930278) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.6021612) q[1];
sx q[1];
rz(-1.9564346) q[1];
sx q[1];
rz(-1.7785317) q[1];
x q[2];
rz(1.8352103) q[3];
sx q[3];
rz(-2.2087277) q[3];
sx q[3];
rz(-0.3141981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.3573542) q[2];
sx q[2];
rz(-2.9604762) q[2];
sx q[2];
rz(2.6757346) q[2];
rz(1.0957796) q[3];
sx q[3];
rz(-2.1491137) q[3];
sx q[3];
rz(-2.5994658) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0453542) q[0];
sx q[0];
rz(-1.5662554) q[0];
sx q[0];
rz(1.5731496) q[0];
rz(-1.0678328) q[1];
sx q[1];
rz(-2.7012431) q[1];
sx q[1];
rz(2.3213097) q[1];
rz(1.3374511) q[2];
sx q[2];
rz(-1.8276311) q[2];
sx q[2];
rz(1.7414321) q[2];
rz(-2.0910083) q[3];
sx q[3];
rz(-2.0512085) q[3];
sx q[3];
rz(-2.7762085) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
