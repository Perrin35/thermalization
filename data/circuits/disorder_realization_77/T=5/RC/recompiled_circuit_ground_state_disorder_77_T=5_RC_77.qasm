OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.90280688) q[0];
sx q[0];
rz(-1.1385318) q[0];
sx q[0];
rz(1.0983374) q[0];
rz(1.1605473) q[1];
sx q[1];
rz(-1.1359954) q[1];
sx q[1];
rz(2.5392037) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9928474) q[0];
sx q[0];
rz(-2.0050328) q[0];
sx q[0];
rz(0.57172267) q[0];
rz(-pi) q[1];
x q[1];
rz(0.7135993) q[2];
sx q[2];
rz(-1.6125792) q[2];
sx q[2];
rz(-3.0609727) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8849395) q[1];
sx q[1];
rz(-0.49044427) q[1];
sx q[1];
rz(1.3262733) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.1040339) q[3];
sx q[3];
rz(-1.9989487) q[3];
sx q[3];
rz(-1.960609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.41523734) q[2];
sx q[2];
rz(-2.0550315) q[2];
sx q[2];
rz(1.497867) q[2];
rz(-0.1114397) q[3];
sx q[3];
rz(-1.2846416) q[3];
sx q[3];
rz(-2.1845412) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9621256) q[0];
sx q[0];
rz(-0.92573708) q[0];
sx q[0];
rz(-0.16394462) q[0];
rz(-0.97490772) q[1];
sx q[1];
rz(-2.4539852) q[1];
sx q[1];
rz(-1.499739) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6234265) q[0];
sx q[0];
rz(-2.2742018) q[0];
sx q[0];
rz(0.33133467) q[0];
x q[1];
rz(2.4663839) q[2];
sx q[2];
rz(-0.9281425) q[2];
sx q[2];
rz(-0.86057885) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.9416391) q[1];
sx q[1];
rz(-2.2806045) q[1];
sx q[1];
rz(1.468424) q[1];
rz(-2.1466931) q[3];
sx q[3];
rz(-2.7585976) q[3];
sx q[3];
rz(2.8378311) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.6770596) q[2];
sx q[2];
rz(-0.50731069) q[2];
sx q[2];
rz(0.23059174) q[2];
rz(-0.683189) q[3];
sx q[3];
rz(-2.4468827) q[3];
sx q[3];
rz(-2.3901239) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36000073) q[0];
sx q[0];
rz(-1.8311904) q[0];
sx q[0];
rz(0.42136425) q[0];
rz(1.6257809) q[1];
sx q[1];
rz(-2.6457364) q[1];
sx q[1];
rz(0.96744195) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66146278) q[0];
sx q[0];
rz(-2.7217743) q[0];
sx q[0];
rz(-0.90645193) q[0];
rz(1.2363628) q[2];
sx q[2];
rz(-0.49116557) q[2];
sx q[2];
rz(-2.9984891) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.142621) q[1];
sx q[1];
rz(-1.843286) q[1];
sx q[1];
rz(2.9914844) q[1];
rz(-pi) q[2];
rz(-1.9032962) q[3];
sx q[3];
rz(-2.1070814) q[3];
sx q[3];
rz(0.09418776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.4030054) q[2];
sx q[2];
rz(-1.0120729) q[2];
sx q[2];
rz(2.9834413) q[2];
rz(-2.3824597) q[3];
sx q[3];
rz(-1.682155) q[3];
sx q[3];
rz(2.7706743) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.247093) q[0];
sx q[0];
rz(-2.5573754) q[0];
sx q[0];
rz(-1.1868813) q[0];
rz(-3.0184556) q[1];
sx q[1];
rz(-1.5712527) q[1];
sx q[1];
rz(1.311696) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9378915) q[0];
sx q[0];
rz(-2.9726056) q[0];
sx q[0];
rz(-2.0800872) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6475296) q[2];
sx q[2];
rz(-1.4686246) q[2];
sx q[2];
rz(0.42021423) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.0801436) q[1];
sx q[1];
rz(-1.176136) q[1];
sx q[1];
rz(-2.6861362) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7100641) q[3];
sx q[3];
rz(-3.1374806) q[3];
sx q[3];
rz(-0.63263661) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.9436403) q[2];
sx q[2];
rz(-1.3200878) q[2];
sx q[2];
rz(2.2894335) q[2];
rz(0.85396829) q[3];
sx q[3];
rz(-1.2205418) q[3];
sx q[3];
rz(1.0796237) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3581486) q[0];
sx q[0];
rz(-1.3842979) q[0];
sx q[0];
rz(0.070601687) q[0];
rz(-2.103503) q[1];
sx q[1];
rz(-1.1425273) q[1];
sx q[1];
rz(0.80931726) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0699969) q[0];
sx q[0];
rz(-0.57868273) q[0];
sx q[0];
rz(-0.13222887) q[0];
x q[1];
rz(2.3951934) q[2];
sx q[2];
rz(-1.047102) q[2];
sx q[2];
rz(-0.46800287) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(3.0496724) q[1];
sx q[1];
rz(-0.61245239) q[1];
sx q[1];
rz(0.46302192) q[1];
x q[2];
rz(2.2693997) q[3];
sx q[3];
rz(-1.707886) q[3];
sx q[3];
rz(1.7952775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.91850963) q[2];
sx q[2];
rz(-0.56879908) q[2];
sx q[2];
rz(-1.1893547) q[2];
rz(-0.14063028) q[3];
sx q[3];
rz(-0.46969241) q[3];
sx q[3];
rz(-2.752059) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4024432) q[0];
sx q[0];
rz(-0.91366714) q[0];
sx q[0];
rz(-1.5341349) q[0];
rz(0.88273478) q[1];
sx q[1];
rz(-0.84461707) q[1];
sx q[1];
rz(2.5260063) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12703757) q[0];
sx q[0];
rz(-1.1786228) q[0];
sx q[0];
rz(1.4711625) q[0];
x q[1];
rz(-2.4886905) q[2];
sx q[2];
rz(-0.98207475) q[2];
sx q[2];
rz(-1.9984286) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.5639685) q[1];
sx q[1];
rz(-1.7610487) q[1];
sx q[1];
rz(-0.14037324) q[1];
x q[2];
rz(2.1911931) q[3];
sx q[3];
rz(-1.5614034) q[3];
sx q[3];
rz(-0.91419807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.1659871) q[2];
sx q[2];
rz(-2.1773982) q[2];
sx q[2];
rz(1.9624422) q[2];
rz(-2.9278582) q[3];
sx q[3];
rz(-1.390018) q[3];
sx q[3];
rz(0.24553044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24446503) q[0];
sx q[0];
rz(-1.9944958) q[0];
sx q[0];
rz(-0.044064673) q[0];
rz(2.7524718) q[1];
sx q[1];
rz(-2.6521284) q[1];
sx q[1];
rz(-0.74180952) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.408106) q[0];
sx q[0];
rz(-1.2950962) q[0];
sx q[0];
rz(2.7635126) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3177708) q[2];
sx q[2];
rz(-1.833311) q[2];
sx q[2];
rz(-1.6561301) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.68315893) q[1];
sx q[1];
rz(-1.6583879) q[1];
sx q[1];
rz(0.32282543) q[1];
x q[2];
rz(0.70176418) q[3];
sx q[3];
rz(-1.7265111) q[3];
sx q[3];
rz(3.0879741) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.6645708) q[2];
sx q[2];
rz(-1.3571813) q[2];
sx q[2];
rz(-0.50194293) q[2];
rz(-1.4255514) q[3];
sx q[3];
rz(-0.52949667) q[3];
sx q[3];
rz(-2.3486923) q[3];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.75941706) q[0];
sx q[0];
rz(-1.734153) q[0];
sx q[0];
rz(-0.39696804) q[0];
rz(1.5396897) q[1];
sx q[1];
rz(-2.868728) q[1];
sx q[1];
rz(-2.4698965) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0881694) q[0];
sx q[0];
rz(-0.4649325) q[0];
sx q[0];
rz(1.2857807) q[0];
rz(-pi) q[1];
rz(2.7717088) q[2];
sx q[2];
rz(-1.1544168) q[2];
sx q[2];
rz(2.2781792) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(3.0782425) q[1];
sx q[1];
rz(-2.1890609) q[1];
sx q[1];
rz(1.9552014) q[1];
rz(2.193161) q[3];
sx q[3];
rz(-2.2978204) q[3];
sx q[3];
rz(-1.2326944) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.4063065) q[2];
sx q[2];
rz(-0.51484171) q[2];
sx q[2];
rz(1.3675281) q[2];
rz(-2.1788518) q[3];
sx q[3];
rz(-1.5509501) q[3];
sx q[3];
rz(-2.4875557) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.830736) q[0];
sx q[0];
rz(-1.5850569) q[0];
sx q[0];
rz(0.32064015) q[0];
rz(2.2036208) q[1];
sx q[1];
rz(-1.9357977) q[1];
sx q[1];
rz(1.7373614) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7787301) q[0];
sx q[0];
rz(-1.4674648) q[0];
sx q[0];
rz(0.10424239) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0530372) q[2];
sx q[2];
rz(-0.27805304) q[2];
sx q[2];
rz(-2.0012358) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.117489) q[1];
sx q[1];
rz(-1.5908594) q[1];
sx q[1];
rz(1.5320918) q[1];
rz(-pi) q[2];
rz(1.8102989) q[3];
sx q[3];
rz(-2.31862) q[3];
sx q[3];
rz(-2.8209958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.9414901) q[2];
sx q[2];
rz(-2.3788171) q[2];
sx q[2];
rz(-2.4007559) q[2];
rz(-2.0790993) q[3];
sx q[3];
rz(-2.0334358) q[3];
sx q[3];
rz(-1.9443996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9366539) q[0];
sx q[0];
rz(-0.42977253) q[0];
sx q[0];
rz(0.75576654) q[0];
rz(1.3142122) q[1];
sx q[1];
rz(-1.6270437) q[1];
sx q[1];
rz(-0.72602138) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1765322) q[0];
sx q[0];
rz(-2.3527761) q[0];
sx q[0];
rz(-1.5223498) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4439031) q[2];
sx q[2];
rz(-2.4524101) q[2];
sx q[2];
rz(0.5529595) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.5948702) q[1];
sx q[1];
rz(-2.7597962) q[1];
sx q[1];
rz(1.1522271) q[1];
rz(-pi) q[2];
rz(-1.2682548) q[3];
sx q[3];
rz(-2.6740251) q[3];
sx q[3];
rz(0.23307652) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.8119592) q[2];
sx q[2];
rz(-2.6128431) q[2];
sx q[2];
rz(0.94919666) q[2];
rz(-0.1720998) q[3];
sx q[3];
rz(-0.97730079) q[3];
sx q[3];
rz(2.7026091) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.028845499) q[0];
sx q[0];
rz(-1.5486568) q[0];
sx q[0];
rz(-1.5966709) q[0];
rz(1.9817837) q[1];
sx q[1];
rz(-1.3216959) q[1];
sx q[1];
rz(-1.6285223) q[1];
rz(2.8904758) q[2];
sx q[2];
rz(-0.56240964) q[2];
sx q[2];
rz(3.0887847) q[2];
rz(-3.0712745) q[3];
sx q[3];
rz(-2.0587772) q[3];
sx q[3];
rz(-1.8785431) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
