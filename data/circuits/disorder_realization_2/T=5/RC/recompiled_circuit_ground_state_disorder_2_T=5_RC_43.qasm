OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.53192294) q[0];
sx q[0];
rz(-1.5250396) q[0];
sx q[0];
rz(-1.8833696) q[0];
rz(-1.582107) q[1];
sx q[1];
rz(-2.6546302) q[1];
sx q[1];
rz(-0.43509405) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5601016) q[0];
sx q[0];
rz(-1.035443) q[0];
sx q[0];
rz(1.0161607) q[0];
rz(-pi) q[1];
rz(-0.53768208) q[2];
sx q[2];
rz(-1.8694365) q[2];
sx q[2];
rz(-0.60288211) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.9261998) q[1];
sx q[1];
rz(-1.6476616) q[1];
sx q[1];
rz(-1.8410626) q[1];
rz(-pi) q[2];
rz(0.70617999) q[3];
sx q[3];
rz(-0.67863388) q[3];
sx q[3];
rz(-2.5265001) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.9293999) q[2];
sx q[2];
rz(-1.2779028) q[2];
sx q[2];
rz(-2.0983992) q[2];
rz(-2.7405401) q[3];
sx q[3];
rz(-1.6845104) q[3];
sx q[3];
rz(-0.57798398) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22724085) q[0];
sx q[0];
rz(-0.14509097) q[0];
sx q[0];
rz(-0.57120728) q[0];
rz(-0.97194833) q[1];
sx q[1];
rz(-2.4010039) q[1];
sx q[1];
rz(1.1245701) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6475424) q[0];
sx q[0];
rz(-0.66141539) q[0];
sx q[0];
rz(-1.1986092) q[0];
rz(2.4907465) q[2];
sx q[2];
rz(-1.2799601) q[2];
sx q[2];
rz(0.42465045) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.31370762) q[1];
sx q[1];
rz(-1.6646302) q[1];
sx q[1];
rz(0.83637823) q[1];
rz(0.092512802) q[3];
sx q[3];
rz(-2.9399151) q[3];
sx q[3];
rz(0.31615931) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.22175114) q[2];
sx q[2];
rz(-1.5849179) q[2];
sx q[2];
rz(2.705503) q[2];
rz(-2.3927205) q[3];
sx q[3];
rz(-2.336453) q[3];
sx q[3];
rz(-0.82912904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0332758) q[0];
sx q[0];
rz(-1.3949787) q[0];
sx q[0];
rz(-0.088031553) q[0];
rz(-2.2875359) q[1];
sx q[1];
rz(-2.4805534) q[1];
sx q[1];
rz(-2.0565654) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2449834) q[0];
sx q[0];
rz(-1.6379781) q[0];
sx q[0];
rz(1.843473) q[0];
rz(1.7872115) q[2];
sx q[2];
rz(-1.7209574) q[2];
sx q[2];
rz(0.79126287) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.8768118) q[1];
sx q[1];
rz(-2.0088043) q[1];
sx q[1];
rz(-2.9609738) q[1];
x q[2];
rz(1.4973699) q[3];
sx q[3];
rz(-1.182397) q[3];
sx q[3];
rz(-1.8457613) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.6701086) q[2];
sx q[2];
rz(-1.4205616) q[2];
sx q[2];
rz(-2.3954083) q[2];
rz(0.21607312) q[3];
sx q[3];
rz(-1.8862628) q[3];
sx q[3];
rz(0.20572534) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
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
rz(-1.7699319) q[0];
sx q[0];
rz(-3.0682204) q[0];
sx q[0];
rz(1.368847) q[0];
rz(2.5313077) q[1];
sx q[1];
rz(-1.3657602) q[1];
sx q[1];
rz(0.38017166) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73715009) q[0];
sx q[0];
rz(-1.8821054) q[0];
sx q[0];
rz(-2.690517) q[0];
rz(0.39003464) q[2];
sx q[2];
rz(-2.6594498) q[2];
sx q[2];
rz(2.6576561) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-3.026555) q[1];
sx q[1];
rz(-1.7174182) q[1];
sx q[1];
rz(2.0762334) q[1];
rz(-1.972547) q[3];
sx q[3];
rz(-1.6309079) q[3];
sx q[3];
rz(1.7463804) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.2569106) q[2];
sx q[2];
rz(-2.1941049) q[2];
sx q[2];
rz(-1.0379637) q[2];
rz(-0.3642309) q[3];
sx q[3];
rz(-2.3528152) q[3];
sx q[3];
rz(-2.4558333) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45201856) q[0];
sx q[0];
rz(-1.4301393) q[0];
sx q[0];
rz(2.3110287) q[0];
rz(-3.0914302) q[1];
sx q[1];
rz(-0.12718931) q[1];
sx q[1];
rz(0.62279472) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68428295) q[0];
sx q[0];
rz(-1.3535392) q[0];
sx q[0];
rz(-1.4909677) q[0];
x q[1];
rz(-0.41983126) q[2];
sx q[2];
rz(-11/(2*pi)) q[2];
sx q[2];
rz(2.3617488) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.078256814) q[1];
sx q[1];
rz(-1.9252658) q[1];
sx q[1];
rz(-0.88651231) q[1];
rz(-pi) q[2];
rz(-0.022449724) q[3];
sx q[3];
rz(-1.5708617) q[3];
sx q[3];
rz(1.8409468) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.097816618) q[2];
sx q[2];
rz(-2.706683) q[2];
sx q[2];
rz(2.3338976) q[2];
rz(-1.5378599) q[3];
sx q[3];
rz(-1.6350428) q[3];
sx q[3];
rz(0.22411552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1767126) q[0];
sx q[0];
rz(-1.8105312) q[0];
sx q[0];
rz(-1.6656732) q[0];
rz(-1.2337947) q[1];
sx q[1];
rz(-1.7314311) q[1];
sx q[1];
rz(-2.9880611) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4013256) q[0];
sx q[0];
rz(-0.10970481) q[0];
sx q[0];
rz(0.033880635) q[0];
x q[1];
rz(-2.3014499) q[2];
sx q[2];
rz(-0.62799997) q[2];
sx q[2];
rz(0.59115228) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.081806) q[1];
sx q[1];
rz(-2.012552) q[1];
sx q[1];
rz(0.55952832) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4755523) q[3];
sx q[3];
rz(-2.207334) q[3];
sx q[3];
rz(-1.7004636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.4393356) q[2];
sx q[2];
rz(-0.80436891) q[2];
sx q[2];
rz(2.5852933) q[2];
rz(0.016911658) q[3];
sx q[3];
rz(-0.97512475) q[3];
sx q[3];
rz(-0.66439605) q[3];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9282114) q[0];
sx q[0];
rz(-3.0688372) q[0];
sx q[0];
rz(0.67759204) q[0];
rz(-1.821359) q[1];
sx q[1];
rz(-1.0878539) q[1];
sx q[1];
rz(0.56603146) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2293733) q[0];
sx q[0];
rz(-1.1347596) q[0];
sx q[0];
rz(-1.1993809) q[0];
rz(-pi) q[1];
rz(-1.8458869) q[2];
sx q[2];
rz(-3.0873211) q[2];
sx q[2];
rz(-1.6942629) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2905419) q[1];
sx q[1];
rz(-0.7194964) q[1];
sx q[1];
rz(1.3403203) q[1];
rz(-pi) q[2];
rz(-1.5291847) q[3];
sx q[3];
rz(-2.0609612) q[3];
sx q[3];
rz(3.0050638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.90938202) q[2];
sx q[2];
rz(-2.7558694) q[2];
sx q[2];
rz(0.56987008) q[2];
rz(-2.0991523) q[3];
sx q[3];
rz(-1.180155) q[3];
sx q[3];
rz(-1.0038143) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7637699) q[0];
sx q[0];
rz(-0.42709392) q[0];
sx q[0];
rz(2.1183993) q[0];
rz(-0.91776735) q[1];
sx q[1];
rz(-0.75420165) q[1];
sx q[1];
rz(0.86407152) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22485944) q[0];
sx q[0];
rz(-1.8141439) q[0];
sx q[0];
rz(3.0152074) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7644453) q[2];
sx q[2];
rz(-1.1349863) q[2];
sx q[2];
rz(-0.47919861) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.9814216) q[1];
sx q[1];
rz(-0.87620253) q[1];
sx q[1];
rz(0.92215529) q[1];
rz(-pi) q[2];
rz(0.72743639) q[3];
sx q[3];
rz(-1.2890649) q[3];
sx q[3];
rz(1.4294595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.0920022) q[2];
sx q[2];
rz(-1.2668173) q[2];
sx q[2];
rz(-2.3769296) q[2];
rz(2.2406254) q[3];
sx q[3];
rz(-2.4668906) q[3];
sx q[3];
rz(2.0425792) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51729274) q[0];
sx q[0];
rz(-2.7462672) q[0];
sx q[0];
rz(-0.98315352) q[0];
rz(2.3120425) q[1];
sx q[1];
rz(-1.364565) q[1];
sx q[1];
rz(0.57873908) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.549141) q[0];
sx q[0];
rz(-1.5835174) q[0];
sx q[0];
rz(-1.1770269) q[0];
x q[1];
rz(-2.5972967) q[2];
sx q[2];
rz(-2.0434663) q[2];
sx q[2];
rz(0.67677697) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.2794046) q[1];
sx q[1];
rz(-1.9721689) q[1];
sx q[1];
rz(-2.1251232) q[1];
x q[2];
rz(-3.0399801) q[3];
sx q[3];
rz(-0.34694537) q[3];
sx q[3];
rz(2.8941151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.65043727) q[2];
sx q[2];
rz(-1.8530242) q[2];
sx q[2];
rz(3.0299661) q[2];
rz(2.1857183) q[3];
sx q[3];
rz(-2.1501232) q[3];
sx q[3];
rz(1.2607695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2535506) q[0];
sx q[0];
rz(-1.0567867) q[0];
sx q[0];
rz(-3.1274617) q[0];
rz(2.0290532) q[1];
sx q[1];
rz(-0.91347778) q[1];
sx q[1];
rz(1.6003476) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3378572) q[0];
sx q[0];
rz(-1.7183665) q[0];
sx q[0];
rz(1.7002559) q[0];
rz(1.4038744) q[2];
sx q[2];
rz(-0.65123616) q[2];
sx q[2];
rz(1.9211574) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.75493357) q[1];
sx q[1];
rz(-1.5112097) q[1];
sx q[1];
rz(-2.7769376) q[1];
rz(-pi) q[2];
x q[2];
rz(1.735481) q[3];
sx q[3];
rz(-1.4030808) q[3];
sx q[3];
rz(-0.27227835) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.17497002) q[2];
sx q[2];
rz(-1.9113767) q[2];
sx q[2];
rz(-0.54565412) q[2];
rz(-3.0555861) q[3];
sx q[3];
rz(-1.6274446) q[3];
sx q[3];
rz(0.29451323) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94588146) q[0];
sx q[0];
rz(-2.0989037) q[0];
sx q[0];
rz(1.209191) q[0];
rz(2.53881) q[1];
sx q[1];
rz(-0.61059112) q[1];
sx q[1];
rz(0.75377725) q[1];
rz(-2.9309737) q[2];
sx q[2];
rz(-0.89952614) q[2];
sx q[2];
rz(2.6256403) q[2];
rz(-0.62282563) q[3];
sx q[3];
rz(-1.6591743) q[3];
sx q[3];
rz(1.7521973) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
