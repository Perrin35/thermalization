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
rz(-0.63646746) q[0];
sx q[0];
rz(-0.31390733) q[0];
sx q[0];
rz(-1.7933581) q[0];
rz(2.1836166) q[1];
sx q[1];
rz(-1.7254683) q[1];
sx q[1];
rz(0.15813601) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30017553) q[0];
sx q[0];
rz(-2.149507) q[0];
sx q[0];
rz(-2.7912223) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8991437) q[2];
sx q[2];
rz(-1.930604) q[2];
sx q[2];
rz(1.4054325) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.6414538) q[1];
sx q[1];
rz(-1.5761426) q[1];
sx q[1];
rz(-3.1405906) q[1];
rz(-pi) q[2];
rz(-1.9495548) q[3];
sx q[3];
rz(-1.1178607) q[3];
sx q[3];
rz(2.4584946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.30369514) q[2];
sx q[2];
rz(-2.1968696) q[2];
sx q[2];
rz(2.5486805) q[2];
rz(-2.8494075) q[3];
sx q[3];
rz(-0.018915011) q[3];
sx q[3];
rz(-1.7570447) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0915461) q[0];
sx q[0];
rz(-2.7757091) q[0];
sx q[0];
rz(0.094060913) q[0];
rz(1.7496109) q[1];
sx q[1];
rz(-1.6110907) q[1];
sx q[1];
rz(-1.403341) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3791524) q[0];
sx q[0];
rz(-1.3032252) q[0];
sx q[0];
rz(1.3570157) q[0];
rz(-pi) q[1];
rz(0.22461598) q[2];
sx q[2];
rz(-1.0550647) q[2];
sx q[2];
rz(1.7040229) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.010621444) q[1];
sx q[1];
rz(-1.567727) q[1];
sx q[1];
rz(-0.9688945) q[1];
rz(0.4085598) q[3];
sx q[3];
rz(-1.6567536) q[3];
sx q[3];
rz(-2.3547821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.2469108) q[2];
sx q[2];
rz(-2.6956788) q[2];
sx q[2];
rz(-1.9692339) q[2];
rz(0.42919484) q[3];
sx q[3];
rz(-2.6515638) q[3];
sx q[3];
rz(-1.4303077) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9046852) q[0];
sx q[0];
rz(-0.58085668) q[0];
sx q[0];
rz(0.3592321) q[0];
rz(-1.563974) q[1];
sx q[1];
rz(-0.79505316) q[1];
sx q[1];
rz(2.1956445) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1157784) q[0];
sx q[0];
rz(-2.9684068) q[0];
sx q[0];
rz(-2.065763) q[0];
x q[1];
rz(-1.456474) q[2];
sx q[2];
rz(-1.5321443) q[2];
sx q[2];
rz(-2.0032515) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-3.0240629) q[1];
sx q[1];
rz(-0.27654031) q[1];
sx q[1];
rz(2.0815954) q[1];
x q[2];
rz(0.69880693) q[3];
sx q[3];
rz(-3.0326789) q[3];
sx q[3];
rz(2.1366936) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.49950162) q[2];
sx q[2];
rz(-1.5311798) q[2];
sx q[2];
rz(-2.0384608) q[2];
rz(2.0842066) q[3];
sx q[3];
rz(-1.5494346) q[3];
sx q[3];
rz(-0.21638432) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10624056) q[0];
sx q[0];
rz(-0.092656605) q[0];
sx q[0];
rz(0.71075359) q[0];
rz(-2.4757929) q[1];
sx q[1];
rz(-3.1328821) q[1];
sx q[1];
rz(-2.8361368) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35985306) q[0];
sx q[0];
rz(-1.6129984) q[0];
sx q[0];
rz(-3.0551987) q[0];
x q[1];
rz(-1.2529073) q[2];
sx q[2];
rz(-1.0895562) q[2];
sx q[2];
rz(2.7900591) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.98887029) q[1];
sx q[1];
rz(-0.45256361) q[1];
sx q[1];
rz(1.3883703) q[1];
x q[2];
rz(-2.2330707) q[3];
sx q[3];
rz(-1.9254596) q[3];
sx q[3];
rz(2.7133872) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.010004369) q[2];
sx q[2];
rz(-1.606521) q[2];
sx q[2];
rz(-0.32001495) q[2];
rz(2.6044676) q[3];
sx q[3];
rz(-2.7607626) q[3];
sx q[3];
rz(0.91184688) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.491275) q[0];
sx q[0];
rz(-2.9170051) q[0];
sx q[0];
rz(3.0564296) q[0];
rz(-2.5501309) q[1];
sx q[1];
rz(-3.1384835) q[1];
sx q[1];
rz(1.9689781) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92511668) q[0];
sx q[0];
rz(-1.5652577) q[0];
sx q[0];
rz(0.001222697) q[0];
rz(0.69967592) q[2];
sx q[2];
rz(-2.778233) q[2];
sx q[2];
rz(2.6154721) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7025801) q[1];
sx q[1];
rz(-0.93573739) q[1];
sx q[1];
rz(-0.61963453) q[1];
x q[2];
rz(2.8211604) q[3];
sx q[3];
rz(-0.51051192) q[3];
sx q[3];
rz(-1.7701469) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.11599596) q[2];
sx q[2];
rz(-1.7946578) q[2];
sx q[2];
rz(1.4541413) q[2];
rz(1.8986374) q[3];
sx q[3];
rz(-2.5386313) q[3];
sx q[3];
rz(2.3292144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1989307) q[0];
sx q[0];
rz(-2.9781065) q[0];
sx q[0];
rz(-1.7401975) q[0];
rz(-2.5247848) q[1];
sx q[1];
rz(-0.016409358) q[1];
sx q[1];
rz(-1.1161463) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45081954) q[0];
sx q[0];
rz(-1.4754533) q[0];
sx q[0];
rz(1.3006163) q[0];
rz(-0.3678488) q[2];
sx q[2];
rz(-0.90613885) q[2];
sx q[2];
rz(-0.14205113) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.20780288) q[1];
sx q[1];
rz(-0.79597616) q[1];
sx q[1];
rz(0.7178623) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7839637) q[3];
sx q[3];
rz(-2.0110594) q[3];
sx q[3];
rz(-2.9235087) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.8650032) q[2];
sx q[2];
rz(-1.9329376) q[2];
sx q[2];
rz(1.2651944) q[2];
rz(-0.64722925) q[3];
sx q[3];
rz(-2.6914458) q[3];
sx q[3];
rz(1.885672) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7728421) q[0];
sx q[0];
rz(-2.5970646) q[0];
sx q[0];
rz(-2.7640589) q[0];
rz(-2.990621) q[1];
sx q[1];
rz(-3.1307463) q[1];
sx q[1];
rz(-0.46802256) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.05968719) q[0];
sx q[0];
rz(-2.0216612) q[0];
sx q[0];
rz(1.5447058) q[0];
x q[1];
rz(0.12590825) q[2];
sx q[2];
rz(-0.85239886) q[2];
sx q[2];
rz(-0.56294051) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.6439754) q[1];
sx q[1];
rz(-2.4598498) q[1];
sx q[1];
rz(-1.6830744) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9472953) q[3];
sx q[3];
rz(-2.0812985) q[3];
sx q[3];
rz(2.0795151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.2630792) q[2];
sx q[2];
rz(-2.897958) q[2];
sx q[2];
rz(-2.0344951) q[2];
rz(2.257972) q[3];
sx q[3];
rz(-0.7928018) q[3];
sx q[3];
rz(0.71225524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20473075) q[0];
sx q[0];
rz(-0.87918133) q[0];
sx q[0];
rz(-2.9469446) q[0];
rz(-1.4600935) q[1];
sx q[1];
rz(-0.016540557) q[1];
sx q[1];
rz(0.3009235) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3022879) q[0];
sx q[0];
rz(-1.1949958) q[0];
sx q[0];
rz(-0.76543937) q[0];
rz(-0.5244523) q[2];
sx q[2];
rz(-1.1202742) q[2];
sx q[2];
rz(-1.9783879) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.1414425) q[1];
sx q[1];
rz(-1.0585691) q[1];
sx q[1];
rz(1.8309831) q[1];
rz(0.17936318) q[3];
sx q[3];
rz(-1.6280884) q[3];
sx q[3];
rz(-3.0843861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.459317) q[2];
sx q[2];
rz(-1.6271017) q[2];
sx q[2];
rz(-2.922399) q[2];
rz(-1.7949665) q[3];
sx q[3];
rz(-3.0248088) q[3];
sx q[3];
rz(0.78484261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-2.858736) q[0];
sx q[0];
rz(-0.63010001) q[0];
sx q[0];
rz(3.1334738) q[0];
rz(2.4722664) q[1];
sx q[1];
rz(-0.012834276) q[1];
sx q[1];
rz(-2.377811) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29944659) q[0];
sx q[0];
rz(-1.9877292) q[0];
sx q[0];
rz(-1.3119158) q[0];
rz(-0.10662756) q[2];
sx q[2];
rz(-0.60938696) q[2];
sx q[2];
rz(-1.0452909) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.67954087) q[1];
sx q[1];
rz(-1.5182919) q[1];
sx q[1];
rz(-2.1173304) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4978746) q[3];
sx q[3];
rz(-2.2016085) q[3];
sx q[3];
rz(0.35177059) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.0236437) q[2];
sx q[2];
rz(-0.014160784) q[2];
sx q[2];
rz(-2.1421049) q[2];
rz(0.1855447) q[3];
sx q[3];
rz(-2.1047635) q[3];
sx q[3];
rz(0.015983494) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.213181) q[0];
sx q[0];
rz(-3.0293334) q[0];
sx q[0];
rz(0.66473329) q[0];
rz(0.96066535) q[1];
sx q[1];
rz(-3.101109) q[1];
sx q[1];
rz(1.5047081) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5211613) q[0];
sx q[0];
rz(-1.3913432) q[0];
sx q[0];
rz(1.8779264) q[0];
rz(-1.6202848) q[2];
sx q[2];
rz(-1.9314543) q[2];
sx q[2];
rz(2.8872763) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9686925) q[1];
sx q[1];
rz(-0.72409276) q[1];
sx q[1];
rz(2.398046) q[1];
rz(-1.0432518) q[3];
sx q[3];
rz(-1.9523373) q[3];
sx q[3];
rz(-2.4455618) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.30161834) q[2];
sx q[2];
rz(-0.0066537298) q[2];
sx q[2];
rz(-1.9807695) q[2];
rz(0.079856722) q[3];
sx q[3];
rz(-3.1294332) q[3];
sx q[3];
rz(1.9399835) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53497159) q[0];
sx q[0];
rz(-1.6666245) q[0];
sx q[0];
rz(1.5635906) q[0];
rz(0.037493575) q[1];
sx q[1];
rz(-0.19352023) q[1];
sx q[1];
rz(0.22307693) q[1];
rz(1.9432595) q[2];
sx q[2];
rz(-1.6732775) q[2];
sx q[2];
rz(0.35107935) q[2];
rz(-1.5885872) q[3];
sx q[3];
rz(-2.4705926) q[3];
sx q[3];
rz(1.8334851) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
