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
rz(0.29851222) q[0];
sx q[0];
rz(-1.8520344) q[0];
sx q[0];
rz(0.02331743) q[0];
rz(0.98195568) q[1];
sx q[1];
rz(-0.73749956) q[1];
sx q[1];
rz(1.4758543) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6029089) q[0];
sx q[0];
rz(-1.6404466) q[0];
sx q[0];
rz(0.04361857) q[0];
rz(-pi) q[1];
rz(-1.2405922) q[2];
sx q[2];
rz(-1.6088952) q[2];
sx q[2];
rz(0.35534258) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.45881328) q[1];
sx q[1];
rz(-0.87535697) q[1];
sx q[1];
rz(-1.0473482) q[1];
rz(0.9528927) q[3];
sx q[3];
rz(-0.77761071) q[3];
sx q[3];
rz(0.67985204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4618571) q[2];
sx q[2];
rz(-1.4718141) q[2];
sx q[2];
rz(-0.30065817) q[2];
rz(0.65827185) q[3];
sx q[3];
rz(-2.7418147) q[3];
sx q[3];
rz(2.5532653) q[3];
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
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41475007) q[0];
sx q[0];
rz(-0.86597935) q[0];
sx q[0];
rz(2.9657189) q[0];
rz(-2.580592) q[1];
sx q[1];
rz(-0.70204061) q[1];
sx q[1];
rz(-2.9452513) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6628274) q[0];
sx q[0];
rz(-1.2365554) q[0];
sx q[0];
rz(-1.1007376) q[0];
x q[1];
rz(-0.77361892) q[2];
sx q[2];
rz(-0.16628312) q[2];
sx q[2];
rz(2.6321049) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4999428) q[1];
sx q[1];
rz(-1.5687935) q[1];
sx q[1];
rz(-1.0553805) q[1];
rz(1.7045872) q[3];
sx q[3];
rz(-1.9113298) q[3];
sx q[3];
rz(1.2029369) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.089513) q[2];
sx q[2];
rz(-2.5255346) q[2];
sx q[2];
rz(1.4698131) q[2];
rz(0.52516627) q[3];
sx q[3];
rz(-1.4273806) q[3];
sx q[3];
rz(-0.66450459) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1216275) q[0];
sx q[0];
rz(-2.2522734) q[0];
sx q[0];
rz(-0.63051939) q[0];
rz(-1.1446674) q[1];
sx q[1];
rz(-1.2455218) q[1];
sx q[1];
rz(-0.05680457) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81679427) q[0];
sx q[0];
rz(-1.1628808) q[0];
sx q[0];
rz(-0.40128951) q[0];
x q[1];
rz(2.2571083) q[2];
sx q[2];
rz(-1.9504042) q[2];
sx q[2];
rz(-1.1432709) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.9020128) q[1];
sx q[1];
rz(-0.48643349) q[1];
sx q[1];
rz(2.835266) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1278387) q[3];
sx q[3];
rz(-0.46681279) q[3];
sx q[3];
rz(1.6396322) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.28222617) q[2];
sx q[2];
rz(-1.6702009) q[2];
sx q[2];
rz(-0.43914208) q[2];
rz(-1.3965083) q[3];
sx q[3];
rz(-1.8789995) q[3];
sx q[3];
rz(0.5784353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0466995) q[0];
sx q[0];
rz(-2.4201604) q[0];
sx q[0];
rz(1.0149581) q[0];
rz(0.062648423) q[1];
sx q[1];
rz(-2.1868314) q[1];
sx q[1];
rz(-2.8573724) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13260176) q[0];
sx q[0];
rz(-1.1418482) q[0];
sx q[0];
rz(-1.1852077) q[0];
rz(-pi) q[1];
rz(2.6437628) q[2];
sx q[2];
rz(-1.4477028) q[2];
sx q[2];
rz(-3.0074938) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.758421) q[1];
sx q[1];
rz(-1.212442) q[1];
sx q[1];
rz(-1.5835632) q[1];
rz(-0.013945097) q[3];
sx q[3];
rz(-1.317465) q[3];
sx q[3];
rz(-1.0903181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.8404428) q[2];
sx q[2];
rz(-1.0089077) q[2];
sx q[2];
rz(0.76048771) q[2];
rz(-2.8216951) q[3];
sx q[3];
rz(-2.0127998) q[3];
sx q[3];
rz(1.0285876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8496721) q[0];
sx q[0];
rz(-0.55067647) q[0];
sx q[0];
rz(-2.736295) q[0];
rz(-2.6536476) q[1];
sx q[1];
rz(-2.0351724) q[1];
sx q[1];
rz(-0.36283666) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80984801) q[0];
sx q[0];
rz(-2.6001626) q[0];
sx q[0];
rz(-2.5141513) q[0];
x q[1];
rz(-1.1929197) q[2];
sx q[2];
rz(-2.123017) q[2];
sx q[2];
rz(0.044980031) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.2606973) q[1];
sx q[1];
rz(-2.915307) q[1];
sx q[1];
rz(-0.21842167) q[1];
rz(-pi) q[2];
x q[2];
rz(0.69546206) q[3];
sx q[3];
rz(-1.8084648) q[3];
sx q[3];
rz(-0.15791751) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.001658) q[2];
sx q[2];
rz(-1.9826865) q[2];
sx q[2];
rz(-2.4190767) q[2];
rz(-2.6288988) q[3];
sx q[3];
rz(-0.86944681) q[3];
sx q[3];
rz(1.2041913) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8743415) q[0];
sx q[0];
rz(-0.63684547) q[0];
sx q[0];
rz(-0.27477086) q[0];
rz(0.35708669) q[1];
sx q[1];
rz(-1.5029224) q[1];
sx q[1];
rz(1.9836609) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6053033) q[0];
sx q[0];
rz(-0.9035631) q[0];
sx q[0];
rz(-0.31179223) q[0];
rz(2.2311796) q[2];
sx q[2];
rz(-1.0779625) q[2];
sx q[2];
rz(3.1255258) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8573945) q[1];
sx q[1];
rz(-1.3669087) q[1];
sx q[1];
rz(-0.60867761) q[1];
rz(2.4729608) q[3];
sx q[3];
rz(-0.9210081) q[3];
sx q[3];
rz(1.1710407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.5623515) q[2];
sx q[2];
rz(-1.7513821) q[2];
sx q[2];
rz(-2.5518899) q[2];
rz(1.8288745) q[3];
sx q[3];
rz(-1.4785942) q[3];
sx q[3];
rz(-0.5635128) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8967459) q[0];
sx q[0];
rz(-2.1488996) q[0];
sx q[0];
rz(-2.8712811) q[0];
rz(0.42701328) q[1];
sx q[1];
rz(-1.3753563) q[1];
sx q[1];
rz(-2.7332773) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6231096) q[0];
sx q[0];
rz(-0.38869959) q[0];
sx q[0];
rz(-0.31950177) q[0];
rz(-pi) q[1];
rz(-2.2683581) q[2];
sx q[2];
rz(-1.6484652) q[2];
sx q[2];
rz(-2.6447062) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.282849) q[1];
sx q[1];
rz(-2.3925684) q[1];
sx q[1];
rz(2.2006459) q[1];
x q[2];
rz(0.116577) q[3];
sx q[3];
rz(-2.4156079) q[3];
sx q[3];
rz(-2.4621012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.2030877) q[2];
sx q[2];
rz(-1.0838584) q[2];
sx q[2];
rz(-0.92424029) q[2];
rz(-2.3067394) q[3];
sx q[3];
rz(-0.82930851) q[3];
sx q[3];
rz(0.99669641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95558178) q[0];
sx q[0];
rz(-0.38014933) q[0];
sx q[0];
rz(-0.36744776) q[0];
rz(0.5683178) q[1];
sx q[1];
rz(-1.808993) q[1];
sx q[1];
rz(0.41499358) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5552705) q[0];
sx q[0];
rz(-1.30997) q[0];
sx q[0];
rz(2.0404979) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.52936036) q[2];
sx q[2];
rz(-1.0870013) q[2];
sx q[2];
rz(2.8279378) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.4492682) q[1];
sx q[1];
rz(-1.0639166) q[1];
sx q[1];
rz(-0.70517069) q[1];
x q[2];
rz(-3.1363963) q[3];
sx q[3];
rz(-1.9165766) q[3];
sx q[3];
rz(2.7246662) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.35855287) q[2];
sx q[2];
rz(-2.2542605) q[2];
sx q[2];
rz(-2.3084194) q[2];
rz(0.35946515) q[3];
sx q[3];
rz(-2.0514252) q[3];
sx q[3];
rz(1.6370714) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1737162) q[0];
sx q[0];
rz(-2.739527) q[0];
sx q[0];
rz(3.1267082) q[0];
rz(1.124565) q[1];
sx q[1];
rz(-2.896307) q[1];
sx q[1];
rz(-0.10717779) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0343362) q[0];
sx q[0];
rz(-1.4217303) q[0];
sx q[0];
rz(-2.3019522) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.066266) q[2];
sx q[2];
rz(-2.3297254) q[2];
sx q[2];
rz(-0.56174413) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.31871492) q[1];
sx q[1];
rz(-1.4221622) q[1];
sx q[1];
rz(-2.6591402) q[1];
rz(-pi) q[2];
x q[2];
rz(0.59085544) q[3];
sx q[3];
rz(-0.31154943) q[3];
sx q[3];
rz(-2.8996866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.5643481) q[2];
sx q[2];
rz(-1.3078835) q[2];
sx q[2];
rz(-2.2819819) q[2];
rz(-2.882242) q[3];
sx q[3];
rz(-1.3602463) q[3];
sx q[3];
rz(0.41659659) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6421826) q[0];
sx q[0];
rz(-3.0152617) q[0];
sx q[0];
rz(1.1580178) q[0];
rz(0.49531373) q[1];
sx q[1];
rz(-1.1664349) q[1];
sx q[1];
rz(0.29702979) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3756145) q[0];
sx q[0];
rz(-0.31904623) q[0];
sx q[0];
rz(-1.5578397) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.530324) q[2];
sx q[2];
rz(-1.2707316) q[2];
sx q[2];
rz(-1.435868) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.7630939) q[1];
sx q[1];
rz(-1.8619672) q[1];
sx q[1];
rz(-1.4062455) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5498497) q[3];
sx q[3];
rz(-2.4789175) q[3];
sx q[3];
rz(1.7998723) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.872252) q[2];
sx q[2];
rz(-2.7597235) q[2];
sx q[2];
rz(2.2807109) q[2];
rz(-0.46368972) q[3];
sx q[3];
rz(-0.90353614) q[3];
sx q[3];
rz(1.1369107) q[3];
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
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9853482) q[0];
sx q[0];
rz(-1.5736268) q[0];
sx q[0];
rz(-1.9851984) q[0];
rz(-1.3337749) q[1];
sx q[1];
rz(-1.540624) q[1];
sx q[1];
rz(-2.435138) q[1];
rz(2.1352519) q[2];
sx q[2];
rz(-0.98802069) q[2];
sx q[2];
rz(1.5715656) q[2];
rz(-0.49336962) q[3];
sx q[3];
rz(-1.9153638) q[3];
sx q[3];
rz(2.3079688) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
