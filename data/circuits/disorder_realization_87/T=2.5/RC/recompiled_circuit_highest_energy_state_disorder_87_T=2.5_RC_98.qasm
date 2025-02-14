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
rz(2.8276853) q[0];
sx q[0];
rz(8.0765434) q[0];
rz(2.1836166) q[1];
sx q[1];
rz(4.557717) q[1];
sx q[1];
rz(9.582914) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8414171) q[0];
sx q[0];
rz(-0.99208562) q[0];
sx q[0];
rz(-2.7912223) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7633171) q[2];
sx q[2];
rz(-1.8774069) q[2];
sx q[2];
rz(3.0956097) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.68542105) q[1];
sx q[1];
rz(-3.1361533) q[1];
sx q[1];
rz(1.3855168) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4918658) q[3];
sx q[3];
rz(-2.5597245) q[3];
sx q[3];
rz(-0.054903809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.30369514) q[2];
sx q[2];
rz(-2.1968696) q[2];
sx q[2];
rz(-0.59291214) q[2];
rz(2.8494075) q[3];
sx q[3];
rz(-3.1226776) q[3];
sx q[3];
rz(-1.7570447) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.050046571) q[0];
sx q[0];
rz(-0.36588359) q[0];
sx q[0];
rz(-3.0475317) q[0];
rz(-1.7496109) q[1];
sx q[1];
rz(-1.530502) q[1];
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
rz(1.196435) q[2];
sx q[2];
rz(-0.5584467) q[2];
sx q[2];
rz(1.2701891) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.5793094) q[1];
sx q[1];
rz(-0.96889773) q[1];
sx q[1];
rz(3.137869) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9280065) q[3];
sx q[3];
rz(-0.41700577) q[3];
sx q[3];
rz(-0.97975376) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.2469108) q[2];
sx q[2];
rz(-2.6956788) q[2];
sx q[2];
rz(-1.1723588) q[2];
rz(-2.7123978) q[3];
sx q[3];
rz(-2.6515638) q[3];
sx q[3];
rz(1.711285) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
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
rz(-1.9046852) q[0];
sx q[0];
rz(-0.58085668) q[0];
sx q[0];
rz(-2.7823606) q[0];
rz(-1.5776186) q[1];
sx q[1];
rz(-2.3465395) q[1];
sx q[1];
rz(-0.94594812) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1079179) q[0];
sx q[0];
rz(-1.4888516) q[0];
sx q[0];
rz(-1.4180523) q[0];
rz(-pi) q[1];
rz(0.038905734) q[2];
sx q[2];
rz(-1.4565598) q[2];
sx q[2];
rz(-2.7047005) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9477152) q[1];
sx q[1];
rz(-1.7046728) q[1];
sx q[1];
rz(1.3280921) q[1];
rz(-pi) q[2];
rz(-1.5005689) q[3];
sx q[3];
rz(-1.4874793) q[3];
sx q[3];
rz(-0.30316363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.49950162) q[2];
sx q[2];
rz(-1.5311798) q[2];
sx q[2];
rz(1.1031319) q[2];
rz(-1.057386) q[3];
sx q[3];
rz(-1.5494346) q[3];
sx q[3];
rz(-0.21638432) q[3];
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
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0353521) q[0];
sx q[0];
rz(-3.048936) q[0];
sx q[0];
rz(-2.4308391) q[0];
rz(2.4757929) q[1];
sx q[1];
rz(-0.0087105287) q[1];
sx q[1];
rz(0.3054558) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75765002) q[0];
sx q[0];
rz(-0.096127495) q[0];
sx q[0];
rz(2.686475) q[0];
rz(-pi) q[1];
rz(-1.8886853) q[2];
sx q[2];
rz(-1.0895562) q[2];
sx q[2];
rz(-2.7900591) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7240746) q[1];
sx q[1];
rz(-1.6502079) q[1];
sx q[1];
rz(2.0168138) q[1];
x q[2];
rz(2.7025619) q[3];
sx q[3];
rz(-0.95618382) q[3];
sx q[3];
rz(0.87814349) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.010004369) q[2];
sx q[2];
rz(-1.5350716) q[2];
sx q[2];
rz(-0.32001495) q[2];
rz(2.6044676) q[3];
sx q[3];
rz(-0.38083005) q[3];
sx q[3];
rz(2.2297458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.491275) q[0];
sx q[0];
rz(-0.22458751) q[0];
sx q[0];
rz(-0.085163072) q[0];
rz(0.59146178) q[1];
sx q[1];
rz(-0.0031091212) q[1];
sx q[1];
rz(1.1726146) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64567287) q[0];
sx q[0];
rz(-1.5695736) q[0];
sx q[0];
rz(-1.5652577) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8109365) q[2];
sx q[2];
rz(-1.2954172) q[2];
sx q[2];
rz(-2.9346043) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.5796645) q[1];
sx q[1];
rz(-2.2854726) q[1];
sx q[1];
rz(-2.2382333) q[1];
x q[2];
rz(0.32043227) q[3];
sx q[3];
rz(-0.51051192) q[3];
sx q[3];
rz(1.7701469) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(3.0255967) q[2];
sx q[2];
rz(-1.7946578) q[2];
sx q[2];
rz(-1.4541413) q[2];
rz(1.2429552) q[3];
sx q[3];
rz(-0.60296139) q[3];
sx q[3];
rz(-0.81237826) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.94266194) q[0];
sx q[0];
rz(-2.9781065) q[0];
sx q[0];
rz(1.7401975) q[0];
rz(0.61680782) q[1];
sx q[1];
rz(-3.1251833) q[1];
sx q[1];
rz(1.1161463) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45081954) q[0];
sx q[0];
rz(-1.6661394) q[0];
sx q[0];
rz(-1.3006163) q[0];
rz(2.7737439) q[2];
sx q[2];
rz(-0.90613885) q[2];
sx q[2];
rz(2.9995415) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.9114218) q[1];
sx q[1];
rz(-1.0814922) q[1];
sx q[1];
rz(2.4858413) q[1];
x q[2];
rz(-0.44915389) q[3];
sx q[3];
rz(-1.7633668) q[3];
sx q[3];
rz(1.4447007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.8650032) q[2];
sx q[2];
rz(-1.9329376) q[2];
sx q[2];
rz(-1.2651944) q[2];
rz(0.64722925) q[3];
sx q[3];
rz(-2.6914458) q[3];
sx q[3];
rz(1.2559206) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7728421) q[0];
sx q[0];
rz(-0.54452801) q[0];
sx q[0];
rz(-0.37753373) q[0];
rz(2.990621) q[1];
sx q[1];
rz(-0.010846373) q[1];
sx q[1];
rz(-0.46802256) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0819055) q[0];
sx q[0];
rz(-1.1199315) q[0];
sx q[0];
rz(-1.5968869) q[0];
rz(0.84845869) q[2];
sx q[2];
rz(-1.6654789) q[2];
sx q[2];
rz(2.2168558) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1605058) q[1];
sx q[1];
rz(-1.5001343) q[1];
sx q[1];
rz(2.2494506) q[1];
rz(0.19429732) q[3];
sx q[3];
rz(-2.0812985) q[3];
sx q[3];
rz(-1.0620775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.2630792) q[2];
sx q[2];
rz(-0.24363467) q[2];
sx q[2];
rz(2.0344951) q[2];
rz(2.257972) q[3];
sx q[3];
rz(-0.7928018) q[3];
sx q[3];
rz(0.71225524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9368619) q[0];
sx q[0];
rz(-0.87918133) q[0];
sx q[0];
rz(-2.9469446) q[0];
rz(-1.4600935) q[1];
sx q[1];
rz(-0.016540557) q[1];
sx q[1];
rz(-2.8406692) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0452354) q[0];
sx q[0];
rz(-0.83556743) q[0];
sx q[0];
rz(2.623933) q[0];
x q[1];
rz(-2.373502) q[2];
sx q[2];
rz(-2.4642307) q[2];
sx q[2];
rz(-1.0528262) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.49773945) q[1];
sx q[1];
rz(-0.56922888) q[1];
sx q[1];
rz(0.42909547) q[1];
x q[2];
rz(-0.31105403) q[3];
sx q[3];
rz(-2.9533953) q[3];
sx q[3];
rz(1.8194906) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.68227565) q[2];
sx q[2];
rz(-1.6271017) q[2];
sx q[2];
rz(-2.922399) q[2];
rz(1.7949665) q[3];
sx q[3];
rz(-3.0248088) q[3];
sx q[3];
rz(-0.78484261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.858736) q[0];
sx q[0];
rz(-0.63010001) q[0];
sx q[0];
rz(-3.1334738) q[0];
rz(-2.4722664) q[1];
sx q[1];
rz(-3.1287584) q[1];
sx q[1];
rz(0.76378167) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8618906) q[0];
sx q[0];
rz(-2.6548626) q[0];
sx q[0];
rz(-2.6175015) q[0];
rz(-pi) q[1];
rz(-1.6449459) q[2];
sx q[2];
rz(-0.96536812) q[2];
sx q[2];
rz(-2.2260967) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.1642832) q[1];
sx q[1];
rz(-0.54879511) q[1];
sx q[1];
rz(-1.4700233) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3108008) q[3];
sx q[3];
rz(-2.0767815) q[3];
sx q[3];
rz(0.80238402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.0236437) q[2];
sx q[2];
rz(-3.1274319) q[2];
sx q[2];
rz(-0.99948779) q[2];
rz(-0.1855447) q[3];
sx q[3];
rz(-2.1047635) q[3];
sx q[3];
rz(3.1256092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.213181) q[0];
sx q[0];
rz(-0.11225926) q[0];
sx q[0];
rz(2.4768594) q[0];
rz(0.96066535) q[1];
sx q[1];
rz(-0.04048368) q[1];
sx q[1];
rz(1.6368846) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6204313) q[0];
sx q[0];
rz(-1.3913432) q[0];
sx q[0];
rz(1.8779264) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.78053) q[2];
sx q[2];
rz(-1.6170986) q[2];
sx q[2];
rz(-1.8425892) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.9686925) q[1];
sx q[1];
rz(-2.4174999) q[1];
sx q[1];
rz(-2.398046) q[1];
rz(-pi) q[2];
x q[2];
rz(0.89791132) q[3];
sx q[3];
rz(-0.6403044) q[3];
sx q[3];
rz(1.698157) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.8399743) q[2];
sx q[2];
rz(-3.1349389) q[2];
sx q[2];
rz(-1.9807695) q[2];
rz(0.079856722) q[3];
sx q[3];
rz(-3.1294332) q[3];
sx q[3];
rz(-1.2016092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6066211) q[0];
sx q[0];
rz(-1.6666245) q[0];
sx q[0];
rz(1.5635906) q[0];
rz(-3.1040991) q[1];
sx q[1];
rz(-0.19352023) q[1];
sx q[1];
rz(0.22307693) q[1];
rz(1.1983331) q[2];
sx q[2];
rz(-1.4683152) q[2];
sx q[2];
rz(-2.7905133) q[2];
rz(0.014122176) q[3];
sx q[3];
rz(-2.2416708) q[3];
sx q[3];
rz(-1.3308217) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
