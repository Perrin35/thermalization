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
rz(4.557717) q[1];
sx q[1];
rz(9.582914) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0682521) q[0];
sx q[0];
rz(-1.86226) q[0];
sx q[0];
rz(-2.1785546) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.70868196) q[2];
sx q[2];
rz(-2.6593609) q[2];
sx q[2];
rz(-0.9672375) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.50013881) q[1];
sx q[1];
rz(-1.5761426) q[1];
sx q[1];
rz(0.0010020573) q[1];
x q[2];
rz(1.9495548) q[3];
sx q[3];
rz(-1.1178607) q[3];
sx q[3];
rz(0.68309802) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.8378975) q[2];
sx q[2];
rz(-0.94472307) q[2];
sx q[2];
rz(0.59291214) q[2];
rz(2.8494075) q[3];
sx q[3];
rz(-3.1226776) q[3];
sx q[3];
rz(1.384548) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0915461) q[0];
sx q[0];
rz(-2.7757091) q[0];
sx q[0];
rz(-0.094060913) q[0];
rz(1.3919818) q[1];
sx q[1];
rz(-1.530502) q[1];
sx q[1];
rz(1.7382517) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3791524) q[0];
sx q[0];
rz(-1.3032252) q[0];
sx q[0];
rz(-1.3570157) q[0];
rz(-1.9451577) q[2];
sx q[2];
rz(-2.583146) q[2];
sx q[2];
rz(-1.2701891) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5793094) q[1];
sx q[1];
rz(-0.96889773) q[1];
sx q[1];
rz(0.0037236969) q[1];
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
x q[1];
rz(2.2469108) q[2];
sx q[2];
rz(-2.6956788) q[2];
sx q[2];
rz(-1.9692339) q[2];
rz(-2.7123978) q[3];
sx q[3];
rz(-0.49002886) q[3];
sx q[3];
rz(-1.711285) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9046852) q[0];
sx q[0];
rz(-0.58085668) q[0];
sx q[0];
rz(-0.3592321) q[0];
rz(-1.563974) q[1];
sx q[1];
rz(-0.79505316) q[1];
sx q[1];
rz(2.1956445) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6170711) q[0];
sx q[0];
rz(-1.4185689) q[0];
sx q[0];
rz(0.082905654) q[0];
rz(3.1026869) q[2];
sx q[2];
rz(-1.4565598) q[2];
sx q[2];
rz(2.7047005) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.1938775) q[1];
sx q[1];
rz(-1.4369198) q[1];
sx q[1];
rz(1.3280921) q[1];
rz(-0.083521997) q[3];
sx q[3];
rz(-1.5008129) q[3];
sx q[3];
rz(1.2617788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.642091) q[2];
sx q[2];
rz(-1.6104128) q[2];
sx q[2];
rz(-2.0384608) q[2];
rz(-1.057386) q[3];
sx q[3];
rz(-1.5494346) q[3];
sx q[3];
rz(-0.21638432) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0353521) q[0];
sx q[0];
rz(-3.048936) q[0];
sx q[0];
rz(-0.71075359) q[0];
rz(2.4757929) q[1];
sx q[1];
rz(-3.1328821) q[1];
sx q[1];
rz(-0.3054558) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9269954) q[0];
sx q[0];
rz(-1.4844795) q[0];
sx q[0];
rz(1.6131562) q[0];
x q[1];
rz(-2.6389559) q[2];
sx q[2];
rz(-1.2900616) q[2];
sx q[2];
rz(-1.0681149) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.1911911) q[1];
sx q[1];
rz(-2.0153075) q[1];
sx q[1];
rz(3.0536122) q[1];
x q[2];
rz(-0.908522) q[3];
sx q[3];
rz(-1.2161331) q[3];
sx q[3];
rz(-0.42820546) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(3.1315883) q[2];
sx q[2];
rz(-1.5350716) q[2];
sx q[2];
rz(2.8215777) q[2];
rz(-0.53712505) q[3];
sx q[3];
rz(-2.7607626) q[3];
sx q[3];
rz(0.91184688) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.491275) q[0];
sx q[0];
rz(-0.22458751) q[0];
sx q[0];
rz(0.085163072) q[0];
rz(-0.59146178) q[1];
sx q[1];
rz(-3.1384835) q[1];
sx q[1];
rz(-1.9689781) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92511668) q[0];
sx q[0];
rz(-1.5763349) q[0];
sx q[0];
rz(3.14037) q[0];
rz(-0.28309254) q[2];
sx q[2];
rz(-1.3398719) q[2];
sx q[2];
rz(-1.4302916) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7025801) q[1];
sx q[1];
rz(-2.2058553) q[1];
sx q[1];
rz(2.5219581) q[1];
rz(-pi) q[2];
rz(-1.745397) q[3];
sx q[3];
rz(-2.0530226) q[3];
sx q[3];
rz(1.7349145) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(3.0255967) q[2];
sx q[2];
rz(-1.7946578) q[2];
sx q[2];
rz(-1.4541413) q[2];
rz(1.8986374) q[3];
sx q[3];
rz(-0.60296139) q[3];
sx q[3];
rz(-2.3292144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94266194) q[0];
sx q[0];
rz(-0.16348612) q[0];
sx q[0];
rz(1.7401975) q[0];
rz(-2.5247848) q[1];
sx q[1];
rz(-0.016409358) q[1];
sx q[1];
rz(2.0254463) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0936174) q[0];
sx q[0];
rz(-1.3018739) q[0];
sx q[0];
rz(3.0426836) q[0];
rz(-2.2693231) q[2];
sx q[2];
rz(-1.8577777) q[2];
sx q[2];
rz(1.6621188) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.20780288) q[1];
sx q[1];
rz(-2.3456165) q[1];
sx q[1];
rz(2.4237304) q[1];
rz(-pi) q[2];
rz(-2.7195144) q[3];
sx q[3];
rz(-2.6554972) q[3];
sx q[3];
rz(2.8897282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8650032) q[2];
sx q[2];
rz(-1.208655) q[2];
sx q[2];
rz(-1.8763982) q[2];
rz(-0.64722925) q[3];
sx q[3];
rz(-0.45014683) q[3];
sx q[3];
rz(1.2559206) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
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
rz(0.3687506) q[0];
sx q[0];
rz(-2.5970646) q[0];
sx q[0];
rz(-0.37753373) q[0];
rz(0.15097161) q[1];
sx q[1];
rz(-0.010846373) q[1];
sx q[1];
rz(-2.6735701) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6418544) q[0];
sx q[0];
rz(-1.5942792) q[0];
sx q[0];
rz(0.45099839) q[0];
x q[1];
rz(0.12590825) q[2];
sx q[2];
rz(-2.2891938) q[2];
sx q[2];
rz(-2.5786521) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.9810869) q[1];
sx q[1];
rz(-1.5001343) q[1];
sx q[1];
rz(2.2494506) q[1];
rz(-pi) q[2];
rz(-1.9028038) q[3];
sx q[3];
rz(-2.598437) q[3];
sx q[3];
rz(2.4623722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.2630792) q[2];
sx q[2];
rz(-2.897958) q[2];
sx q[2];
rz(1.1070975) q[2];
rz(2.257972) q[3];
sx q[3];
rz(-2.3487909) q[3];
sx q[3];
rz(-0.71225524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
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
rz(2.9368619) q[0];
sx q[0];
rz(-2.2624113) q[0];
sx q[0];
rz(0.19464807) q[0];
rz(-1.6814992) q[1];
sx q[1];
rz(-3.1250521) q[1];
sx q[1];
rz(-2.8406692) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3022879) q[0];
sx q[0];
rz(-1.1949958) q[0];
sx q[0];
rz(2.3761533) q[0];
rz(2.6171404) q[2];
sx q[2];
rz(-2.0213184) q[2];
sx q[2];
rz(-1.1632048) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.7006939) q[1];
sx q[1];
rz(-1.3446302) q[1];
sx q[1];
rz(-2.6146099) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8305386) q[3];
sx q[3];
rz(-0.18819735) q[3];
sx q[3];
rz(-1.8194906) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.459317) q[2];
sx q[2];
rz(-1.6271017) q[2];
sx q[2];
rz(0.21919361) q[2];
rz(-1.7949665) q[3];
sx q[3];
rz(-3.0248088) q[3];
sx q[3];
rz(0.78484261) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.858736) q[0];
sx q[0];
rz(-2.5114926) q[0];
sx q[0];
rz(-0.0081188763) q[0];
rz(0.66932622) q[1];
sx q[1];
rz(-3.1287584) q[1];
sx q[1];
rz(0.76378167) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29944659) q[0];
sx q[0];
rz(-1.9877292) q[0];
sx q[0];
rz(-1.3119158) q[0];
x q[1];
rz(-0.10662756) q[2];
sx q[2];
rz(-2.5322057) q[2];
sx q[2];
rz(-2.0963017) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.4620518) q[1];
sx q[1];
rz(-1.6233007) q[1];
sx q[1];
rz(2.1173304) q[1];
rz(-pi) q[2];
rz(0.88293255) q[3];
sx q[3];
rz(-2.2730423) q[3];
sx q[3];
rz(-1.8851999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.0236437) q[2];
sx q[2];
rz(-3.1274319) q[2];
sx q[2];
rz(2.1421049) q[2];
rz(-2.956048) q[3];
sx q[3];
rz(-2.1047635) q[3];
sx q[3];
rz(0.015983494) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.213181) q[0];
sx q[0];
rz(-3.0293334) q[0];
sx q[0];
rz(-0.66473329) q[0];
rz(2.1809273) q[1];
sx q[1];
rz(-3.101109) q[1];
sx q[1];
rz(-1.5047081) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5793631) q[0];
sx q[0];
rz(-0.35427552) q[0];
sx q[0];
rz(-1.0303524) q[0];
rz(-pi) q[1];
rz(0.36106266) q[2];
sx q[2];
rz(-1.5244941) q[2];
sx q[2];
rz(1.8425892) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.3469228) q[1];
sx q[1];
rz(-1.1058013) q[1];
sx q[1];
rz(-0.57703206) q[1];
x q[2];
rz(0.43470498) q[3];
sx q[3];
rz(-1.0846595) q[3];
sx q[3];
rz(2.4804573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.30161834) q[2];
sx q[2];
rz(-0.0066537298) q[2];
sx q[2];
rz(1.1608231) q[2];
rz(-3.0617359) q[3];
sx q[3];
rz(-0.012159497) q[3];
sx q[3];
rz(1.2016092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
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
rz(2.6066211) q[0];
sx q[0];
rz(-1.6666245) q[0];
sx q[0];
rz(1.5635906) q[0];
rz(-0.037493575) q[1];
sx q[1];
rz(-2.9480724) q[1];
sx q[1];
rz(-2.9185157) q[1];
rz(-3.0316261) q[2];
sx q[2];
rz(-1.2003821) q[2];
sx q[2];
rz(-1.1797689) q[2];
rz(-2.2417193) q[3];
sx q[3];
rz(-1.5818578) q[3];
sx q[3];
rz(0.24875438) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
