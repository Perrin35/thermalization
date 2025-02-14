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
rz(-0.95797602) q[1];
sx q[1];
rz(-1.4161243) q[1];
sx q[1];
rz(2.9834566) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0733405) q[0];
sx q[0];
rz(-1.86226) q[0];
sx q[0];
rz(2.1785546) q[0];
x q[1];
rz(0.70868196) q[2];
sx q[2];
rz(-0.48223178) q[2];
sx q[2];
rz(-0.9672375) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.68542105) q[1];
sx q[1];
rz(-0.0054393689) q[1];
sx q[1];
rz(-1.7560759) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6590822) q[3];
sx q[3];
rz(-1.231862) q[3];
sx q[3];
rz(-1.0601251) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.30369514) q[2];
sx q[2];
rz(-0.94472307) q[2];
sx q[2];
rz(-2.5486805) q[2];
rz(2.8494075) q[3];
sx q[3];
rz(-0.018915011) q[3];
sx q[3];
rz(-1.384548) q[3];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.050046571) q[0];
sx q[0];
rz(-0.36588359) q[0];
sx q[0];
rz(0.094060913) q[0];
rz(1.7496109) q[1];
sx q[1];
rz(-1.530502) q[1];
sx q[1];
rz(-1.7382517) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8926131) q[0];
sx q[0];
rz(-1.7768581) q[0];
sx q[0];
rz(0.27350105) q[0];
rz(1.044079) q[2];
sx q[2];
rz(-1.3758012) q[2];
sx q[2];
rz(-3.1205683) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.5622833) q[1];
sx q[1];
rz(-0.96889773) q[1];
sx q[1];
rz(-3.137869) q[1];
rz(-0.4085598) q[3];
sx q[3];
rz(-1.6567536) q[3];
sx q[3];
rz(2.3547821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.89468181) q[2];
sx q[2];
rz(-2.6956788) q[2];
sx q[2];
rz(1.9692339) q[2];
rz(2.7123978) q[3];
sx q[3];
rz(-0.49002886) q[3];
sx q[3];
rz(-1.4303077) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2369075) q[0];
sx q[0];
rz(-2.560736) q[0];
sx q[0];
rz(-0.3592321) q[0];
rz(1.563974) q[1];
sx q[1];
rz(-0.79505316) q[1];
sx q[1];
rz(-2.1956445) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52452156) q[0];
sx q[0];
rz(-1.7230238) q[0];
sx q[0];
rz(0.082905654) q[0];
rz(-pi) q[1];
x q[1];
rz(0.038905734) q[2];
sx q[2];
rz(-1.4565598) q[2];
sx q[2];
rz(0.43689219) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.11752973) q[1];
sx q[1];
rz(-2.8650523) q[1];
sx q[1];
rz(1.0599972) q[1];
x q[2];
rz(-0.69880693) q[3];
sx q[3];
rz(-0.10891373) q[3];
sx q[3];
rz(-1.004899) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.49950162) q[2];
sx q[2];
rz(-1.6104128) q[2];
sx q[2];
rz(-2.0384608) q[2];
rz(2.0842066) q[3];
sx q[3];
rz(-1.5494346) q[3];
sx q[3];
rz(2.9252083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10624056) q[0];
sx q[0];
rz(-0.092656605) q[0];
sx q[0];
rz(2.4308391) q[0];
rz(0.6657998) q[1];
sx q[1];
rz(-3.1328821) q[1];
sx q[1];
rz(0.3054558) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9269954) q[0];
sx q[0];
rz(-1.6571132) q[0];
sx q[0];
rz(1.6131562) q[0];
x q[1];
rz(0.50263672) q[2];
sx q[2];
rz(-1.2900616) q[2];
sx q[2];
rz(-1.0681149) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.1911911) q[1];
sx q[1];
rz(-2.0153075) q[1];
sx q[1];
rz(-3.0536122) q[1];
rz(-pi) q[2];
rz(2.2330707) q[3];
sx q[3];
rz(-1.9254596) q[3];
sx q[3];
rz(-2.7133872) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.1315883) q[2];
sx q[2];
rz(-1.5350716) q[2];
sx q[2];
rz(2.8215777) q[2];
rz(0.53712505) q[3];
sx q[3];
rz(-2.7607626) q[3];
sx q[3];
rz(-0.91184688) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.491275) q[0];
sx q[0];
rz(-2.9170051) q[0];
sx q[0];
rz(-0.085163072) q[0];
rz(0.59146178) q[1];
sx q[1];
rz(-3.1384835) q[1];
sx q[1];
rz(-1.1726146) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4959198) q[0];
sx q[0];
rz(-1.572019) q[0];
sx q[0];
rz(1.576335) q[0];
rz(-pi) q[1];
rz(1.8109365) q[2];
sx q[2];
rz(-1.2954172) q[2];
sx q[2];
rz(-2.9346043) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.53211305) q[1];
sx q[1];
rz(-1.0843186) q[1];
sx q[1];
rz(-0.83513135) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6530442) q[3];
sx q[3];
rz(-1.7253157) q[3];
sx q[3];
rz(0.082495436) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.11599596) q[2];
sx q[2];
rz(-1.7946578) q[2];
sx q[2];
rz(-1.6874514) q[2];
rz(1.8986374) q[3];
sx q[3];
rz(-0.60296139) q[3];
sx q[3];
rz(-2.3292144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
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
rz(0.94266194) q[0];
sx q[0];
rz(-0.16348612) q[0];
sx q[0];
rz(1.7401975) q[0];
rz(-2.5247848) q[1];
sx q[1];
rz(-0.016409358) q[1];
sx q[1];
rz(-1.1161463) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4510613) q[0];
sx q[0];
rz(-2.8554718) q[0];
sx q[0];
rz(1.2267434) q[0];
rz(-0.87226954) q[2];
sx q[2];
rz(-1.8577777) q[2];
sx q[2];
rz(1.4794738) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.9114218) q[1];
sx q[1];
rz(-2.0601005) q[1];
sx q[1];
rz(-2.4858413) q[1];
x q[2];
rz(0.42207828) q[3];
sx q[3];
rz(-2.6554972) q[3];
sx q[3];
rz(-0.25186447) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.8650032) q[2];
sx q[2];
rz(-1.208655) q[2];
sx q[2];
rz(1.2651944) q[2];
rz(-0.64722925) q[3];
sx q[3];
rz(-2.6914458) q[3];
sx q[3];
rz(-1.2559206) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3687506) q[0];
sx q[0];
rz(-0.54452801) q[0];
sx q[0];
rz(0.37753373) q[0];
rz(-0.15097161) q[1];
sx q[1];
rz(-0.010846373) q[1];
sx q[1];
rz(2.6735701) q[1];
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
rz(-1.5968869) q[0];
rz(-pi) q[1];
rz(-1.4281316) q[2];
sx q[2];
rz(-2.4141867) q[2];
sx q[2];
rz(-0.75293604) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.7881778) q[1];
sx q[1];
rz(-0.89415293) q[1];
sx q[1];
rz(0.090678399) q[1];
rz(-pi) q[2];
rz(1.2387889) q[3];
sx q[3];
rz(-2.598437) q[3];
sx q[3];
rz(2.4623722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.2630792) q[2];
sx q[2];
rz(-2.897958) q[2];
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
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20473075) q[0];
sx q[0];
rz(-0.87918133) q[0];
sx q[0];
rz(2.9469446) q[0];
rz(-1.6814992) q[1];
sx q[1];
rz(-0.016540557) q[1];
sx q[1];
rz(2.8406692) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3022879) q[0];
sx q[0];
rz(-1.9465969) q[0];
sx q[0];
rz(0.76543937) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.5244523) q[2];
sx q[2];
rz(-2.0213184) q[2];
sx q[2];
rz(1.9783879) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.6438532) q[1];
sx q[1];
rz(-0.56922888) q[1];
sx q[1];
rz(-2.7124972) q[1];
x q[2];
rz(2.8305386) q[3];
sx q[3];
rz(-0.18819735) q[3];
sx q[3];
rz(1.3221021) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.68227565) q[2];
sx q[2];
rz(-1.514491) q[2];
sx q[2];
rz(2.922399) q[2];
rz(-1.7949665) q[3];
sx q[3];
rz(-0.11678385) q[3];
sx q[3];
rz(2.35675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.858736) q[0];
sx q[0];
rz(-2.5114926) q[0];
sx q[0];
rz(-3.1334738) q[0];
rz(-0.66932622) q[1];
sx q[1];
rz(-0.012834276) q[1];
sx q[1];
rz(0.76378167) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8618906) q[0];
sx q[0];
rz(-0.48673004) q[0];
sx q[0];
rz(-0.52409117) q[0];
rz(-pi) q[1];
rz(2.5348762) q[2];
sx q[2];
rz(-1.5098443) q[2];
sx q[2];
rz(-2.5285442) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4620518) q[1];
sx q[1];
rz(-1.6233007) q[1];
sx q[1];
rz(2.1173304) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3108008) q[3];
sx q[3];
rz(-1.0648111) q[3];
sx q[3];
rz(2.3392086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.0236437) q[2];
sx q[2];
rz(-3.1274319) q[2];
sx q[2];
rz(-0.99948779) q[2];
rz(0.1855447) q[3];
sx q[3];
rz(-2.1047635) q[3];
sx q[3];
rz(-3.1256092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9284116) q[0];
sx q[0];
rz(-0.11225926) q[0];
sx q[0];
rz(-0.66473329) q[0];
rz(-0.96066535) q[1];
sx q[1];
rz(-0.04048368) q[1];
sx q[1];
rz(-1.6368846) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5211613) q[0];
sx q[0];
rz(-1.7502494) q[0];
sx q[0];
rz(1.8779264) q[0];
rz(-pi) q[1];
rz(0.13041741) q[2];
sx q[2];
rz(-2.7777024) q[2];
sx q[2];
rz(-0.39375776) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.79466984) q[1];
sx q[1];
rz(-2.0357913) q[1];
sx q[1];
rz(-2.5645606) q[1];
rz(2.2436813) q[3];
sx q[3];
rz(-0.6403044) q[3];
sx q[3];
rz(-1.698157) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.8399743) q[2];
sx q[2];
rz(-3.1349389) q[2];
sx q[2];
rz(1.1608231) q[2];
rz(-0.079856722) q[3];
sx q[3];
rz(-3.1294332) q[3];
sx q[3];
rz(-1.9399835) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53497159) q[0];
sx q[0];
rz(-1.4749682) q[0];
sx q[0];
rz(-1.5780021) q[0];
rz(3.1040991) q[1];
sx q[1];
rz(-2.9480724) q[1];
sx q[1];
rz(-2.9185157) q[1];
rz(-1.846215) q[2];
sx q[2];
rz(-0.38566914) q[2];
sx q[2];
rz(1.6657823) q[2];
rz(1.5885872) q[3];
sx q[3];
rz(-0.67100009) q[3];
sx q[3];
rz(-1.3081076) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
