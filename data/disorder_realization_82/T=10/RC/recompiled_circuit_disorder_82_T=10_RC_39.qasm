OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.3368971) q[0];
sx q[0];
rz(-2.1043632) q[0];
sx q[0];
rz(-0.35559911) q[0];
rz(-0.30272499) q[1];
sx q[1];
rz(-2.0974789) q[1];
sx q[1];
rz(1.8619327) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1188388) q[0];
sx q[0];
rz(-1.4316598) q[0];
sx q[0];
rz(-1.820178) q[0];
x q[1];
rz(-1.7945292) q[2];
sx q[2];
rz(-1.9828084) q[2];
sx q[2];
rz(1.3734693) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.0919839) q[1];
sx q[1];
rz(-2.1503452) q[1];
sx q[1];
rz(2.0642573) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7765331) q[3];
sx q[3];
rz(-1.7997026) q[3];
sx q[3];
rz(-2.8347901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.1352284) q[2];
sx q[2];
rz(-2.2029115) q[2];
sx q[2];
rz(-0.65650666) q[2];
rz(-2.4025829) q[3];
sx q[3];
rz(-2.6813172) q[3];
sx q[3];
rz(2.7242993) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1332557) q[0];
sx q[0];
rz(-0.79611859) q[0];
sx q[0];
rz(0.59536368) q[0];
rz(0.061925109) q[1];
sx q[1];
rz(-1.2281111) q[1];
sx q[1];
rz(2.6541236) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9092642) q[0];
sx q[0];
rz(-0.10350138) q[0];
sx q[0];
rz(-2.2444025) q[0];
x q[1];
rz(-2.3357453) q[2];
sx q[2];
rz(-0.57493756) q[2];
sx q[2];
rz(-0.88428674) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-3.0812379) q[1];
sx q[1];
rz(-1.758467) q[1];
sx q[1];
rz(-0.16138046) q[1];
x q[2];
rz(-1.8879714) q[3];
sx q[3];
rz(-1.653169) q[3];
sx q[3];
rz(1.9272643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.1797103) q[2];
sx q[2];
rz(-1.2625182) q[2];
sx q[2];
rz(2.7462192) q[2];
rz(-1.0428492) q[3];
sx q[3];
rz(-2.6104749) q[3];
sx q[3];
rz(3.1029207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9449126) q[0];
sx q[0];
rz(-1.1059462) q[0];
sx q[0];
rz(-2.9512067) q[0];
rz(-0.12292513) q[1];
sx q[1];
rz(-0.38750896) q[1];
sx q[1];
rz(-0.22274676) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5771772) q[0];
sx q[0];
rz(-2.5860997) q[0];
sx q[0];
rz(2.3781858) q[0];
x q[1];
rz(-2.9163755) q[2];
sx q[2];
rz(-1.2080964) q[2];
sx q[2];
rz(-0.97937102) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.45914868) q[1];
sx q[1];
rz(-2.5249081) q[1];
sx q[1];
rz(3.0241806) q[1];
x q[2];
rz(0.95007105) q[3];
sx q[3];
rz(-2.7347703) q[3];
sx q[3];
rz(-1.0328968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.8824076) q[2];
sx q[2];
rz(-1.2732482) q[2];
sx q[2];
rz(1.1941236) q[2];
rz(-2.0866701) q[3];
sx q[3];
rz(-1.6566365) q[3];
sx q[3];
rz(0.96364337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7524183) q[0];
sx q[0];
rz(-0.27802935) q[0];
sx q[0];
rz(-1.7383204) q[0];
rz(1.6482884) q[1];
sx q[1];
rz(-1.5416668) q[1];
sx q[1];
rz(0.9202252) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5266787) q[0];
sx q[0];
rz(-0.81256142) q[0];
sx q[0];
rz(-0.99292361) q[0];
rz(2.780519) q[2];
sx q[2];
rz(-1.4898584) q[2];
sx q[2];
rz(2.0780448) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.11776609) q[1];
sx q[1];
rz(-2.1823366) q[1];
sx q[1];
rz(2.141342) q[1];
rz(-pi) q[2];
rz(2.782015) q[3];
sx q[3];
rz(-1.3782856) q[3];
sx q[3];
rz(1.0958835) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.197864) q[2];
sx q[2];
rz(-1.7284164) q[2];
sx q[2];
rz(-1.8614004) q[2];
rz(-0.81930339) q[3];
sx q[3];
rz(-1.3179444) q[3];
sx q[3];
rz(-1.7839446) q[3];
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
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6827877) q[0];
sx q[0];
rz(-1.3764494) q[0];
sx q[0];
rz(-1.4047594) q[0];
rz(0.76830307) q[1];
sx q[1];
rz(-2.4914425) q[1];
sx q[1];
rz(0.4531025) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26950726) q[0];
sx q[0];
rz(-0.082490248) q[0];
sx q[0];
rz(-2.0399658) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4278973) q[2];
sx q[2];
rz(-1.60534) q[2];
sx q[2];
rz(-0.75682109) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.52757712) q[1];
sx q[1];
rz(-1.3506538) q[1];
sx q[1];
rz(-0.084053587) q[1];
rz(0.971332) q[3];
sx q[3];
rz(-1.9595651) q[3];
sx q[3];
rz(0.99861162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.12864628) q[2];
sx q[2];
rz(-2.3036028) q[2];
sx q[2];
rz(1.6298693) q[2];
rz(-2.417918) q[3];
sx q[3];
rz(-1.8975763) q[3];
sx q[3];
rz(2.2912912) q[3];
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
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1081651) q[0];
sx q[0];
rz(-1.4135679) q[0];
sx q[0];
rz(-0.23813716) q[0];
rz(-2.761633) q[1];
sx q[1];
rz(-2.0894876) q[1];
sx q[1];
rz(-1.628081) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3753525) q[0];
sx q[0];
rz(-1.1002812) q[0];
sx q[0];
rz(1.9348295) q[0];
rz(0.061821826) q[2];
sx q[2];
rz(-1.1218058) q[2];
sx q[2];
rz(-2.2668554) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.43135333) q[1];
sx q[1];
rz(-0.87190404) q[1];
sx q[1];
rz(2.5297574) q[1];
rz(2.7558277) q[3];
sx q[3];
rz(-1.9915951) q[3];
sx q[3];
rz(-2.5667131) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.05904077) q[2];
sx q[2];
rz(-0.5849134) q[2];
sx q[2];
rz(1.9980105) q[2];
rz(3.011076) q[3];
sx q[3];
rz(-1.427622) q[3];
sx q[3];
rz(2.9706764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0156353) q[0];
sx q[0];
rz(-2.1645249) q[0];
sx q[0];
rz(-0.39757279) q[0];
rz(1.3145087) q[1];
sx q[1];
rz(-0.54388261) q[1];
sx q[1];
rz(-1.9546753) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8966973) q[0];
sx q[0];
rz(-0.58422409) q[0];
sx q[0];
rz(-1.0446192) q[0];
x q[1];
rz(2.6830964) q[2];
sx q[2];
rz(-0.065932238) q[2];
sx q[2];
rz(2.3022848) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6656875) q[1];
sx q[1];
rz(-0.74531065) q[1];
sx q[1];
rz(1.3728632) q[1];
rz(-pi) q[2];
rz(0.39077057) q[3];
sx q[3];
rz(-0.6676538) q[3];
sx q[3];
rz(0.30927502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.879803) q[2];
sx q[2];
rz(-1.3061085) q[2];
sx q[2];
rz(1.3558033) q[2];
rz(1.5073744) q[3];
sx q[3];
rz(-0.58871388) q[3];
sx q[3];
rz(2.142895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3223406) q[0];
sx q[0];
rz(-1.1871908) q[0];
sx q[0];
rz(2.2858455) q[0];
rz(0.021727173) q[1];
sx q[1];
rz(-1.0236434) q[1];
sx q[1];
rz(1.0303248) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54661575) q[0];
sx q[0];
rz(-0.82406509) q[0];
sx q[0];
rz(1.0975518) q[0];
rz(1.4872929) q[2];
sx q[2];
rz(-2.6132772) q[2];
sx q[2];
rz(2.5022142) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.028138782) q[1];
sx q[1];
rz(-2.8061562) q[1];
sx q[1];
rz(2.0466652) q[1];
x q[2];
rz(0.52328531) q[3];
sx q[3];
rz(-1.0417632) q[3];
sx q[3];
rz(-0.39213359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.29256233) q[2];
sx q[2];
rz(-1.8886671) q[2];
sx q[2];
rz(-0.070177468) q[2];
rz(0.82693806) q[3];
sx q[3];
rz(-1.6707872) q[3];
sx q[3];
rz(2.5551445) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4940015) q[0];
sx q[0];
rz(-1.7140472) q[0];
sx q[0];
rz(-0.31627396) q[0];
rz(2.0896185) q[1];
sx q[1];
rz(-0.72967044) q[1];
sx q[1];
rz(1.1351599) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0907744) q[0];
sx q[0];
rz(-0.43681112) q[0];
sx q[0];
rz(0.45849053) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.029149292) q[2];
sx q[2];
rz(-1.2505194) q[2];
sx q[2];
rz(-2.0707891) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3865349) q[1];
sx q[1];
rz(-0.48663501) q[1];
sx q[1];
rz(1.4450032) q[1];
rz(-pi) q[2];
x q[2];
rz(3.1177403) q[3];
sx q[3];
rz(-2.1937222) q[3];
sx q[3];
rz(2.9305487) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.51745522) q[2];
sx q[2];
rz(-0.74666658) q[2];
sx q[2];
rz(1.5448145) q[2];
rz(2.4638713) q[3];
sx q[3];
rz(-0.8422519) q[3];
sx q[3];
rz(2.8519582) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90821663) q[0];
sx q[0];
rz(-2.4514618) q[0];
sx q[0];
rz(2.7897575) q[0];
rz(0.31967638) q[1];
sx q[1];
rz(-2.7668178) q[1];
sx q[1];
rz(2.9454254) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88403945) q[0];
sx q[0];
rz(-1.0172052) q[0];
sx q[0];
rz(-2.2474225) q[0];
rz(-1.9253795) q[2];
sx q[2];
rz(-2.9973929) q[2];
sx q[2];
rz(-1.9229638) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.0386104) q[1];
sx q[1];
rz(-2.3048232) q[1];
sx q[1];
rz(-0.02070133) q[1];
x q[2];
rz(-0.84368002) q[3];
sx q[3];
rz(-1.0495249) q[3];
sx q[3];
rz(2.664444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.7508042) q[2];
sx q[2];
rz(-1.5020341) q[2];
sx q[2];
rz(-0.081136726) q[2];
rz(1.1817415) q[3];
sx q[3];
rz(-2.2406082) q[3];
sx q[3];
rz(-2.0666163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1186196) q[0];
sx q[0];
rz(-1.58283) q[0];
sx q[0];
rz(1.3118623) q[0];
rz(-1.8267869) q[1];
sx q[1];
rz(-2.343315) q[1];
sx q[1];
rz(1.5409484) q[1];
rz(-1.6191471) q[2];
sx q[2];
rz(-1.804525) q[2];
sx q[2];
rz(-1.3934025) q[2];
rz(-1.7668016) q[3];
sx q[3];
rz(-2.9075665) q[3];
sx q[3];
rz(-2.8961765) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
