OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.3172265) q[0];
sx q[0];
rz(-2.0269725) q[0];
sx q[0];
rz(0.00014076509) q[0];
rz(-1.8074942) q[1];
sx q[1];
rz(-0.9642095) q[1];
sx q[1];
rz(-1.1934086) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6407335) q[0];
sx q[0];
rz(-1.2841604) q[0];
sx q[0];
rz(-2.9564234) q[0];
rz(-pi) q[1];
x q[1];
rz(0.46618669) q[2];
sx q[2];
rz(-0.59980118) q[2];
sx q[2];
rz(-0.28238645) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.57230091) q[1];
sx q[1];
rz(-0.83218677) q[1];
sx q[1];
rz(-0.65170793) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7883349) q[3];
sx q[3];
rz(-1.6780403) q[3];
sx q[3];
rz(-1.4835964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.6821735) q[2];
sx q[2];
rz(-3.1176304) q[2];
sx q[2];
rz(1.2288644) q[2];
rz(-1.4131644) q[3];
sx q[3];
rz(-2.0404405) q[3];
sx q[3];
rz(1.4878954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5380149) q[0];
sx q[0];
rz(-1.5025654) q[0];
sx q[0];
rz(2.1287825) q[0];
rz(3.1139328) q[1];
sx q[1];
rz(-2.467997) q[1];
sx q[1];
rz(2.0181296) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9155974) q[0];
sx q[0];
rz(-1.5132656) q[0];
sx q[0];
rz(-2.0321839) q[0];
rz(-pi) q[1];
rz(-0.063617184) q[2];
sx q[2];
rz(-2.3607495) q[2];
sx q[2];
rz(2.7240679) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0151129) q[1];
sx q[1];
rz(-1.0547332) q[1];
sx q[1];
rz(-0.59241398) q[1];
x q[2];
rz(0.68617679) q[3];
sx q[3];
rz(-2.3585329) q[3];
sx q[3];
rz(-2.1427597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3479487) q[2];
sx q[2];
rz(-1.0898033) q[2];
sx q[2];
rz(2.222555) q[2];
rz(-2.4675026) q[3];
sx q[3];
rz(-0.6522817) q[3];
sx q[3];
rz(-1.526171) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27750257) q[0];
sx q[0];
rz(-0.16177495) q[0];
sx q[0];
rz(1.2751689) q[0];
rz(-0.69349849) q[1];
sx q[1];
rz(-1.2561412) q[1];
sx q[1];
rz(2.0085874) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61921652) q[0];
sx q[0];
rz(-1.4285061) q[0];
sx q[0];
rz(-0.0033358047) q[0];
rz(-pi) q[1];
rz(2.0793545) q[2];
sx q[2];
rz(-0.8478176) q[2];
sx q[2];
rz(1.6194956) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.36724597) q[1];
sx q[1];
rz(-2.5173752) q[1];
sx q[1];
rz(1.063785) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5127701) q[3];
sx q[3];
rz(-2.0938794) q[3];
sx q[3];
rz(-2.3716225) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.2514078) q[2];
sx q[2];
rz(-0.79139411) q[2];
sx q[2];
rz(1.2934925) q[2];
rz(0.039316468) q[3];
sx q[3];
rz(-1.9226363) q[3];
sx q[3];
rz(-1.2600651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8816198) q[0];
sx q[0];
rz(-0.078475229) q[0];
sx q[0];
rz(1.9807293) q[0];
rz(-0.89598957) q[1];
sx q[1];
rz(-1.4410102) q[1];
sx q[1];
rz(-0.13555759) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4281222) q[0];
sx q[0];
rz(-2.2566416) q[0];
sx q[0];
rz(1.1629521) q[0];
rz(-pi) q[1];
rz(-1.3877669) q[2];
sx q[2];
rz(-0.39003885) q[2];
sx q[2];
rz(2.4412465) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.902066) q[1];
sx q[1];
rz(-0.24589989) q[1];
sx q[1];
rz(-0.86218254) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0566063) q[3];
sx q[3];
rz(-2.8286472) q[3];
sx q[3];
rz(-1.7900975) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.23665145) q[2];
sx q[2];
rz(-2.1950978) q[2];
sx q[2];
rz(0.87990749) q[2];
rz(3.0974292) q[3];
sx q[3];
rz(-1.6396089) q[3];
sx q[3];
rz(-2.8529609) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0376461) q[0];
sx q[0];
rz(-2.7665311) q[0];
sx q[0];
rz(-1.0132382) q[0];
rz(-3.0918616) q[1];
sx q[1];
rz(-2.2278992) q[1];
sx q[1];
rz(-1.0838881) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5686544) q[0];
sx q[0];
rz(-1.3875811) q[0];
sx q[0];
rz(1.818548) q[0];
rz(-1.8439699) q[2];
sx q[2];
rz(-1.8130842) q[2];
sx q[2];
rz(-0.66852409) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0005972) q[1];
sx q[1];
rz(-2.6037569) q[1];
sx q[1];
rz(1.3586033) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6277222) q[3];
sx q[3];
rz(-1.6146982) q[3];
sx q[3];
rz(1.9572789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.9087387) q[2];
sx q[2];
rz(-2.8149657) q[2];
sx q[2];
rz(-2.8971635) q[2];
rz(-2.7092253) q[3];
sx q[3];
rz(-1.3997388) q[3];
sx q[3];
rz(2.6385245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2844834) q[0];
sx q[0];
rz(-1.720022) q[0];
sx q[0];
rz(-0.094141468) q[0];
rz(-2.969818) q[1];
sx q[1];
rz(-1.1356907) q[1];
sx q[1];
rz(-2.24618) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6796414) q[0];
sx q[0];
rz(-1.9129176) q[0];
sx q[0];
rz(-1.6109986) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8199811) q[2];
sx q[2];
rz(-0.65882896) q[2];
sx q[2];
rz(-2.4668601) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.10455924) q[1];
sx q[1];
rz(-1.4022786) q[1];
sx q[1];
rz(-0.044635459) q[1];
rz(-pi) q[2];
x q[2];
rz(3.1392519) q[3];
sx q[3];
rz(-1.4961092) q[3];
sx q[3];
rz(0.60126388) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.0078997) q[2];
sx q[2];
rz(-2.7329625) q[2];
sx q[2];
rz(0.80319476) q[2];
rz(-1.1903654) q[3];
sx q[3];
rz(-1.9093711) q[3];
sx q[3];
rz(-0.41263321) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.068709277) q[0];
sx q[0];
rz(-0.16462737) q[0];
sx q[0];
rz(2.6224526) q[0];
rz(-2.5601162) q[1];
sx q[1];
rz(-2.0362208) q[1];
sx q[1];
rz(1.2566459) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2557247) q[0];
sx q[0];
rz(-2.8868544) q[0];
sx q[0];
rz(1.4101009) q[0];
rz(-1.16876) q[2];
sx q[2];
rz(-2.6425344) q[2];
sx q[2];
rz(-1.7390342) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7278442) q[1];
sx q[1];
rz(-1.2894948) q[1];
sx q[1];
rz(-1.6378228) q[1];
rz(-pi) q[2];
rz(2.4323746) q[3];
sx q[3];
rz(-1.8172482) q[3];
sx q[3];
rz(2.646288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1533623) q[2];
sx q[2];
rz(-2.1116657) q[2];
sx q[2];
rz(-1.3640277) q[2];
rz(2.2310232) q[3];
sx q[3];
rz(-1.1547337) q[3];
sx q[3];
rz(1.5301269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0664739) q[0];
sx q[0];
rz(-0.56448889) q[0];
sx q[0];
rz(-2.8334154) q[0];
rz(0.072487436) q[1];
sx q[1];
rz(-1.0132353) q[1];
sx q[1];
rz(2.7546308) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2012126) q[0];
sx q[0];
rz(-1.7965172) q[0];
sx q[0];
rz(-0.97198995) q[0];
rz(-0.87256356) q[2];
sx q[2];
rz(-1.1485032) q[2];
sx q[2];
rz(2.1625105) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.19883979) q[1];
sx q[1];
rz(-1.6457335) q[1];
sx q[1];
rz(-2.4692175) q[1];
rz(-pi) q[2];
rz(2.6450063) q[3];
sx q[3];
rz(-1.7644617) q[3];
sx q[3];
rz(2.1878939) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5179634) q[2];
sx q[2];
rz(-1.364418) q[2];
sx q[2];
rz(-2.7015838) q[2];
rz(0.7157588) q[3];
sx q[3];
rz(-1.7093168) q[3];
sx q[3];
rz(-1.0796775) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0004262) q[0];
sx q[0];
rz(-2.3957802) q[0];
sx q[0];
rz(2.0429042) q[0];
rz(0.72775841) q[1];
sx q[1];
rz(-2.7658503) q[1];
sx q[1];
rz(-3.0922906) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4816149) q[0];
sx q[0];
rz(-1.0064631) q[0];
sx q[0];
rz(0.83156395) q[0];
rz(-2.1575035) q[2];
sx q[2];
rz(-2.3684635) q[2];
sx q[2];
rz(0.099260515) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2652492) q[1];
sx q[1];
rz(-1.1117522) q[1];
sx q[1];
rz(-0.52258073) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6155869) q[3];
sx q[3];
rz(-0.81323871) q[3];
sx q[3];
rz(-0.40532743) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.07842841) q[2];
sx q[2];
rz(-2.5421263) q[2];
sx q[2];
rz(2.4196529) q[2];
rz(2.1980964) q[3];
sx q[3];
rz(-2.3908581) q[3];
sx q[3];
rz(0.25434428) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14770517) q[0];
sx q[0];
rz(-1.9839956) q[0];
sx q[0];
rz(-1.0797427) q[0];
rz(1.059277) q[1];
sx q[1];
rz(-0.22288999) q[1];
sx q[1];
rz(1.7396897) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7584383) q[0];
sx q[0];
rz(-0.24807319) q[0];
sx q[0];
rz(1.0323348) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9981668) q[2];
sx q[2];
rz(-2.9529245) q[2];
sx q[2];
rz(0.27486899) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.460443) q[1];
sx q[1];
rz(-2.3773758) q[1];
sx q[1];
rz(-1.1230254) q[1];
rz(-pi) q[2];
rz(1.363477) q[3];
sx q[3];
rz(-0.94982409) q[3];
sx q[3];
rz(-0.45351115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.6587276) q[2];
sx q[2];
rz(-1.3663224) q[2];
sx q[2];
rz(1.6213017) q[2];
rz(-0.55082095) q[3];
sx q[3];
rz(-2.3362624) q[3];
sx q[3];
rz(2.4441392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.993492) q[0];
sx q[0];
rz(-1.8363331) q[0];
sx q[0];
rz(1.6114417) q[0];
rz(0.91611721) q[1];
sx q[1];
rz(-0.59090186) q[1];
sx q[1];
rz(-0.59060243) q[1];
rz(-2.3475636) q[2];
sx q[2];
rz(-2.2326438) q[2];
sx q[2];
rz(2.2868962) q[2];
rz(-1.0536853) q[3];
sx q[3];
rz(-1.5862982) q[3];
sx q[3];
rz(1.3261212) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
