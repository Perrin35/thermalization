OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.91650668) q[0];
sx q[0];
rz(-2.9986311) q[0];
sx q[0];
rz(-1.6530871) q[0];
rz(1.1324963) q[1];
sx q[1];
rz(-2.7096665) q[1];
sx q[1];
rz(-0.63634029) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36447517) q[0];
sx q[0];
rz(-2.416673) q[0];
sx q[0];
rz(-0.97282797) q[0];
rz(0.77677988) q[2];
sx q[2];
rz(-2.619639) q[2];
sx q[2];
rz(1.7372434) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.591037) q[1];
sx q[1];
rz(-1.6948957) q[1];
sx q[1];
rz(-1.1731338) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2872055) q[3];
sx q[3];
rz(-2.1463089) q[3];
sx q[3];
rz(-2.7137386) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1119614) q[2];
sx q[2];
rz(-2.0165069) q[2];
sx q[2];
rz(-2.4911528) q[2];
rz(-1.316831) q[3];
sx q[3];
rz(-1.3464758) q[3];
sx q[3];
rz(0.10183798) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0874852) q[0];
sx q[0];
rz(-2.5547145) q[0];
sx q[0];
rz(0.45463872) q[0];
rz(3.1184323) q[1];
sx q[1];
rz(-1.4440447) q[1];
sx q[1];
rz(0.32618162) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9232193) q[0];
sx q[0];
rz(-1.0207286) q[0];
sx q[0];
rz(2.0327507) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1670879) q[2];
sx q[2];
rz(-0.49879227) q[2];
sx q[2];
rz(1.6287664) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4224976) q[1];
sx q[1];
rz(-2.0764253) q[1];
sx q[1];
rz(-1.0984332) q[1];
x q[2];
rz(2.2315027) q[3];
sx q[3];
rz(-1.6636563) q[3];
sx q[3];
rz(1.6805436) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0252016) q[2];
sx q[2];
rz(-1.2025183) q[2];
sx q[2];
rz(-2.1873059) q[2];
rz(1.2497692) q[3];
sx q[3];
rz(-0.3229177) q[3];
sx q[3];
rz(3.1402816) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-1.8702451) q[0];
sx q[0];
rz(-1.988669) q[0];
sx q[0];
rz(1.1767607) q[0];
rz(-0.36508834) q[1];
sx q[1];
rz(-1.3026214) q[1];
sx q[1];
rz(-1.5286068) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5919892) q[0];
sx q[0];
rz(-3.1168584) q[0];
sx q[0];
rz(1.0036841) q[0];
x q[1];
rz(-0.079054376) q[2];
sx q[2];
rz(-2.0131265) q[2];
sx q[2];
rz(-2.4238264) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.6512201) q[1];
sx q[1];
rz(-1.4191322) q[1];
sx q[1];
rz(0.030812736) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5131475) q[3];
sx q[3];
rz(-1.5022657) q[3];
sx q[3];
rz(1.203361) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.24885808) q[2];
sx q[2];
rz(-2.5881519) q[2];
sx q[2];
rz(-1.152732) q[2];
rz(1.8966127) q[3];
sx q[3];
rz(-2.1608519) q[3];
sx q[3];
rz(-0.28321701) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(0.98274851) q[0];
sx q[0];
rz(-11*pi/12) q[0];
sx q[0];
rz(0.68956476) q[0];
rz(-0.025029643) q[1];
sx q[1];
rz(-1.6233147) q[1];
sx q[1];
rz(1.7207346) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7574437) q[0];
sx q[0];
rz(-1.4673646) q[0];
sx q[0];
rz(-2.5162656) q[0];
x q[1];
rz(1.1406607) q[2];
sx q[2];
rz(-1.4073512) q[2];
sx q[2];
rz(0.30980936) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.017548718) q[1];
sx q[1];
rz(-0.50395233) q[1];
sx q[1];
rz(-0.057226463) q[1];
x q[2];
rz(2.8042364) q[3];
sx q[3];
rz(-1.3135664) q[3];
sx q[3];
rz(-1.6599865) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.932852) q[2];
sx q[2];
rz(-0.17540652) q[2];
sx q[2];
rz(1.1669195) q[2];
rz(2.6973727) q[3];
sx q[3];
rz(-1.2362213) q[3];
sx q[3];
rz(-2.8319632) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2883478) q[0];
sx q[0];
rz(-1.7720368) q[0];
sx q[0];
rz(0.41611588) q[0];
rz(0.5101282) q[1];
sx q[1];
rz(-2.3704539) q[1];
sx q[1];
rz(-0.77521926) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32019553) q[0];
sx q[0];
rz(-0.96548432) q[0];
sx q[0];
rz(1.8138252) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.40521605) q[2];
sx q[2];
rz(-1.7381258) q[2];
sx q[2];
rz(0.93003231) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.21501274) q[1];
sx q[1];
rz(-1.6832507) q[1];
sx q[1];
rz(1.6968898) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6121077) q[3];
sx q[3];
rz(-2.4765402) q[3];
sx q[3];
rz(-0.82529008) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.2148718) q[2];
sx q[2];
rz(-2.1472609) q[2];
sx q[2];
rz(-2.1057687) q[2];
rz(-0.18946798) q[3];
sx q[3];
rz(-2.6697956) q[3];
sx q[3];
rz(-2.6389879) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8871317) q[0];
sx q[0];
rz(-0.04763617) q[0];
sx q[0];
rz(0.3558085) q[0];
rz(1.4954781) q[1];
sx q[1];
rz(-2.1864086) q[1];
sx q[1];
rz(1.1411508) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1352859) q[0];
sx q[0];
rz(-1.7809488) q[0];
sx q[0];
rz(-2.2457613) q[0];
rz(-pi) q[1];
rz(-0.009013225) q[2];
sx q[2];
rz(-2.1823332) q[2];
sx q[2];
rz(2.1076941) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.0547949) q[1];
sx q[1];
rz(-1.7112477) q[1];
sx q[1];
rz(1.676692) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0047046) q[3];
sx q[3];
rz(-1.5597789) q[3];
sx q[3];
rz(-2.9316255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.63153875) q[2];
sx q[2];
rz(-0.71439356) q[2];
sx q[2];
rz(1.4368524) q[2];
rz(1.0531744) q[3];
sx q[3];
rz(-1.5318233) q[3];
sx q[3];
rz(-2.6385782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18008867) q[0];
sx q[0];
rz(-0.29619521) q[0];
sx q[0];
rz(-0.12251138) q[0];
rz(0.65613121) q[1];
sx q[1];
rz(-1.8729112) q[1];
sx q[1];
rz(-1.8001385) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.487566) q[0];
sx q[0];
rz(-1.5254505) q[0];
sx q[0];
rz(1.5349755) q[0];
x q[1];
rz(-0.28320988) q[2];
sx q[2];
rz(-1.4669751) q[2];
sx q[2];
rz(1.4701172) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.0602051) q[1];
sx q[1];
rz(-1.3045746) q[1];
sx q[1];
rz(2.2039444) q[1];
rz(-2.9060494) q[3];
sx q[3];
rz(-0.45484124) q[3];
sx q[3];
rz(1.7419635) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.573367) q[2];
sx q[2];
rz(-1.2218916) q[2];
sx q[2];
rz(-1.4637671) q[2];
rz(0.73417869) q[3];
sx q[3];
rz(-0.28407431) q[3];
sx q[3];
rz(-1.3045788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3463335) q[0];
sx q[0];
rz(-1.3672375) q[0];
sx q[0];
rz(-2.6829868) q[0];
rz(-2.1268225) q[1];
sx q[1];
rz(-1.1943123) q[1];
sx q[1];
rz(1.0626622) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39358562) q[0];
sx q[0];
rz(-0.17608041) q[0];
sx q[0];
rz(0.3349456) q[0];
rz(-2.6133461) q[2];
sx q[2];
rz(-2.2604803) q[2];
sx q[2];
rz(1.6843585) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.52673756) q[1];
sx q[1];
rz(-1.7018082) q[1];
sx q[1];
rz(3.0796771) q[1];
x q[2];
rz(2.0902781) q[3];
sx q[3];
rz(-1.2881345) q[3];
sx q[3];
rz(-1.9611028) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.8899272) q[2];
sx q[2];
rz(-1.3883611) q[2];
sx q[2];
rz(-1.3512705) q[2];
rz(1.2190367) q[3];
sx q[3];
rz(-2.7128897) q[3];
sx q[3];
rz(-0.8333227) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0625921) q[0];
sx q[0];
rz(-1.3840249) q[0];
sx q[0];
rz(-0.90743995) q[0];
rz(-2.9145248) q[1];
sx q[1];
rz(-1.0457958) q[1];
sx q[1];
rz(-0.89920941) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2405906) q[0];
sx q[0];
rz(-2.0792897) q[0];
sx q[0];
rz(2.2700538) q[0];
x q[1];
rz(1.6532517) q[2];
sx q[2];
rz(-1.4926148) q[2];
sx q[2];
rz(2.9417335) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.9252852) q[1];
sx q[1];
rz(-1.3853964) q[1];
sx q[1];
rz(-2.0865296) q[1];
rz(3.090631) q[3];
sx q[3];
rz(-0.82333857) q[3];
sx q[3];
rz(-2.2281102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.9565309) q[2];
sx q[2];
rz(-1.9756292) q[2];
sx q[2];
rz(0.60066191) q[2];
rz(-1.0968084) q[3];
sx q[3];
rz(-0.59022248) q[3];
sx q[3];
rz(-2.2685331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3908865) q[0];
sx q[0];
rz(-2.0426671) q[0];
sx q[0];
rz(0.98439687) q[0];
rz(2.6495433) q[1];
sx q[1];
rz(-1.4130519) q[1];
sx q[1];
rz(1.1955998) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9162684) q[0];
sx q[0];
rz(-1.9231503) q[0];
sx q[0];
rz(-0.63017643) q[0];
x q[1];
rz(-0.97862794) q[2];
sx q[2];
rz(-2.371863) q[2];
sx q[2];
rz(-2.9362188) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8859133) q[1];
sx q[1];
rz(-1.1885841) q[1];
sx q[1];
rz(-0.087319386) q[1];
rz(-pi) q[2];
rz(-1.1117842) q[3];
sx q[3];
rz(-2.2892366) q[3];
sx q[3];
rz(-1.8000923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.0857346) q[2];
sx q[2];
rz(-1.5208289) q[2];
sx q[2];
rz(0.98304191) q[2];
rz(-2.8899657) q[3];
sx q[3];
rz(-2.7139137) q[3];
sx q[3];
rz(1.0035286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6428103) q[0];
sx q[0];
rz(-2.1693873) q[0];
sx q[0];
rz(0.65709773) q[0];
rz(1.2596754) q[1];
sx q[1];
rz(-1.4114264) q[1];
sx q[1];
rz(-2.2687601) q[1];
rz(-0.85763422) q[2];
sx q[2];
rz(-1.7339755) q[2];
sx q[2];
rz(-0.82790464) q[2];
rz(-1.2498517) q[3];
sx q[3];
rz(-1.8886011) q[3];
sx q[3];
rz(-0.26816751) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
