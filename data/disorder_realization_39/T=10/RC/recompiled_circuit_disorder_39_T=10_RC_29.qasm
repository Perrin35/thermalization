OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.0948148) q[0];
sx q[0];
rz(4.2098213) q[0];
sx q[0];
rz(9.8888483) q[0];
rz(1.9595454) q[1];
sx q[1];
rz(6.2160677) q[1];
sx q[1];
rz(10.709229) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0811631) q[0];
sx q[0];
rz(-1.2788749) q[0];
sx q[0];
rz(1.1638327) q[0];
x q[1];
rz(1.5377827) q[2];
sx q[2];
rz(-2.0602977) q[2];
sx q[2];
rz(-2.2061493) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.1844993) q[1];
sx q[1];
rz(-2.543499) q[1];
sx q[1];
rz(-0.68006541) q[1];
rz(-pi) q[2];
rz(-0.10847096) q[3];
sx q[3];
rz(-2.9900108) q[3];
sx q[3];
rz(-0.82419318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.7444732) q[2];
sx q[2];
rz(-2.2019272) q[2];
sx q[2];
rz(-2.822067) q[2];
rz(-0.5784936) q[3];
sx q[3];
rz(-2.6632023) q[3];
sx q[3];
rz(2.4676676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7330866) q[0];
sx q[0];
rz(-1.5698313) q[0];
sx q[0];
rz(2.615036) q[0];
rz(2.5780442) q[1];
sx q[1];
rz(-1.6105517) q[1];
sx q[1];
rz(0.79663509) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63882534) q[0];
sx q[0];
rz(-2.4135547) q[0];
sx q[0];
rz(-1.0685705) q[0];
rz(2.8393306) q[2];
sx q[2];
rz(-2.5455591) q[2];
sx q[2];
rz(-2.9994534) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.7318774) q[1];
sx q[1];
rz(-2.5341923) q[1];
sx q[1];
rz(3.1190447) q[1];
x q[2];
rz(-0.91005743) q[3];
sx q[3];
rz(-2.4168192) q[3];
sx q[3];
rz(-1.7006601) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.18156302) q[2];
sx q[2];
rz(-1.3788362) q[2];
sx q[2];
rz(-0.78655085) q[2];
rz(-0.49318796) q[3];
sx q[3];
rz(-1.9599873) q[3];
sx q[3];
rz(2.6942159) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86984533) q[0];
sx q[0];
rz(-2.2646876) q[0];
sx q[0];
rz(-1.7012117) q[0];
rz(2.4213743) q[1];
sx q[1];
rz(-0.59599382) q[1];
sx q[1];
rz(2.4386141) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5876578) q[0];
sx q[0];
rz(-0.30447391) q[0];
sx q[0];
rz(-1.9703883) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.643232) q[2];
sx q[2];
rz(-1.8866072) q[2];
sx q[2];
rz(0.74893803) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.30333334) q[1];
sx q[1];
rz(-1.6175744) q[1];
sx q[1];
rz(-0.10951885) q[1];
rz(2.1941575) q[3];
sx q[3];
rz(-1.7121592) q[3];
sx q[3];
rz(-2.3879104) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.8748223) q[2];
sx q[2];
rz(-1.6569933) q[2];
sx q[2];
rz(2.2300569) q[2];
rz(-0.91397816) q[3];
sx q[3];
rz(-1.4303047) q[3];
sx q[3];
rz(2.3142557) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5660969) q[0];
sx q[0];
rz(-2.5056582) q[0];
sx q[0];
rz(-2.3655868) q[0];
rz(-1.2674747) q[1];
sx q[1];
rz(-1.1088561) q[1];
sx q[1];
rz(-2.5783096) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7960498) q[0];
sx q[0];
rz(-2.4420218) q[0];
sx q[0];
rz(-1.1388586) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7411555) q[2];
sx q[2];
rz(-2.4106328) q[2];
sx q[2];
rz(-1.2618582) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.1174417) q[1];
sx q[1];
rz(-1.5772595) q[1];
sx q[1];
rz(0.27177377) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.71743439) q[3];
sx q[3];
rz(-1.7955901) q[3];
sx q[3];
rz(-0.7625398) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6528066) q[2];
sx q[2];
rz(-2.2061009) q[2];
sx q[2];
rz(0.164786) q[2];
rz(-0.22848836) q[3];
sx q[3];
rz(-2.861664) q[3];
sx q[3];
rz(-0.94648186) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27424681) q[0];
sx q[0];
rz(-2.2773401) q[0];
sx q[0];
rz(2.2891323) q[0];
rz(-0.35119855) q[1];
sx q[1];
rz(-0.57792592) q[1];
sx q[1];
rz(1.9794827) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0162109) q[0];
sx q[0];
rz(-2.1124766) q[0];
sx q[0];
rz(-2.0340232) q[0];
rz(-pi) q[1];
rz(-0.39761333) q[2];
sx q[2];
rz(-1.0828472) q[2];
sx q[2];
rz(0.60148009) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.1197549) q[1];
sx q[1];
rz(-0.86958414) q[1];
sx q[1];
rz(-0.55418684) q[1];
rz(-pi) q[2];
rz(2.4200053) q[3];
sx q[3];
rz(-1.5820832) q[3];
sx q[3];
rz(-1.23502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.6953485) q[2];
sx q[2];
rz(-2.7480795) q[2];
sx q[2];
rz(1.2472786) q[2];
rz(-0.013109664) q[3];
sx q[3];
rz(-1.7316827) q[3];
sx q[3];
rz(-2.9124027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2728249) q[0];
sx q[0];
rz(-1.2478963) q[0];
sx q[0];
rz(-0.75063467) q[0];
rz(2.0320832) q[1];
sx q[1];
rz(-2.7710319) q[1];
sx q[1];
rz(1.3060588) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7795777) q[0];
sx q[0];
rz(-2.0834196) q[0];
sx q[0];
rz(0.4431475) q[0];
rz(1.7734217) q[2];
sx q[2];
rz(-1.1154004) q[2];
sx q[2];
rz(2.1855598) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.27481025) q[1];
sx q[1];
rz(-0.85046235) q[1];
sx q[1];
rz(-1.5809098) q[1];
x q[2];
rz(-1.4218876) q[3];
sx q[3];
rz(-2.1544666) q[3];
sx q[3];
rz(-2.7910809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.3874454) q[2];
sx q[2];
rz(-1.0070609) q[2];
sx q[2];
rz(-0.60069096) q[2];
rz(2.1389652) q[3];
sx q[3];
rz(-2.6032344) q[3];
sx q[3];
rz(-0.35282648) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7464741) q[0];
sx q[0];
rz(-1.4557319) q[0];
sx q[0];
rz(-1.0466928) q[0];
rz(-1.612161) q[1];
sx q[1];
rz(-1.4271095) q[1];
sx q[1];
rz(-0.41710645) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5865267) q[0];
sx q[0];
rz(-1.546372) q[0];
sx q[0];
rz(1.7909554) q[0];
rz(-pi) q[1];
rz(-1.276937) q[2];
sx q[2];
rz(-1.8655348) q[2];
sx q[2];
rz(-0.35640946) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.6228094) q[1];
sx q[1];
rz(-2.2768524) q[1];
sx q[1];
rz(-0.6219567) q[1];
rz(-pi) q[2];
rz(-2.4884175) q[3];
sx q[3];
rz(-1.8110868) q[3];
sx q[3];
rz(-2.0865292) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.1356915) q[2];
sx q[2];
rz(-2.6689172) q[2];
sx q[2];
rz(1.3674412) q[2];
rz(2.588429) q[3];
sx q[3];
rz(-1.7382712) q[3];
sx q[3];
rz(-0.52545351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.816514) q[0];
sx q[0];
rz(-1.3229609) q[0];
sx q[0];
rz(1.1897855) q[0];
rz(-1.7143543) q[1];
sx q[1];
rz(-0.27806565) q[1];
sx q[1];
rz(3.0292125) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.053514078) q[0];
sx q[0];
rz(-1.3980165) q[0];
sx q[0];
rz(-2.0515576) q[0];
rz(-2.8078812) q[2];
sx q[2];
rz(-1.6534272) q[2];
sx q[2];
rz(-2.4090648) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.97265128) q[1];
sx q[1];
rz(-2.2691138) q[1];
sx q[1];
rz(2.6886743) q[1];
rz(-1.1376082) q[3];
sx q[3];
rz(-2.8578651) q[3];
sx q[3];
rz(-1.9445436) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.0671976) q[2];
sx q[2];
rz(-1.9125568) q[2];
sx q[2];
rz(1.7929662) q[2];
rz(1.9366692) q[3];
sx q[3];
rz(-1.4068853) q[3];
sx q[3];
rz(0.19781923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7444721) q[0];
sx q[0];
rz(-1.67698) q[0];
sx q[0];
rz(3.1066185) q[0];
rz(0.84683013) q[1];
sx q[1];
rz(-1.433082) q[1];
sx q[1];
rz(-0.91167489) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9104886) q[0];
sx q[0];
rz(-1.3470874) q[0];
sx q[0];
rz(1.847812) q[0];
rz(-pi) q[1];
rz(-1.4439092) q[2];
sx q[2];
rz(-1.03089) q[2];
sx q[2];
rz(1.1586231) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.24278495) q[1];
sx q[1];
rz(-0.50711942) q[1];
sx q[1];
rz(-0.42906638) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4089912) q[3];
sx q[3];
rz(-2.4676975) q[3];
sx q[3];
rz(-3.1042838) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.4336865) q[2];
sx q[2];
rz(-0.5138548) q[2];
sx q[2];
rz(0.24924499) q[2];
rz(-2.3748659) q[3];
sx q[3];
rz(-1.8079575) q[3];
sx q[3];
rz(-2.7419817) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7837759) q[0];
sx q[0];
rz(-2.0932842) q[0];
sx q[0];
rz(-2.7695079) q[0];
rz(0.58139873) q[1];
sx q[1];
rz(-0.77722469) q[1];
sx q[1];
rz(-1.7262329) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.088179528) q[0];
sx q[0];
rz(-1.3706494) q[0];
sx q[0];
rz(-2.8739909) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.66843372) q[2];
sx q[2];
rz(-2.4032776) q[2];
sx q[2];
rz(2.7065606) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.2409754) q[1];
sx q[1];
rz(-1.6458578) q[1];
sx q[1];
rz(1.5013298) q[1];
rz(-0.38378999) q[3];
sx q[3];
rz(-1.0597611) q[3];
sx q[3];
rz(-2.8904861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.8156478) q[2];
sx q[2];
rz(-3.0013913) q[2];
sx q[2];
rz(2.9369205) q[2];
rz(1.7278016) q[3];
sx q[3];
rz(-1.1433583) q[3];
sx q[3];
rz(2.0457101) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-2.6407912) q[0];
sx q[0];
rz(-1.4688107) q[0];
sx q[0];
rz(-2.2039913) q[0];
rz(-1.5851371) q[1];
sx q[1];
rz(-1.382348) q[1];
sx q[1];
rz(1.4437645) q[1];
rz(1.8150868) q[2];
sx q[2];
rz(-0.46637022) q[2];
sx q[2];
rz(-2.7844219) q[2];
rz(3.0588991) q[3];
sx q[3];
rz(-1.0001928) q[3];
sx q[3];
rz(2.115888) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
