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
rz(2.3020701) q[0];
sx q[0];
rz(-1.5538202) q[0];
sx q[0];
rz(2.177218) q[0];
rz(0.40274611) q[1];
sx q[1];
rz(-1.6778814) q[1];
sx q[1];
rz(-0.95199624) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.096199252) q[0];
sx q[0];
rz(-1.081776) q[0];
sx q[0];
rz(0.10639356) q[0];
x q[1];
rz(3.0430284) q[2];
sx q[2];
rz(-1.9910887) q[2];
sx q[2];
rz(-1.9914876) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.1237028) q[1];
sx q[1];
rz(-2.5933939) q[1];
sx q[1];
rz(2.3831559) q[1];
rz(-pi) q[2];
x q[2];
rz(0.49523109) q[3];
sx q[3];
rz(-2.3217045) q[3];
sx q[3];
rz(0.26671975) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.52557785) q[2];
sx q[2];
rz(-1.792045) q[2];
sx q[2];
rz(-2.7363321) q[2];
rz(-2.0112093) q[3];
sx q[3];
rz(-0.78961343) q[3];
sx q[3];
rz(1.214504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2363215) q[0];
sx q[0];
rz(-3.1089678) q[0];
sx q[0];
rz(-0.76500934) q[0];
rz(-0.24461034) q[1];
sx q[1];
rz(-1.2449934) q[1];
sx q[1];
rz(2.5027221) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8804358) q[0];
sx q[0];
rz(-0.64468274) q[0];
sx q[0];
rz(0.042244367) q[0];
x q[1];
rz(0.62221726) q[2];
sx q[2];
rz(-1.2583311) q[2];
sx q[2];
rz(1.086233) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.32989001) q[1];
sx q[1];
rz(-2.5479758) q[1];
sx q[1];
rz(2.4827186) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3078626) q[3];
sx q[3];
rz(-1.7639188) q[3];
sx q[3];
rz(1.5964791) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.2834393) q[2];
sx q[2];
rz(-2.6946113) q[2];
sx q[2];
rz(-0.24432527) q[2];
rz(-1.1089995) q[3];
sx q[3];
rz(-1.2764443) q[3];
sx q[3];
rz(-0.23787704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37794381) q[0];
sx q[0];
rz(-2.0574594) q[0];
sx q[0];
rz(-1.2373244) q[0];
rz(-1.7297277) q[1];
sx q[1];
rz(-0.4487764) q[1];
sx q[1];
rz(-1.3317187) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6645637) q[0];
sx q[0];
rz(-2.3630099) q[0];
sx q[0];
rz(0.12523366) q[0];
x q[1];
rz(-0.21816605) q[2];
sx q[2];
rz(-1.2274172) q[2];
sx q[2];
rz(-0.81264466) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.0983372) q[1];
sx q[1];
rz(-1.5374633) q[1];
sx q[1];
rz(3.0860677) q[1];
x q[2];
rz(2.3490155) q[3];
sx q[3];
rz(-0.78781414) q[3];
sx q[3];
rz(-0.28258309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.1639013) q[2];
sx q[2];
rz(-2.0745664) q[2];
sx q[2];
rz(0.65995556) q[2];
rz(-1.6948949) q[3];
sx q[3];
rz(-1.7132672) q[3];
sx q[3];
rz(-1.1032633) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(2.8834943) q[0];
sx q[0];
rz(-0.75089184) q[0];
sx q[0];
rz(2.0080436) q[0];
rz(-0.51620382) q[1];
sx q[1];
rz(-1.8323003) q[1];
sx q[1];
rz(-2.5925327) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7507197) q[0];
sx q[0];
rz(-1.735512) q[0];
sx q[0];
rz(-2.9956942) q[0];
rz(-1.8070776) q[2];
sx q[2];
rz(-2.5388814) q[2];
sx q[2];
rz(-0.38379764) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.258865) q[1];
sx q[1];
rz(-2.8587864) q[1];
sx q[1];
rz(2.905719) q[1];
rz(-1.9677343) q[3];
sx q[3];
rz(-1.8324495) q[3];
sx q[3];
rz(1.73571) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.05222008) q[2];
sx q[2];
rz(-2.1286271) q[2];
sx q[2];
rz(0.29903665) q[2];
rz(-0.61797577) q[3];
sx q[3];
rz(-2.0682014) q[3];
sx q[3];
rz(1.7834024) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
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
rz(-0.95737902) q[0];
sx q[0];
rz(-0.03802499) q[0];
sx q[0];
rz(0.46603611) q[0];
rz(2.274463) q[1];
sx q[1];
rz(-2.6301818) q[1];
sx q[1];
rz(-2.3066511) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4049618) q[0];
sx q[0];
rz(-1.5190655) q[0];
sx q[0];
rz(-1.5993656) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4136426) q[2];
sx q[2];
rz(-0.11874983) q[2];
sx q[2];
rz(-0.70357067) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7606733) q[1];
sx q[1];
rz(-1.5173418) q[1];
sx q[1];
rz(-1.0905114) q[1];
x q[2];
rz(-0.94423202) q[3];
sx q[3];
rz(-2.196518) q[3];
sx q[3];
rz(-2.0056796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.9338108) q[2];
sx q[2];
rz(-2.2489579) q[2];
sx q[2];
rz(3.1100698) q[2];
rz(-0.75782123) q[3];
sx q[3];
rz(-2.0407245) q[3];
sx q[3];
rz(2.5789564) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9071478) q[0];
sx q[0];
rz(-1.3749159) q[0];
sx q[0];
rz(0.031540792) q[0];
rz(0.08983287) q[1];
sx q[1];
rz(-2.641771) q[1];
sx q[1];
rz(-2.8057742) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9851097) q[0];
sx q[0];
rz(-1.8976189) q[0];
sx q[0];
rz(0.37460297) q[0];
x q[1];
rz(0.33863314) q[2];
sx q[2];
rz(-1.0446608) q[2];
sx q[2];
rz(-0.76802415) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.3396158) q[1];
sx q[1];
rz(-1.3286137) q[1];
sx q[1];
rz(-2.1276317) q[1];
x q[2];
rz(-2.8434136) q[3];
sx q[3];
rz(-2.0070873) q[3];
sx q[3];
rz(2.0661708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.1699367) q[2];
sx q[2];
rz(-1.1145096) q[2];
sx q[2];
rz(-1.2706832) q[2];
rz(3.0830248) q[3];
sx q[3];
rz(-1.6978369) q[3];
sx q[3];
rz(0.34846714) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.612065) q[0];
sx q[0];
rz(-0.53851524) q[0];
sx q[0];
rz(-2.8756323) q[0];
rz(2.3174441) q[1];
sx q[1];
rz(-1.8310841) q[1];
sx q[1];
rz(1.6110274) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8088246) q[0];
sx q[0];
rz(-2.7746189) q[0];
sx q[0];
rz(0.76353188) q[0];
x q[1];
rz(-1.62362) q[2];
sx q[2];
rz(-0.84047707) q[2];
sx q[2];
rz(-2.2234302) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.33324533) q[1];
sx q[1];
rz(-0.63831224) q[1];
sx q[1];
rz(2.6898333) q[1];
rz(-pi) q[2];
rz(-2.7830151) q[3];
sx q[3];
rz(-2.3365575) q[3];
sx q[3];
rz(-0.81881675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.019235762) q[2];
sx q[2];
rz(-1.1160206) q[2];
sx q[2];
rz(2.1717211) q[2];
rz(2.9979749) q[3];
sx q[3];
rz(-0.93846455) q[3];
sx q[3];
rz(-2.5041049) q[3];
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
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73760167) q[0];
sx q[0];
rz(-2.8856475) q[0];
sx q[0];
rz(3.0373489) q[0];
rz(2.3999816) q[1];
sx q[1];
rz(-0.18709083) q[1];
sx q[1];
rz(0.43824497) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61442626) q[0];
sx q[0];
rz(-2.2402555) q[0];
sx q[0];
rz(1.7015083) q[0];
rz(-pi) q[1];
rz(-0.23542685) q[2];
sx q[2];
rz(-0.98611438) q[2];
sx q[2];
rz(0.92008051) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.150369) q[1];
sx q[1];
rz(-2.5349778) q[1];
sx q[1];
rz(-1.6305417) q[1];
x q[2];
rz(-0.53536587) q[3];
sx q[3];
rz(-1.6965877) q[3];
sx q[3];
rz(2.4411502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7354342) q[2];
sx q[2];
rz(-2.8696852) q[2];
sx q[2];
rz(-0.3212277) q[2];
rz(2.1494703) q[3];
sx q[3];
rz(-2.0480053) q[3];
sx q[3];
rz(1.4120302) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.853249) q[0];
sx q[0];
rz(-0.11985954) q[0];
sx q[0];
rz(-2.993592) q[0];
rz(-2.0026228) q[1];
sx q[1];
rz(-0.6593467) q[1];
sx q[1];
rz(-2.151087) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5873224) q[0];
sx q[0];
rz(-2.0519216) q[0];
sx q[0];
rz(1.4863798) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1075145) q[2];
sx q[2];
rz(-2.2474237) q[2];
sx q[2];
rz(-2.7156242) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.79117486) q[1];
sx q[1];
rz(-1.5759849) q[1];
sx q[1];
rz(-1.6020011) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.93659525) q[3];
sx q[3];
rz(-2.5504125) q[3];
sx q[3];
rz(-1.9612966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.972435) q[2];
sx q[2];
rz(-0.7908228) q[2];
sx q[2];
rz(0.54231826) q[2];
rz(1.6769241) q[3];
sx q[3];
rz(-1.2277536) q[3];
sx q[3];
rz(0.98384583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9024314) q[0];
sx q[0];
rz(-0.79817525) q[0];
sx q[0];
rz(2.8620128) q[0];
rz(-0.78508776) q[1];
sx q[1];
rz(-0.79477349) q[1];
sx q[1];
rz(0.40207544) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.94984834) q[0];
sx q[0];
rz(-1.254891) q[0];
sx q[0];
rz(2.6788968) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.80569141) q[2];
sx q[2];
rz(-2.1955954) q[2];
sx q[2];
rz(-0.8936409) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.86568923) q[1];
sx q[1];
rz(-2.4498057) q[1];
sx q[1];
rz(-1.2671641) q[1];
rz(-pi) q[2];
rz(-1.9832595) q[3];
sx q[3];
rz(-1.8329046) q[3];
sx q[3];
rz(0.35221186) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.1606007) q[2];
sx q[2];
rz(-2.4805785) q[2];
sx q[2];
rz(-0.92657363) q[2];
rz(2.5911234) q[3];
sx q[3];
rz(-2.2678943) q[3];
sx q[3];
rz(-1.2800823) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4704623) q[0];
sx q[0];
rz(-1.067938) q[0];
sx q[0];
rz(-2.9965042) q[0];
rz(3.0628224) q[1];
sx q[1];
rz(-1.4046184) q[1];
sx q[1];
rz(1.2974993) q[1];
rz(1.5552022) q[2];
sx q[2];
rz(-0.77257665) q[2];
sx q[2];
rz(2.0333124) q[2];
rz(2.0918431) q[3];
sx q[3];
rz(-1.6853942) q[3];
sx q[3];
rz(-2.7484721) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
