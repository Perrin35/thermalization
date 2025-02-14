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
rz(4.7293651) q[0];
sx q[0];
rz(11.601996) q[0];
rz(-2.7388465) q[1];
sx q[1];
rz(-1.4637113) q[1];
sx q[1];
rz(0.95199624) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0453934) q[0];
sx q[0];
rz(-2.0598167) q[0];
sx q[0];
rz(-3.0351991) q[0];
rz(-pi) q[1];
rz(0.098564221) q[2];
sx q[2];
rz(-1.1505039) q[2];
sx q[2];
rz(-1.9914876) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.0178899) q[1];
sx q[1];
rz(-0.54819878) q[1];
sx q[1];
rz(-2.3831559) q[1];
x q[2];
rz(2.3856569) q[3];
sx q[3];
rz(-1.2159676) q[3];
sx q[3];
rz(-0.95099041) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6160148) q[2];
sx q[2];
rz(-1.3495477) q[2];
sx q[2];
rz(-0.40526059) q[2];
rz(-2.0112093) q[3];
sx q[3];
rz(-0.78961343) q[3];
sx q[3];
rz(-1.9270886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2363215) q[0];
sx q[0];
rz(-3.1089678) q[0];
sx q[0];
rz(2.3765833) q[0];
rz(-2.8969823) q[1];
sx q[1];
rz(-1.2449934) q[1];
sx q[1];
rz(-2.5027221) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9332704) q[0];
sx q[0];
rz(-2.2148085) q[0];
sx q[0];
rz(-1.6025376) q[0];
x q[1];
rz(-0.62221726) q[2];
sx q[2];
rz(-1.2583311) q[2];
sx q[2];
rz(-1.086233) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.42145211) q[1];
sx q[1];
rz(-1.1126592) q[1];
sx q[1];
rz(1.1790183) q[1];
rz(-pi) q[2];
x q[2];
rz(0.83373006) q[3];
sx q[3];
rz(-1.7639188) q[3];
sx q[3];
rz(1.5451136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.8581533) q[2];
sx q[2];
rz(-2.6946113) q[2];
sx q[2];
rz(-0.24432527) q[2];
rz(1.1089995) q[3];
sx q[3];
rz(-1.2764443) q[3];
sx q[3];
rz(0.23787704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37794381) q[0];
sx q[0];
rz(-1.0841333) q[0];
sx q[0];
rz(-1.2373244) q[0];
rz(1.411865) q[1];
sx q[1];
rz(-0.4487764) q[1];
sx q[1];
rz(-1.3317187) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.004382172) q[0];
sx q[0];
rz(-1.4829652) q[0];
sx q[0];
rz(2.3669404) q[0];
x q[1];
rz(2.1151416) q[2];
sx q[2];
rz(-2.7371001) q[2];
sx q[2];
rz(-2.9112688) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.6159042) q[1];
sx q[1];
rz(-1.6262904) q[1];
sx q[1];
rz(-1.5374119) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5272433) q[3];
sx q[3];
rz(-2.099937) q[3];
sx q[3];
rz(1.9093311) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1639013) q[2];
sx q[2];
rz(-1.0670263) q[2];
sx q[2];
rz(-2.4816371) q[2];
rz(1.4466977) q[3];
sx q[3];
rz(-1.7132672) q[3];
sx q[3];
rz(-1.1032633) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8834943) q[0];
sx q[0];
rz(-2.3907008) q[0];
sx q[0];
rz(1.1335491) q[0];
rz(0.51620382) q[1];
sx q[1];
rz(-1.8323003) q[1];
sx q[1];
rz(2.5925327) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20401317) q[0];
sx q[0];
rz(-1.7147062) q[0];
sx q[0];
rz(1.4043441) q[0];
rz(-2.1604161) q[2];
sx q[2];
rz(-1.7038888) q[2];
sx q[2];
rz(0.99118628) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.8827277) q[1];
sx q[1];
rz(-2.8587864) q[1];
sx q[1];
rz(-2.905719) q[1];
rz(-pi) q[2];
rz(1.9677343) q[3];
sx q[3];
rz(-1.8324495) q[3];
sx q[3];
rz(1.4058827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-3.0893726) q[2];
sx q[2];
rz(-2.1286271) q[2];
sx q[2];
rz(-0.29903665) q[2];
rz(0.61797577) q[3];
sx q[3];
rz(-1.0733913) q[3];
sx q[3];
rz(1.7834024) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1842136) q[0];
sx q[0];
rz(-3.1035677) q[0];
sx q[0];
rz(-2.6755565) q[0];
rz(-2.274463) q[1];
sx q[1];
rz(-0.51141089) q[1];
sx q[1];
rz(-2.3066511) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83564312) q[0];
sx q[0];
rz(-1.5993274) q[0];
sx q[0];
rz(3.0898407) q[0];
x q[1];
rz(-1.4915799) q[2];
sx q[2];
rz(-1.4822373) q[2];
sx q[2];
rz(-1.7065602) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.38091936) q[1];
sx q[1];
rz(-1.5173418) q[1];
sx q[1];
rz(1.0905114) q[1];
x q[2];
rz(2.4132009) q[3];
sx q[3];
rz(-1.0755223) q[3];
sx q[3];
rz(0.033897696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.20778188) q[2];
sx q[2];
rz(-2.2489579) q[2];
sx q[2];
rz(-0.031522838) q[2];
rz(-2.3837714) q[3];
sx q[3];
rz(-1.1008681) q[3];
sx q[3];
rz(-0.56263629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9071478) q[0];
sx q[0];
rz(-1.3749159) q[0];
sx q[0];
rz(3.1100519) q[0];
rz(-0.08983287) q[1];
sx q[1];
rz(-0.49982163) q[1];
sx q[1];
rz(-2.8057742) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8528362) q[0];
sx q[0];
rz(-1.9246708) q[0];
sx q[0];
rz(1.920098) q[0];
rz(-2.8029595) q[2];
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
rz(-1.5427067) q[1];
sx q[1];
rz(-2.5395091) q[1];
sx q[1];
rz(-2.00804) q[1];
x q[2];
rz(-2.1330413) q[3];
sx q[3];
rz(-2.6186071) q[3];
sx q[3];
rz(2.6949835) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.1699367) q[2];
sx q[2];
rz(-2.0270831) q[2];
sx q[2];
rz(-1.8709095) q[2];
rz(-0.058567889) q[3];
sx q[3];
rz(-1.4437557) q[3];
sx q[3];
rz(-0.34846714) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5295277) q[0];
sx q[0];
rz(-2.6030774) q[0];
sx q[0];
rz(2.8756323) q[0];
rz(-2.3174441) q[1];
sx q[1];
rz(-1.3105086) q[1];
sx q[1];
rz(1.6110274) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3327681) q[0];
sx q[0];
rz(-2.7746189) q[0];
sx q[0];
rz(-0.76353188) q[0];
rz(-pi) q[1];
x q[1];
rz(0.73101298) q[2];
sx q[2];
rz(-1.6101398) q[2];
sx q[2];
rz(0.61737663) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.33324533) q[1];
sx q[1];
rz(-2.5032804) q[1];
sx q[1];
rz(2.6898333) q[1];
x q[2];
rz(-1.2208185) q[3];
sx q[3];
rz(-0.82982291) q[3];
sx q[3];
rz(-0.32311026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(3.1223569) q[2];
sx q[2];
rz(-1.1160206) q[2];
sx q[2];
rz(-0.96987152) q[2];
rz(2.9979749) q[3];
sx q[3];
rz(-2.2031281) q[3];
sx q[3];
rz(2.5041049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73760167) q[0];
sx q[0];
rz(-0.25594512) q[0];
sx q[0];
rz(3.0373489) q[0];
rz(0.741611) q[1];
sx q[1];
rz(-0.18709083) q[1];
sx q[1];
rz(2.7033477) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1038228) q[0];
sx q[0];
rz(-1.4684104) q[0];
sx q[0];
rz(-0.67363221) q[0];
x q[1];
rz(-0.23542685) q[2];
sx q[2];
rz(-2.1554783) q[2];
sx q[2];
rz(-0.92008051) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2230473) q[1];
sx q[1];
rz(-2.1761736) q[1];
sx q[1];
rz(3.1001841) q[1];
rz(-pi) q[2];
rz(2.6062268) q[3];
sx q[3];
rz(-1.4450049) q[3];
sx q[3];
rz(-2.4411502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7354342) q[2];
sx q[2];
rz(-2.8696852) q[2];
sx q[2];
rz(0.3212277) q[2];
rz(-2.1494703) q[3];
sx q[3];
rz(-2.0480053) q[3];
sx q[3];
rz(1.7295624) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.853249) q[0];
sx q[0];
rz(-0.11985954) q[0];
sx q[0];
rz(2.993592) q[0];
rz(1.1389698) q[1];
sx q[1];
rz(-0.6593467) q[1];
sx q[1];
rz(0.99050561) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97738701) q[0];
sx q[0];
rz(-1.6456104) q[0];
sx q[0];
rz(-0.48259026) q[0];
x q[1];
rz(2.0340781) q[2];
sx q[2];
rz(-2.2474237) q[2];
sx q[2];
rz(2.7156242) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.79117486) q[1];
sx q[1];
rz(-1.5759849) q[1];
sx q[1];
rz(-1.5395916) q[1];
rz(2.7630291) q[3];
sx q[3];
rz(-2.0364015) q[3];
sx q[3];
rz(-1.2363124) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.972435) q[2];
sx q[2];
rz(-2.3507698) q[2];
sx q[2];
rz(2.5992744) q[2];
rz(-1.4646685) q[3];
sx q[3];
rz(-1.9138391) q[3];
sx q[3];
rz(-0.98384583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9024314) q[0];
sx q[0];
rz(-0.79817525) q[0];
sx q[0];
rz(-2.8620128) q[0];
rz(-2.3565049) q[1];
sx q[1];
rz(-2.3468192) q[1];
sx q[1];
rz(0.40207544) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94984834) q[0];
sx q[0];
rz(-1.254891) q[0];
sx q[0];
rz(-2.6788968) q[0];
x q[1];
rz(-0.78530757) q[2];
sx q[2];
rz(-0.97427893) q[2];
sx q[2];
rz(1.9526838) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.86568923) q[1];
sx q[1];
rz(-2.4498057) q[1];
sx q[1];
rz(-1.8744286) q[1];
rz(-pi) q[2];
rz(2.8567186) q[3];
sx q[3];
rz(-1.9683629) q[3];
sx q[3];
rz(-1.8101102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1606007) q[2];
sx q[2];
rz(-2.4805785) q[2];
sx q[2];
rz(-0.92657363) q[2];
rz(-2.5911234) q[3];
sx q[3];
rz(-0.87369839) q[3];
sx q[3];
rz(-1.2800823) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4704623) q[0];
sx q[0];
rz(-1.067938) q[0];
sx q[0];
rz(-2.9965042) q[0];
rz(-0.078770272) q[1];
sx q[1];
rz(-1.4046184) q[1];
sx q[1];
rz(1.2974993) q[1];
rz(-0.79828046) q[2];
sx q[2];
rz(-1.5816806) q[2];
sx q[2];
rz(0.47368373) q[2];
rz(-2.0918431) q[3];
sx q[3];
rz(-1.4561984) q[3];
sx q[3];
rz(0.39312058) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
