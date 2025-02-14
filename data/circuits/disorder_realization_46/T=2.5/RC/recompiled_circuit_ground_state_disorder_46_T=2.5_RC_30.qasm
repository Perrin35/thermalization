OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.0134861) q[0];
sx q[0];
rz(-0.82733893) q[0];
sx q[0];
rz(1.7537355) q[0];
rz(-2.2881621) q[1];
sx q[1];
rz(-2.7397459) q[1];
sx q[1];
rz(-2.0274577) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8134817) q[0];
sx q[0];
rz(-0.8719647) q[0];
sx q[0];
rz(-0.33150406) q[0];
rz(1.3783437) q[2];
sx q[2];
rz(-2.2666605) q[2];
sx q[2];
rz(-3.1163094) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.6889894) q[1];
sx q[1];
rz(-1.5905979) q[1];
sx q[1];
rz(-0.0048586998) q[1];
rz(-pi) q[2];
rz(2.6403342) q[3];
sx q[3];
rz(-1.8717531) q[3];
sx q[3];
rz(0.00084998849) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.9806597) q[2];
sx q[2];
rz(-2.7200343) q[2];
sx q[2];
rz(1.4483615) q[2];
rz(2.8990922) q[3];
sx q[3];
rz(-0.91553965) q[3];
sx q[3];
rz(1.8814794) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5777957) q[0];
sx q[0];
rz(-1.8255434) q[0];
sx q[0];
rz(2.1412361) q[0];
rz(-0.73161221) q[1];
sx q[1];
rz(-1.7121406) q[1];
sx q[1];
rz(-1.9614722) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6574616) q[0];
sx q[0];
rz(-0.13741446) q[0];
sx q[0];
rz(-2.4262587) q[0];
x q[1];
rz(-0.033887788) q[2];
sx q[2];
rz(-1.292997) q[2];
sx q[2];
rz(0.44521618) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.6669162) q[1];
sx q[1];
rz(-1.70494) q[1];
sx q[1];
rz(0.0018906126) q[1];
x q[2];
rz(-0.7970771) q[3];
sx q[3];
rz(-2.7555572) q[3];
sx q[3];
rz(-2.4222056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0483027) q[2];
sx q[2];
rz(-0.90128171) q[2];
sx q[2];
rz(-1.0789336) q[2];
rz(0.087873936) q[3];
sx q[3];
rz(-1.352997) q[3];
sx q[3];
rz(2.8972076) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74697772) q[0];
sx q[0];
rz(-1.1629539) q[0];
sx q[0];
rz(-1.278247) q[0];
rz(-2.6784015) q[1];
sx q[1];
rz(-1.3563503) q[1];
sx q[1];
rz(2.3203704) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.074479178) q[0];
sx q[0];
rz(-0.68148208) q[0];
sx q[0];
rz(-1.2740178) q[0];
rz(-pi) q[1];
rz(-2.98396) q[2];
sx q[2];
rz(-2.7660884) q[2];
sx q[2];
rz(-0.65345308) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.48135172) q[1];
sx q[1];
rz(-1.4235745) q[1];
sx q[1];
rz(2.6805356) q[1];
rz(-pi) q[2];
rz(-2.3825112) q[3];
sx q[3];
rz(-1.0831175) q[3];
sx q[3];
rz(-1.9239192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.897573) q[2];
sx q[2];
rz(-1.579318) q[2];
sx q[2];
rz(-1.0025586) q[2];
rz(-1.9133441) q[3];
sx q[3];
rz(-1.4017665) q[3];
sx q[3];
rz(1.7206515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(2.2380075) q[0];
sx q[0];
rz(-2.0915732) q[0];
sx q[0];
rz(-0.28582698) q[0];
rz(-0.26126513) q[1];
sx q[1];
rz(-1.8087872) q[1];
sx q[1];
rz(-0.030698311) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8809533) q[0];
sx q[0];
rz(-0.27361037) q[0];
sx q[0];
rz(-0.27283313) q[0];
x q[1];
rz(1.6692363) q[2];
sx q[2];
rz(-0.67924309) q[2];
sx q[2];
rz(-1.498865) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.0568119) q[1];
sx q[1];
rz(-1.1263945) q[1];
sx q[1];
rz(2.8612575) q[1];
x q[2];
rz(1.3251036) q[3];
sx q[3];
rz(-2.8188955) q[3];
sx q[3];
rz(0.58354577) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.3280481) q[2];
sx q[2];
rz(-1.8009461) q[2];
sx q[2];
rz(-0.0030585232) q[2];
rz(-2.3004153) q[3];
sx q[3];
rz(-2.5523461) q[3];
sx q[3];
rz(-0.35471788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
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
rz(2.4641814) q[0];
sx q[0];
rz(-0.84742904) q[0];
sx q[0];
rz(-1.4671951) q[0];
rz(1.5122308) q[1];
sx q[1];
rz(-0.96577516) q[1];
sx q[1];
rz(-0.95219749) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95785054) q[0];
sx q[0];
rz(-1.6410367) q[0];
sx q[0];
rz(1.5635256) q[0];
x q[1];
rz(1.1811851) q[2];
sx q[2];
rz(-2.2305373) q[2];
sx q[2];
rz(-2.8174741) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.044301) q[1];
sx q[1];
rz(-0.11137577) q[1];
sx q[1];
rz(-1.9457413) q[1];
rz(-pi) q[2];
rz(-0.59479338) q[3];
sx q[3];
rz(-2.1457971) q[3];
sx q[3];
rz(1.9146384) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.3516922) q[2];
sx q[2];
rz(-1.5396427) q[2];
sx q[2];
rz(0.44463012) q[2];
rz(2.6897258) q[3];
sx q[3];
rz(-0.63932747) q[3];
sx q[3];
rz(3.0795081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.708798) q[0];
sx q[0];
rz(-2.4169156) q[0];
sx q[0];
rz(-0.92025796) q[0];
rz(-2.1814749) q[1];
sx q[1];
rz(-1.3939539) q[1];
sx q[1];
rz(2.6422909) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0408142) q[0];
sx q[0];
rz(-1.4199323) q[0];
sx q[0];
rz(-0.34025451) q[0];
x q[1];
rz(0.22412207) q[2];
sx q[2];
rz(-2.4680228) q[2];
sx q[2];
rz(-2.3526255) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4319246) q[1];
sx q[1];
rz(-1.0076771) q[1];
sx q[1];
rz(1.9456882) q[1];
rz(-pi) q[2];
rz(0.74120993) q[3];
sx q[3];
rz(-1.831372) q[3];
sx q[3];
rz(-1.8678566) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.3906735) q[2];
sx q[2];
rz(-1.5385188) q[2];
sx q[2];
rz(-1.0323367) q[2];
rz(2.27683) q[3];
sx q[3];
rz(-1.4356177) q[3];
sx q[3];
rz(2.5800956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
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
rz(0.432935) q[0];
sx q[0];
rz(-1.9517887) q[0];
sx q[0];
rz(1.7565961) q[0];
rz(0.9224836) q[1];
sx q[1];
rz(-1.205516) q[1];
sx q[1];
rz(2.9965957) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7940878) q[0];
sx q[0];
rz(-1.4055168) q[0];
sx q[0];
rz(-1.7519622) q[0];
rz(-3.112193) q[2];
sx q[2];
rz(-2.2839632) q[2];
sx q[2];
rz(-2.5273163) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.089894421) q[1];
sx q[1];
rz(-1.6015918) q[1];
sx q[1];
rz(1.6440065) q[1];
rz(-pi) q[2];
rz(-2.5707158) q[3];
sx q[3];
rz(-2.2789189) q[3];
sx q[3];
rz(-2.8857846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.9292235) q[2];
sx q[2];
rz(-2.9085458) q[2];
sx q[2];
rz(1.2082427) q[2];
rz(-2.3182747) q[3];
sx q[3];
rz(-1.1813141) q[3];
sx q[3];
rz(1.7821144) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.061148297) q[0];
sx q[0];
rz(-2.3539982) q[0];
sx q[0];
rz(-2.0620692) q[0];
rz(-2.1615255) q[1];
sx q[1];
rz(-1.020224) q[1];
sx q[1];
rz(-3.0316839) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3138374) q[0];
sx q[0];
rz(-1.2775363) q[0];
sx q[0];
rz(0.37773962) q[0];
x q[1];
rz(1.5797516) q[2];
sx q[2];
rz(-1.8491462) q[2];
sx q[2];
rz(-3.0305733) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.3802783) q[1];
sx q[1];
rz(-1.8683234) q[1];
sx q[1];
rz(-2.6975836) q[1];
rz(1.0411383) q[3];
sx q[3];
rz(-1.3762489) q[3];
sx q[3];
rz(1.1916849) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.93426934) q[2];
sx q[2];
rz(-2.0685652) q[2];
sx q[2];
rz(0.70460021) q[2];
rz(-2.8568824) q[3];
sx q[3];
rz(-0.73211089) q[3];
sx q[3];
rz(0.43340161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-2.1339742) q[0];
sx q[0];
rz(-1.824279) q[0];
sx q[0];
rz(-2.2229069) q[0];
rz(-0.48565117) q[1];
sx q[1];
rz(-2.698027) q[1];
sx q[1];
rz(-1.1007016) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0328631) q[0];
sx q[0];
rz(-0.19685611) q[0];
sx q[0];
rz(-2.3545676) q[0];
x q[1];
rz(1.8728016) q[2];
sx q[2];
rz(-1.7922316) q[2];
sx q[2];
rz(-1.0332274) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4816669) q[1];
sx q[1];
rz(-1.9710014) q[1];
sx q[1];
rz(-3.0277962) q[1];
rz(-0.086035919) q[3];
sx q[3];
rz(-0.63736308) q[3];
sx q[3];
rz(-1.6364545) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.6653768) q[2];
sx q[2];
rz(-2.2874338) q[2];
sx q[2];
rz(-2.3821135) q[2];
rz(0.94799834) q[3];
sx q[3];
rz(-1.4576603) q[3];
sx q[3];
rz(-1.6215526) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
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
rz(-0.043902472) q[0];
sx q[0];
rz(-0.49711415) q[0];
sx q[0];
rz(0.33101606) q[0];
rz(-0.76200062) q[1];
sx q[1];
rz(-2.8490729) q[1];
sx q[1];
rz(-2.5133572) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8593438) q[0];
sx q[0];
rz(-2.4971136) q[0];
sx q[0];
rz(-2.6656239) q[0];
rz(2.8339001) q[2];
sx q[2];
rz(-1.3465704) q[2];
sx q[2];
rz(-0.7219519) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.64092161) q[1];
sx q[1];
rz(-2.1643049) q[1];
sx q[1];
rz(-1.9954084) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9167494) q[3];
sx q[3];
rz(-1.1802434) q[3];
sx q[3];
rz(-1.2929163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.021775333) q[2];
sx q[2];
rz(-2.6675318) q[2];
sx q[2];
rz(2.9277756) q[2];
rz(-2.1553433) q[3];
sx q[3];
rz(-1.8503559) q[3];
sx q[3];
rz(0.10828644) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4719791) q[0];
sx q[0];
rz(-1.4431974) q[0];
sx q[0];
rz(1.6786014) q[0];
rz(-1.9646473) q[1];
sx q[1];
rz(-1.6564449) q[1];
sx q[1];
rz(0.0080531837) q[1];
rz(-1.0604924) q[2];
sx q[2];
rz(-2.7955187) q[2];
sx q[2];
rz(-2.1375755) q[2];
rz(0.10569345) q[3];
sx q[3];
rz(-2.4873545) q[3];
sx q[3];
rz(2.8771709) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
