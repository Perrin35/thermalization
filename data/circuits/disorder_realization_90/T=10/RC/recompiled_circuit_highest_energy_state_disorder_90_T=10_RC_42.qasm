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
rz(0.69819063) q[0];
sx q[0];
rz(-0.31261045) q[0];
sx q[0];
rz(-0.99553776) q[0];
rz(-8.3182316) q[1];
sx q[1];
rz(3.9099524) q[1];
sx q[1];
rz(16.04019) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6512079) q[0];
sx q[0];
rz(-0.92726196) q[0];
sx q[0];
rz(-1.0585375) q[0];
rz(-1.3008437) q[2];
sx q[2];
rz(-2.6630262) q[2];
sx q[2];
rz(-1.1978483) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.7417833) q[1];
sx q[1];
rz(-0.64983778) q[1];
sx q[1];
rz(-2.3122961) q[1];
rz(-pi) q[2];
x q[2];
rz(0.99523441) q[3];
sx q[3];
rz(-2.0166631) q[3];
sx q[3];
rz(0.80328548) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.535061) q[2];
sx q[2];
rz(-2.1273095) q[2];
sx q[2];
rz(-0.052058546) q[2];
rz(1.082513) q[3];
sx q[3];
rz(-0.94729298) q[3];
sx q[3];
rz(1.989971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.624619) q[0];
sx q[0];
rz(-3.1095412) q[0];
sx q[0];
rz(-2.8042941) q[0];
rz(-0.15268046) q[1];
sx q[1];
rz(-0.76714194) q[1];
sx q[1];
rz(1.618636) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6628159) q[0];
sx q[0];
rz(-0.95774507) q[0];
sx q[0];
rz(1.7752035) q[0];
rz(-pi) q[1];
rz(2.5039429) q[2];
sx q[2];
rz(-1.7343688) q[2];
sx q[2];
rz(2.7113544) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.7584658) q[1];
sx q[1];
rz(-2.7039006) q[1];
sx q[1];
rz(-2.0630702) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4159774) q[3];
sx q[3];
rz(-0.972675) q[3];
sx q[3];
rz(0.036916669) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.0308257) q[2];
sx q[2];
rz(-0.36588565) q[2];
sx q[2];
rz(-2.7552674) q[2];
rz(-1.857916) q[3];
sx q[3];
rz(-1.9648896) q[3];
sx q[3];
rz(-1.9234689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.06685) q[0];
sx q[0];
rz(-1.0800986) q[0];
sx q[0];
rz(2.1225488) q[0];
rz(-2.3661803) q[1];
sx q[1];
rz(-1.8653899) q[1];
sx q[1];
rz(2.6554241) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81091034) q[0];
sx q[0];
rz(-0.031477246) q[0];
sx q[0];
rz(-0.093491836) q[0];
rz(-0.48392754) q[2];
sx q[2];
rz(-0.84887767) q[2];
sx q[2];
rz(-2.1579822) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(3.0224494) q[1];
sx q[1];
rz(-1.0876942) q[1];
sx q[1];
rz(1.8181605) q[1];
rz(-pi) q[2];
rz(2.4047732) q[3];
sx q[3];
rz(-0.42356682) q[3];
sx q[3];
rz(-1.119348) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7458618) q[2];
sx q[2];
rz(-1.1660601) q[2];
sx q[2];
rz(-0.80500785) q[2];
rz(-2.4896367) q[3];
sx q[3];
rz(-2.0396353) q[3];
sx q[3];
rz(2.2074047) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
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
rz(2.4790799) q[0];
sx q[0];
rz(-1.8204239) q[0];
sx q[0];
rz(-1.8030193) q[0];
rz(1.592912) q[1];
sx q[1];
rz(-0.70392307) q[1];
sx q[1];
rz(-0.42565638) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6451745) q[0];
sx q[0];
rz(-1.4789464) q[0];
sx q[0];
rz(1.5625815) q[0];
rz(-pi) q[1];
rz(1.4623619) q[2];
sx q[2];
rz(-0.30955704) q[2];
sx q[2];
rz(-2.5029687) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7158) q[1];
sx q[1];
rz(-1.1669901) q[1];
sx q[1];
rz(0.76622643) q[1];
x q[2];
rz(-2.4243746) q[3];
sx q[3];
rz(-1.7857496) q[3];
sx q[3];
rz(1.4357097) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.272133) q[2];
sx q[2];
rz(-1.3674066) q[2];
sx q[2];
rz(0.41932219) q[2];
rz(-1.9648633) q[3];
sx q[3];
rz(-0.71507016) q[3];
sx q[3];
rz(0.56593219) q[3];
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
x q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6856573) q[0];
sx q[0];
rz(-2.3929907) q[0];
sx q[0];
rz(-1.5022044) q[0];
rz(1.4337076) q[1];
sx q[1];
rz(-1.2805484) q[1];
sx q[1];
rz(-2.7388403) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1486738) q[0];
sx q[0];
rz(-0.67280992) q[0];
sx q[0];
rz(1.9901754) q[0];
rz(-pi) q[1];
rz(3.0482015) q[2];
sx q[2];
rz(-0.56677188) q[2];
sx q[2];
rz(0.10223481) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.6357543) q[1];
sx q[1];
rz(-0.16316667) q[1];
sx q[1];
rz(-2.5058305) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4536269) q[3];
sx q[3];
rz(-0.90725431) q[3];
sx q[3];
rz(-2.508916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.5002284) q[2];
sx q[2];
rz(-1.0470752) q[2];
sx q[2];
rz(0.25263986) q[2];
rz(-0.80444515) q[3];
sx q[3];
rz(-0.52144709) q[3];
sx q[3];
rz(0.49853244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8135524) q[0];
sx q[0];
rz(-3.033162) q[0];
sx q[0];
rz(0.88388467) q[0];
rz(-2.6416595) q[1];
sx q[1];
rz(-1.4686613) q[1];
sx q[1];
rz(2.6271741) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.894569) q[0];
sx q[0];
rz(-1.2366364) q[0];
sx q[0];
rz(2.1804125) q[0];
x q[1];
rz(-1.8158004) q[2];
sx q[2];
rz(-1.0520237) q[2];
sx q[2];
rz(0.10344783) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.0658256) q[1];
sx q[1];
rz(-1.3525241) q[1];
sx q[1];
rz(1.9820007) q[1];
rz(-pi) q[2];
rz(0.026192709) q[3];
sx q[3];
rz(-2.3520654) q[3];
sx q[3];
rz(1.7408371) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.9875235) q[2];
sx q[2];
rz(-1.46393) q[2];
sx q[2];
rz(3.0164111) q[2];
rz(0.82031885) q[3];
sx q[3];
rz(-1.9671974) q[3];
sx q[3];
rz(-0.34631795) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5888551) q[0];
sx q[0];
rz(-2.4947566) q[0];
sx q[0];
rz(0.11189017) q[0];
rz(2.0018068) q[1];
sx q[1];
rz(-2.9264989) q[1];
sx q[1];
rz(1.8591759) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8837785) q[0];
sx q[0];
rz(-2.1436267) q[0];
sx q[0];
rz(0.32578592) q[0];
rz(-pi) q[1];
rz(-1.7159903) q[2];
sx q[2];
rz(-1.8000126) q[2];
sx q[2];
rz(2.9289233) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.7805182) q[1];
sx q[1];
rz(-2.2782234) q[1];
sx q[1];
rz(-0.49496469) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5130643) q[3];
sx q[3];
rz(-2.2542076) q[3];
sx q[3];
rz(0.47799712) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.68975) q[2];
sx q[2];
rz(-2.8359783) q[2];
sx q[2];
rz(-1.8243194) q[2];
rz(3.1210323) q[3];
sx q[3];
rz(-2.8097184) q[3];
sx q[3];
rz(2.1338972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9356215) q[0];
sx q[0];
rz(-0.99140778) q[0];
sx q[0];
rz(2.8427065) q[0];
rz(1.0609974) q[1];
sx q[1];
rz(-1.3875049) q[1];
sx q[1];
rz(1.908173) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4087112) q[0];
sx q[0];
rz(-1.2378843) q[0];
sx q[0];
rz(-0.10245086) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.685254) q[2];
sx q[2];
rz(-1.5298163) q[2];
sx q[2];
rz(1.6270527) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.021344443) q[1];
sx q[1];
rz(-1.8665736) q[1];
sx q[1];
rz(-2.2463754) q[1];
rz(2.4806907) q[3];
sx q[3];
rz(-1.9993625) q[3];
sx q[3];
rz(-2.1235583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.3491106) q[2];
sx q[2];
rz(-1.0976378) q[2];
sx q[2];
rz(2.9343572) q[2];
rz(-2.4810897) q[3];
sx q[3];
rz(-2.1410172) q[3];
sx q[3];
rz(1.8628666) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5527363) q[0];
sx q[0];
rz(-1.739946) q[0];
sx q[0];
rz(1.2218342) q[0];
rz(2.6264722) q[1];
sx q[1];
rz(-0.11038596) q[1];
sx q[1];
rz(0.55330223) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3344655) q[0];
sx q[0];
rz(-0.57274023) q[0];
sx q[0];
rz(-1.80869) q[0];
x q[1];
rz(-3.0764941) q[2];
sx q[2];
rz(-1.1506478) q[2];
sx q[2];
rz(-1.9646027) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.8043704) q[1];
sx q[1];
rz(-2.9644199) q[1];
sx q[1];
rz(1.5371662) q[1];
rz(-pi) q[2];
rz(2.4724602) q[3];
sx q[3];
rz(-1.5520943) q[3];
sx q[3];
rz(1.6220055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.8052266) q[2];
sx q[2];
rz(-2.0311425) q[2];
sx q[2];
rz(2.6160348) q[2];
rz(-1.4793652) q[3];
sx q[3];
rz(-0.95366228) q[3];
sx q[3];
rz(0.76013887) q[3];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0568327) q[0];
sx q[0];
rz(-2.3469717) q[0];
sx q[0];
rz(2.7123465) q[0];
rz(1.6191354) q[1];
sx q[1];
rz(-1.8712964) q[1];
sx q[1];
rz(0.92990184) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7741755) q[0];
sx q[0];
rz(-0.40433592) q[0];
sx q[0];
rz(1.0003759) q[0];
rz(-pi) q[1];
rz(1.11082) q[2];
sx q[2];
rz(-2.8975652) q[2];
sx q[2];
rz(0.27966732) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.279544) q[1];
sx q[1];
rz(-2.8144208) q[1];
sx q[1];
rz(-2.0930834) q[1];
rz(1.6305805) q[3];
sx q[3];
rz(-2.1230704) q[3];
sx q[3];
rz(-0.58638257) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.4055736) q[2];
sx q[2];
rz(-2.6655727) q[2];
sx q[2];
rz(-0.19301566) q[2];
rz(3.0613101) q[3];
sx q[3];
rz(-0.86997) q[3];
sx q[3];
rz(-1.75753) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7667465) q[0];
sx q[0];
rz(-1.9369047) q[0];
sx q[0];
rz(1.9932224) q[0];
rz(-0.97024067) q[1];
sx q[1];
rz(-1.2087676) q[1];
sx q[1];
rz(-1.2420775) q[1];
rz(-2.316723) q[2];
sx q[2];
rz(-1.7879267) q[2];
sx q[2];
rz(0.59404165) q[2];
rz(2.8945782) q[3];
sx q[3];
rz(-3.044496) q[3];
sx q[3];
rz(-2.501957) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
