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
rz(0.99474466) q[0];
sx q[0];
rz(3.7288546) q[0];
sx q[0];
rz(10.050339) q[0];
rz(-2.5211531) q[1];
sx q[1];
rz(-0.99045366) q[1];
sx q[1];
rz(-0.78392309) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44255689) q[0];
sx q[0];
rz(-2.7074506) q[0];
sx q[0];
rz(-3.1331557) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6033635) q[2];
sx q[2];
rz(-2.9799649) q[2];
sx q[2];
rz(-1.5618351) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.82134889) q[1];
sx q[1];
rz(-0.87961266) q[1];
sx q[1];
rz(0.76637474) q[1];
rz(0.089139912) q[3];
sx q[3];
rz(-0.77971346) q[3];
sx q[3];
rz(-0.69621315) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.1938532) q[2];
sx q[2];
rz(-2.096602) q[2];
sx q[2];
rz(2.4742773) q[2];
rz(0.31467485) q[3];
sx q[3];
rz(-2.5421725) q[3];
sx q[3];
rz(2.2144894) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2941403) q[0];
sx q[0];
rz(-0.84646928) q[0];
sx q[0];
rz(-0.2997998) q[0];
rz(1.0046129) q[1];
sx q[1];
rz(-2.424898) q[1];
sx q[1];
rz(-2.2915548) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93545234) q[0];
sx q[0];
rz(-2.5295534) q[0];
sx q[0];
rz(-1.403359) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9710567) q[2];
sx q[2];
rz(-1.4922304) q[2];
sx q[2];
rz(-0.54135347) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.5932754) q[1];
sx q[1];
rz(-1.2139582) q[1];
sx q[1];
rz(0.22290454) q[1];
rz(-pi) q[2];
rz(-1.5570693) q[3];
sx q[3];
rz(-1.1662332) q[3];
sx q[3];
rz(1.8006397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.5726418) q[2];
sx q[2];
rz(-2.6093542) q[2];
sx q[2];
rz(-2.4786095) q[2];
rz(0.7524544) q[3];
sx q[3];
rz(-1.821725) q[3];
sx q[3];
rz(0.19645709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3982584) q[0];
sx q[0];
rz(-2.7404116) q[0];
sx q[0];
rz(2.4617526) q[0];
rz(2.4196449) q[1];
sx q[1];
rz(-0.55689055) q[1];
sx q[1];
rz(-0.78071761) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3328898) q[0];
sx q[0];
rz(-1.5631544) q[0];
sx q[0];
rz(1.1701258) q[0];
rz(-3.137792) q[2];
sx q[2];
rz(-1.6569306) q[2];
sx q[2];
rz(1.116556) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.36759427) q[1];
sx q[1];
rz(-1.9450608) q[1];
sx q[1];
rz(2.4623929) q[1];
x q[2];
rz(-2.9227681) q[3];
sx q[3];
rz(-0.84753643) q[3];
sx q[3];
rz(1.3506324) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.22689936) q[2];
sx q[2];
rz(-1.2279899) q[2];
sx q[2];
rz(2.4719888) q[2];
rz(0.90034825) q[3];
sx q[3];
rz(-1.0122296) q[3];
sx q[3];
rz(-2.3646234) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8825619) q[0];
sx q[0];
rz(-2.7029523) q[0];
sx q[0];
rz(-0.90748179) q[0];
rz(-2.5471845) q[1];
sx q[1];
rz(-0.92955697) q[1];
sx q[1];
rz(-1.1297191) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27760273) q[0];
sx q[0];
rz(-0.66501319) q[0];
sx q[0];
rz(-1.8825085) q[0];
x q[1];
rz(1.3200892) q[2];
sx q[2];
rz(-2.5959229) q[2];
sx q[2];
rz(-2.9755993) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.2767375) q[1];
sx q[1];
rz(-1.841479) q[1];
sx q[1];
rz(-0.60941831) q[1];
rz(-pi) q[2];
rz(1.8640737) q[3];
sx q[3];
rz(-1.3012816) q[3];
sx q[3];
rz(-0.5130918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-3.0323459) q[2];
sx q[2];
rz(-1.3999908) q[2];
sx q[2];
rz(2.4692811) q[2];
rz(-0.0191056) q[3];
sx q[3];
rz(-0.34135434) q[3];
sx q[3];
rz(3.0066971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30447176) q[0];
sx q[0];
rz(-2.3470375) q[0];
sx q[0];
rz(0.16803148) q[0];
rz(-1.2270323) q[1];
sx q[1];
rz(-1.2201759) q[1];
sx q[1];
rz(0.29774818) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15140238) q[0];
sx q[0];
rz(-2.08648) q[0];
sx q[0];
rz(-0.97008743) q[0];
rz(-0.57712675) q[2];
sx q[2];
rz(-1.071605) q[2];
sx q[2];
rz(-1.9124203) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.9492901) q[1];
sx q[1];
rz(-1.5956892) q[1];
sx q[1];
rz(-0.20583238) q[1];
rz(1.4120123) q[3];
sx q[3];
rz(-1.5701236) q[3];
sx q[3];
rz(-1.0345881) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.9510522) q[2];
sx q[2];
rz(-1.0995883) q[2];
sx q[2];
rz(0.61868787) q[2];
rz(0.80858532) q[3];
sx q[3];
rz(-0.4777258) q[3];
sx q[3];
rz(-1.8019069) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.488778) q[0];
sx q[0];
rz(-1.6365016) q[0];
sx q[0];
rz(-0.46522796) q[0];
rz(0.48509994) q[1];
sx q[1];
rz(-0.41050375) q[1];
sx q[1];
rz(3.0533275) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53118491) q[0];
sx q[0];
rz(-2.0636798) q[0];
sx q[0];
rz(2.6238074) q[0];
rz(3.1287304) q[2];
sx q[2];
rz(-2.2996827) q[2];
sx q[2];
rz(2.5358605) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.084296249) q[1];
sx q[1];
rz(-2.3688593) q[1];
sx q[1];
rz(2.904416) q[1];
rz(2.8507455) q[3];
sx q[3];
rz(-0.29971443) q[3];
sx q[3];
rz(-2.7951473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.39468592) q[2];
sx q[2];
rz(-2.2862819) q[2];
sx q[2];
rz(-0.49760231) q[2];
rz(-0.48007128) q[3];
sx q[3];
rz(-0.92104715) q[3];
sx q[3];
rz(3.0729496) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8467872) q[0];
sx q[0];
rz(-0.46600309) q[0];
sx q[0];
rz(-1.71126) q[0];
rz(-1.826674) q[1];
sx q[1];
rz(-2.7480835) q[1];
sx q[1];
rz(1.5550782) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4710812) q[0];
sx q[0];
rz(-1.1789315) q[0];
sx q[0];
rz(-3.0125822) q[0];
x q[1];
rz(-2.2847963) q[2];
sx q[2];
rz(-2.0371549) q[2];
sx q[2];
rz(0.12596065) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.3433229) q[1];
sx q[1];
rz(-2.1499691) q[1];
sx q[1];
rz(-0.12076853) q[1];
rz(-pi) q[2];
rz(-1.4659381) q[3];
sx q[3];
rz(-2.5075081) q[3];
sx q[3];
rz(0.88632562) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.2985208) q[2];
sx q[2];
rz(-0.71121794) q[2];
sx q[2];
rz(0.92740518) q[2];
rz(-1.1585506) q[3];
sx q[3];
rz(-2.4535593) q[3];
sx q[3];
rz(-0.096175171) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7189099) q[0];
sx q[0];
rz(-2.2490608) q[0];
sx q[0];
rz(0.476015) q[0];
rz(0.99126518) q[1];
sx q[1];
rz(-2.3920993) q[1];
sx q[1];
rz(3.057726) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0821486) q[0];
sx q[0];
rz(-2.5184439) q[0];
sx q[0];
rz(-0.92654689) q[0];
rz(-0.4617347) q[2];
sx q[2];
rz(-1.2424011) q[2];
sx q[2];
rz(0.4730313) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.8955651) q[1];
sx q[1];
rz(-1.152134) q[1];
sx q[1];
rz(3.0484867) q[1];
rz(-pi) q[2];
rz(-0.08633498) q[3];
sx q[3];
rz(-1.9342967) q[3];
sx q[3];
rz(-0.32237651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.91741651) q[2];
sx q[2];
rz(-2.8755964) q[2];
sx q[2];
rz(-2.4931397) q[2];
rz(2.2367541) q[3];
sx q[3];
rz(-1.6831393) q[3];
sx q[3];
rz(0.31074935) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34173486) q[0];
sx q[0];
rz(-0.71001995) q[0];
sx q[0];
rz(2.572686) q[0];
rz(-2.3078602) q[1];
sx q[1];
rz(-1.2833475) q[1];
sx q[1];
rz(2.1368829) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89812714) q[0];
sx q[0];
rz(-1.7633738) q[0];
sx q[0];
rz(-0.63611998) q[0];
rz(-0.23022454) q[2];
sx q[2];
rz(-1.0448928) q[2];
sx q[2];
rz(-2.3221225) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9819599) q[1];
sx q[1];
rz(-1.6430322) q[1];
sx q[1];
rz(2.2700929) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.58752693) q[3];
sx q[3];
rz(-1.5024795) q[3];
sx q[3];
rz(0.63623601) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.4363165) q[2];
sx q[2];
rz(-0.87203163) q[2];
sx q[2];
rz(0.3479859) q[2];
rz(0.82585382) q[3];
sx q[3];
rz(-2.9233942) q[3];
sx q[3];
rz(0.60514778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0572877) q[0];
sx q[0];
rz(-1.0658406) q[0];
sx q[0];
rz(2.9344946) q[0];
rz(-0.53889489) q[1];
sx q[1];
rz(-2.5205044) q[1];
sx q[1];
rz(2.5394687) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83610632) q[0];
sx q[0];
rz(-2.0417111) q[0];
sx q[0];
rz(-2.6013646) q[0];
x q[1];
rz(-1.1634356) q[2];
sx q[2];
rz(-2.1199787) q[2];
sx q[2];
rz(-0.34005766) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.6844897) q[1];
sx q[1];
rz(-2.7278363) q[1];
sx q[1];
rz(1.3828418) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8982417) q[3];
sx q[3];
rz(-2.2909938) q[3];
sx q[3];
rz(-0.44046381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3567317) q[2];
sx q[2];
rz(-0.95180231) q[2];
sx q[2];
rz(-1.3447364) q[2];
rz(2.6209659) q[3];
sx q[3];
rz(-0.27025637) q[3];
sx q[3];
rz(2.3394913) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6945334) q[0];
sx q[0];
rz(-2.4166528) q[0];
sx q[0];
rz(1.7550533) q[0];
rz(-2.3583892) q[1];
sx q[1];
rz(-1.422311) q[1];
sx q[1];
rz(1.5625988) q[1];
rz(0.65880792) q[2];
sx q[2];
rz(-0.42429069) q[2];
sx q[2];
rz(-1.9784603) q[2];
rz(0.24012031) q[3];
sx q[3];
rz(-1.5821531) q[3];
sx q[3];
rz(1.5021363) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
