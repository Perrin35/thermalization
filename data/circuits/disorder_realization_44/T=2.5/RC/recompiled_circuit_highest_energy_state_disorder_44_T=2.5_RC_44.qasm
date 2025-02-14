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
rz(-0.62701464) q[0];
sx q[0];
rz(3.7778683) q[0];
sx q[0];
rz(10.502622) q[0];
rz(1.0765422) q[1];
sx q[1];
rz(-0.62907469) q[1];
sx q[1];
rz(2.10973) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2159696) q[0];
sx q[0];
rz(-1.5711938) q[0];
sx q[0];
rz(0.0041603869) q[0];
rz(0.35030318) q[2];
sx q[2];
rz(-1.3418578) q[2];
sx q[2];
rz(0.13267429) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.061432211) q[1];
sx q[1];
rz(-2.5974063) q[1];
sx q[1];
rz(0.70744343) q[1];
rz(-pi) q[2];
rz(-2.9239465) q[3];
sx q[3];
rz(-2.7983694) q[3];
sx q[3];
rz(1.9453334) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.91118139) q[2];
sx q[2];
rz(-1.1831256) q[2];
sx q[2];
rz(-2.6739056) q[2];
rz(-2.6905401) q[3];
sx q[3];
rz(-0.39877287) q[3];
sx q[3];
rz(-0.27802813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69538799) q[0];
sx q[0];
rz(-0.22286335) q[0];
sx q[0];
rz(0.15431246) q[0];
rz(-0.34126869) q[1];
sx q[1];
rz(-2.7879265) q[1];
sx q[1];
rz(2.7214859) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.882928) q[0];
sx q[0];
rz(-0.9416343) q[0];
sx q[0];
rz(1.8101519) q[0];
x q[1];
rz(-0.046816512) q[2];
sx q[2];
rz(-2.280378) q[2];
sx q[2];
rz(3.0474718) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.5678011) q[1];
sx q[1];
rz(-0.47079362) q[1];
sx q[1];
rz(-0.0030229645) q[1];
rz(-pi) q[2];
rz(1.4147725) q[3];
sx q[3];
rz(-0.82956639) q[3];
sx q[3];
rz(1.6680731) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.5306065) q[2];
sx q[2];
rz(-0.89908081) q[2];
sx q[2];
rz(-0.77303028) q[2];
rz(2.2981339) q[3];
sx q[3];
rz(-1.5638899) q[3];
sx q[3];
rz(2.5696866) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17324363) q[0];
sx q[0];
rz(-2.1143715) q[0];
sx q[0];
rz(-0.20208836) q[0];
rz(-2.659722) q[1];
sx q[1];
rz(-2.9402969) q[1];
sx q[1];
rz(0.906382) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8608241) q[0];
sx q[0];
rz(-2.021179) q[0];
sx q[0];
rz(0.76669873) q[0];
rz(-1.1163122) q[2];
sx q[2];
rz(-1.2589728) q[2];
sx q[2];
rz(0.21909595) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.14531586) q[1];
sx q[1];
rz(-1.4106803) q[1];
sx q[1];
rz(2.5969863) q[1];
x q[2];
rz(1.2754457) q[3];
sx q[3];
rz(-1.7315627) q[3];
sx q[3];
rz(-2.8637342) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0704982) q[2];
sx q[2];
rz(-1.1134032) q[2];
sx q[2];
rz(-2.5615198) q[2];
rz(1.5813367) q[3];
sx q[3];
rz(-1.3908849) q[3];
sx q[3];
rz(1.4238547) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73109126) q[0];
sx q[0];
rz(-3.0803362) q[0];
sx q[0];
rz(-0.047792338) q[0];
rz(1.2740678) q[1];
sx q[1];
rz(-1.2893226) q[1];
sx q[1];
rz(3.0963669) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52818283) q[0];
sx q[0];
rz(-2.4827213) q[0];
sx q[0];
rz(-1.3263561) q[0];
rz(-pi) q[1];
rz(0.32856648) q[2];
sx q[2];
rz(-0.96146482) q[2];
sx q[2];
rz(-2.6549465) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8020683) q[1];
sx q[1];
rz(-1.5763723) q[1];
sx q[1];
rz(0.0030203621) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1158695) q[3];
sx q[3];
rz(-1.1909606) q[3];
sx q[3];
rz(-1.7120509) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.73146003) q[2];
sx q[2];
rz(-1.8803909) q[2];
sx q[2];
rz(0.88978466) q[2];
rz(-0.54640031) q[3];
sx q[3];
rz(-2.140464) q[3];
sx q[3];
rz(-0.31118292) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6379717) q[0];
sx q[0];
rz(-1.6872971) q[0];
sx q[0];
rz(1.4171492) q[0];
rz(-2.0218938) q[1];
sx q[1];
rz(-1.7599301) q[1];
sx q[1];
rz(-1.959257) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6851824) q[0];
sx q[0];
rz(-1.6401924) q[0];
sx q[0];
rz(2.7539537) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.117302) q[2];
sx q[2];
rz(-1.1536479) q[2];
sx q[2];
rz(0.81062775) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.4880138) q[1];
sx q[1];
rz(-1.6199448) q[1];
sx q[1];
rz(0.7213416) q[1];
x q[2];
rz(0.74966689) q[3];
sx q[3];
rz(-1.9056068) q[3];
sx q[3];
rz(2.1129169) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.22415367) q[2];
sx q[2];
rz(-2.3771693) q[2];
sx q[2];
rz(-2.7177641) q[2];
rz(2.0297) q[3];
sx q[3];
rz(-2.6424776) q[3];
sx q[3];
rz(-2.4901938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79913419) q[0];
sx q[0];
rz(-3.0693711) q[0];
sx q[0];
rz(2.1891731) q[0];
rz(0.77596387) q[1];
sx q[1];
rz(-2.2064078) q[1];
sx q[1];
rz(1.5752569) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1492406) q[0];
sx q[0];
rz(-1.5910205) q[0];
sx q[0];
rz(0.038073564) q[0];
rz(-pi) q[1];
x q[1];
rz(0.03869073) q[2];
sx q[2];
rz(-0.96444791) q[2];
sx q[2];
rz(2.4232466) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.9989062) q[1];
sx q[1];
rz(-2.493788) q[1];
sx q[1];
rz(2.7728989) q[1];
x q[2];
rz(1.0799584) q[3];
sx q[3];
rz(-1.7340525) q[3];
sx q[3];
rz(-1.8051683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.7427407) q[2];
sx q[2];
rz(-0.80790085) q[2];
sx q[2];
rz(-0.74637949) q[2];
rz(-1.0910723) q[3];
sx q[3];
rz(-1.4336136) q[3];
sx q[3];
rz(2.5892042) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93254507) q[0];
sx q[0];
rz(-0.046367558) q[0];
sx q[0];
rz(2.5982502) q[0];
rz(-0.53687334) q[1];
sx q[1];
rz(-2.8165292) q[1];
sx q[1];
rz(-3.0184025) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7167599) q[0];
sx q[0];
rz(-0.93532978) q[0];
sx q[0];
rz(3.1034971) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.67233147) q[2];
sx q[2];
rz(-0.62158442) q[2];
sx q[2];
rz(0.42447916) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(3.1063042) q[1];
sx q[1];
rz(-1.0266487) q[1];
sx q[1];
rz(2.418251) q[1];
rz(-pi) q[2];
rz(-0.79027883) q[3];
sx q[3];
rz(-2.1690282) q[3];
sx q[3];
rz(1.7733396) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.3049551) q[2];
sx q[2];
rz(-1.7540437) q[2];
sx q[2];
rz(-2.6564964) q[2];
rz(1.1525611) q[3];
sx q[3];
rz(-1.3644812) q[3];
sx q[3];
rz(3.0936354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3845859) q[0];
sx q[0];
rz(-2.204019) q[0];
sx q[0];
rz(-0.49569976) q[0];
rz(-3.0777625) q[1];
sx q[1];
rz(-0.7494691) q[1];
sx q[1];
rz(-0.071050342) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4249504) q[0];
sx q[0];
rz(-2.0293689) q[0];
sx q[0];
rz(1.3032662) q[0];
rz(-pi) q[1];
rz(2.536473) q[2];
sx q[2];
rz(-1.4265043) q[2];
sx q[2];
rz(2.6168106) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.7361765) q[1];
sx q[1];
rz(-1.725493) q[1];
sx q[1];
rz(1.3265893) q[1];
rz(-pi) q[2];
rz(-1.6556173) q[3];
sx q[3];
rz(-1.124998) q[3];
sx q[3];
rz(0.63157192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.0157328) q[2];
sx q[2];
rz(-1.4733529) q[2];
sx q[2];
rz(2.2802172) q[2];
rz(-1.3043978) q[3];
sx q[3];
rz(-0.574489) q[3];
sx q[3];
rz(1.5931574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(2.9651589) q[0];
sx q[0];
rz(-2.6506944) q[0];
sx q[0];
rz(-1.7023671) q[0];
rz(-3.0610415) q[1];
sx q[1];
rz(-1.4713902) q[1];
sx q[1];
rz(-1.7505987) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4347067) q[0];
sx q[0];
rz(-1.1483743) q[0];
sx q[0];
rz(-2.1377853) q[0];
rz(-pi) q[1];
rz(-1.8755781) q[2];
sx q[2];
rz(-1.551318) q[2];
sx q[2];
rz(1.2661042) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.7841107) q[1];
sx q[1];
rz(-2.3598089) q[1];
sx q[1];
rz(0.89872054) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7709522) q[3];
sx q[3];
rz(-1.9270867) q[3];
sx q[3];
rz(0.72539893) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.78802687) q[2];
sx q[2];
rz(-1.2949508) q[2];
sx q[2];
rz(-3.0143152) q[2];
rz(2.2150529) q[3];
sx q[3];
rz(-0.24181952) q[3];
sx q[3];
rz(-1.4827137) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38630286) q[0];
sx q[0];
rz(-0.64978623) q[0];
sx q[0];
rz(-1.5007098) q[0];
rz(0.81037784) q[1];
sx q[1];
rz(-2.2893298) q[1];
sx q[1];
rz(-2.3694029) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0836853) q[0];
sx q[0];
rz(-2.2026688) q[0];
sx q[0];
rz(2.4993012) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7842845) q[2];
sx q[2];
rz(-1.8597892) q[2];
sx q[2];
rz(-1.1505145) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.82622426) q[1];
sx q[1];
rz(-1.5591219) q[1];
sx q[1];
rz(1.0213901) q[1];
rz(-2.2507319) q[3];
sx q[3];
rz(-1.5308342) q[3];
sx q[3];
rz(3.0639632) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.7605674) q[2];
sx q[2];
rz(-0.80485359) q[2];
sx q[2];
rz(0.62925657) q[2];
rz(-2.4106846) q[3];
sx q[3];
rz(-2.2980502) q[3];
sx q[3];
rz(1.6985016) q[3];
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
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44169852) q[0];
sx q[0];
rz(-0.48342539) q[0];
sx q[0];
rz(-0.40406686) q[0];
rz(2.7370257) q[1];
sx q[1];
rz(-1.7351983) q[1];
sx q[1];
rz(2.4529967) q[1];
rz(-0.28180939) q[2];
sx q[2];
rz(-1.5055613) q[2];
sx q[2];
rz(0.67508634) q[2];
rz(2.6647207) q[3];
sx q[3];
rz(-0.58860368) q[3];
sx q[3];
rz(0.5640201) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
