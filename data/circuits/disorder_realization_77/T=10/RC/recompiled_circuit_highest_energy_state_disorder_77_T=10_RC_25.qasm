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
rz(-0.51802975) q[0];
sx q[0];
rz(-2.0002444) q[0];
sx q[0];
rz(3.1014693) q[0];
rz(-2.8934381) q[1];
sx q[1];
rz(-2.0791972) q[1];
sx q[1];
rz(0.084029347) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18672189) q[0];
sx q[0];
rz(-1.6419517) q[0];
sx q[0];
rz(1.6297718) q[0];
rz(-pi) q[1];
rz(-1.749) q[2];
sx q[2];
rz(-1.9629109) q[2];
sx q[2];
rz(-0.41022476) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.89584589) q[1];
sx q[1];
rz(-1.7406643) q[1];
sx q[1];
rz(-2.6636811) q[1];
rz(-pi) q[2];
rz(-1.9840475) q[3];
sx q[3];
rz(-1.0344369) q[3];
sx q[3];
rz(2.9265917) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.4485126) q[2];
sx q[2];
rz(-1.7757519) q[2];
sx q[2];
rz(-2.8903294) q[2];
rz(-1.1723899) q[3];
sx q[3];
rz(-1.1029714) q[3];
sx q[3];
rz(-0.049886543) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4548816) q[0];
sx q[0];
rz(-0.48191163) q[0];
sx q[0];
rz(-1.3739817) q[0];
rz(-2.8107367) q[1];
sx q[1];
rz(-1.7287946) q[1];
sx q[1];
rz(-0.27423283) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2811063) q[0];
sx q[0];
rz(-1.8869068) q[0];
sx q[0];
rz(-1.4962026) q[0];
rz(-pi) q[1];
rz(-1.577967) q[2];
sx q[2];
rz(-1.0375377) q[2];
sx q[2];
rz(2.8523142) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3365605) q[1];
sx q[1];
rz(-0.39440409) q[1];
sx q[1];
rz(0.59275643) q[1];
rz(-1.7555439) q[3];
sx q[3];
rz(-0.42506252) q[3];
sx q[3];
rz(-0.63068542) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8809044) q[2];
sx q[2];
rz(-1.3284677) q[2];
sx q[2];
rz(1.0880967) q[2];
rz(1.043383) q[3];
sx q[3];
rz(-2.9731396) q[3];
sx q[3];
rz(2.8823631) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.008721) q[0];
sx q[0];
rz(-2.6368124) q[0];
sx q[0];
rz(2.0399427) q[0];
rz(-0.77611008) q[1];
sx q[1];
rz(-1.848369) q[1];
sx q[1];
rz(0.64220846) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5284536) q[0];
sx q[0];
rz(-2.0106843) q[0];
sx q[0];
rz(-0.47275193) q[0];
rz(-2.4218161) q[2];
sx q[2];
rz(-1.7584561) q[2];
sx q[2];
rz(-1.8328779) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.8485198) q[1];
sx q[1];
rz(-2.7024033) q[1];
sx q[1];
rz(1.3249012) q[1];
rz(-2.1884584) q[3];
sx q[3];
rz(-1.0079621) q[3];
sx q[3];
rz(1.1783534) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.8316101) q[2];
sx q[2];
rz(-1.1534561) q[2];
sx q[2];
rz(-0.93008053) q[2];
rz(-0.81405226) q[3];
sx q[3];
rz(-1.9508773) q[3];
sx q[3];
rz(1.9109776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.880068) q[0];
sx q[0];
rz(-1.9517169) q[0];
sx q[0];
rz(0.50450182) q[0];
rz(-2.2241459) q[1];
sx q[1];
rz(-0.73926273) q[1];
sx q[1];
rz(2.9333072) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84346164) q[0];
sx q[0];
rz(-1.4220474) q[0];
sx q[0];
rz(0.24555969) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4217626) q[2];
sx q[2];
rz(-1.7832758) q[2];
sx q[2];
rz(1.8393387) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.54305824) q[1];
sx q[1];
rz(-2.0428791) q[1];
sx q[1];
rz(-0.92453875) q[1];
rz(-pi) q[2];
rz(-1.3127906) q[3];
sx q[3];
rz(-1.0258706) q[3];
sx q[3];
rz(-1.8817924) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.31671277) q[2];
sx q[2];
rz(-1.0944288) q[2];
sx q[2];
rz(1.1701976) q[2];
rz(2.1757388) q[3];
sx q[3];
rz(-1.0713157) q[3];
sx q[3];
rz(-0.14537183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.589094) q[0];
sx q[0];
rz(-0.67411244) q[0];
sx q[0];
rz(-0.66459769) q[0];
rz(-2.873114) q[1];
sx q[1];
rz(-1.6842027) q[1];
sx q[1];
rz(2.1427515) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.474668) q[0];
sx q[0];
rz(-2.5380236) q[0];
sx q[0];
rz(-1.3003028) q[0];
x q[1];
rz(3.0078857) q[2];
sx q[2];
rz(-1.1503476) q[2];
sx q[2];
rz(-2.9131817) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.70390648) q[1];
sx q[1];
rz(-1.6900628) q[1];
sx q[1];
rz(1.7410623) q[1];
rz(-pi) q[2];
x q[2];
rz(3.1206297) q[3];
sx q[3];
rz(-0.87910324) q[3];
sx q[3];
rz(1.616459) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.9773679) q[2];
sx q[2];
rz(-0.94027844) q[2];
sx q[2];
rz(0.11717907) q[2];
rz(3.0818648) q[3];
sx q[3];
rz(-1.9220587) q[3];
sx q[3];
rz(-2.3265649) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3672459) q[0];
sx q[0];
rz(-0.4137488) q[0];
sx q[0];
rz(1.2600979) q[0];
rz(-0.11481181) q[1];
sx q[1];
rz(-1.7962619) q[1];
sx q[1];
rz(0.095349163) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0336825) q[0];
sx q[0];
rz(-2.9386407) q[0];
sx q[0];
rz(0.072921948) q[0];
x q[1];
rz(-1.2749407) q[2];
sx q[2];
rz(-0.5230147) q[2];
sx q[2];
rz(2.366334) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1667249) q[1];
sx q[1];
rz(-2.0151867) q[1];
sx q[1];
rz(-2.6728898) q[1];
x q[2];
rz(-1.5524149) q[3];
sx q[3];
rz(-0.52747969) q[3];
sx q[3];
rz(-0.35273409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.054333869) q[2];
sx q[2];
rz(-0.91780535) q[2];
sx q[2];
rz(2.4799662) q[2];
rz(-0.76964393) q[3];
sx q[3];
rz(-1.4504284) q[3];
sx q[3];
rz(2.2746287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
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
rz(1.0041644) q[0];
sx q[0];
rz(-0.34808174) q[0];
sx q[0];
rz(2.3837756) q[0];
rz(0.025253145) q[1];
sx q[1];
rz(-0.39552894) q[1];
sx q[1];
rz(1.1753561) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7987418) q[0];
sx q[0];
rz(-1.9254058) q[0];
sx q[0];
rz(-0.50574982) q[0];
rz(-pi) q[1];
rz(2.25039) q[2];
sx q[2];
rz(-1.1783349) q[2];
sx q[2];
rz(0.061937112) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.7832406) q[1];
sx q[1];
rz(-0.84642998) q[1];
sx q[1];
rz(-1.0388166) q[1];
x q[2];
rz(0.050042764) q[3];
sx q[3];
rz(-0.42946766) q[3];
sx q[3];
rz(-0.95145254) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.5777099) q[2];
sx q[2];
rz(-1.5241728) q[2];
sx q[2];
rz(-1.9263402) q[2];
rz(2.3828714) q[3];
sx q[3];
rz(-1.7251549) q[3];
sx q[3];
rz(-1.5569713) q[3];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6587081) q[0];
sx q[0];
rz(-1.8526798) q[0];
sx q[0];
rz(-2.7523852) q[0];
rz(0.088118531) q[1];
sx q[1];
rz(-1.4382818) q[1];
sx q[1];
rz(1.5315936) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54231468) q[0];
sx q[0];
rz(-1.8070343) q[0];
sx q[0];
rz(-1.3482984) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0247718) q[2];
sx q[2];
rz(-2.1810227) q[2];
sx q[2];
rz(-0.35200859) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.64634215) q[1];
sx q[1];
rz(-1.6227437) q[1];
sx q[1];
rz(2.6990128) q[1];
rz(-pi) q[2];
rz(2.76042) q[3];
sx q[3];
rz(-2.2444674) q[3];
sx q[3];
rz(-0.98667373) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.4072998) q[2];
sx q[2];
rz(-0.6137085) q[2];
sx q[2];
rz(1.8180004) q[2];
rz(-1.7981516) q[3];
sx q[3];
rz(-0.7898134) q[3];
sx q[3];
rz(0.33671236) q[3];
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
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6650498) q[0];
sx q[0];
rz(-2.2344868) q[0];
sx q[0];
rz(-2.6764349) q[0];
rz(0.05050412) q[1];
sx q[1];
rz(-2.2513794) q[1];
sx q[1];
rz(0.49096289) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57632724) q[0];
sx q[0];
rz(-2.6748228) q[0];
sx q[0];
rz(-0.76961036) q[0];
rz(-pi) q[1];
rz(0.30740909) q[2];
sx q[2];
rz(-1.4239862) q[2];
sx q[2];
rz(-1.6830144) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.0081353) q[1];
sx q[1];
rz(-2.244859) q[1];
sx q[1];
rz(0.73456711) q[1];
x q[2];
rz(2.773953) q[3];
sx q[3];
rz(-0.30045569) q[3];
sx q[3];
rz(-2.5641172) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.65840536) q[2];
sx q[2];
rz(-2.5453973) q[2];
sx q[2];
rz(-0.8636221) q[2];
rz(-0.18381271) q[3];
sx q[3];
rz(-2.176216) q[3];
sx q[3];
rz(-0.98627728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7154536) q[0];
sx q[0];
rz(-1.5197536) q[0];
sx q[0];
rz(-0.62582985) q[0];
rz(-2.4848056) q[1];
sx q[1];
rz(-2.1757809) q[1];
sx q[1];
rz(1.2084557) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0330808) q[0];
sx q[0];
rz(-1.7038561) q[0];
sx q[0];
rz(-3.0390374) q[0];
rz(-pi) q[1];
rz(2.7371128) q[2];
sx q[2];
rz(-0.70044971) q[2];
sx q[2];
rz(0.03534783) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.085026) q[1];
sx q[1];
rz(-1.657227) q[1];
sx q[1];
rz(1.5734476) q[1];
rz(-pi) q[2];
rz(1.3937065) q[3];
sx q[3];
rz(-1.0064126) q[3];
sx q[3];
rz(1.0799112) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.5586231) q[2];
sx q[2];
rz(-1.308459) q[2];
sx q[2];
rz(2.6226131) q[2];
rz(2.6988622) q[3];
sx q[3];
rz(-1.818592) q[3];
sx q[3];
rz(1.7980661) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9291572) q[0];
sx q[0];
rz(-0.99016187) q[0];
sx q[0];
rz(-1.1708175) q[0];
rz(2.9875372) q[1];
sx q[1];
rz(-0.93692056) q[1];
sx q[1];
rz(-0.9137203) q[1];
rz(-0.072617037) q[2];
sx q[2];
rz(-0.5177064) q[2];
sx q[2];
rz(0.032042423) q[2];
rz(2.2440425) q[3];
sx q[3];
rz(-1.8445476) q[3];
sx q[3];
rz(-2.6060819) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
