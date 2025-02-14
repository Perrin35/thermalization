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
rz(-1.929317) q[0];
sx q[0];
rz(-1.1317929) q[0];
sx q[0];
rz(-1.3728859) q[0];
rz(2.0139439) q[1];
sx q[1];
rz(-1.7665266) q[1];
sx q[1];
rz(2.3553203) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4117811) q[0];
sx q[0];
rz(-1.775911) q[0];
sx q[0];
rz(1.1954444) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.588618) q[2];
sx q[2];
rz(-2.1510501) q[2];
sx q[2];
rz(-2.0748823) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8189803) q[1];
sx q[1];
rz(-2.8507345) q[1];
sx q[1];
rz(1.4664654) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0065828) q[3];
sx q[3];
rz(-0.67970961) q[3];
sx q[3];
rz(-0.93384472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.11975153) q[2];
sx q[2];
rz(-1.4234875) q[2];
sx q[2];
rz(-1.2500259) q[2];
rz(0.73578468) q[3];
sx q[3];
rz(-2.9671228) q[3];
sx q[3];
rz(1.5031987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64330548) q[0];
sx q[0];
rz(-1.1730288) q[0];
sx q[0];
rz(-0.071320891) q[0];
rz(1.1530863) q[1];
sx q[1];
rz(-2.3801897) q[1];
sx q[1];
rz(0.52871314) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.081093719) q[0];
sx q[0];
rz(-1.1046263) q[0];
sx q[0];
rz(-2.7866298) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3350851) q[2];
sx q[2];
rz(-1.7545752) q[2];
sx q[2];
rz(1.9102131) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.9829927) q[1];
sx q[1];
rz(-2.0865177) q[1];
sx q[1];
rz(-0.92293585) q[1];
rz(-0.22457476) q[3];
sx q[3];
rz(-1.16515) q[3];
sx q[3];
rz(-0.39164513) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.66136393) q[2];
sx q[2];
rz(-2.7429548) q[2];
sx q[2];
rz(3.1129692) q[2];
rz(-0.35401595) q[3];
sx q[3];
rz(-2.386644) q[3];
sx q[3];
rz(-0.55907512) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4383168) q[0];
sx q[0];
rz(-1.0030614) q[0];
sx q[0];
rz(2.2502374) q[0];
rz(-0.50093961) q[1];
sx q[1];
rz(-1.5762065) q[1];
sx q[1];
rz(3.0827789) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56486216) q[0];
sx q[0];
rz(-1.65671) q[0];
sx q[0];
rz(1.6187173) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1040014) q[2];
sx q[2];
rz(-2.5639471) q[2];
sx q[2];
rz(1.6922598) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.14172983) q[1];
sx q[1];
rz(-2.0577621) q[1];
sx q[1];
rz(0.75090639) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9887588) q[3];
sx q[3];
rz(-0.18251576) q[3];
sx q[3];
rz(-2.07978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.7848876) q[2];
sx q[2];
rz(-2.5829743) q[2];
sx q[2];
rz(-2.9050262) q[2];
rz(-2.2208354) q[3];
sx q[3];
rz(-1.0509138) q[3];
sx q[3];
rz(-1.4111655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9193566) q[0];
sx q[0];
rz(-1.8574497) q[0];
sx q[0];
rz(0.87523571) q[0];
rz(-1.5454166) q[1];
sx q[1];
rz(-2.2794006) q[1];
sx q[1];
rz(0.045305591) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52462477) q[0];
sx q[0];
rz(-3.0201206) q[0];
sx q[0];
rz(-2.139702) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1854109) q[2];
sx q[2];
rz(-0.76757694) q[2];
sx q[2];
rz(-2.0337348) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.44573024) q[1];
sx q[1];
rz(-0.67477422) q[1];
sx q[1];
rz(1.4314907) q[1];
rz(-pi) q[2];
rz(2.153715) q[3];
sx q[3];
rz(-1.4081892) q[3];
sx q[3];
rz(0.041253003) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.8594325) q[2];
sx q[2];
rz(-2.1789447) q[2];
sx q[2];
rz(-1.1939987) q[2];
rz(2.2512839) q[3];
sx q[3];
rz(-1.9186391) q[3];
sx q[3];
rz(-2.1931026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0328261) q[0];
sx q[0];
rz(-2.325401) q[0];
sx q[0];
rz(-0.68242514) q[0];
rz(1.2884864) q[1];
sx q[1];
rz(-1.5623743) q[1];
sx q[1];
rz(1.8103745) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85859834) q[0];
sx q[0];
rz(-2.4147075) q[0];
sx q[0];
rz(-1.5461393) q[0];
rz(-2.8814828) q[2];
sx q[2];
rz(-1.1887728) q[2];
sx q[2];
rz(-1.1191238) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.23003787) q[1];
sx q[1];
rz(-1.7111519) q[1];
sx q[1];
rz(-0.36169238) q[1];
rz(1.0430452) q[3];
sx q[3];
rz(-1.966779) q[3];
sx q[3];
rz(0.92786232) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0266626) q[2];
sx q[2];
rz(-2.5422091) q[2];
sx q[2];
rz(0.66693532) q[2];
rz(2.5866348) q[3];
sx q[3];
rz(-0.68735492) q[3];
sx q[3];
rz(2.8426898) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32894593) q[0];
sx q[0];
rz(-1.8829367) q[0];
sx q[0];
rz(2.4378648) q[0];
rz(-0.16920371) q[1];
sx q[1];
rz(-1.5237944) q[1];
sx q[1];
rz(-2.9827548) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.279351) q[0];
sx q[0];
rz(-0.18687525) q[0];
sx q[0];
rz(-0.86743768) q[0];
rz(2.9877404) q[2];
sx q[2];
rz(-1.6922608) q[2];
sx q[2];
rz(2.6328994) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.3212686) q[1];
sx q[1];
rz(-1.0918198) q[1];
sx q[1];
rz(1.7198573) q[1];
rz(-pi) q[2];
rz(1.3424642) q[3];
sx q[3];
rz(-0.94551802) q[3];
sx q[3];
rz(-0.6354374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.3682897) q[2];
sx q[2];
rz(-0.94179073) q[2];
sx q[2];
rz(2.2840195) q[2];
rz(1.1693906) q[3];
sx q[3];
rz(-1.8604859) q[3];
sx q[3];
rz(2.9331971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17539772) q[0];
sx q[0];
rz(-2.3747787) q[0];
sx q[0];
rz(-0.71863693) q[0];
rz(-0.28193685) q[1];
sx q[1];
rz(-1.3749296) q[1];
sx q[1];
rz(0.66863543) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40067264) q[0];
sx q[0];
rz(-2.1869279) q[0];
sx q[0];
rz(0.4192062) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1988772) q[2];
sx q[2];
rz(-1.2758024) q[2];
sx q[2];
rz(1.6571972) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.1437415) q[1];
sx q[1];
rz(-2.6392548) q[1];
sx q[1];
rz(0.84059244) q[1];
rz(-pi) q[2];
rz(-1.578978) q[3];
sx q[3];
rz(-0.93957179) q[3];
sx q[3];
rz(-0.62731987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.64688993) q[2];
sx q[2];
rz(-1.9560445) q[2];
sx q[2];
rz(1.0364214) q[2];
rz(-0.63452619) q[3];
sx q[3];
rz(-2.1536638) q[3];
sx q[3];
rz(-0.79276597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-0.46102872) q[0];
sx q[0];
rz(-3.0240318) q[0];
sx q[0];
rz(2.4155937) q[0];
rz(1.1622102) q[1];
sx q[1];
rz(-1.9644968) q[1];
sx q[1];
rz(-1.9451709) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0045753) q[0];
sx q[0];
rz(-1.5992303) q[0];
sx q[0];
rz(1.397055) q[0];
rz(0.047646626) q[2];
sx q[2];
rz(-1.2638014) q[2];
sx q[2];
rz(0.36587151) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.4774306) q[1];
sx q[1];
rz(-0.71738418) q[1];
sx q[1];
rz(-0.75023164) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4713424) q[3];
sx q[3];
rz(-1.847451) q[3];
sx q[3];
rz(1.9346227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.90765816) q[2];
sx q[2];
rz(-1.5479167) q[2];
sx q[2];
rz(0.93351239) q[2];
rz(-2.524611) q[3];
sx q[3];
rz(-2.139293) q[3];
sx q[3];
rz(0.84701076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9762978) q[0];
sx q[0];
rz(-3.0369861) q[0];
sx q[0];
rz(-3.0883375) q[0];
rz(2.8540197) q[1];
sx q[1];
rz(-0.92929274) q[1];
sx q[1];
rz(-0.25921777) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0236146) q[0];
sx q[0];
rz(-3.0833457) q[0];
sx q[0];
rz(-1.144676) q[0];
x q[1];
rz(-2.0854112) q[2];
sx q[2];
rz(-2.8764203) q[2];
sx q[2];
rz(2.3725703) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.99154918) q[1];
sx q[1];
rz(-2.3880312) q[1];
sx q[1];
rz(1.6607453) q[1];
rz(-pi) q[2];
x q[2];
rz(0.76595347) q[3];
sx q[3];
rz(-1.4608835) q[3];
sx q[3];
rz(0.069086941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.5549434) q[2];
sx q[2];
rz(-1.2536851) q[2];
sx q[2];
rz(0.1813691) q[2];
rz(-0.3950611) q[3];
sx q[3];
rz(-2.919988) q[3];
sx q[3];
rz(0.980353) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7419389) q[0];
sx q[0];
rz(-1.5927277) q[0];
sx q[0];
rz(2.7594866) q[0];
rz(2.8451879) q[1];
sx q[1];
rz(-2.1370685) q[1];
sx q[1];
rz(1.872725) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.60703245) q[0];
sx q[0];
rz(-1.9968613) q[0];
sx q[0];
rz(-2.0924278) q[0];
rz(-pi) q[1];
rz(-1.4051132) q[2];
sx q[2];
rz(-0.68000092) q[2];
sx q[2];
rz(2.9939637) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.30291828) q[1];
sx q[1];
rz(-1.5494635) q[1];
sx q[1];
rz(-1.231074) q[1];
x q[2];
rz(-1.4488698) q[3];
sx q[3];
rz(-1.1200179) q[3];
sx q[3];
rz(-1.4565695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.2962013) q[2];
sx q[2];
rz(-2.5808344) q[2];
sx q[2];
rz(2.7517547) q[2];
rz(1.2975533) q[3];
sx q[3];
rz(-1.1417979) q[3];
sx q[3];
rz(0.84387422) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1525477) q[0];
sx q[0];
rz(-1.4023517) q[0];
sx q[0];
rz(-1.5040816) q[0];
rz(-1.7851495) q[1];
sx q[1];
rz(-2.103613) q[1];
sx q[1];
rz(1.6115859) q[1];
rz(2.6806954) q[2];
sx q[2];
rz(-0.97958889) q[2];
sx q[2];
rz(-0.19340672) q[2];
rz(0.33792321) q[3];
sx q[3];
rz(-2.3281964) q[3];
sx q[3];
rz(-2.7184814) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
