OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.7656443) q[0];
sx q[0];
rz(-0.41416895) q[0];
sx q[0];
rz(-0.8015269) q[0];
rz(-0.0030567788) q[1];
sx q[1];
rz(-0.86441511) q[1];
sx q[1];
rz(-0.094223082) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72529257) q[0];
sx q[0];
rz(-2.2803118) q[0];
sx q[0];
rz(-2.6363591) q[0];
rz(-1.8434486) q[2];
sx q[2];
rz(-1.5412207) q[2];
sx q[2];
rz(2.3675413) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.73477106) q[1];
sx q[1];
rz(-1.1048006) q[1];
sx q[1];
rz(0.32225208) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5729957) q[3];
sx q[3];
rz(-1.4181976) q[3];
sx q[3];
rz(-2.8650177) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.6039383) q[2];
sx q[2];
rz(-2.5014169) q[2];
sx q[2];
rz(-1.630265) q[2];
rz(0.18621914) q[3];
sx q[3];
rz(-0.43034601) q[3];
sx q[3];
rz(2.490624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
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
rz(0.52361012) q[0];
sx q[0];
rz(-1.1815434) q[0];
sx q[0];
rz(0.13667662) q[0];
rz(0.094820529) q[1];
sx q[1];
rz(-2.6780728) q[1];
sx q[1];
rz(-0.76900855) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.626132) q[0];
sx q[0];
rz(-0.40053408) q[0];
sx q[0];
rz(-2.787091) q[0];
rz(0.693635) q[2];
sx q[2];
rz(-1.2918378) q[2];
sx q[2];
rz(0.50194937) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.0817854) q[1];
sx q[1];
rz(-1.4084653) q[1];
sx q[1];
rz(-0.74682208) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3202928) q[3];
sx q[3];
rz(-0.13938306) q[3];
sx q[3];
rz(-1.2047307) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.7424221) q[2];
sx q[2];
rz(-2.6593282) q[2];
sx q[2];
rz(-1.7102309) q[2];
rz(0.4906022) q[3];
sx q[3];
rz(-0.49121818) q[3];
sx q[3];
rz(-0.77559364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(-1.8785777) q[0];
sx q[0];
rz(-2.7624625) q[0];
sx q[0];
rz(1.0947134) q[0];
rz(1.8393983) q[1];
sx q[1];
rz(-2.1523988) q[1];
sx q[1];
rz(-3.1375695) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72633171) q[0];
sx q[0];
rz(-1.7021966) q[0];
sx q[0];
rz(0.8432998) q[0];
rz(-2.3962069) q[2];
sx q[2];
rz(-0.89085397) q[2];
sx q[2];
rz(-0.58733515) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.0672863) q[1];
sx q[1];
rz(-1.3137215) q[1];
sx q[1];
rz(-1.6021614) q[1];
rz(-0.28272776) q[3];
sx q[3];
rz(-2.4077031) q[3];
sx q[3];
rz(-2.3343488) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.5912938) q[2];
sx q[2];
rz(-2.55547) q[2];
sx q[2];
rz(-2.6934521) q[2];
rz(0.48629931) q[3];
sx q[3];
rz(-2.2794162) q[3];
sx q[3];
rz(-2.015131) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9124741) q[0];
sx q[0];
rz(-1.7121226) q[0];
sx q[0];
rz(0.66977704) q[0];
rz(1.9122596) q[1];
sx q[1];
rz(-2.6826617) q[1];
sx q[1];
rz(2.5965447) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16825039) q[0];
sx q[0];
rz(-1.1118044) q[0];
sx q[0];
rz(-0.58990546) q[0];
rz(-pi) q[1];
rz(0.9093693) q[2];
sx q[2];
rz(-2.0880359) q[2];
sx q[2];
rz(1.6632207) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.4008444) q[1];
sx q[1];
rz(-1.5103589) q[1];
sx q[1];
rz(1.3558186) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3432076) q[3];
sx q[3];
rz(-1.2252062) q[3];
sx q[3];
rz(1.8963199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5547319) q[2];
sx q[2];
rz(-1.502159) q[2];
sx q[2];
rz(-2.9736605) q[2];
rz(1.2083017) q[3];
sx q[3];
rz(-0.30253634) q[3];
sx q[3];
rz(-1.4471853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5912882) q[0];
sx q[0];
rz(-1.2606324) q[0];
sx q[0];
rz(-0.0038797832) q[0];
rz(2.8344391) q[1];
sx q[1];
rz(-2.3852564) q[1];
sx q[1];
rz(2.602813) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8545495) q[0];
sx q[0];
rz(-1.7998621) q[0];
sx q[0];
rz(1.4255217) q[0];
rz(1.4319375) q[2];
sx q[2];
rz(-1.9966239) q[2];
sx q[2];
rz(2.6151997) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.6314124) q[1];
sx q[1];
rz(-0.69563991) q[1];
sx q[1];
rz(0.56711491) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6461824) q[3];
sx q[3];
rz(-0.79548478) q[3];
sx q[3];
rz(1.0496248) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.6946081) q[2];
sx q[2];
rz(-1.0993404) q[2];
sx q[2];
rz(1.851932) q[2];
rz(1.6192216) q[3];
sx q[3];
rz(-0.74519849) q[3];
sx q[3];
rz(-1.6600018) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30037844) q[0];
sx q[0];
rz(-2.2024246) q[0];
sx q[0];
rz(-0.10264957) q[0];
rz(-0.73219055) q[1];
sx q[1];
rz(-0.79473549) q[1];
sx q[1];
rz(-0.57714677) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6711376) q[0];
sx q[0];
rz(-1.5671434) q[0];
sx q[0];
rz(1.5994344) q[0];
rz(-2.5546165) q[2];
sx q[2];
rz(-1.4329301) q[2];
sx q[2];
rz(1.1816813) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.6070093) q[1];
sx q[1];
rz(-1.680169) q[1];
sx q[1];
rz(1.4916199) q[1];
rz(-pi) q[2];
x q[2];
rz(0.34586819) q[3];
sx q[3];
rz(-0.78639275) q[3];
sx q[3];
rz(-2.6707471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.010043667) q[2];
sx q[2];
rz(-2.4250344) q[2];
sx q[2];
rz(-1.6183759) q[2];
rz(-3.0736382) q[3];
sx q[3];
rz(-2.8165635) q[3];
sx q[3];
rz(-2.3006191) q[3];
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
x q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1286569) q[0];
sx q[0];
rz(-1.034863) q[0];
sx q[0];
rz(-0.77142429) q[0];
rz(-0.5386638) q[1];
sx q[1];
rz(-1.795105) q[1];
sx q[1];
rz(1.0661941) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26771256) q[0];
sx q[0];
rz(-1.4940049) q[0];
sx q[0];
rz(1.8241054) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.23287878) q[2];
sx q[2];
rz(-1.0693197) q[2];
sx q[2];
rz(-1.5053144) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9061588) q[1];
sx q[1];
rz(-1.3536436) q[1];
sx q[1];
rz(1.557334) q[1];
x q[2];
rz(-1.3049502) q[3];
sx q[3];
rz(-1.2472594) q[3];
sx q[3];
rz(-0.27184799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.6923339) q[2];
sx q[2];
rz(-2.3350495) q[2];
sx q[2];
rz(0.44245693) q[2];
rz(-0.19065204) q[3];
sx q[3];
rz(-0.72398829) q[3];
sx q[3];
rz(0.75907069) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8964748) q[0];
sx q[0];
rz(-2.9560095) q[0];
sx q[0];
rz(1.3099439) q[0];
rz(-0.39778057) q[1];
sx q[1];
rz(-2.3186627) q[1];
sx q[1];
rz(1.3936874) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3041954) q[0];
sx q[0];
rz(-1.0365573) q[0];
sx q[0];
rz(-2.6287352) q[0];
rz(-pi) q[1];
rz(1.2692503) q[2];
sx q[2];
rz(-1.3450977) q[2];
sx q[2];
rz(-1.3653473) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.9164394) q[1];
sx q[1];
rz(-1.8691366) q[1];
sx q[1];
rz(2.6544996) q[1];
x q[2];
rz(-3.0354064) q[3];
sx q[3];
rz(-1.7683889) q[3];
sx q[3];
rz(-1.9492409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.028367793) q[2];
sx q[2];
rz(-0.71030474) q[2];
sx q[2];
rz(-3.0681211) q[2];
rz(2.4469589) q[3];
sx q[3];
rz(-1.2725384) q[3];
sx q[3];
rz(-2.1139483) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6756814) q[0];
sx q[0];
rz(-0.2008734) q[0];
sx q[0];
rz(-0.95941108) q[0];
rz(-0.77964669) q[1];
sx q[1];
rz(-2.1387073) q[1];
sx q[1];
rz(2.8693105) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.773943) q[0];
sx q[0];
rz(-3.127357) q[0];
sx q[0];
rz(-1.1971426) q[0];
x q[1];
rz(2.831976) q[2];
sx q[2];
rz(-2.2520492) q[2];
sx q[2];
rz(-2.3818784) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.6710558) q[1];
sx q[1];
rz(-1.9771641) q[1];
sx q[1];
rz(-2.0996835) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.168448) q[3];
sx q[3];
rz(-2.2242507) q[3];
sx q[3];
rz(-2.0367095) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4043364) q[2];
sx q[2];
rz(-1.5055089) q[2];
sx q[2];
rz(2.9340202) q[2];
rz(-0.10457822) q[3];
sx q[3];
rz(-0.50555491) q[3];
sx q[3];
rz(1.6149717) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4614457) q[0];
sx q[0];
rz(-1.9141645) q[0];
sx q[0];
rz(-2.0848059) q[0];
rz(-2.5959004) q[1];
sx q[1];
rz(-0.97815198) q[1];
sx q[1];
rz(-2.812885) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.064081505) q[0];
sx q[0];
rz(-1.9010549) q[0];
sx q[0];
rz(1.0868549) q[0];
rz(-3.0691807) q[2];
sx q[2];
rz(-1.6254404) q[2];
sx q[2];
rz(1.3018198) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.199177) q[1];
sx q[1];
rz(-2.7991382) q[1];
sx q[1];
rz(1.5210995) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2966818) q[3];
sx q[3];
rz(-2.3233827) q[3];
sx q[3];
rz(2.3880098) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1591961) q[2];
sx q[2];
rz(-3.0089162) q[2];
sx q[2];
rz(-1.1245493) q[2];
rz(-0.076796181) q[3];
sx q[3];
rz(-0.43282893) q[3];
sx q[3];
rz(2.2176149) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15976739) q[0];
sx q[0];
rz(-1.2564909) q[0];
sx q[0];
rz(1.5782574) q[0];
rz(-0.81623296) q[1];
sx q[1];
rz(-1.7385794) q[1];
sx q[1];
rz(2.0508918) q[1];
rz(-0.53173595) q[2];
sx q[2];
rz(-0.81960631) q[2];
sx q[2];
rz(-1.3255957) q[2];
rz(-2.7362551) q[3];
sx q[3];
rz(-1.5673076) q[3];
sx q[3];
rz(0.81893541) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
