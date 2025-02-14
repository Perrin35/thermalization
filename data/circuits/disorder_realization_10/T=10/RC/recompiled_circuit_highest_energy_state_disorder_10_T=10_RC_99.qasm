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
rz(-2.7918474) q[0];
sx q[0];
rz(-2.1748769) q[0];
sx q[0];
rz(2.7803687) q[0];
rz(-0.53520441) q[1];
sx q[1];
rz(-1.1517795) q[1];
sx q[1];
rz(1.8395543) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63473782) q[0];
sx q[0];
rz(-1.8950721) q[0];
sx q[0];
rz(1.8659331) q[0];
rz(-pi) q[1];
rz(0.7187219) q[2];
sx q[2];
rz(-1.5978004) q[2];
sx q[2];
rz(-0.57098963) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.9141317) q[1];
sx q[1];
rz(-1.9771076) q[1];
sx q[1];
rz(1.3106457) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8252772) q[3];
sx q[3];
rz(-0.22238734) q[3];
sx q[3];
rz(-0.60222821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5754622) q[2];
sx q[2];
rz(-2.8464816) q[2];
sx q[2];
rz(-2.5556514) q[2];
rz(1.5015886) q[3];
sx q[3];
rz(-1.8668886) q[3];
sx q[3];
rz(-1.9408102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5457299) q[0];
sx q[0];
rz(-1.0965309) q[0];
sx q[0];
rz(2.0134266) q[0];
rz(-2.3919171) q[1];
sx q[1];
rz(-1.3355037) q[1];
sx q[1];
rz(-1.658176) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5779232) q[0];
sx q[0];
rz(-2.2240608) q[0];
sx q[0];
rz(2.1640414) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.5107611) q[2];
sx q[2];
rz(-2.4756749) q[2];
sx q[2];
rz(-0.65424742) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-3.0771328) q[1];
sx q[1];
rz(-1.049548) q[1];
sx q[1];
rz(-1.8690994) q[1];
x q[2];
rz(0.97256487) q[3];
sx q[3];
rz(-0.56275193) q[3];
sx q[3];
rz(0.59845729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.7538606) q[2];
sx q[2];
rz(-2.3591177) q[2];
sx q[2];
rz(-2.1369047) q[2];
rz(0.0082155148) q[3];
sx q[3];
rz(-1.3722082) q[3];
sx q[3];
rz(-0.43732244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8179271) q[0];
sx q[0];
rz(-1.8765457) q[0];
sx q[0];
rz(-1.0311968) q[0];
rz(2.9464856) q[1];
sx q[1];
rz(-2.52774) q[1];
sx q[1];
rz(0.88417792) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4219727) q[0];
sx q[0];
rz(-0.0032711903) q[0];
sx q[0];
rz(-0.070912675) q[0];
rz(-pi) q[1];
rz(0.43135204) q[2];
sx q[2];
rz(-2.6425053) q[2];
sx q[2];
rz(-1.8519397) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.2078731) q[1];
sx q[1];
rz(-1.6657511) q[1];
sx q[1];
rz(1.1594561) q[1];
rz(-pi) q[2];
rz(2.1572491) q[3];
sx q[3];
rz(-1.801681) q[3];
sx q[3];
rz(-1.6572052) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.5174133) q[2];
sx q[2];
rz(-1.8884582) q[2];
sx q[2];
rz(-0.45428983) q[2];
rz(-1.0034466) q[3];
sx q[3];
rz(-0.61383057) q[3];
sx q[3];
rz(-1.412089) q[3];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91785947) q[0];
sx q[0];
rz(-0.35686019) q[0];
sx q[0];
rz(-0.80765635) q[0];
rz(1.3937021) q[1];
sx q[1];
rz(-1.423577) q[1];
sx q[1];
rz(1.7234939) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7634228) q[0];
sx q[0];
rz(-0.45158197) q[0];
sx q[0];
rz(-1.9096791) q[0];
x q[1];
rz(1.3803064) q[2];
sx q[2];
rz(-2.2675626) q[2];
sx q[2];
rz(-0.89542365) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.56546375) q[1];
sx q[1];
rz(-1.9525146) q[1];
sx q[1];
rz(1.5576499) q[1];
rz(0.40748458) q[3];
sx q[3];
rz(-2.0563246) q[3];
sx q[3];
rz(0.21903506) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7915446) q[2];
sx q[2];
rz(-1.9106661) q[2];
sx q[2];
rz(-3.0080504) q[2];
rz(-2.7327025) q[3];
sx q[3];
rz(-0.10433993) q[3];
sx q[3];
rz(-1.0094118) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.908476) q[0];
sx q[0];
rz(-2.3148843) q[0];
sx q[0];
rz(-2.6772461) q[0];
rz(-2.6404479) q[1];
sx q[1];
rz(-0.25096384) q[1];
sx q[1];
rz(0.73890013) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8221677) q[0];
sx q[0];
rz(-2.86045) q[0];
sx q[0];
rz(-2.7168363) q[0];
x q[1];
rz(-0.22245714) q[2];
sx q[2];
rz(-0.55402606) q[2];
sx q[2];
rz(1.2740327) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.0245783) q[1];
sx q[1];
rz(-2.0110333) q[1];
sx q[1];
rz(1.3280895) q[1];
rz(0.64844267) q[3];
sx q[3];
rz(-1.1243226) q[3];
sx q[3];
rz(2.1677599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.240694) q[2];
sx q[2];
rz(-1.7186761) q[2];
sx q[2];
rz(0.23441976) q[2];
rz(-2.8068986) q[3];
sx q[3];
rz(-0.60939279) q[3];
sx q[3];
rz(1.7083683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1471106) q[0];
sx q[0];
rz(-2.2587977) q[0];
sx q[0];
rz(2.2301646) q[0];
rz(-1.8270739) q[1];
sx q[1];
rz(-0.85845566) q[1];
sx q[1];
rz(1.7427157) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1778622) q[0];
sx q[0];
rz(-1.1060004) q[0];
sx q[0];
rz(-2.5739118) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6009459) q[2];
sx q[2];
rz(-1.4316403) q[2];
sx q[2];
rz(0.88983941) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.8533852) q[1];
sx q[1];
rz(-2.7010272) q[1];
sx q[1];
rz(1.3390932) q[1];
x q[2];
rz(0.26246385) q[3];
sx q[3];
rz(-1.0499716) q[3];
sx q[3];
rz(0.14674047) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.2167902) q[2];
sx q[2];
rz(-2.3239467) q[2];
sx q[2];
rz(1.7693046) q[2];
rz(1.4745845) q[3];
sx q[3];
rz(-1.0627012) q[3];
sx q[3];
rz(2.8180928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.737616) q[0];
sx q[0];
rz(-2.1196899) q[0];
sx q[0];
rz(-1.2603941) q[0];
rz(2.8003287) q[1];
sx q[1];
rz(-2.212846) q[1];
sx q[1];
rz(2.0492679) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72427801) q[0];
sx q[0];
rz(-1.5515741) q[0];
sx q[0];
rz(2.6249983) q[0];
rz(-pi) q[1];
rz(1.3767042) q[2];
sx q[2];
rz(-1.9067099) q[2];
sx q[2];
rz(-2.9011992) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.2637529) q[1];
sx q[1];
rz(-1.8108551) q[1];
sx q[1];
rz(-0.93024436) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.41938765) q[3];
sx q[3];
rz(-1.8598781) q[3];
sx q[3];
rz(-1.737864) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.1616538) q[2];
sx q[2];
rz(-1.3875763) q[2];
sx q[2];
rz(-2.2765344) q[2];
rz(3.0480399) q[3];
sx q[3];
rz(-1.4253987) q[3];
sx q[3];
rz(3.0430999) q[3];
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
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.426067) q[0];
sx q[0];
rz(-2.1676097) q[0];
sx q[0];
rz(1.6851715) q[0];
rz(2.9086225) q[1];
sx q[1];
rz(-1.3825682) q[1];
sx q[1];
rz(1.2196541) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9838966) q[0];
sx q[0];
rz(-0.86787631) q[0];
sx q[0];
rz(-1.9546175) q[0];
x q[1];
rz(-0.43841534) q[2];
sx q[2];
rz(-1.9389956) q[2];
sx q[2];
rz(-0.52939289) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-3.0387242) q[1];
sx q[1];
rz(-1.3967434) q[1];
sx q[1];
rz(2.7636488) q[1];
rz(-2.0601012) q[3];
sx q[3];
rz(-0.93554745) q[3];
sx q[3];
rz(2.7528544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.098027078) q[2];
sx q[2];
rz(-1.4412014) q[2];
sx q[2];
rz(-0.11858693) q[2];
rz(-2.3218527) q[3];
sx q[3];
rz(-0.44410646) q[3];
sx q[3];
rz(0.31258252) q[3];
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
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0176625) q[0];
sx q[0];
rz(-0.95159641) q[0];
sx q[0];
rz(0.96624017) q[0];
rz(0.36695925) q[1];
sx q[1];
rz(-2.1789813) q[1];
sx q[1];
rz(-2.5746094) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9241739) q[0];
sx q[0];
rz(-1.3710306) q[0];
sx q[0];
rz(0.49446797) q[0];
rz(-0.77436297) q[2];
sx q[2];
rz(-1.0972112) q[2];
sx q[2];
rz(0.87606664) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2046656) q[1];
sx q[1];
rz(-1.7742397) q[1];
sx q[1];
rz(2.4406747) q[1];
rz(2.7893355) q[3];
sx q[3];
rz(-2.2116419) q[3];
sx q[3];
rz(-0.0059222277) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.6734267) q[2];
sx q[2];
rz(-1.4665073) q[2];
sx q[2];
rz(1.0183938) q[2];
rz(-2.9112725) q[3];
sx q[3];
rz(-0.46263328) q[3];
sx q[3];
rz(-3.0751626) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51909834) q[0];
sx q[0];
rz(-1.007217) q[0];
sx q[0];
rz(1.8100716) q[0];
rz(-1.9271756) q[1];
sx q[1];
rz(-1.6969095) q[1];
sx q[1];
rz(-0.4164947) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.331401) q[0];
sx q[0];
rz(-1.7749359) q[0];
sx q[0];
rz(2.2603574) q[0];
rz(-pi) q[1];
rz(1.3168066) q[2];
sx q[2];
rz(-1.87623) q[2];
sx q[2];
rz(-1.0930201) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.70717274) q[1];
sx q[1];
rz(-1.0796629) q[1];
sx q[1];
rz(1.0855326) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3630886) q[3];
sx q[3];
rz(-0.45942222) q[3];
sx q[3];
rz(0.28693084) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.3930964) q[2];
sx q[2];
rz(-1.8668207) q[2];
sx q[2];
rz(1.3416802) q[2];
rz(1.120535) q[3];
sx q[3];
rz(-1.7399961) q[3];
sx q[3];
rz(0.11693461) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31711598) q[0];
sx q[0];
rz(-1.308029) q[0];
sx q[0];
rz(-1.8608004) q[0];
rz(-0.9723797) q[1];
sx q[1];
rz(-2.0449816) q[1];
sx q[1];
rz(2.4349946) q[1];
rz(-0.35928753) q[2];
sx q[2];
rz(-2.4854598) q[2];
sx q[2];
rz(1.796915) q[2];
rz(-1.9119143) q[3];
sx q[3];
rz(-3.0451265) q[3];
sx q[3];
rz(1.9652386) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
