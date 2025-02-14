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
rz(0.25536594) q[0];
sx q[0];
rz(2.3998883) q[0];
sx q[0];
rz(11.381082) q[0];
rz(-1.2019295) q[1];
sx q[1];
rz(-2.1639106) q[1];
sx q[1];
rz(-0.77562195) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14549832) q[0];
sx q[0];
rz(-0.87255635) q[0];
sx q[0];
rz(-0.32693873) q[0];
rz(0.70574944) q[2];
sx q[2];
rz(-2.305491) q[2];
sx q[2];
rz(2.9388464) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.66167261) q[1];
sx q[1];
rz(-1.993517) q[1];
sx q[1];
rz(-1.8289966) q[1];
rz(-pi) q[2];
rz(1.5585639) q[3];
sx q[3];
rz(-2.1656007) q[3];
sx q[3];
rz(-0.4629688) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.97950196) q[2];
sx q[2];
rz(-0.68715954) q[2];
sx q[2];
rz(-3.0882351) q[2];
rz(-0.067367628) q[3];
sx q[3];
rz(-2.8991883) q[3];
sx q[3];
rz(-1.7072201) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.963221) q[0];
sx q[0];
rz(-2.0942978) q[0];
sx q[0];
rz(-2.6317327) q[0];
rz(1.1222104) q[1];
sx q[1];
rz(-2.0041859) q[1];
sx q[1];
rz(3.1013464) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48692856) q[0];
sx q[0];
rz(-1.860284) q[0];
sx q[0];
rz(-2.3280294) q[0];
rz(1.613125) q[2];
sx q[2];
rz(-2.5045103) q[2];
sx q[2];
rz(2.5893164) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.27850917) q[1];
sx q[1];
rz(-2.1467904) q[1];
sx q[1];
rz(1.7830677) q[1];
rz(-pi) q[2];
rz(-1.5533617) q[3];
sx q[3];
rz(-2.2477345) q[3];
sx q[3];
rz(-1.5172322) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.3736854) q[2];
sx q[2];
rz(-1.3482956) q[2];
sx q[2];
rz(0.56758869) q[2];
rz(-2.9034985) q[3];
sx q[3];
rz(-1.6425491) q[3];
sx q[3];
rz(-2.8897918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0955065) q[0];
sx q[0];
rz(-2.6987785) q[0];
sx q[0];
rz(-2.1777022) q[0];
rz(-1.5267728) q[1];
sx q[1];
rz(-1.9183728) q[1];
sx q[1];
rz(-2.8098106) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.362865) q[0];
sx q[0];
rz(-0.84072564) q[0];
sx q[0];
rz(-1.8018468) q[0];
rz(-pi) q[1];
x q[1];
rz(1.502723) q[2];
sx q[2];
rz(-0.34095665) q[2];
sx q[2];
rz(2.6457977) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0283734) q[1];
sx q[1];
rz(-3.0455596) q[1];
sx q[1];
rz(1.5043831) q[1];
x q[2];
rz(1.863722) q[3];
sx q[3];
rz(-2.6832242) q[3];
sx q[3];
rz(-1.019695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8277467) q[2];
sx q[2];
rz(-1.670994) q[2];
sx q[2];
rz(-0.016679114) q[2];
rz(-1.553675) q[3];
sx q[3];
rz(-1.1976306) q[3];
sx q[3];
rz(2.8488819) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68506771) q[0];
sx q[0];
rz(-0.86251384) q[0];
sx q[0];
rz(2.9271794) q[0];
rz(0.21687493) q[1];
sx q[1];
rz(-0.82019359) q[1];
sx q[1];
rz(2.0490501) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32347816) q[0];
sx q[0];
rz(-1.7351377) q[0];
sx q[0];
rz(-1.3413317) q[0];
x q[1];
rz(0.09437807) q[2];
sx q[2];
rz(-0.98446956) q[2];
sx q[2];
rz(-0.5273158) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.7937524) q[1];
sx q[1];
rz(-1.5303601) q[1];
sx q[1];
rz(-0.83918013) q[1];
rz(-pi) q[2];
rz(0.86534604) q[3];
sx q[3];
rz(-2.1095037) q[3];
sx q[3];
rz(-2.5821843) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.07) q[2];
sx q[2];
rz(-0.53539521) q[2];
sx q[2];
rz(0.52616057) q[2];
rz(-1.5633265) q[3];
sx q[3];
rz(-1.9901146) q[3];
sx q[3];
rz(2.9621942) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(0.18463317) q[0];
sx q[0];
rz(-1.3728091) q[0];
sx q[0];
rz(2.2233295) q[0];
rz(-2.1271162) q[1];
sx q[1];
rz(-2.4160073) q[1];
sx q[1];
rz(-2.0065506) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6776575) q[0];
sx q[0];
rz(-2.8615263) q[0];
sx q[0];
rz(2.8051462) q[0];
rz(-pi) q[1];
rz(-0.19164284) q[2];
sx q[2];
rz(-1.1517081) q[2];
sx q[2];
rz(-0.32391325) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.0088494) q[1];
sx q[1];
rz(-2.4970876) q[1];
sx q[1];
rz(-2.3185541) q[1];
rz(-pi) q[2];
rz(-0.5253864) q[3];
sx q[3];
rz(-0.30934096) q[3];
sx q[3];
rz(1.8745195) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.501118) q[2];
sx q[2];
rz(-1.9116118) q[2];
sx q[2];
rz(-0.099700363) q[2];
rz(0.68263188) q[3];
sx q[3];
rz(-1.3411025) q[3];
sx q[3];
rz(-2.7734267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30609983) q[0];
sx q[0];
rz(-1.2149128) q[0];
sx q[0];
rz(3.1021571) q[0];
rz(0.71190747) q[1];
sx q[1];
rz(-2.6075677) q[1];
sx q[1];
rz(1.3077516) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1354254) q[0];
sx q[0];
rz(-2.8949304) q[0];
sx q[0];
rz(-1.7390854) q[0];
x q[1];
rz(0.4591934) q[2];
sx q[2];
rz(-1.9373496) q[2];
sx q[2];
rz(-1.9286148) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.73346123) q[1];
sx q[1];
rz(-2.0234589) q[1];
sx q[1];
rz(-0.83475097) q[1];
rz(-pi) q[2];
rz(-0.51898308) q[3];
sx q[3];
rz(-0.95967889) q[3];
sx q[3];
rz(-0.42762363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.51440614) q[2];
sx q[2];
rz(-1.6986948) q[2];
sx q[2];
rz(0.31748873) q[2];
rz(-2.9182331) q[3];
sx q[3];
rz(-0.76627365) q[3];
sx q[3];
rz(-1.6754231) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8873972) q[0];
sx q[0];
rz(-1.6228209) q[0];
sx q[0];
rz(1.9018824) q[0];
rz(-2.8511035) q[1];
sx q[1];
rz(-1.5668642) q[1];
sx q[1];
rz(1.3715502) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1216095) q[0];
sx q[0];
rz(-1.0399315) q[0];
sx q[0];
rz(1.68501) q[0];
x q[1];
rz(2.7206037) q[2];
sx q[2];
rz(-1.1596646) q[2];
sx q[2];
rz(-0.56356424) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.54372245) q[1];
sx q[1];
rz(-0.67670435) q[1];
sx q[1];
rz(1.1110439) q[1];
rz(-pi) q[2];
rz(0.79381277) q[3];
sx q[3];
rz(-1.6598637) q[3];
sx q[3];
rz(1.7410914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.74662465) q[2];
sx q[2];
rz(-2.0396502) q[2];
sx q[2];
rz(1.9071707) q[2];
rz(-1.5359991) q[3];
sx q[3];
rz(-0.77080828) q[3];
sx q[3];
rz(-2.3710282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(-1.4602373) q[0];
sx q[0];
rz(-1.3139895) q[0];
sx q[0];
rz(2.6621057) q[0];
rz(1.208249) q[1];
sx q[1];
rz(-1.5346601) q[1];
sx q[1];
rz(-1.23034) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63060627) q[0];
sx q[0];
rz(-2.2842151) q[0];
sx q[0];
rz(-3.0632714) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6012457) q[2];
sx q[2];
rz(-1.4518988) q[2];
sx q[2];
rz(0.24248294) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8248511) q[1];
sx q[1];
rz(-1.6607641) q[1];
sx q[1];
rz(1.2895256) q[1];
rz(-pi) q[2];
x q[2];
rz(0.14982381) q[3];
sx q[3];
rz(-1.3725159) q[3];
sx q[3];
rz(1.0396569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.6464107) q[2];
sx q[2];
rz(-2.9718706) q[2];
sx q[2];
rz(1.0037615) q[2];
rz(2.0206001) q[3];
sx q[3];
rz(-1.5209773) q[3];
sx q[3];
rz(-2.3641172) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0178888) q[0];
sx q[0];
rz(-2.5214891) q[0];
sx q[0];
rz(1.7271127) q[0];
rz(-1.428528) q[1];
sx q[1];
rz(-2.2500549) q[1];
sx q[1];
rz(-3.0895272) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49567006) q[0];
sx q[0];
rz(-2.033253) q[0];
sx q[0];
rz(-2.511419) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.33627908) q[2];
sx q[2];
rz(-1.1999201) q[2];
sx q[2];
rz(-0.337643) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.3754943) q[1];
sx q[1];
rz(-1.3693989) q[1];
sx q[1];
rz(-0.010013568) q[1];
rz(-pi) q[2];
rz(-1.8776953) q[3];
sx q[3];
rz(-0.53987316) q[3];
sx q[3];
rz(2.8835874) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.97466126) q[2];
sx q[2];
rz(-0.31376803) q[2];
sx q[2];
rz(-1.3775919) q[2];
rz(-2.9396131) q[3];
sx q[3];
rz(-1.6986676) q[3];
sx q[3];
rz(2.0418237) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5438775) q[0];
sx q[0];
rz(-1.7622204) q[0];
sx q[0];
rz(-2.5482225) q[0];
rz(-1.6186391) q[1];
sx q[1];
rz(-2.7227089) q[1];
sx q[1];
rz(0.12817344) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0759715) q[0];
sx q[0];
rz(-0.21766454) q[0];
sx q[0];
rz(2.429921) q[0];
rz(-1.173557) q[2];
sx q[2];
rz(-2.4284869) q[2];
sx q[2];
rz(-1.8548059) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4306698) q[1];
sx q[1];
rz(-1.4928977) q[1];
sx q[1];
rz(1.0820612) q[1];
x q[2];
rz(1.0331328) q[3];
sx q[3];
rz(-2.3735316) q[3];
sx q[3];
rz(1.4466172) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0494978) q[2];
sx q[2];
rz(-2.8785661) q[2];
sx q[2];
rz(2.5035739) q[2];
rz(-1.1370031) q[3];
sx q[3];
rz(-0.59320366) q[3];
sx q[3];
rz(-3.1021281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20208246) q[0];
sx q[0];
rz(-1.408229) q[0];
sx q[0];
rz(-0.73233488) q[0];
rz(-2.9843075) q[1];
sx q[1];
rz(-1.597225) q[1];
sx q[1];
rz(1.0953915) q[1];
rz(-1.1284053) q[2];
sx q[2];
rz(-0.55991713) q[2];
sx q[2];
rz(-1.9922064) q[2];
rz(0.19098115) q[3];
sx q[3];
rz(-2.2195787) q[3];
sx q[3];
rz(1.7968404) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
