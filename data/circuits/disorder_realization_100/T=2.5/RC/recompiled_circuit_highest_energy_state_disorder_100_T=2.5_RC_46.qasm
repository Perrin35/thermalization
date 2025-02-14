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
rz(2.9313791) q[0];
sx q[0];
rz(5.9271521) q[0];
sx q[0];
rz(7.9071101) q[0];
rz(1.8310945) q[1];
sx q[1];
rz(-2.2902391) q[1];
sx q[1];
rz(-0.89827615) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2668928) q[0];
sx q[0];
rz(-1.9277097) q[0];
sx q[0];
rz(0.92706417) q[0];
rz(0.29062985) q[2];
sx q[2];
rz(-2.0465474) q[2];
sx q[2];
rz(1.0889183) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.15816244) q[1];
sx q[1];
rz(-1.753141) q[1];
sx q[1];
rz(-1.4901864) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6226852) q[3];
sx q[3];
rz(-0.16434419) q[3];
sx q[3];
rz(2.5353894) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.5400759) q[2];
sx q[2];
rz(-2.8105152) q[2];
sx q[2];
rz(2.4430821) q[2];
rz(-1.4861594) q[3];
sx q[3];
rz(-0.58659068) q[3];
sx q[3];
rz(1.9894039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0450714) q[0];
sx q[0];
rz(-0.69913816) q[0];
sx q[0];
rz(-0.29749468) q[0];
rz(2.7104764) q[1];
sx q[1];
rz(-0.86687207) q[1];
sx q[1];
rz(-1.8284304) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0348059) q[0];
sx q[0];
rz(-1.2277831) q[0];
sx q[0];
rz(2.4136132) q[0];
rz(-1.1651016) q[2];
sx q[2];
rz(-1.5303474) q[2];
sx q[2];
rz(-1.2834398) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6152108) q[1];
sx q[1];
rz(-1.5632249) q[1];
sx q[1];
rz(0.79736276) q[1];
x q[2];
rz(-1.6974988) q[3];
sx q[3];
rz(-1.8167348) q[3];
sx q[3];
rz(1.9501231) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.3678652) q[2];
sx q[2];
rz(-0.44168681) q[2];
sx q[2];
rz(2.3750677) q[2];
rz(2.4567228) q[3];
sx q[3];
rz(-1.6133512) q[3];
sx q[3];
rz(-2.6510356) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(1.3016475) q[0];
sx q[0];
rz(-1.3293581) q[0];
sx q[0];
rz(-2.8369821) q[0];
rz(-2.1412762) q[1];
sx q[1];
rz(-2.0392923) q[1];
sx q[1];
rz(0.89835483) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2294126) q[0];
sx q[0];
rz(-2.1026975) q[0];
sx q[0];
rz(3.0454265) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8091891) q[2];
sx q[2];
rz(-0.89611182) q[2];
sx q[2];
rz(1.5522267) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5542986) q[1];
sx q[1];
rz(-2.0600494) q[1];
sx q[1];
rz(0.6502519) q[1];
rz(-pi) q[2];
rz(-1.0482084) q[3];
sx q[3];
rz(-1.4915803) q[3];
sx q[3];
rz(2.0780711) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.1343214) q[2];
sx q[2];
rz(-1.7223225) q[2];
sx q[2];
rz(2.814494) q[2];
rz(1.5056115) q[3];
sx q[3];
rz(-1.6749195) q[3];
sx q[3];
rz(1.8065642) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1148249) q[0];
sx q[0];
rz(-0.46634316) q[0];
sx q[0];
rz(-0.53792167) q[0];
rz(2.3665358) q[1];
sx q[1];
rz(-1.8102976) q[1];
sx q[1];
rz(-2.7937826) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8224768) q[0];
sx q[0];
rz(-0.70981423) q[0];
sx q[0];
rz(-0.7546203) q[0];
rz(1.4653019) q[2];
sx q[2];
rz(-1.6755793) q[2];
sx q[2];
rz(-1.3380761) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.2217933) q[1];
sx q[1];
rz(-1.9927653) q[1];
sx q[1];
rz(-0.74733644) q[1];
rz(-pi) q[2];
rz(-2.9628539) q[3];
sx q[3];
rz(-1.1832269) q[3];
sx q[3];
rz(-0.77690682) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9378918) q[2];
sx q[2];
rz(-2.4163279) q[2];
sx q[2];
rz(2.4791278) q[2];
rz(1.0558111) q[3];
sx q[3];
rz(-2.8231088) q[3];
sx q[3];
rz(-1.2663579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1556959) q[0];
sx q[0];
rz(-2.6031384) q[0];
sx q[0];
rz(1.023531) q[0];
rz(2.0899978) q[1];
sx q[1];
rz(-1.6513499) q[1];
sx q[1];
rz(-3.1255299) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4996289) q[0];
sx q[0];
rz(-1.0283386) q[0];
sx q[0];
rz(-2.9835975) q[0];
rz(2.4098924) q[2];
sx q[2];
rz(-0.83544399) q[2];
sx q[2];
rz(-1.4994061) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.2019883) q[1];
sx q[1];
rz(-1.3496205) q[1];
sx q[1];
rz(-2.6084234) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.315867) q[3];
sx q[3];
rz(-2.8615017) q[3];
sx q[3];
rz(0.051163604) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.84231275) q[2];
sx q[2];
rz(-0.88025847) q[2];
sx q[2];
rz(0.24097815) q[2];
rz(2.0321417) q[3];
sx q[3];
rz(-1.2883319) q[3];
sx q[3];
rz(0.39458767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1269162) q[0];
sx q[0];
rz(-0.83634818) q[0];
sx q[0];
rz(0.29689223) q[0];
rz(-1.9706005) q[1];
sx q[1];
rz(-1.3208656) q[1];
sx q[1];
rz(2.6062633) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7479582) q[0];
sx q[0];
rz(-1.7553333) q[0];
sx q[0];
rz(-0.18963295) q[0];
x q[1];
rz(-0.81670834) q[2];
sx q[2];
rz(-1.3544519) q[2];
sx q[2];
rz(2.3522481) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.6397093) q[1];
sx q[1];
rz(-2.0514666) q[1];
sx q[1];
rz(-0.23940803) q[1];
rz(-pi) q[2];
rz(-0.13958474) q[3];
sx q[3];
rz(-1.727316) q[3];
sx q[3];
rz(0.28790441) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.31763306) q[2];
sx q[2];
rz(-0.77482688) q[2];
sx q[2];
rz(0.63014692) q[2];
rz(1.0050425) q[3];
sx q[3];
rz(-0.21868394) q[3];
sx q[3];
rz(-0.69766831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9610577) q[0];
sx q[0];
rz(-3.0735885) q[0];
sx q[0];
rz(-1.2766174) q[0];
rz(-2.0969773) q[1];
sx q[1];
rz(-1.4199384) q[1];
sx q[1];
rz(-0.23342625) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0208527) q[0];
sx q[0];
rz(-2.4197612) q[0];
sx q[0];
rz(1.634725) q[0];
rz(-2.7415375) q[2];
sx q[2];
rz(-2.4597485) q[2];
sx q[2];
rz(-1.2712196) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.0835159) q[1];
sx q[1];
rz(-2.3928071) q[1];
sx q[1];
rz(0.056549055) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0393321) q[3];
sx q[3];
rz(-0.58057154) q[3];
sx q[3];
rz(-1.5619265) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.6650247) q[2];
sx q[2];
rz(-1.2580405) q[2];
sx q[2];
rz(-2.8705719) q[2];
rz(2.792231) q[3];
sx q[3];
rz(-0.91785279) q[3];
sx q[3];
rz(2.5640986) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2006328) q[0];
sx q[0];
rz(-0.9372434) q[0];
sx q[0];
rz(-1.3295133) q[0];
rz(-2.8726574) q[1];
sx q[1];
rz(-2.4766998) q[1];
sx q[1];
rz(1.7609133) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5912) q[0];
sx q[0];
rz(-1.368528) q[0];
sx q[0];
rz(1.8254542) q[0];
rz(-0.024928851) q[2];
sx q[2];
rz(-0.59621323) q[2];
sx q[2];
rz(-0.91170694) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.24522752) q[1];
sx q[1];
rz(-0.62313634) q[1];
sx q[1];
rz(-1.8521897) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4333357) q[3];
sx q[3];
rz(-2.1270184) q[3];
sx q[3];
rz(-1.463365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.55714503) q[2];
sx q[2];
rz(-0.86981213) q[2];
sx q[2];
rz(0.73428854) q[2];
rz(-2.0860489) q[3];
sx q[3];
rz(-1.5922092) q[3];
sx q[3];
rz(-0.9744823) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4771117) q[0];
sx q[0];
rz(-0.25150126) q[0];
sx q[0];
rz(2.2736881) q[0];
rz(-1.8382629) q[1];
sx q[1];
rz(-1.5497327) q[1];
sx q[1];
rz(-2.910639) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.445914) q[0];
sx q[0];
rz(-1.983108) q[0];
sx q[0];
rz(-0.59230174) q[0];
rz(-pi) q[1];
rz(0.43416331) q[2];
sx q[2];
rz(-1.6953529) q[2];
sx q[2];
rz(1.9366154) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4296023) q[1];
sx q[1];
rz(-1.3531405) q[1];
sx q[1];
rz(0.46516402) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4695314) q[3];
sx q[3];
rz(-1.5717744) q[3];
sx q[3];
rz(-3.1100215) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.93297282) q[2];
sx q[2];
rz(-1.695637) q[2];
sx q[2];
rz(0.43819532) q[2];
rz(-2.4402601) q[3];
sx q[3];
rz(-2.0739136) q[3];
sx q[3];
rz(0.35593885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0375131) q[0];
sx q[0];
rz(-2.6689745) q[0];
sx q[0];
rz(-1.0802826) q[0];
rz(1.0516306) q[1];
sx q[1];
rz(-2.0104505) q[1];
sx q[1];
rz(-1.0952605) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0290463) q[0];
sx q[0];
rz(-1.284082) q[0];
sx q[0];
rz(-0.53347828) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.3682179) q[2];
sx q[2];
rz(-1.5309586) q[2];
sx q[2];
rz(-1.7263655) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.124123) q[1];
sx q[1];
rz(-2.4061205) q[1];
sx q[1];
rz(1.7361438) q[1];
x q[2];
rz(0.67948158) q[3];
sx q[3];
rz(-0.71820161) q[3];
sx q[3];
rz(-1.3905202) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.1209391) q[2];
sx q[2];
rz(-1.8755269) q[2];
sx q[2];
rz(2.2361163) q[2];
rz(2.4261273) q[3];
sx q[3];
rz(-0.72532907) q[3];
sx q[3];
rz(-0.65413094) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41659551) q[0];
sx q[0];
rz(-1.409197) q[0];
sx q[0];
rz(-2.5247164) q[0];
rz(-1.1299409) q[1];
sx q[1];
rz(-1.8604953) q[1];
sx q[1];
rz(-1.4166191) q[1];
rz(-0.62192179) q[2];
sx q[2];
rz(-0.94246943) q[2];
sx q[2];
rz(3.111919) q[2];
rz(1.3601639) q[3];
sx q[3];
rz(-1.4223301) q[3];
sx q[3];
rz(-2.8098047) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
