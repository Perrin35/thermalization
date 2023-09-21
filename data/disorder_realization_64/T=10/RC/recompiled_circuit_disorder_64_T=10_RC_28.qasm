OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.7862608) q[0];
sx q[0];
rz(-0.064602764) q[0];
sx q[0];
rz(0.021615418) q[0];
rz(-0.99524438) q[1];
sx q[1];
rz(-1.3270451) q[1];
sx q[1];
rz(1.8099161) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3992213) q[0];
sx q[0];
rz(-1.0065777) q[0];
sx q[0];
rz(0.39682927) q[0];
x q[1];
rz(1.7982499) q[2];
sx q[2];
rz(-1.6709575) q[2];
sx q[2];
rz(1.4129461) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0845619) q[1];
sx q[1];
rz(-1.7006526) q[1];
sx q[1];
rz(-2.3741541) q[1];
rz(0.94727912) q[3];
sx q[3];
rz(-1.2926658) q[3];
sx q[3];
rz(-2.5840685) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.33064476) q[2];
sx q[2];
rz(-1.4928879) q[2];
sx q[2];
rz(0.60418207) q[2];
rz(-2.1172681) q[3];
sx q[3];
rz(-1.9842792) q[3];
sx q[3];
rz(1.1272875) q[3];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3246831) q[0];
sx q[0];
rz(-0.01318251) q[0];
sx q[0];
rz(1.0634134) q[0];
rz(2.2564607) q[1];
sx q[1];
rz(-1.5567895) q[1];
sx q[1];
rz(3.1399472) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9379399) q[0];
sx q[0];
rz(-1.5814039) q[0];
sx q[0];
rz(-1.525735) q[0];
x q[1];
rz(1.8284945) q[2];
sx q[2];
rz(-0.94444599) q[2];
sx q[2];
rz(-1.4305654) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.26047036) q[1];
sx q[1];
rz(-0.80658856) q[1];
sx q[1];
rz(1.9116198) q[1];
rz(0.1699154) q[3];
sx q[3];
rz(-1.3193469) q[3];
sx q[3];
rz(-2.2268695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.985618) q[2];
sx q[2];
rz(-1.6321471) q[2];
sx q[2];
rz(0.70811159) q[2];
rz(-1.0937141) q[3];
sx q[3];
rz(-1.8023068) q[3];
sx q[3];
rz(-2.1480907) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22096069) q[0];
sx q[0];
rz(-1.7376124) q[0];
sx q[0];
rz(-2.3143342) q[0];
rz(3.1365085) q[1];
sx q[1];
rz(-1.2132443) q[1];
sx q[1];
rz(1.089383) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2978965) q[0];
sx q[0];
rz(-1.9934405) q[0];
sx q[0];
rz(-2.1556913) q[0];
rz(-2.5554352) q[2];
sx q[2];
rz(-0.68550368) q[2];
sx q[2];
rz(-3.0229178) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.5685801) q[1];
sx q[1];
rz(-1.9553361) q[1];
sx q[1];
rz(-3.0494855) q[1];
x q[2];
rz(0.79077625) q[3];
sx q[3];
rz(-1.9867992) q[3];
sx q[3];
rz(-2.5115867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.0777145) q[2];
sx q[2];
rz(-0.87746799) q[2];
sx q[2];
rz(-2.2568978) q[2];
rz(-1.1832773) q[3];
sx q[3];
rz(-1.8235455) q[3];
sx q[3];
rz(-1.4900835) q[3];
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
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5692212) q[0];
sx q[0];
rz(-2.5771993) q[0];
sx q[0];
rz(-0.72682056) q[0];
rz(0.80365333) q[1];
sx q[1];
rz(-1.0875965) q[1];
sx q[1];
rz(-2.7817536) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.026304631) q[0];
sx q[0];
rz(-1.7548772) q[0];
sx q[0];
rz(0.91362761) q[0];
rz(-pi) q[1];
rz(-0.14881046) q[2];
sx q[2];
rz(-0.65650046) q[2];
sx q[2];
rz(-1.7800926) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6872014) q[1];
sx q[1];
rz(-1.9481716) q[1];
sx q[1];
rz(-0.45339938) q[1];
rz(-pi) q[2];
rz(-0.8538586) q[3];
sx q[3];
rz(-2.6918908) q[3];
sx q[3];
rz(-0.85292294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.56746733) q[2];
sx q[2];
rz(-1.987477) q[2];
sx q[2];
rz(0.7652258) q[2];
rz(2.3848173) q[3];
sx q[3];
rz(-2.5496343) q[3];
sx q[3];
rz(2.9005907) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1446447) q[0];
sx q[0];
rz(-0.47676555) q[0];
sx q[0];
rz(-1.416052) q[0];
rz(-2.7744746) q[1];
sx q[1];
rz(-1.7558302) q[1];
sx q[1];
rz(2.1062772) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0400378) q[0];
sx q[0];
rz(-0.77417513) q[0];
sx q[0];
rz(-0.9206307) q[0];
rz(-pi) q[1];
rz(3.0683124) q[2];
sx q[2];
rz(-0.49256941) q[2];
sx q[2];
rz(0.57893334) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4673748) q[1];
sx q[1];
rz(-1.6661223) q[1];
sx q[1];
rz(1.6256888) q[1];
x q[2];
rz(1.3607929) q[3];
sx q[3];
rz(-1.8404507) q[3];
sx q[3];
rz(-2.8402013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.7845903) q[2];
sx q[2];
rz(-1.4214397) q[2];
sx q[2];
rz(2.3727097) q[2];
rz(2.8055577) q[3];
sx q[3];
rz(-0.78032812) q[3];
sx q[3];
rz(1.6736354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96200213) q[0];
sx q[0];
rz(-0.63350326) q[0];
sx q[0];
rz(1.996421) q[0];
rz(1.1046474) q[1];
sx q[1];
rz(-1.9112174) q[1];
sx q[1];
rz(-0.2072269) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3375895) q[0];
sx q[0];
rz(-2.1875256) q[0];
sx q[0];
rz(-1.7417275) q[0];
x q[1];
rz(2.5099498) q[2];
sx q[2];
rz(-2.4372059) q[2];
sx q[2];
rz(2.7001675) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.0761557) q[1];
sx q[1];
rz(-1.8199311) q[1];
sx q[1];
rz(0.25030987) q[1];
x q[2];
rz(2.3239922) q[3];
sx q[3];
rz(-1.451965) q[3];
sx q[3];
rz(0.51018754) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.8032288) q[2];
sx q[2];
rz(-2.2613566) q[2];
sx q[2];
rz(-2.4198789) q[2];
rz(-1.8668113) q[3];
sx q[3];
rz(-1.5694247) q[3];
sx q[3];
rz(-2.8815564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8969144) q[0];
sx q[0];
rz(-1.5748064) q[0];
sx q[0];
rz(-2.4095643) q[0];
rz(-0.0094982068) q[1];
sx q[1];
rz(-2.5932725) q[1];
sx q[1];
rz(-2.8498555) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1563819) q[0];
sx q[0];
rz(-1.1522326) q[0];
sx q[0];
rz(1.6267938) q[0];
x q[1];
rz(-3.0351699) q[2];
sx q[2];
rz(-1.2692045) q[2];
sx q[2];
rz(-2.0472722) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.3953494) q[1];
sx q[1];
rz(-0.33848539) q[1];
sx q[1];
rz(-2.6836718) q[1];
x q[2];
rz(-1.2613867) q[3];
sx q[3];
rz(-2.072812) q[3];
sx q[3];
rz(-1.7878143) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.26414028) q[2];
sx q[2];
rz(-1.9984657) q[2];
sx q[2];
rz(0.83731246) q[2];
rz(-1.1710179) q[3];
sx q[3];
rz(-1.5191017) q[3];
sx q[3];
rz(3.0795735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0983122) q[0];
sx q[0];
rz(-0.68798143) q[0];
sx q[0];
rz(-0.6638546) q[0];
rz(0.10617667) q[1];
sx q[1];
rz(-0.60634923) q[1];
sx q[1];
rz(0.95867872) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6129235) q[0];
sx q[0];
rz(-1.0597502) q[0];
sx q[0];
rz(0.59707609) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1088976) q[2];
sx q[2];
rz(-1.3287462) q[2];
sx q[2];
rz(0.48697105) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.2221335) q[1];
sx q[1];
rz(-2.6047915) q[1];
sx q[1];
rz(-1.3608576) q[1];
rz(1.6931375) q[3];
sx q[3];
rz(-2.1515176) q[3];
sx q[3];
rz(-1.8078705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.133193) q[2];
sx q[2];
rz(-1.9063176) q[2];
sx q[2];
rz(-0.91910249) q[2];
rz(1.7637926) q[3];
sx q[3];
rz(-1.2104687) q[3];
sx q[3];
rz(-1.9581883) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8063426) q[0];
sx q[0];
rz(-0.85334539) q[0];
sx q[0];
rz(-0.43689716) q[0];
rz(2.4412952) q[1];
sx q[1];
rz(-1.4236139) q[1];
sx q[1];
rz(0.46554309) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3532928) q[0];
sx q[0];
rz(-2.4294937) q[0];
sx q[0];
rz(-0.13029356) q[0];
x q[1];
rz(-3.0131857) q[2];
sx q[2];
rz(-1.1616716) q[2];
sx q[2];
rz(0.94305925) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.000396) q[1];
sx q[1];
rz(-1.9138412) q[1];
sx q[1];
rz(-2.6617906) q[1];
rz(-pi) q[2];
rz(-2.8903928) q[3];
sx q[3];
rz(-1.2217055) q[3];
sx q[3];
rz(0.72124764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.41032252) q[2];
sx q[2];
rz(-0.53933829) q[2];
sx q[2];
rz(1.4979866) q[2];
rz(0.20478976) q[3];
sx q[3];
rz(-2.1361165) q[3];
sx q[3];
rz(0.27206102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64514226) q[0];
sx q[0];
rz(-2.5800939) q[0];
sx q[0];
rz(1.1219332) q[0];
rz(2.3902068) q[1];
sx q[1];
rz(-1.5263999) q[1];
sx q[1];
rz(-2.5591992) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84401417) q[0];
sx q[0];
rz(-1.3676924) q[0];
sx q[0];
rz(2.2556979) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0803797) q[2];
sx q[2];
rz(-1.4368125) q[2];
sx q[2];
rz(1.6756563) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.6813587) q[1];
sx q[1];
rz(-1.453062) q[1];
sx q[1];
rz(1.387499) q[1];
x q[2];
rz(2.9243745) q[3];
sx q[3];
rz(-1.5981042) q[3];
sx q[3];
rz(-2.581493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.1311243) q[2];
sx q[2];
rz(-0.98103395) q[2];
sx q[2];
rz(-2.3790512) q[2];
rz(1.7307581) q[3];
sx q[3];
rz(-2.2236731) q[3];
sx q[3];
rz(-1.5738752) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
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
rz(2.0653771) q[0];
sx q[0];
rz(-0.98631728) q[0];
sx q[0];
rz(1.4022723) q[0];
rz(1.3394042) q[1];
sx q[1];
rz(-1.4504455) q[1];
sx q[1];
rz(1.6557678) q[1];
rz(1.273524) q[2];
sx q[2];
rz(-0.77902972) q[2];
sx q[2];
rz(-3.0697889) q[2];
rz(-1.0675666) q[3];
sx q[3];
rz(-1.9964841) q[3];
sx q[3];
rz(-1.869429) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
