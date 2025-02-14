OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.5247076) q[0];
sx q[0];
rz(-2.4029713) q[0];
sx q[0];
rz(3.024616) q[0];
rz(-2.1891131) q[1];
sx q[1];
rz(-1.4701966) q[1];
sx q[1];
rz(2.267946) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3461733) q[0];
sx q[0];
rz(-1.367698) q[0];
sx q[0];
rz(-0.83195306) q[0];
rz(0.9175542) q[2];
sx q[2];
rz(-2.5981256) q[2];
sx q[2];
rz(0.57910165) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3368766) q[1];
sx q[1];
rz(-2.118944) q[1];
sx q[1];
rz(3.0744661) q[1];
rz(2.1527613) q[3];
sx q[3];
rz(-2.6675794) q[3];
sx q[3];
rz(1.0128761) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.17254193) q[2];
sx q[2];
rz(-1.3323063) q[2];
sx q[2];
rz(-1.8633899) q[2];
rz(1.5365907) q[3];
sx q[3];
rz(-1.9032109) q[3];
sx q[3];
rz(-2.8240375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69739598) q[0];
sx q[0];
rz(-0.2187271) q[0];
sx q[0];
rz(-2.2636407) q[0];
rz(-2.3776993) q[1];
sx q[1];
rz(-2.0593675) q[1];
sx q[1];
rz(-1.1306995) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5661086) q[0];
sx q[0];
rz(-0.16666767) q[0];
sx q[0];
rz(-1.1995633) q[0];
x q[1];
rz(2.0085336) q[2];
sx q[2];
rz(-1.4302876) q[2];
sx q[2];
rz(2.2214691) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.4731684) q[1];
sx q[1];
rz(-2.3239909) q[1];
sx q[1];
rz(-0.017488285) q[1];
x q[2];
rz(-1.6783707) q[3];
sx q[3];
rz(-2.1940596) q[3];
sx q[3];
rz(-3.0855765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.45103249) q[2];
sx q[2];
rz(-1.6310383) q[2];
sx q[2];
rz(0.65563273) q[2];
rz(1.2304652) q[3];
sx q[3];
rz(-0.96223193) q[3];
sx q[3];
rz(2.3088096) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3967628) q[0];
sx q[0];
rz(-1.1711045) q[0];
sx q[0];
rz(-0.50286621) q[0];
rz(-2.9961131) q[1];
sx q[1];
rz(-1.611404) q[1];
sx q[1];
rz(-1.1605638) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7877951) q[0];
sx q[0];
rz(-2.702454) q[0];
sx q[0];
rz(2.3697858) q[0];
rz(-pi) q[1];
rz(1.4394748) q[2];
sx q[2];
rz(-1.223837) q[2];
sx q[2];
rz(1.0441213) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.5215456) q[1];
sx q[1];
rz(-1.1898387) q[1];
sx q[1];
rz(0.89922616) q[1];
rz(-1.3011906) q[3];
sx q[3];
rz(-1.9072215) q[3];
sx q[3];
rz(2.2158282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.1955371) q[2];
sx q[2];
rz(-2.1058407) q[2];
sx q[2];
rz(-0.76988402) q[2];
rz(0.59512538) q[3];
sx q[3];
rz(-2.6502521) q[3];
sx q[3];
rz(0.46014053) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50764099) q[0];
sx q[0];
rz(-2.0946298) q[0];
sx q[0];
rz(-1.3557583) q[0];
rz(0.62375623) q[1];
sx q[1];
rz(-1.535894) q[1];
sx q[1];
rz(0.49130586) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3771916) q[0];
sx q[0];
rz(-1.0391512) q[0];
sx q[0];
rz(1.7488983) q[0];
rz(-pi) q[1];
rz(-1.7801966) q[2];
sx q[2];
rz(-1.2970957) q[2];
sx q[2];
rz(1.4865246) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.32403261) q[1];
sx q[1];
rz(-2.3563317) q[1];
sx q[1];
rz(1.1263442) q[1];
x q[2];
rz(-0.66670615) q[3];
sx q[3];
rz(-1.9940064) q[3];
sx q[3];
rz(0.70021399) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.19646159) q[2];
sx q[2];
rz(-1.6170231) q[2];
sx q[2];
rz(2.5932236) q[2];
rz(-1.8114926) q[3];
sx q[3];
rz(-2.8434704) q[3];
sx q[3];
rz(2.8032081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2746628) q[0];
sx q[0];
rz(-2.6882956) q[0];
sx q[0];
rz(-2.0931661) q[0];
rz(-3.0323845) q[1];
sx q[1];
rz(-1.5914773) q[1];
sx q[1];
rz(-1.106326) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8620548) q[0];
sx q[0];
rz(-1.3509271) q[0];
sx q[0];
rz(-0.1105652) q[0];
x q[1];
rz(-0.335131) q[2];
sx q[2];
rz(-0.8742399) q[2];
sx q[2];
rz(2.1256688) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.2611003) q[1];
sx q[1];
rz(-1.893781) q[1];
sx q[1];
rz(2.9862981) q[1];
rz(-pi) q[2];
x q[2];
rz(0.89792395) q[3];
sx q[3];
rz(-0.37386383) q[3];
sx q[3];
rz(-1.7927235) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.76935524) q[2];
sx q[2];
rz(-1.8753588) q[2];
sx q[2];
rz(2.9388169) q[2];
rz(0.77962223) q[3];
sx q[3];
rz(-2.0519665) q[3];
sx q[3];
rz(0.37105086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62696537) q[0];
sx q[0];
rz(-0.045128673) q[0];
sx q[0];
rz(0.87920642) q[0];
rz(0.59576398) q[1];
sx q[1];
rz(-0.87409449) q[1];
sx q[1];
rz(1.8858058) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0107675) q[0];
sx q[0];
rz(-1.6629476) q[0];
sx q[0];
rz(1.6238201) q[0];
x q[1];
rz(-2.5686092) q[2];
sx q[2];
rz(-1.8610916) q[2];
sx q[2];
rz(3.109755) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.75460282) q[1];
sx q[1];
rz(-2.5856254) q[1];
sx q[1];
rz(2.703859) q[1];
rz(-2.7967795) q[3];
sx q[3];
rz(-2.4969144) q[3];
sx q[3];
rz(-0.0078474069) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6229728) q[2];
sx q[2];
rz(-2.8039248) q[2];
sx q[2];
rz(1.6597623) q[2];
rz(1.6984113) q[3];
sx q[3];
rz(-1.7549763) q[3];
sx q[3];
rz(-2.0229009) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5959394) q[0];
sx q[0];
rz(-0.82939363) q[0];
sx q[0];
rz(0.19350061) q[0];
rz(-0.045348383) q[1];
sx q[1];
rz(-1.1646611) q[1];
sx q[1];
rz(-2.4605816) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4880002) q[0];
sx q[0];
rz(-0.18463273) q[0];
sx q[0];
rz(-2.6499278) q[0];
x q[1];
rz(1.1521336) q[2];
sx q[2];
rz(-0.56850009) q[2];
sx q[2];
rz(2.1529039) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.6521218) q[1];
sx q[1];
rz(-0.014111405) q[1];
sx q[1];
rz(-1.1443787) q[1];
rz(-pi) q[2];
rz(1.4459386) q[3];
sx q[3];
rz(-1.9546851) q[3];
sx q[3];
rz(-0.53016183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.8650774) q[2];
sx q[2];
rz(-0.64537185) q[2];
sx q[2];
rz(-1.3481677) q[2];
rz(1.7776141) q[3];
sx q[3];
rz(-1.3185225) q[3];
sx q[3];
rz(-1.1905504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
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
rz(2.6337223) q[0];
sx q[0];
rz(-0.81955925) q[0];
sx q[0];
rz(-0.67935294) q[0];
rz(-0.23718111) q[1];
sx q[1];
rz(-0.91865426) q[1];
sx q[1];
rz(1.0666581) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0268008) q[0];
sx q[0];
rz(-1.566883) q[0];
sx q[0];
rz(2.6312909) q[0];
x q[1];
rz(-0.094631976) q[2];
sx q[2];
rz(-0.96853791) q[2];
sx q[2];
rz(-2.8631543) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.6427463) q[1];
sx q[1];
rz(-0.53697911) q[1];
sx q[1];
rz(1.3259423) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8130624) q[3];
sx q[3];
rz(-1.9764998) q[3];
sx q[3];
rz(2.2689087) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.25881585) q[2];
sx q[2];
rz(-1.9100185) q[2];
sx q[2];
rz(0.41137496) q[2];
rz(1.8607633) q[3];
sx q[3];
rz(-0.6461834) q[3];
sx q[3];
rz(2.083174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48788747) q[0];
sx q[0];
rz(-1.7937086) q[0];
sx q[0];
rz(2.2109798) q[0];
rz(2.6764684) q[1];
sx q[1];
rz(-0.5828751) q[1];
sx q[1];
rz(1.7238269) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0865308) q[0];
sx q[0];
rz(-1.3320005) q[0];
sx q[0];
rz(0.60542528) q[0];
x q[1];
rz(-0.64196824) q[2];
sx q[2];
rz(-2.1479508) q[2];
sx q[2];
rz(-2.746144) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.237933) q[1];
sx q[1];
rz(-1.2301908) q[1];
sx q[1];
rz(-1.8800859) q[1];
x q[2];
rz(-2.1499436) q[3];
sx q[3];
rz(-1.4738899) q[3];
sx q[3];
rz(-2.4976418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.96665367) q[2];
sx q[2];
rz(-0.83028364) q[2];
sx q[2];
rz(-1.9261599) q[2];
rz(-2.4231966) q[3];
sx q[3];
rz(-1.4677488) q[3];
sx q[3];
rz(0.64232701) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4833118) q[0];
sx q[0];
rz(-2.719306) q[0];
sx q[0];
rz(1.8102113) q[0];
rz(2.4347958) q[1];
sx q[1];
rz(-1.7168609) q[1];
sx q[1];
rz(-1.4250863) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15062103) q[0];
sx q[0];
rz(-2.2036834) q[0];
sx q[0];
rz(-2.1314243) q[0];
rz(-0.064878929) q[2];
sx q[2];
rz(-0.91514313) q[2];
sx q[2];
rz(2.3196089) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4689961) q[1];
sx q[1];
rz(-2.0083659) q[1];
sx q[1];
rz(0.5794748) q[1];
rz(-pi) q[2];
rz(-0.8031527) q[3];
sx q[3];
rz(-1.8997822) q[3];
sx q[3];
rz(-2.064049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.66296545) q[2];
sx q[2];
rz(-2.4586283) q[2];
sx q[2];
rz(1.1753725) q[2];
rz(3.1183682) q[3];
sx q[3];
rz(-0.5718137) q[3];
sx q[3];
rz(2.8779252) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5671253) q[0];
sx q[0];
rz(-0.99526417) q[0];
sx q[0];
rz(-1.9010726) q[0];
rz(-0.65470882) q[1];
sx q[1];
rz(-1.2798825) q[1];
sx q[1];
rz(1.9725694) q[1];
rz(2.2562182) q[2];
sx q[2];
rz(-0.17687951) q[2];
sx q[2];
rz(2.542873) q[2];
rz(0.67778604) q[3];
sx q[3];
rz(-0.98193632) q[3];
sx q[3];
rz(2.2241398) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
