OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.6168851) q[0];
sx q[0];
rz(-0.73862139) q[0];
sx q[0];
rz(0.11697669) q[0];
rz(0.95247954) q[1];
sx q[1];
rz(-1.6713961) q[1];
sx q[1];
rz(-2.267946) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7352493) q[0];
sx q[0];
rz(-0.85052089) q[0];
sx q[0];
rz(2.8699037) q[0];
rz(1.123465) q[2];
sx q[2];
rz(-1.8904933) q[2];
sx q[2];
rz(1.5699707) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.80471604) q[1];
sx q[1];
rz(-1.0226486) q[1];
sx q[1];
rz(3.0744661) q[1];
x q[2];
rz(-2.8667371) q[3];
sx q[3];
rz(-1.9620202) q[3];
sx q[3];
rz(-1.6496342) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.9690507) q[2];
sx q[2];
rz(-1.8092864) q[2];
sx q[2];
rz(1.8633899) q[2];
rz(-1.5365907) q[3];
sx q[3];
rz(-1.2383818) q[3];
sx q[3];
rz(-2.8240375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69739598) q[0];
sx q[0];
rz(-2.9228656) q[0];
sx q[0];
rz(-2.2636407) q[0];
rz(0.76389337) q[1];
sx q[1];
rz(-2.0593675) q[1];
sx q[1];
rz(2.0108932) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7797403) q[0];
sx q[0];
rz(-1.6310143) q[0];
sx q[0];
rz(-1.7262993) q[0];
x q[1];
rz(1.133059) q[2];
sx q[2];
rz(-1.4302876) q[2];
sx q[2];
rz(0.92012355) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.66842428) q[1];
sx q[1];
rz(-2.3239909) q[1];
sx q[1];
rz(0.017488285) q[1];
x q[2];
rz(2.5155792) q[3];
sx q[3];
rz(-1.6580868) q[3];
sx q[3];
rz(-1.5638626) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6905602) q[2];
sx q[2];
rz(-1.5105543) q[2];
sx q[2];
rz(2.4859599) q[2];
rz(-1.2304652) q[3];
sx q[3];
rz(-0.96223193) q[3];
sx q[3];
rz(0.83278304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7448298) q[0];
sx q[0];
rz(-1.1711045) q[0];
sx q[0];
rz(0.50286621) q[0];
rz(-2.9961131) q[1];
sx q[1];
rz(-1.5301887) q[1];
sx q[1];
rz(1.1605638) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9391286) q[0];
sx q[0];
rz(-1.8718429) q[0];
sx q[0];
rz(-0.32471628) q[0];
rz(-0.34744743) q[2];
sx q[2];
rz(-0.37003741) q[2];
sx q[2];
rz(1.7269788) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7532901) q[1];
sx q[1];
rz(-0.75725746) q[1];
sx q[1];
rz(0.99885916) q[1];
x q[2];
rz(-2.4907001) q[3];
sx q[3];
rz(-0.42789352) q[3];
sx q[3];
rz(-0.22881697) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.9460556) q[2];
sx q[2];
rz(-2.1058407) q[2];
sx q[2];
rz(2.3717086) q[2];
rz(-2.5464673) q[3];
sx q[3];
rz(-0.49134058) q[3];
sx q[3];
rz(2.6814521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50764099) q[0];
sx q[0];
rz(-2.0946298) q[0];
sx q[0];
rz(-1.3557583) q[0];
rz(2.5178364) q[1];
sx q[1];
rz(-1.6056986) q[1];
sx q[1];
rz(-2.6502868) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7183863) q[0];
sx q[0];
rz(-2.583637) q[0];
sx q[0];
rz(-2.8490081) q[0];
rz(-0.27950269) q[2];
sx q[2];
rz(-1.3692999) q[2];
sx q[2];
rz(-3.1147) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.2196893) q[1];
sx q[1];
rz(-1.8796726) q[1];
sx q[1];
rz(2.3050453) q[1];
rz(-pi) q[2];
rz(-0.66670615) q[3];
sx q[3];
rz(-1.1475862) q[3];
sx q[3];
rz(2.4413787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.19646159) q[2];
sx q[2];
rz(-1.6170231) q[2];
sx q[2];
rz(0.54836908) q[2];
rz(1.8114926) q[3];
sx q[3];
rz(-0.29812223) q[3];
sx q[3];
rz(2.8032081) q[3];
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
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86692989) q[0];
sx q[0];
rz(-0.45329705) q[0];
sx q[0];
rz(-1.0484265) q[0];
rz(-3.0323845) q[1];
sx q[1];
rz(-1.5914773) q[1];
sx q[1];
rz(-1.106326) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8620548) q[0];
sx q[0];
rz(-1.3509271) q[0];
sx q[0];
rz(0.1105652) q[0];
rz(2.8064617) q[2];
sx q[2];
rz(-0.8742399) q[2];
sx q[2];
rz(-1.0159238) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.2611003) q[1];
sx q[1];
rz(-1.893781) q[1];
sx q[1];
rz(-2.9862981) q[1];
rz(2.2436687) q[3];
sx q[3];
rz(-2.7677288) q[3];
sx q[3];
rz(-1.7927235) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.76935524) q[2];
sx q[2];
rz(-1.8753588) q[2];
sx q[2];
rz(-2.9388169) q[2];
rz(0.77962223) q[3];
sx q[3];
rz(-2.0519665) q[3];
sx q[3];
rz(-2.7705418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62696537) q[0];
sx q[0];
rz(-3.096464) q[0];
sx q[0];
rz(-2.2623862) q[0];
rz(-0.59576398) q[1];
sx q[1];
rz(-2.2674982) q[1];
sx q[1];
rz(-1.2557868) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1308252) q[0];
sx q[0];
rz(-1.4786451) q[0];
sx q[0];
rz(-1.6238201) q[0];
rz(1.2292172) q[2];
sx q[2];
rz(-2.1170127) q[2];
sx q[2];
rz(1.721576) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.8906011) q[1];
sx q[1];
rz(-1.0724147) q[1];
sx q[1];
rz(-1.313264) q[1];
rz(2.5257931) q[3];
sx q[3];
rz(-1.7753505) q[3];
sx q[3];
rz(1.2990862) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.51861989) q[2];
sx q[2];
rz(-0.33766782) q[2];
sx q[2];
rz(-1.4818304) q[2];
rz(1.4431813) q[3];
sx q[3];
rz(-1.3866164) q[3];
sx q[3];
rz(1.1186918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5959394) q[0];
sx q[0];
rz(-2.312199) q[0];
sx q[0];
rz(-0.19350061) q[0];
rz(3.0962443) q[1];
sx q[1];
rz(-1.1646611) q[1];
sx q[1];
rz(0.68101105) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.986833) q[0];
sx q[0];
rz(-1.4082419) q[0];
sx q[0];
rz(1.6587371) q[0];
rz(1.0424642) q[2];
sx q[2];
rz(-1.7914504) q[2];
sx q[2];
rz(-0.94086066) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.4894708) q[1];
sx q[1];
rz(-0.014111405) q[1];
sx q[1];
rz(1.997214) q[1];
rz(-0.29905921) q[3];
sx q[3];
rz(-2.7388696) q[3];
sx q[3];
rz(-2.2880768) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.8650774) q[2];
sx q[2];
rz(-2.4962208) q[2];
sx q[2];
rz(-1.7934249) q[2];
rz(-1.7776141) q[3];
sx q[3];
rz(-1.8230702) q[3];
sx q[3];
rz(1.9510423) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50787038) q[0];
sx q[0];
rz(-2.3220334) q[0];
sx q[0];
rz(-0.67935294) q[0];
rz(-2.9044115) q[1];
sx q[1];
rz(-0.91865426) q[1];
sx q[1];
rz(2.0749345) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1147918) q[0];
sx q[0];
rz(-1.5747096) q[0];
sx q[0];
rz(2.6312909) q[0];
rz(-pi) q[1];
rz(2.1751498) q[2];
sx q[2];
rz(-1.4928515) q[2];
sx q[2];
rz(1.9029531) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.8581526) q[1];
sx q[1];
rz(-1.4464708) q[1];
sx q[1];
rz(-2.0945626) q[1];
rz(-pi) q[2];
rz(2.2150463) q[3];
sx q[3];
rz(-2.6253581) q[3];
sx q[3];
rz(-2.981319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.8827768) q[2];
sx q[2];
rz(-1.9100185) q[2];
sx q[2];
rz(0.41137496) q[2];
rz(-1.8607633) q[3];
sx q[3];
rz(-0.6461834) q[3];
sx q[3];
rz(-2.083174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48788747) q[0];
sx q[0];
rz(-1.7937086) q[0];
sx q[0];
rz(0.93061289) q[0];
rz(0.46512428) q[1];
sx q[1];
rz(-2.5587176) q[1];
sx q[1];
rz(1.7238269) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8450369) q[0];
sx q[0];
rz(-0.64528685) q[0];
sx q[0];
rz(-0.40420239) q[0];
rz(-2.2533855) q[2];
sx q[2];
rz(-2.0964551) q[2];
sx q[2];
rz(-1.5789248) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4747367) q[1];
sx q[1];
rz(-0.45598331) q[1];
sx q[1];
rz(0.70959301) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7465789) q[3];
sx q[3];
rz(-0.58628288) q[3];
sx q[3];
rz(2.0678584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.174939) q[2];
sx q[2];
rz(-2.311309) q[2];
sx q[2];
rz(1.2154328) q[2];
rz(0.71839607) q[3];
sx q[3];
rz(-1.4677488) q[3];
sx q[3];
rz(-2.4992656) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4833118) q[0];
sx q[0];
rz(-2.719306) q[0];
sx q[0];
rz(1.8102113) q[0];
rz(-0.70679682) q[1];
sx q[1];
rz(-1.4247318) q[1];
sx q[1];
rz(-1.7165064) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15062103) q[0];
sx q[0];
rz(-0.93790927) q[0];
sx q[0];
rz(-2.1314243) q[0];
rz(-pi) q[1];
rz(0.91412506) q[2];
sx q[2];
rz(-1.5193835) q[2];
sx q[2];
rz(2.35319) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.37230587) q[1];
sx q[1];
rz(-2.0897749) q[1];
sx q[1];
rz(-2.080588) q[1];
rz(-pi) q[2];
x q[2];
rz(0.8031527) q[3];
sx q[3];
rz(-1.2418104) q[3];
sx q[3];
rz(-2.064049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.4786272) q[2];
sx q[2];
rz(-2.4586283) q[2];
sx q[2];
rz(-1.9662201) q[2];
rz(3.1183682) q[3];
sx q[3];
rz(-0.5718137) q[3];
sx q[3];
rz(-0.2636675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
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
rz(-2.2562182) q[2];
sx q[2];
rz(-2.9647131) q[2];
sx q[2];
rz(-0.59871968) q[2];
rz(0.8169218) q[3];
sx q[3];
rz(-0.86614843) q[3];
sx q[3];
rz(-1.8843005) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
