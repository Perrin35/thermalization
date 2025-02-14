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
rz(-2.8415866) q[0];
sx q[0];
rz(-2.0553148) q[0];
sx q[0];
rz(-2.3334184) q[0];
rz(1.9975245) q[1];
sx q[1];
rz(-2.3699528) q[1];
sx q[1];
rz(1.5789403) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.232036) q[0];
sx q[0];
rz(-2.252451) q[0];
sx q[0];
rz(-0.85083346) q[0];
x q[1];
rz(-1.3035266) q[2];
sx q[2];
rz(-1.7198945) q[2];
sx q[2];
rz(-1.2318512) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.4889604) q[1];
sx q[1];
rz(-0.51981407) q[1];
sx q[1];
rz(-1.7955154) q[1];
rz(-2.4993131) q[3];
sx q[3];
rz(-0.056996973) q[3];
sx q[3];
rz(1.1617203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.943104) q[2];
sx q[2];
rz(-2.9297332) q[2];
sx q[2];
rz(2.1073821) q[2];
rz(-0.69283038) q[3];
sx q[3];
rz(-2.0878017) q[3];
sx q[3];
rz(-2.0809765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8387872) q[0];
sx q[0];
rz(-3.1133339) q[0];
sx q[0];
rz(0.57648188) q[0];
rz(3.1209962) q[1];
sx q[1];
rz(-2.7437904) q[1];
sx q[1];
rz(-2.1131262) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45555025) q[0];
sx q[0];
rz(-2.8032113) q[0];
sx q[0];
rz(0.69940059) q[0];
x q[1];
rz(-2.8345077) q[2];
sx q[2];
rz(-1.9283617) q[2];
sx q[2];
rz(1.7890695) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8081144) q[1];
sx q[1];
rz(-1.9030182) q[1];
sx q[1];
rz(1.9971041) q[1];
x q[2];
rz(1.1787492) q[3];
sx q[3];
rz(-2.1806751) q[3];
sx q[3];
rz(1.4489685) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.91745201) q[2];
sx q[2];
rz(-0.96106207) q[2];
sx q[2];
rz(-1.8769598) q[2];
rz(-2.1680016) q[3];
sx q[3];
rz(-1.6907938) q[3];
sx q[3];
rz(-0.26315954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6384386) q[0];
sx q[0];
rz(-1.8373024) q[0];
sx q[0];
rz(2.2880182) q[0];
rz(2.0897934) q[1];
sx q[1];
rz(-0.97426668) q[1];
sx q[1];
rz(0.015017088) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5511709) q[0];
sx q[0];
rz(-0.71066868) q[0];
sx q[0];
rz(1.5842313) q[0];
rz(0.21785801) q[2];
sx q[2];
rz(-3.0431722) q[2];
sx q[2];
rz(-2.276536) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.4045118) q[1];
sx q[1];
rz(-2.28302) q[1];
sx q[1];
rz(1.668307) q[1];
x q[2];
rz(-1.6347359) q[3];
sx q[3];
rz(-1.2815164) q[3];
sx q[3];
rz(0.48529406) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.1198279) q[2];
sx q[2];
rz(-1.652635) q[2];
sx q[2];
rz(-2.8601698) q[2];
rz(1.4540539) q[3];
sx q[3];
rz(-1.931793) q[3];
sx q[3];
rz(-2.37229) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5948831) q[0];
sx q[0];
rz(-1.895772) q[0];
sx q[0];
rz(-2.967714) q[0];
rz(-1.8100544) q[1];
sx q[1];
rz(-1.7200108) q[1];
sx q[1];
rz(-1.4283659) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7998909) q[0];
sx q[0];
rz(-0.19289324) q[0];
sx q[0];
rz(1.6454205) q[0];
x q[1];
rz(-2.8160353) q[2];
sx q[2];
rz(-2.7635305) q[2];
sx q[2];
rz(-0.69272536) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.54879872) q[1];
sx q[1];
rz(-0.6585291) q[1];
sx q[1];
rz(-3.0786773) q[1];
x q[2];
rz(1.9589728) q[3];
sx q[3];
rz(-1.2283162) q[3];
sx q[3];
rz(-1.4701209) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.4398769) q[2];
sx q[2];
rz(-2.3150847) q[2];
sx q[2];
rz(-0.41713777) q[2];
rz(2.9719628) q[3];
sx q[3];
rz(-3.0483584) q[3];
sx q[3];
rz(2.4197742) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4193264) q[0];
sx q[0];
rz(-0.37501431) q[0];
sx q[0];
rz(0.30935031) q[0];
rz(0.7695235) q[1];
sx q[1];
rz(-1.0898432) q[1];
sx q[1];
rz(1.6228898) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2941336) q[0];
sx q[0];
rz(-1.3193697) q[0];
sx q[0];
rz(2.2537116) q[0];
rz(-0.73878397) q[2];
sx q[2];
rz(-2.1257943) q[2];
sx q[2];
rz(-0.273663) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.57598439) q[1];
sx q[1];
rz(-0.19153015) q[1];
sx q[1];
rz(1.4148786) q[1];
rz(-pi) q[2];
rz(0.098747323) q[3];
sx q[3];
rz(-1.704146) q[3];
sx q[3];
rz(-0.46252201) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.3994483) q[2];
sx q[2];
rz(-2.1592906) q[2];
sx q[2];
rz(2.1461416) q[2];
rz(-2.4230912) q[3];
sx q[3];
rz(-1.813443) q[3];
sx q[3];
rz(2.3464581) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8566078) q[0];
sx q[0];
rz(-1.4549078) q[0];
sx q[0];
rz(-1.6108151) q[0];
rz(1.6948304) q[1];
sx q[1];
rz(-1.6981354) q[1];
sx q[1];
rz(1.7036499) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79057825) q[0];
sx q[0];
rz(-1.5503128) q[0];
sx q[0];
rz(-1.9786506) q[0];
rz(-pi) q[1];
rz(-1.1968568) q[2];
sx q[2];
rz(-1.0122006) q[2];
sx q[2];
rz(-0.93511673) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.2338381) q[1];
sx q[1];
rz(-1.4871039) q[1];
sx q[1];
rz(0.5096883) q[1];
rz(-1.5737278) q[3];
sx q[3];
rz(-0.81122196) q[3];
sx q[3];
rz(2.7377759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.97856727) q[2];
sx q[2];
rz(-1.7820396) q[2];
sx q[2];
rz(1.2320409) q[2];
rz(1.1466522) q[3];
sx q[3];
rz(-1.8012643) q[3];
sx q[3];
rz(2.8310217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6335886) q[0];
sx q[0];
rz(-2.3242943) q[0];
sx q[0];
rz(2.001413) q[0];
rz(-2.2445402) q[1];
sx q[1];
rz(-1.8042118) q[1];
sx q[1];
rz(-3.1111029) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.545118) q[0];
sx q[0];
rz(-1.6870572) q[0];
sx q[0];
rz(1.3482987) q[0];
rz(0.42564904) q[2];
sx q[2];
rz(-2.4535745) q[2];
sx q[2];
rz(1.4378215) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.5329689) q[1];
sx q[1];
rz(-2.2734959) q[1];
sx q[1];
rz(-2.7547902) q[1];
rz(-0.45400158) q[3];
sx q[3];
rz(-2.2488065) q[3];
sx q[3];
rz(-1.0832204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.9474779) q[2];
sx q[2];
rz(-1.0741445) q[2];
sx q[2];
rz(-1.2987632) q[2];
rz(1.4878368) q[3];
sx q[3];
rz(-2.1066809) q[3];
sx q[3];
rz(1.8208767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3461935) q[0];
sx q[0];
rz(-2.8454056) q[0];
sx q[0];
rz(-1.7561703) q[0];
rz(1.9253383) q[1];
sx q[1];
rz(-1.6938035) q[1];
sx q[1];
rz(-2.6522955) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2953257) q[0];
sx q[0];
rz(-1.5510484) q[0];
sx q[0];
rz(0.67157816) q[0];
rz(-pi) q[1];
rz(1.4665746) q[2];
sx q[2];
rz(-2.8036593) q[2];
sx q[2];
rz(-2.3399885) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.41114901) q[1];
sx q[1];
rz(-0.48161067) q[1];
sx q[1];
rz(1.3661307) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9527444) q[3];
sx q[3];
rz(-2.8738465) q[3];
sx q[3];
rz(1.6588039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.74799246) q[2];
sx q[2];
rz(-1.1425428) q[2];
sx q[2];
rz(0.46621123) q[2];
rz(1.5318058) q[3];
sx q[3];
rz(-1.5038871) q[3];
sx q[3];
rz(0.32061779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6521249) q[0];
sx q[0];
rz(-0.89972275) q[0];
sx q[0];
rz(-0.049064431) q[0];
rz(-0.24066726) q[1];
sx q[1];
rz(-1.0180232) q[1];
sx q[1];
rz(-0.022445591) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81240679) q[0];
sx q[0];
rz(-2.5617449) q[0];
sx q[0];
rz(-1.5270698) q[0];
rz(-pi) q[1];
rz(-2.4149553) q[2];
sx q[2];
rz(-0.49984806) q[2];
sx q[2];
rz(-2.8195087) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.6091303) q[1];
sx q[1];
rz(-2.3098619) q[1];
sx q[1];
rz(0.66166454) q[1];
rz(2.4565036) q[3];
sx q[3];
rz(-1.2108608) q[3];
sx q[3];
rz(-1.0459096) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.13751328) q[2];
sx q[2];
rz(-1.235032) q[2];
sx q[2];
rz(-0.9683041) q[2];
rz(-2.3051895) q[3];
sx q[3];
rz(-2.8439549) q[3];
sx q[3];
rz(1.3692726) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4150998) q[0];
sx q[0];
rz(-1.2905755) q[0];
sx q[0];
rz(0.37877628) q[0];
rz(-3.1355766) q[1];
sx q[1];
rz(-0.52979398) q[1];
sx q[1];
rz(2.5078497) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0519052) q[0];
sx q[0];
rz(-1.9080592) q[0];
sx q[0];
rz(-0.71292065) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.815755) q[2];
sx q[2];
rz(-0.58734054) q[2];
sx q[2];
rz(-2.7224685) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.56713518) q[1];
sx q[1];
rz(-2.0525041) q[1];
sx q[1];
rz(-2.5888799) q[1];
rz(-0.87234906) q[3];
sx q[3];
rz(-2.3309338) q[3];
sx q[3];
rz(-1.3083506) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.70574957) q[2];
sx q[2];
rz(-1.6195932) q[2];
sx q[2];
rz(-0.61326927) q[2];
rz(-0.67522007) q[3];
sx q[3];
rz(-2.7628511) q[3];
sx q[3];
rz(-2.3238382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2262065) q[0];
sx q[0];
rz(-1.362726) q[0];
sx q[0];
rz(1.2363731) q[0];
rz(2.9551103) q[1];
sx q[1];
rz(-1.7541371) q[1];
sx q[1];
rz(-1.6288155) q[1];
rz(1.5507422) q[2];
sx q[2];
rz(-2.0303844) q[2];
sx q[2];
rz(-2.5302237) q[2];
rz(1.4433031) q[3];
sx q[3];
rz(-2.5155407) q[3];
sx q[3];
rz(-1.2609353) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
