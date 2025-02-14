OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.4021969) q[0];
sx q[0];
rz(-0.88391179) q[0];
sx q[0];
rz(-2.28595) q[0];
rz(-0.17172509) q[1];
sx q[1];
rz(6.167616) q[1];
sx q[1];
rz(12.003916) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5006258) q[0];
sx q[0];
rz(-1.62407) q[0];
sx q[0];
rz(1.2059653) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9859574) q[2];
sx q[2];
rz(-1.068972) q[2];
sx q[2];
rz(1.4841472) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.89002883) q[1];
sx q[1];
rz(-1.152413) q[1];
sx q[1];
rz(-2.0189925) q[1];
rz(0.39333087) q[3];
sx q[3];
rz(-1.9489904) q[3];
sx q[3];
rz(2.8761169) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.90613753) q[2];
sx q[2];
rz(-1.1780585) q[2];
sx q[2];
rz(2.3423024) q[2];
rz(0.47131395) q[3];
sx q[3];
rz(-0.91336942) q[3];
sx q[3];
rz(0.93588626) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1021295) q[0];
sx q[0];
rz(-1.3251745) q[0];
sx q[0];
rz(-1.7934196) q[0];
rz(-1.3735636) q[1];
sx q[1];
rz(-1.1496239) q[1];
sx q[1];
rz(-1.9893533) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3593529) q[0];
sx q[0];
rz(-0.77762654) q[0];
sx q[0];
rz(1.9944784) q[0];
rz(-2.466408) q[2];
sx q[2];
rz(-1.0567046) q[2];
sx q[2];
rz(-2.6509283) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.79346953) q[1];
sx q[1];
rz(-1.7382227) q[1];
sx q[1];
rz(2.4504285) q[1];
rz(2.3186734) q[3];
sx q[3];
rz(-2.4429818) q[3];
sx q[3];
rz(0.60958344) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.79933244) q[2];
sx q[2];
rz(-0.41877425) q[2];
sx q[2];
rz(0.93117923) q[2];
rz(3.0350507) q[3];
sx q[3];
rz(-1.1139161) q[3];
sx q[3];
rz(-0.60025269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.057864144) q[0];
sx q[0];
rz(-1.3902384) q[0];
sx q[0];
rz(-2.6237543) q[0];
rz(-2.2528516) q[1];
sx q[1];
rz(-2.4338212) q[1];
sx q[1];
rz(-0.4471561) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.06181051) q[0];
sx q[0];
rz(-1.5070276) q[0];
sx q[0];
rz(2.879854) q[0];
rz(-pi) q[1];
rz(2.1250399) q[2];
sx q[2];
rz(-0.92281658) q[2];
sx q[2];
rz(-0.78314834) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.7501156) q[1];
sx q[1];
rz(-0.98408723) q[1];
sx q[1];
rz(-0.16280414) q[1];
x q[2];
rz(-1.2918043) q[3];
sx q[3];
rz(-1.0615665) q[3];
sx q[3];
rz(1.694547) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.37457028) q[2];
sx q[2];
rz(-2.3775103) q[2];
sx q[2];
rz(2.6089597) q[2];
rz(0.24752188) q[3];
sx q[3];
rz(-2.4041924) q[3];
sx q[3];
rz(-1.0386764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5175051) q[0];
sx q[0];
rz(-0.085260304) q[0];
sx q[0];
rz(3.0349773) q[0];
rz(-0.34635776) q[1];
sx q[1];
rz(-2.296591) q[1];
sx q[1];
rz(-1.9812298) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.896987) q[0];
sx q[0];
rz(-1.3972613) q[0];
sx q[0];
rz(0.24994295) q[0];
rz(2.3290655) q[2];
sx q[2];
rz(-0.90355325) q[2];
sx q[2];
rz(0.086712547) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.44768279) q[1];
sx q[1];
rz(-0.47397787) q[1];
sx q[1];
rz(-2.276027) q[1];
rz(-pi) q[2];
rz(-1.0456234) q[3];
sx q[3];
rz(-1.1743288) q[3];
sx q[3];
rz(0.1016271) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0059263) q[2];
sx q[2];
rz(-2.8432507) q[2];
sx q[2];
rz(-0.076210991) q[2];
rz(-2.5557319) q[3];
sx q[3];
rz(-1.9867089) q[3];
sx q[3];
rz(-1.3425672) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(3.0214486) q[0];
sx q[0];
rz(-1.8360538) q[0];
sx q[0];
rz(-2.7914877) q[0];
rz(2.1856951) q[1];
sx q[1];
rz(-1.8616385) q[1];
sx q[1];
rz(1.9116481) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1317539) q[0];
sx q[0];
rz(-1.3309877) q[0];
sx q[0];
rz(1.9480223) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.036418865) q[2];
sx q[2];
rz(-1.7360188) q[2];
sx q[2];
rz(-0.18319229) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.38795162) q[1];
sx q[1];
rz(-1.3820096) q[1];
sx q[1];
rz(-1.9061879) q[1];
x q[2];
rz(0.98389174) q[3];
sx q[3];
rz(-1.7376889) q[3];
sx q[3];
rz(-1.2323762) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2436287) q[2];
sx q[2];
rz(-2.882759) q[2];
sx q[2];
rz(-1.492307) q[2];
rz(0.40488511) q[3];
sx q[3];
rz(-0.76571524) q[3];
sx q[3];
rz(-2.305472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.006007) q[0];
sx q[0];
rz(-2.0642991) q[0];
sx q[0];
rz(0.58240044) q[0];
rz(-2.4549585) q[1];
sx q[1];
rz(-1.4056987) q[1];
sx q[1];
rz(-1.2581717) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4674038) q[0];
sx q[0];
rz(-1.4454525) q[0];
sx q[0];
rz(2.0095429) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1297087) q[2];
sx q[2];
rz(-2.5311573) q[2];
sx q[2];
rz(-1.287078) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.4542089) q[1];
sx q[1];
rz(-1.5011678) q[1];
sx q[1];
rz(-1.6268262) q[1];
x q[2];
rz(0.21922501) q[3];
sx q[3];
rz(-0.98837438) q[3];
sx q[3];
rz(-1.5150013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.76576343) q[2];
sx q[2];
rz(-1.148369) q[2];
sx q[2];
rz(2.1235535) q[2];
rz(-0.24614075) q[3];
sx q[3];
rz(-1.7459511) q[3];
sx q[3];
rz(-1.6233981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70596424) q[0];
sx q[0];
rz(-0.80393296) q[0];
sx q[0];
rz(-0.58018082) q[0];
rz(2.9980581) q[1];
sx q[1];
rz(-2.6519471) q[1];
sx q[1];
rz(-2.9339583) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26941368) q[0];
sx q[0];
rz(-2.2996443) q[0];
sx q[0];
rz(-0.4704041) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7933664) q[2];
sx q[2];
rz(-2.1219606) q[2];
sx q[2];
rz(-1.4081692) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.8676198) q[1];
sx q[1];
rz(-1.311353) q[1];
sx q[1];
rz(-2.2745489) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8945699) q[3];
sx q[3];
rz(-0.37370703) q[3];
sx q[3];
rz(-1.5242239) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.84270728) q[2];
sx q[2];
rz(-1.2879813) q[2];
sx q[2];
rz(0.60619727) q[2];
rz(-0.68743622) q[3];
sx q[3];
rz(-2.0076553) q[3];
sx q[3];
rz(-2.4351951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6523478) q[0];
sx q[0];
rz(-3.1135961) q[0];
sx q[0];
rz(-2.0835173) q[0];
rz(-3.0319013) q[1];
sx q[1];
rz(-2.0249764) q[1];
sx q[1];
rz(1.441997) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88272754) q[0];
sx q[0];
rz(-1.2611212) q[0];
sx q[0];
rz(-2.726) q[0];
rz(-pi) q[1];
rz(1.2953561) q[2];
sx q[2];
rz(-2.0552962) q[2];
sx q[2];
rz(-0.4415919) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7322526) q[1];
sx q[1];
rz(-2.2870018) q[1];
sx q[1];
rz(2.2316037) q[1];
x q[2];
rz(0.85956802) q[3];
sx q[3];
rz(-2.7894434) q[3];
sx q[3];
rz(-1.5052049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.480392) q[2];
sx q[2];
rz(-0.19442393) q[2];
sx q[2];
rz(1.9332168) q[2];
rz(0.66323534) q[3];
sx q[3];
rz(-1.4570718) q[3];
sx q[3];
rz(-1.2164345) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-0.97063589) q[0];
sx q[0];
rz(-0.96681505) q[0];
sx q[0];
rz(0.092967689) q[0];
rz(1.8611106) q[1];
sx q[1];
rz(-0.72851506) q[1];
sx q[1];
rz(-0.13519898) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2743535) q[0];
sx q[0];
rz(-1.6238191) q[0];
sx q[0];
rz(1.5252602) q[0];
x q[1];
rz(-2.8835758) q[2];
sx q[2];
rz(-1.7191559) q[2];
sx q[2];
rz(2.6515863) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.9433121) q[1];
sx q[1];
rz(-1.3168646) q[1];
sx q[1];
rz(-1.5928245) q[1];
rz(0.074196176) q[3];
sx q[3];
rz(-1.418494) q[3];
sx q[3];
rz(-1.1049113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.70904237) q[2];
sx q[2];
rz(-1.21864) q[2];
sx q[2];
rz(0.77152983) q[2];
rz(2.6500474) q[3];
sx q[3];
rz(-1.2794269) q[3];
sx q[3];
rz(-2.5468723) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58525697) q[0];
sx q[0];
rz(-2.519643) q[0];
sx q[0];
rz(2.0167895) q[0];
rz(1.8428165) q[1];
sx q[1];
rz(-0.61538428) q[1];
sx q[1];
rz(-2.3769456) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1021834) q[0];
sx q[0];
rz(-1.4847857) q[0];
sx q[0];
rz(0.10552222) q[0];
rz(-pi) q[1];
rz(2.2510193) q[2];
sx q[2];
rz(-1.169551) q[2];
sx q[2];
rz(0.97637343) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.1435755) q[1];
sx q[1];
rz(-1.7735266) q[1];
sx q[1];
rz(3.0122117) q[1];
rz(-2.5872676) q[3];
sx q[3];
rz(-1.8968685) q[3];
sx q[3];
rz(2.8248685) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.3760959) q[2];
sx q[2];
rz(-1.8509879) q[2];
sx q[2];
rz(2.9343904) q[2];
rz(0.97992212) q[3];
sx q[3];
rz(-1.1332952) q[3];
sx q[3];
rz(0.19237147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1205263) q[0];
sx q[0];
rz(-1.6854032) q[0];
sx q[0];
rz(2.2717463) q[0];
rz(0.5008685) q[1];
sx q[1];
rz(-2.905838) q[1];
sx q[1];
rz(0.9314608) q[1];
rz(-0.13774112) q[2];
sx q[2];
rz(-1.6888035) q[2];
sx q[2];
rz(2.9464108) q[2];
rz(-0.59400995) q[3];
sx q[3];
rz(-2.8566711) q[3];
sx q[3];
rz(2.2950238) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
