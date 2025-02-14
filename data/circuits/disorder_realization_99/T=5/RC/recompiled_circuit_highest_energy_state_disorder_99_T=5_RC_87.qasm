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
rz(1.4003657) q[0];
sx q[0];
rz(-0.23973149) q[0];
sx q[0];
rz(0.32455197) q[0];
rz(-0.95207721) q[1];
sx q[1];
rz(-0.2258741) q[1];
sx q[1];
rz(1.3083375) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1150779) q[0];
sx q[0];
rz(-2.240772) q[0];
sx q[0];
rz(-1.7448533) q[0];
rz(-0.24271528) q[2];
sx q[2];
rz(-2.6449892) q[2];
sx q[2];
rz(-2.9951823) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.89186016) q[1];
sx q[1];
rz(-1.8292055) q[1];
sx q[1];
rz(-1.0904466) q[1];
rz(-pi) q[2];
rz(-1.4300554) q[3];
sx q[3];
rz(-2.2651197) q[3];
sx q[3];
rz(-1.3545881) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7294881) q[2];
sx q[2];
rz(-1.2594014) q[2];
sx q[2];
rz(-0.35120249) q[2];
rz(-1.2900194) q[3];
sx q[3];
rz(-1.131564) q[3];
sx q[3];
rz(-1.9180408) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7555162) q[0];
sx q[0];
rz(-1.2048683) q[0];
sx q[0];
rz(0.93038857) q[0];
rz(-1.8366086) q[1];
sx q[1];
rz(-0.84140673) q[1];
sx q[1];
rz(-2.5321541) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0770438) q[0];
sx q[0];
rz(-1.4676369) q[0];
sx q[0];
rz(-1.1431689) q[0];
x q[1];
rz(-3.0514293) q[2];
sx q[2];
rz(-1.8675641) q[2];
sx q[2];
rz(-2.0903843) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7623065) q[1];
sx q[1];
rz(-1.9348782) q[1];
sx q[1];
rz(-1.9696139) q[1];
rz(-0.16709631) q[3];
sx q[3];
rz(-1.0067938) q[3];
sx q[3];
rz(0.43844863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.096752) q[2];
sx q[2];
rz(-0.98016206) q[2];
sx q[2];
rz(3.0668104) q[2];
rz(-0.52037248) q[3];
sx q[3];
rz(-0.64295355) q[3];
sx q[3];
rz(0.61409942) q[3];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5388913) q[0];
sx q[0];
rz(-0.7170054) q[0];
sx q[0];
rz(-2.8394748) q[0];
rz(-2.865454) q[1];
sx q[1];
rz(-1.3172251) q[1];
sx q[1];
rz(1.709323) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3870704) q[0];
sx q[0];
rz(-1.9517448) q[0];
sx q[0];
rz(-2.1552298) q[0];
rz(-pi) q[1];
rz(1.309607) q[2];
sx q[2];
rz(-0.99203324) q[2];
sx q[2];
rz(-0.28489339) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.773743) q[1];
sx q[1];
rz(-0.58102312) q[1];
sx q[1];
rz(-1.2817205) q[1];
rz(1.7114221) q[3];
sx q[3];
rz(-2.6717477) q[3];
sx q[3];
rz(1.8610561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.3064731) q[2];
sx q[2];
rz(-0.46572954) q[2];
sx q[2];
rz(-1.845537) q[2];
rz(-0.45201388) q[3];
sx q[3];
rz(-2.0343503) q[3];
sx q[3];
rz(0.38255102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(0.52962676) q[0];
sx q[0];
rz(-1.1671678) q[0];
sx q[0];
rz(1.8091328) q[0];
rz(-2.0924856) q[1];
sx q[1];
rz(-2.0499947) q[1];
sx q[1];
rz(-1.5528991) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4047239) q[0];
sx q[0];
rz(-1.8717614) q[0];
sx q[0];
rz(-2.7851581) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.013601549) q[2];
sx q[2];
rz(-2.0715393) q[2];
sx q[2];
rz(2.406943) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.2237593) q[1];
sx q[1];
rz(-0.67549113) q[1];
sx q[1];
rz(1.89066) q[1];
rz(2.2838628) q[3];
sx q[3];
rz(-1.1489043) q[3];
sx q[3];
rz(-0.5798012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.66854746) q[2];
sx q[2];
rz(-1.0205525) q[2];
sx q[2];
rz(0.2505396) q[2];
rz(2.1446832) q[3];
sx q[3];
rz(-2.0018115) q[3];
sx q[3];
rz(0.38058773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32960358) q[0];
sx q[0];
rz(-2.6530837) q[0];
sx q[0];
rz(1.3478152) q[0];
rz(-2.7373121) q[1];
sx q[1];
rz(-0.75899044) q[1];
sx q[1];
rz(-1.1291198) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30325748) q[0];
sx q[0];
rz(-1.1764488) q[0];
sx q[0];
rz(-2.3189937) q[0];
rz(0.26910946) q[2];
sx q[2];
rz(-0.92293149) q[2];
sx q[2];
rz(-1.8108167) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3394062) q[1];
sx q[1];
rz(-1.4660343) q[1];
sx q[1];
rz(0.89749194) q[1];
x q[2];
rz(1.4604767) q[3];
sx q[3];
rz(-1.8838143) q[3];
sx q[3];
rz(-0.68389326) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.79444844) q[2];
sx q[2];
rz(-2.2114387) q[2];
sx q[2];
rz(-1.8974737) q[2];
rz(2.0813023) q[3];
sx q[3];
rz(-0.62842193) q[3];
sx q[3];
rz(-2.8384143) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-2.4853972) q[0];
sx q[0];
rz(-3.0363016) q[0];
sx q[0];
rz(2.7144077) q[0];
rz(0.034491388) q[1];
sx q[1];
rz(-1.5692312) q[1];
sx q[1];
rz(-0.010206612) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1447574) q[0];
sx q[0];
rz(-3.138665) q[0];
sx q[0];
rz(-2.5064431) q[0];
rz(-0.16049196) q[2];
sx q[2];
rz(-2.1969323) q[2];
sx q[2];
rz(-0.4753091) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.2040703) q[1];
sx q[1];
rz(-2.2559705) q[1];
sx q[1];
rz(-3.0584945) q[1];
rz(-pi) q[2];
x q[2];
rz(0.7179435) q[3];
sx q[3];
rz(-0.69980757) q[3];
sx q[3];
rz(1.0216624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.5222142) q[2];
sx q[2];
rz(-2.7284315) q[2];
sx q[2];
rz(-2.348032) q[2];
rz(-0.86269745) q[3];
sx q[3];
rz(-1.5669275) q[3];
sx q[3];
rz(-2.4376552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70235395) q[0];
sx q[0];
rz(-2.2965501) q[0];
sx q[0];
rz(1.1335565) q[0];
rz(1.0776862) q[1];
sx q[1];
rz(-0.51529854) q[1];
sx q[1];
rz(0.30212197) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39224276) q[0];
sx q[0];
rz(-1.2686522) q[0];
sx q[0];
rz(2.7878615) q[0];
rz(-pi) q[1];
rz(-0.59692817) q[2];
sx q[2];
rz(-1.9608627) q[2];
sx q[2];
rz(-0.82118195) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.4430284) q[1];
sx q[1];
rz(-2.5686567) q[1];
sx q[1];
rz(-2.7354476) q[1];
x q[2];
rz(-0.98432912) q[3];
sx q[3];
rz(-2.1696089) q[3];
sx q[3];
rz(2.5564155) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.2332396) q[2];
sx q[2];
rz(-0.88963228) q[2];
sx q[2];
rz(1.7944149) q[2];
rz(-1.8917278) q[3];
sx q[3];
rz(-1.9439387) q[3];
sx q[3];
rz(0.10716042) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5477448) q[0];
sx q[0];
rz(-0.43314728) q[0];
sx q[0];
rz(3.0331842) q[0];
rz(-2.0938865) q[1];
sx q[1];
rz(-2.3240418) q[1];
sx q[1];
rz(2.122706) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86025809) q[0];
sx q[0];
rz(-2.0906013) q[0];
sx q[0];
rz(0.61183527) q[0];
rz(-pi) q[1];
rz(-1.564304) q[2];
sx q[2];
rz(-2.3320856) q[2];
sx q[2];
rz(-2.6463375) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.093204) q[1];
sx q[1];
rz(-1.1133514) q[1];
sx q[1];
rz(2.0162575) q[1];
rz(-pi) q[2];
rz(-0.15411994) q[3];
sx q[3];
rz(-0.74609038) q[3];
sx q[3];
rz(-2.8733159) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.63742796) q[2];
sx q[2];
rz(-0.97680682) q[2];
sx q[2];
rz(0.43761474) q[2];
rz(0.49312433) q[3];
sx q[3];
rz(-2.2212641) q[3];
sx q[3];
rz(-1.5776618) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2757932) q[0];
sx q[0];
rz(-1.2340622) q[0];
sx q[0];
rz(0.27780521) q[0];
rz(-1.5996784) q[1];
sx q[1];
rz(-0.96822396) q[1];
sx q[1];
rz(-0.42617282) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0539249) q[0];
sx q[0];
rz(-3.1172981) q[0];
sx q[0];
rz(2.3870275) q[0];
rz(-pi) q[1];
x q[1];
rz(2.369719) q[2];
sx q[2];
rz(-1.0091595) q[2];
sx q[2];
rz(1.2489212) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.26697576) q[1];
sx q[1];
rz(-2.8966581) q[1];
sx q[1];
rz(2.0744978) q[1];
x q[2];
rz(0.73240273) q[3];
sx q[3];
rz(-0.97627813) q[3];
sx q[3];
rz(1.5624832) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.6466732) q[2];
sx q[2];
rz(-2.9743331) q[2];
sx q[2];
rz(-1.0529998) q[2];
rz(-2.7827175) q[3];
sx q[3];
rz(-1.5263298) q[3];
sx q[3];
rz(-1.0927965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36295715) q[0];
sx q[0];
rz(-0.21720049) q[0];
sx q[0];
rz(2.2208075) q[0];
rz(2.6329363) q[1];
sx q[1];
rz(-2.3743036) q[1];
sx q[1];
rz(-2.7210534) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32146183) q[0];
sx q[0];
rz(-2.252914) q[0];
sx q[0];
rz(0.90773037) q[0];
rz(2.8996572) q[2];
sx q[2];
rz(-0.95429776) q[2];
sx q[2];
rz(2.7665507) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.9982165) q[1];
sx q[1];
rz(-2.3399379) q[1];
sx q[1];
rz(-2.2364535) q[1];
x q[2];
rz(2.2617662) q[3];
sx q[3];
rz(-1.4925957) q[3];
sx q[3];
rz(-1.9281338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.73021182) q[2];
sx q[2];
rz(-2.3513942) q[2];
sx q[2];
rz(0.31663695) q[2];
rz(0.5591048) q[3];
sx q[3];
rz(-2.5059301) q[3];
sx q[3];
rz(2.3234698) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6560787) q[0];
sx q[0];
rz(-2.06388) q[0];
sx q[0];
rz(-2.0583454) q[0];
rz(-0.47646933) q[1];
sx q[1];
rz(-2.10119) q[1];
sx q[1];
rz(1.4269921) q[1];
rz(-0.086924788) q[2];
sx q[2];
rz(-0.82155052) q[2];
sx q[2];
rz(-0.62053298) q[2];
rz(-1.9540167) q[3];
sx q[3];
rz(-1.4837617) q[3];
sx q[3];
rz(-1.0401534) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
