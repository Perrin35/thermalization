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
rz(0.90647107) q[0];
sx q[0];
rz(-1.4426761) q[0];
sx q[0];
rz(0.28191167) q[0];
rz(0.52892041) q[1];
sx q[1];
rz(4.79098) q[1];
sx q[1];
rz(11.00287) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0391386) q[0];
sx q[0];
rz(-0.41101563) q[0];
sx q[0];
rz(-1.1059815) q[0];
rz(0.27484244) q[2];
sx q[2];
rz(-1.5656398) q[2];
sx q[2];
rz(-0.33493638) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3696182) q[1];
sx q[1];
rz(-1.4450412) q[1];
sx q[1];
rz(-2.5096276) q[1];
rz(2.7759477) q[3];
sx q[3];
rz(-2.4198101) q[3];
sx q[3];
rz(0.66058285) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.4549183) q[2];
sx q[2];
rz(-2.8439971) q[2];
sx q[2];
rz(-2.5285524) q[2];
rz(-0.4736627) q[3];
sx q[3];
rz(-1.9434171) q[3];
sx q[3];
rz(1.7090428) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
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
rz(-1.4050201) q[0];
sx q[0];
rz(-2.9614145) q[0];
sx q[0];
rz(-2.3905684) q[0];
rz(-0.48149064) q[1];
sx q[1];
rz(-1.0844237) q[1];
sx q[1];
rz(-0.96985936) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50502993) q[0];
sx q[0];
rz(-1.2004932) q[0];
sx q[0];
rz(0.67815336) q[0];
rz(-pi) q[1];
rz(-1.5300691) q[2];
sx q[2];
rz(-1.3551522) q[2];
sx q[2];
rz(3.0037896) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.6578601) q[1];
sx q[1];
rz(-0.93440234) q[1];
sx q[1];
rz(-0.32245335) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0131575) q[3];
sx q[3];
rz(-1.0201038) q[3];
sx q[3];
rz(0.1503508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.91886175) q[2];
sx q[2];
rz(-1.6657882) q[2];
sx q[2];
rz(0.53885031) q[2];
rz(-3.1239964) q[3];
sx q[3];
rz(-0.16418695) q[3];
sx q[3];
rz(-2.2331451) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4410412) q[0];
sx q[0];
rz(-2.7563162) q[0];
sx q[0];
rz(3.0803296) q[0];
rz(-1.5047269) q[1];
sx q[1];
rz(-2.442339) q[1];
sx q[1];
rz(1.1963074) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17249566) q[0];
sx q[0];
rz(-2.3845256) q[0];
sx q[0];
rz(-1.2189381) q[0];
x q[1];
rz(-0.26232403) q[2];
sx q[2];
rz(-1.6379426) q[2];
sx q[2];
rz(0.1131499) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.34431008) q[1];
sx q[1];
rz(-2.4786665) q[1];
sx q[1];
rz(-1.0271038) q[1];
rz(-pi) q[2];
rz(-1.9798093) q[3];
sx q[3];
rz(-0.73506415) q[3];
sx q[3];
rz(-2.4215557) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.68941826) q[2];
sx q[2];
rz(-1.8345366) q[2];
sx q[2];
rz(0.43251953) q[2];
rz(-2.050926) q[3];
sx q[3];
rz(-2.5011823) q[3];
sx q[3];
rz(-2.045491) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5890305) q[0];
sx q[0];
rz(-0.93695372) q[0];
sx q[0];
rz(-0.20137782) q[0];
rz(-2.6248113) q[1];
sx q[1];
rz(-0.38025451) q[1];
sx q[1];
rz(-2.1929599) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3654914) q[0];
sx q[0];
rz(-2.1946215) q[0];
sx q[0];
rz(-2.7073949) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2791512) q[2];
sx q[2];
rz(-1.7956327) q[2];
sx q[2];
rz(2.1612957) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3740179) q[1];
sx q[1];
rz(-1.5841055) q[1];
sx q[1];
rz(-2.637251) q[1];
x q[2];
rz(0.52093769) q[3];
sx q[3];
rz(-1.4864941) q[3];
sx q[3];
rz(-1.4985633) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.2739233) q[2];
sx q[2];
rz(-0.40931585) q[2];
sx q[2];
rz(2.5841827) q[2];
rz(-0.080816001) q[3];
sx q[3];
rz(-1.9010474) q[3];
sx q[3];
rz(1.2062937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7252561) q[0];
sx q[0];
rz(-1.4755604) q[0];
sx q[0];
rz(-2.1112554) q[0];
rz(-2.985785) q[1];
sx q[1];
rz(-2.0741597) q[1];
sx q[1];
rz(-0.20419289) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1847398) q[0];
sx q[0];
rz(-2.889688) q[0];
sx q[0];
rz(0.27916081) q[0];
rz(-pi) q[1];
x q[1];
rz(0.5459909) q[2];
sx q[2];
rz(-0.24790774) q[2];
sx q[2];
rz(-1.1746097) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0818644) q[1];
sx q[1];
rz(-1.7769741) q[1];
sx q[1];
rz(-2.6463406) q[1];
rz(-2.3685826) q[3];
sx q[3];
rz(-1.3979619) q[3];
sx q[3];
rz(-2.5484249) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.2121409) q[2];
sx q[2];
rz(-1.4931623) q[2];
sx q[2];
rz(0.41352752) q[2];
rz(-2.2996969) q[3];
sx q[3];
rz(-1.7588408) q[3];
sx q[3];
rz(2.1690185) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.095116422) q[0];
sx q[0];
rz(-2.0501417) q[0];
sx q[0];
rz(0.0055775642) q[0];
rz(-1.7976409) q[1];
sx q[1];
rz(-0.19897142) q[1];
sx q[1];
rz(0.70095789) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0049135) q[0];
sx q[0];
rz(-2.4304996) q[0];
sx q[0];
rz(0.85791608) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1326601) q[2];
sx q[2];
rz(-0.72266662) q[2];
sx q[2];
rz(0.36107963) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.720916) q[1];
sx q[1];
rz(-1.6117967) q[1];
sx q[1];
rz(0.27537217) q[1];
rz(-pi) q[2];
x q[2];
rz(0.72356059) q[3];
sx q[3];
rz(-0.71092194) q[3];
sx q[3];
rz(-1.5953181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.31412101) q[2];
sx q[2];
rz(-0.33385971) q[2];
sx q[2];
rz(1.1673048) q[2];
rz(0.88268924) q[3];
sx q[3];
rz(-2.3736931) q[3];
sx q[3];
rz(0.32514969) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(-2.2348787) q[0];
sx q[0];
rz(-1.119708) q[0];
sx q[0];
rz(-3.0562905) q[0];
rz(2.3872497) q[1];
sx q[1];
rz(-0.5032379) q[1];
sx q[1];
rz(-1.0275966) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9011544) q[0];
sx q[0];
rz(-1.5605956) q[0];
sx q[0];
rz(-0.21252327) q[0];
x q[1];
rz(-1.5761887) q[2];
sx q[2];
rz(-0.9741592) q[2];
sx q[2];
rz(-0.83534681) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.8366521) q[1];
sx q[1];
rz(-1.8936367) q[1];
sx q[1];
rz(-3.0129827) q[1];
rz(2.1291162) q[3];
sx q[3];
rz(-1.1577679) q[3];
sx q[3];
rz(-2.922204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.9895642) q[2];
sx q[2];
rz(-0.98735183) q[2];
sx q[2];
rz(-0.41934553) q[2];
rz(-0.63465214) q[3];
sx q[3];
rz(-0.80497634) q[3];
sx q[3];
rz(-2.0161207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80970508) q[0];
sx q[0];
rz(-0.87600791) q[0];
sx q[0];
rz(-0.69212717) q[0];
rz(-0.57811111) q[1];
sx q[1];
rz(-2.7222241) q[1];
sx q[1];
rz(1.3828166) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3457926) q[0];
sx q[0];
rz(-1.5032856) q[0];
sx q[0];
rz(0.019492143) q[0];
x q[1];
rz(-0.67202576) q[2];
sx q[2];
rz(-2.1161072) q[2];
sx q[2];
rz(1.3915075) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.2410215) q[1];
sx q[1];
rz(-1.8041148) q[1];
sx q[1];
rz(2.1395348) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0973849) q[3];
sx q[3];
rz(-2.3447737) q[3];
sx q[3];
rz(2.6897893) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3708923) q[2];
sx q[2];
rz(-2.6110677) q[2];
sx q[2];
rz(1.4250866) q[2];
rz(2.6254081) q[3];
sx q[3];
rz(-0.97847146) q[3];
sx q[3];
rz(-1.2396575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4952963) q[0];
sx q[0];
rz(-1.7072059) q[0];
sx q[0];
rz(0.37149757) q[0];
rz(0.83483541) q[1];
sx q[1];
rz(-2.3284262) q[1];
sx q[1];
rz(0.017597839) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71401063) q[0];
sx q[0];
rz(-1.0508063) q[0];
sx q[0];
rz(-2.8218357) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8210635) q[2];
sx q[2];
rz(-0.8508209) q[2];
sx q[2];
rz(-1.7257476) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.0068215) q[1];
sx q[1];
rz(-2.4398514) q[1];
sx q[1];
rz(-2.546258) q[1];
rz(-pi) q[2];
rz(0.44638388) q[3];
sx q[3];
rz(-2.4842815) q[3];
sx q[3];
rz(0.039358649) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.1060433) q[2];
sx q[2];
rz(-0.86842662) q[2];
sx q[2];
rz(0.10031984) q[2];
rz(2.912168) q[3];
sx q[3];
rz(-1.0474297) q[3];
sx q[3];
rz(2.6917698) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68596524) q[0];
sx q[0];
rz(-2.8534511) q[0];
sx q[0];
rz(2.567754) q[0];
rz(-2.0876743) q[1];
sx q[1];
rz(-2.0951447) q[1];
sx q[1];
rz(0.67646772) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8343) q[0];
sx q[0];
rz(-0.3737475) q[0];
sx q[0];
rz(2.4227002) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8847938) q[2];
sx q[2];
rz(-1.1780103) q[2];
sx q[2];
rz(-1.8671672) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.0645731) q[1];
sx q[1];
rz(-1.6592245) q[1];
sx q[1];
rz(0.52845533) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1031688) q[3];
sx q[3];
rz(-1.7422973) q[3];
sx q[3];
rz(0.72491437) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.68710589) q[2];
sx q[2];
rz(-2.3837619) q[2];
sx q[2];
rz(-1.494361) q[2];
rz(-2.2729661) q[3];
sx q[3];
rz(-2.4354911) q[3];
sx q[3];
rz(0.47006616) q[3];
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
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0710707) q[0];
sx q[0];
rz(-1.1916397) q[0];
sx q[0];
rz(1.2217039) q[0];
rz(-1.3337878) q[1];
sx q[1];
rz(-1.8232657) q[1];
sx q[1];
rz(2.5008536) q[1];
rz(1.6619353) q[2];
sx q[2];
rz(-2.3341134) q[2];
sx q[2];
rz(-2.3011617) q[2];
rz(-0.044338772) q[3];
sx q[3];
rz(-1.3953184) q[3];
sx q[3];
rz(0.75147006) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
