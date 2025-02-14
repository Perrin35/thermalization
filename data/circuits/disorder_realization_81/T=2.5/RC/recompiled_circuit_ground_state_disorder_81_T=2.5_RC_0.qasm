OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.17570198) q[0];
sx q[0];
rz(-1.269729) q[0];
sx q[0];
rz(1.9667392) q[0];
rz(2.0570316) q[1];
sx q[1];
rz(-1.6697872) q[1];
sx q[1];
rz(-0.56301277) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36328735) q[0];
sx q[0];
rz(-1.3749003) q[0];
sx q[0];
rz(-1.5165853) q[0];
rz(-pi) q[1];
rz(-2.1839574) q[2];
sx q[2];
rz(-2.4102825) q[2];
sx q[2];
rz(-2.4201408) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.9809499) q[1];
sx q[1];
rz(-1.8997571) q[1];
sx q[1];
rz(1.739722) q[1];
x q[2];
rz(1.2123763) q[3];
sx q[3];
rz(-1.2120768) q[3];
sx q[3];
rz(-1.7061074) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.53434831) q[2];
sx q[2];
rz(-1.6515942) q[2];
sx q[2];
rz(1.472817) q[2];
rz(1.8913174) q[3];
sx q[3];
rz(-1.6142802) q[3];
sx q[3];
rz(2.1691624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2431353) q[0];
sx q[0];
rz(-0.055483015) q[0];
sx q[0];
rz(0.45795101) q[0];
rz(3.0398439) q[1];
sx q[1];
rz(-1.190217) q[1];
sx q[1];
rz(-0.32624689) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3889623) q[0];
sx q[0];
rz(-1.9174308) q[0];
sx q[0];
rz(1.1244357) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7918705) q[2];
sx q[2];
rz(-1.2794297) q[2];
sx q[2];
rz(1.2779026) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.5811903) q[1];
sx q[1];
rz(-0.96424864) q[1];
sx q[1];
rz(0.62967695) q[1];
rz(1.2539034) q[3];
sx q[3];
rz(-0.85318789) q[3];
sx q[3];
rz(2.9760391) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-3.132167) q[2];
sx q[2];
rz(-2.1475809) q[2];
sx q[2];
rz(-1.2635292) q[2];
rz(3.1215014) q[3];
sx q[3];
rz(-0.76954904) q[3];
sx q[3];
rz(1.2300389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.024071368) q[0];
sx q[0];
rz(-3.0776403) q[0];
sx q[0];
rz(2.5653895) q[0];
rz(1.4441215) q[1];
sx q[1];
rz(-1.8918119) q[1];
sx q[1];
rz(-1.8796657) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44718868) q[0];
sx q[0];
rz(-1.2617636) q[0];
sx q[0];
rz(3.1080327) q[0];
rz(-pi) q[1];
rz(2.0327225) q[2];
sx q[2];
rz(-1.0903768) q[2];
sx q[2];
rz(-1.9146718) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.8518119) q[1];
sx q[1];
rz(-1.8153615) q[1];
sx q[1];
rz(-0.62965392) q[1];
rz(-2.1274873) q[3];
sx q[3];
rz(-0.93645778) q[3];
sx q[3];
rz(1.2403099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.5835517) q[2];
sx q[2];
rz(-0.59528196) q[2];
sx q[2];
rz(-2.9884647) q[2];
rz(2.2554452) q[3];
sx q[3];
rz(-1.359442) q[3];
sx q[3];
rz(-1.2019633) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-1.2578289) q[0];
sx q[0];
rz(-1.9572636) q[0];
sx q[0];
rz(0.15876874) q[0];
rz(-1.0348882) q[1];
sx q[1];
rz(-2.1885469) q[1];
sx q[1];
rz(0.75611702) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0949815) q[0];
sx q[0];
rz(-1.621052) q[0];
sx q[0];
rz(-0.057970164) q[0];
x q[1];
rz(-3.1140585) q[2];
sx q[2];
rz(-1.3176271) q[2];
sx q[2];
rz(2.262676) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3269006) q[1];
sx q[1];
rz(-1.3291825) q[1];
sx q[1];
rz(1.7820226) q[1];
rz(-0.71847312) q[3];
sx q[3];
rz(-2.8333896) q[3];
sx q[3];
rz(-1.1402477) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.5990344) q[2];
sx q[2];
rz(-0.90347806) q[2];
sx q[2];
rz(-2.2912045) q[2];
rz(-1.6709857) q[3];
sx q[3];
rz(-1.3266404) q[3];
sx q[3];
rz(-2.3175122) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65664148) q[0];
sx q[0];
rz(-1.739772) q[0];
sx q[0];
rz(0.17855074) q[0];
rz(-1.8932331) q[1];
sx q[1];
rz(-1.6105885) q[1];
sx q[1];
rz(2.607883) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7631012) q[0];
sx q[0];
rz(-0.39787492) q[0];
sx q[0];
rz(-0.48872013) q[0];
rz(-pi) q[1];
rz(0.49939807) q[2];
sx q[2];
rz(-0.93448567) q[2];
sx q[2];
rz(-2.9884574) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.7099784) q[1];
sx q[1];
rz(-1.5182765) q[1];
sx q[1];
rz(2.2389328) q[1];
rz(1.0664135) q[3];
sx q[3];
rz(-1.6565859) q[3];
sx q[3];
rz(-2.498621) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.3137714) q[2];
sx q[2];
rz(-1.2510108) q[2];
sx q[2];
rz(-2.0884464) q[2];
rz(-0.18038067) q[3];
sx q[3];
rz(-0.80612055) q[3];
sx q[3];
rz(0.36261305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(2.5798222) q[0];
sx q[0];
rz(-1.0012015) q[0];
sx q[0];
rz(1.2286105) q[0];
rz(-2.6749532) q[1];
sx q[1];
rz(-1.9231611) q[1];
sx q[1];
rz(3.012588) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0848964) q[0];
sx q[0];
rz(-1.7040623) q[0];
sx q[0];
rz(3.0396118) q[0];
rz(-pi) q[1];
rz(1.751288) q[2];
sx q[2];
rz(-2.4036088) q[2];
sx q[2];
rz(0.55115095) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.1946757) q[1];
sx q[1];
rz(-1.8652152) q[1];
sx q[1];
rz(0.63398449) q[1];
x q[2];
rz(2.0444851) q[3];
sx q[3];
rz(-2.8420181) q[3];
sx q[3];
rz(-2.4921474) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.0064319) q[2];
sx q[2];
rz(-0.60144037) q[2];
sx q[2];
rz(2.637291) q[2];
rz(1.0935316) q[3];
sx q[3];
rz(-1.3239599) q[3];
sx q[3];
rz(2.1849911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4920376) q[0];
sx q[0];
rz(-1.858359) q[0];
sx q[0];
rz(-1.5564224) q[0];
rz(0.57836142) q[1];
sx q[1];
rz(-1.882694) q[1];
sx q[1];
rz(-0.10231054) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1083173) q[0];
sx q[0];
rz(-2.1731097) q[0];
sx q[0];
rz(0.48745103) q[0];
rz(-pi) q[1];
rz(2.7447276) q[2];
sx q[2];
rz(-1.7831274) q[2];
sx q[2];
rz(3.120829) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.7150517) q[1];
sx q[1];
rz(-2.7965961) q[1];
sx q[1];
rz(0.47899063) q[1];
rz(-pi) q[2];
x q[2];
rz(0.7500521) q[3];
sx q[3];
rz(-2.1382209) q[3];
sx q[3];
rz(-0.22818434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.1954701) q[2];
sx q[2];
rz(-1.8778233) q[2];
sx q[2];
rz(-0.61402399) q[2];
rz(0.92769235) q[3];
sx q[3];
rz(-1.5988908) q[3];
sx q[3];
rz(2.4641002) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55352655) q[0];
sx q[0];
rz(-0.29356846) q[0];
sx q[0];
rz(1.9422096) q[0];
rz(1.9623494) q[1];
sx q[1];
rz(-1.0284547) q[1];
sx q[1];
rz(0.95338043) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48483585) q[0];
sx q[0];
rz(-1.0547045) q[0];
sx q[0];
rz(0.94658312) q[0];
rz(-pi) q[1];
x q[1];
rz(0.10886635) q[2];
sx q[2];
rz(-1.035454) q[2];
sx q[2];
rz(-1.9809198) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9645365) q[1];
sx q[1];
rz(-0.74443084) q[1];
sx q[1];
rz(-1.7339279) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6497506) q[3];
sx q[3];
rz(-1.0903666) q[3];
sx q[3];
rz(1.7460268) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.74932468) q[2];
sx q[2];
rz(-3.0159123) q[2];
sx q[2];
rz(0.18088642) q[2];
rz(2.0285897) q[3];
sx q[3];
rz(-1.0180611) q[3];
sx q[3];
rz(0.38016144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(1.0893294) q[0];
sx q[0];
rz(-2.2398529) q[0];
sx q[0];
rz(2.3774636) q[0];
rz(-2.1185421) q[1];
sx q[1];
rz(-1.6108395) q[1];
sx q[1];
rz(1.4952362) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21432971) q[0];
sx q[0];
rz(-1.8678878) q[0];
sx q[0];
rz(-1.4100177) q[0];
x q[1];
rz(-2.1713397) q[2];
sx q[2];
rz(-1.8489328) q[2];
sx q[2];
rz(0.32626611) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.047840441) q[1];
sx q[1];
rz(-1.4994267) q[1];
sx q[1];
rz(-2.755295) q[1];
x q[2];
rz(-0.71919433) q[3];
sx q[3];
rz(-0.56108863) q[3];
sx q[3];
rz(-0.51855519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.19227795) q[2];
sx q[2];
rz(-2.1239026) q[2];
sx q[2];
rz(-2.6940572) q[2];
rz(-3.0053075) q[3];
sx q[3];
rz(-0.49301967) q[3];
sx q[3];
rz(0.59806699) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3560781) q[0];
sx q[0];
rz(-1.0250174) q[0];
sx q[0];
rz(0.4726952) q[0];
rz(2.7418432) q[1];
sx q[1];
rz(-2.4374297) q[1];
sx q[1];
rz(-2.3230816) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8560573) q[0];
sx q[0];
rz(-1.7489079) q[0];
sx q[0];
rz(-0.31383236) q[0];
rz(-pi) q[1];
rz(-2.632896) q[2];
sx q[2];
rz(-0.74328586) q[2];
sx q[2];
rz(-2.7270779) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(3.1295197) q[1];
sx q[1];
rz(-0.32271472) q[1];
sx q[1];
rz(1.5606776) q[1];
rz(-1.3977526) q[3];
sx q[3];
rz(-1.2283843) q[3];
sx q[3];
rz(2.4416123) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.5588536) q[2];
sx q[2];
rz(-2.0595198) q[2];
sx q[2];
rz(1.8633128) q[2];
rz(1.1996972) q[3];
sx q[3];
rz(-1.4408828) q[3];
sx q[3];
rz(-0.27967683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3661135) q[0];
sx q[0];
rz(-1.7178602) q[0];
sx q[0];
rz(-1.3445509) q[0];
rz(-1.4115502) q[1];
sx q[1];
rz(-0.81137864) q[1];
sx q[1];
rz(2.484533) q[1];
rz(1.5700339) q[2];
sx q[2];
rz(-1.4298202) q[2];
sx q[2];
rz(-0.50411171) q[2];
rz(2.2349002) q[3];
sx q[3];
rz(-1.8421052) q[3];
sx q[3];
rz(2.6143034) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
