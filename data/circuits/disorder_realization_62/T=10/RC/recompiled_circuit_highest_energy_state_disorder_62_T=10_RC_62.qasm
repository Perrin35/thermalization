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
rz(-0.57617968) q[0];
sx q[0];
rz(-1.4644858) q[0];
sx q[0];
rz(-0.65281868) q[0];
rz(-0.44261143) q[1];
sx q[1];
rz(-2.0191329) q[1];
sx q[1];
rz(5.6607487) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63875893) q[0];
sx q[0];
rz(-1.2593102) q[0];
sx q[0];
rz(-3.019096) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9046225) q[2];
sx q[2];
rz(-0.99572748) q[2];
sx q[2];
rz(-1.0525557) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.6316526) q[1];
sx q[1];
rz(-1.3871683) q[1];
sx q[1];
rz(-1.0650403) q[1];
rz(-pi) q[2];
x q[2];
rz(0.15776026) q[3];
sx q[3];
rz(-2.1737242) q[3];
sx q[3];
rz(-1.2186528) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.5114674) q[2];
sx q[2];
rz(-1.3518535) q[2];
sx q[2];
rz(-2.5959065) q[2];
rz(2.4557579) q[3];
sx q[3];
rz(-2.4617709) q[3];
sx q[3];
rz(2.1849476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0551374) q[0];
sx q[0];
rz(-0.70239037) q[0];
sx q[0];
rz(2.0461244) q[0];
rz(0.99758863) q[1];
sx q[1];
rz(-1.9065403) q[1];
sx q[1];
rz(-1.9445317) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67162469) q[0];
sx q[0];
rz(-1.3177682) q[0];
sx q[0];
rz(-2.1718911) q[0];
rz(-pi) q[1];
rz(2.3648974) q[2];
sx q[2];
rz(-1.1001081) q[2];
sx q[2];
rz(1.7913417) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8833739) q[1];
sx q[1];
rz(-2.0257468) q[1];
sx q[1];
rz(-0.31724288) q[1];
rz(-pi) q[2];
x q[2];
rz(0.32917804) q[3];
sx q[3];
rz(-1.7478895) q[3];
sx q[3];
rz(2.8383423) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.41584388) q[2];
sx q[2];
rz(-1.3943322) q[2];
sx q[2];
rz(1.7572629) q[2];
rz(1.3437126) q[3];
sx q[3];
rz(-1.5735156) q[3];
sx q[3];
rz(1.2725007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
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
rz(2.5196853) q[0];
sx q[0];
rz(-2.3540731) q[0];
sx q[0];
rz(-2.7144879) q[0];
rz(-1.9097795) q[1];
sx q[1];
rz(-1.6424664) q[1];
sx q[1];
rz(-0.61418358) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6525567) q[0];
sx q[0];
rz(-1.8985735) q[0];
sx q[0];
rz(-1.2851738) q[0];
x q[1];
rz(-1.7730447) q[2];
sx q[2];
rz(-2.3073811) q[2];
sx q[2];
rz(1.5327765) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.8607448) q[1];
sx q[1];
rz(-0.13361803) q[1];
sx q[1];
rz(1.0722776) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7216629) q[3];
sx q[3];
rz(-2.5961317) q[3];
sx q[3];
rz(2.2614516) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.4517639) q[2];
sx q[2];
rz(-0.58388766) q[2];
sx q[2];
rz(2.8705987) q[2];
rz(0.48218918) q[3];
sx q[3];
rz(-0.8420344) q[3];
sx q[3];
rz(-1.2838001) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2778306) q[0];
sx q[0];
rz(-1.3348802) q[0];
sx q[0];
rz(-0.41896391) q[0];
rz(-0.22467443) q[1];
sx q[1];
rz(-1.9326262) q[1];
sx q[1];
rz(3.0640501) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52151742) q[0];
sx q[0];
rz(-1.8635611) q[0];
sx q[0];
rz(1.0241246) q[0];
rz(-2.783804) q[2];
sx q[2];
rz(-2.3364106) q[2];
sx q[2];
rz(-1.5186269) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.3540707) q[1];
sx q[1];
rz(-1.7131299) q[1];
sx q[1];
rz(-2.700129) q[1];
x q[2];
rz(2.5285012) q[3];
sx q[3];
rz(-0.35260751) q[3];
sx q[3];
rz(-2.4168454) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-3.0858687) q[2];
sx q[2];
rz(-2.8102977) q[2];
sx q[2];
rz(-1.2223988) q[2];
rz(-0.27082768) q[3];
sx q[3];
rz(-1.2374249) q[3];
sx q[3];
rz(-2.6868611) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5643352) q[0];
sx q[0];
rz(-0.99522796) q[0];
sx q[0];
rz(1.7876392) q[0];
rz(2.1063781) q[1];
sx q[1];
rz(-1.926492) q[1];
sx q[1];
rz(-1.6114906) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28883067) q[0];
sx q[0];
rz(-2.0977328) q[0];
sx q[0];
rz(2.2599392) q[0];
rz(-pi) q[1];
x q[1];
rz(0.3283073) q[2];
sx q[2];
rz(-1.7264778) q[2];
sx q[2];
rz(0.98634431) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.51401201) q[1];
sx q[1];
rz(-2.7360533) q[1];
sx q[1];
rz(-3.0058014) q[1];
x q[2];
rz(3.0102481) q[3];
sx q[3];
rz(-0.41413158) q[3];
sx q[3];
rz(1.1826178) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.3349541) q[2];
sx q[2];
rz(-0.97794509) q[2];
sx q[2];
rz(0.063892603) q[2];
rz(2.4334) q[3];
sx q[3];
rz(-1.6876561) q[3];
sx q[3];
rz(-1.8074869) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42234364) q[0];
sx q[0];
rz(-3.1184734) q[0];
sx q[0];
rz(-2.0656021) q[0];
rz(-0.05323449) q[1];
sx q[1];
rz(-2.2884171) q[1];
sx q[1];
rz(-2.7630973) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.196825) q[0];
sx q[0];
rz(-1.7009957) q[0];
sx q[0];
rz(2.8848335) q[0];
rz(-pi) q[1];
rz(-1.4827864) q[2];
sx q[2];
rz(-1.529379) q[2];
sx q[2];
rz(-2.8633022) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.82847257) q[1];
sx q[1];
rz(-2.8250029) q[1];
sx q[1];
rz(2.5324138) q[1];
rz(-pi) q[2];
x q[2];
rz(0.58121292) q[3];
sx q[3];
rz(-0.85678116) q[3];
sx q[3];
rz(0.17457822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.1899015) q[2];
sx q[2];
rz(-1.4611763) q[2];
sx q[2];
rz(-1.0339197) q[2];
rz(-0.90384358) q[3];
sx q[3];
rz(-2.240182) q[3];
sx q[3];
rz(0.044053642) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1164301) q[0];
sx q[0];
rz(-1.4465541) q[0];
sx q[0];
rz(-0.59301162) q[0];
rz(-3.016839) q[1];
sx q[1];
rz(-1.1589103) q[1];
sx q[1];
rz(0.4932901) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.497555) q[0];
sx q[0];
rz(-1.293566) q[0];
sx q[0];
rz(1.2406741) q[0];
rz(-1.5957082) q[2];
sx q[2];
rz(-1.3427715) q[2];
sx q[2];
rz(-1.4268503) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.3354937) q[1];
sx q[1];
rz(-0.68365288) q[1];
sx q[1];
rz(2.4324576) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4629355) q[3];
sx q[3];
rz(-1.8450122) q[3];
sx q[3];
rz(-1.3317684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.3844246) q[2];
sx q[2];
rz(-1.1606471) q[2];
sx q[2];
rz(-1.9942795) q[2];
rz(0.82289639) q[3];
sx q[3];
rz(-1.9673248) q[3];
sx q[3];
rz(-2.4790922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61854521) q[0];
sx q[0];
rz(-1.86919) q[0];
sx q[0];
rz(-2.7230895) q[0];
rz(-0.11323994) q[1];
sx q[1];
rz(-1.6721882) q[1];
sx q[1];
rz(-1.07771) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0200054) q[0];
sx q[0];
rz(-1.5527871) q[0];
sx q[0];
rz(1.2822582) q[0];
rz(-pi) q[1];
rz(1.8968717) q[2];
sx q[2];
rz(-0.80279826) q[2];
sx q[2];
rz(1.3981512) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.4227161) q[1];
sx q[1];
rz(-0.36564974) q[1];
sx q[1];
rz(1.4332194) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5111132) q[3];
sx q[3];
rz(-2.2264997) q[3];
sx q[3];
rz(-2.7969517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.8354127) q[2];
sx q[2];
rz(-2.4330008) q[2];
sx q[2];
rz(-1.4554679) q[2];
rz(-2.9742187) q[3];
sx q[3];
rz(-0.51515976) q[3];
sx q[3];
rz(1.0719871) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6587332) q[0];
sx q[0];
rz(-1.8287683) q[0];
sx q[0];
rz(2.7400548) q[0];
rz(2.7032779) q[1];
sx q[1];
rz(-0.79151789) q[1];
sx q[1];
rz(1.4567136) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4493713) q[0];
sx q[0];
rz(-1.6047932) q[0];
sx q[0];
rz(-0.0095937455) q[0];
x q[1];
rz(-2.8598665) q[2];
sx q[2];
rz(-2.4664219) q[2];
sx q[2];
rz(0.074568579) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.46692013) q[1];
sx q[1];
rz(-2.2385257) q[1];
sx q[1];
rz(0.93514644) q[1];
rz(-1.0965804) q[3];
sx q[3];
rz(-1.1062396) q[3];
sx q[3];
rz(-1.43917) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.7598286) q[2];
sx q[2];
rz(-2.9324053) q[2];
sx q[2];
rz(-0.71167243) q[2];
rz(0.034218637) q[3];
sx q[3];
rz(-0.70845571) q[3];
sx q[3];
rz(-1.9862407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3592247) q[0];
sx q[0];
rz(-1.6509667) q[0];
sx q[0];
rz(2.3427298) q[0];
rz(1.0460188) q[1];
sx q[1];
rz(-1.9452399) q[1];
sx q[1];
rz(2.9050713) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49451429) q[0];
sx q[0];
rz(-3.0768659) q[0];
sx q[0];
rz(1.9231803) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2138519) q[2];
sx q[2];
rz(-1.0272214) q[2];
sx q[2];
rz(-1.7801931) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0482208) q[1];
sx q[1];
rz(-1.908675) q[1];
sx q[1];
rz(0.45856234) q[1];
rz(-pi) q[2];
x q[2];
rz(0.13670758) q[3];
sx q[3];
rz(-2.4259704) q[3];
sx q[3];
rz(-0.8807883) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.025042621) q[2];
sx q[2];
rz(-2.5849085) q[2];
sx q[2];
rz(2.4994948) q[2];
rz(2.5721278) q[3];
sx q[3];
rz(-2.6537708) q[3];
sx q[3];
rz(-0.08629442) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4702598) q[0];
sx q[0];
rz(-1.8076121) q[0];
sx q[0];
rz(-2.1847834) q[0];
rz(-1.9500465) q[1];
sx q[1];
rz(-1.5793431) q[1];
sx q[1];
rz(1.6021077) q[1];
rz(0.79193476) q[2];
sx q[2];
rz(-1.4832693) q[2];
sx q[2];
rz(-2.5434189) q[2];
rz(-0.91611923) q[3];
sx q[3];
rz(-2.6226433) q[3];
sx q[3];
rz(-0.15437689) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
