OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.9590149) q[0];
sx q[0];
rz(-2.9380517) q[0];
sx q[0];
rz(-1.8939053) q[0];
rz(-1.7893451) q[1];
sx q[1];
rz(-2.4440553) q[1];
sx q[1];
rz(0.69812671) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9331935) q[0];
sx q[0];
rz(-1.0675479) q[0];
sx q[0];
rz(-1.1017975) q[0];
rz(-pi) q[1];
rz(3.1280367) q[2];
sx q[2];
rz(-0.45326158) q[2];
sx q[2];
rz(2.1238985) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.9351077) q[1];
sx q[1];
rz(-0.84151959) q[1];
sx q[1];
rz(-2.7323099) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.1269253) q[3];
sx q[3];
rz(-1.3157744) q[3];
sx q[3];
rz(-1.5262814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.84833604) q[2];
sx q[2];
rz(-1.9170599) q[2];
sx q[2];
rz(2.2621034) q[2];
rz(0.73421156) q[3];
sx q[3];
rz(-2.0101533) q[3];
sx q[3];
rz(-0.82936275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5010928) q[0];
sx q[0];
rz(-1.9556029) q[0];
sx q[0];
rz(-2.1511141) q[0];
rz(0.99539202) q[1];
sx q[1];
rz(-1.4442911) q[1];
sx q[1];
rz(-0.46517864) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.746691) q[0];
sx q[0];
rz(-1.1550265) q[0];
sx q[0];
rz(1.3711934) q[0];
rz(-pi) q[1];
x q[1];
rz(0.37613531) q[2];
sx q[2];
rz(-2.1636164) q[2];
sx q[2];
rz(1.7629634) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.097660573) q[1];
sx q[1];
rz(-2.2273835) q[1];
sx q[1];
rz(2.2022579) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4898759) q[3];
sx q[3];
rz(-2.8164125) q[3];
sx q[3];
rz(-2.1212494) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5140932) q[2];
sx q[2];
rz(-1.9072615) q[2];
sx q[2];
rz(-2.9631183) q[2];
rz(0.56525362) q[3];
sx q[3];
rz(-1.3924799) q[3];
sx q[3];
rz(-0.55725151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31767118) q[0];
sx q[0];
rz(-1.3452106) q[0];
sx q[0];
rz(2.4406261) q[0];
rz(2.7340381) q[1];
sx q[1];
rz(-1.2113672) q[1];
sx q[1];
rz(1.0754546) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70301818) q[0];
sx q[0];
rz(-1.0319627) q[0];
sx q[0];
rz(-1.9608905) q[0];
x q[1];
rz(-1.2811866) q[2];
sx q[2];
rz(-0.67669213) q[2];
sx q[2];
rz(-1.8829568) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.13710626) q[1];
sx q[1];
rz(-1.1710694) q[1];
sx q[1];
rz(0.14182472) q[1];
rz(-pi) q[2];
rz(-0.050450669) q[3];
sx q[3];
rz(-1.8352494) q[3];
sx q[3];
rz(-1.3095625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.0740697) q[2];
sx q[2];
rz(-1.3347551) q[2];
sx q[2];
rz(1.4074116) q[2];
rz(2.6635026) q[3];
sx q[3];
rz(-2.7985088) q[3];
sx q[3];
rz(-0.34496719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(2.5673229) q[0];
sx q[0];
rz(-1.3936717) q[0];
sx q[0];
rz(-2.0950914) q[0];
rz(0.25431713) q[1];
sx q[1];
rz(-0.81552243) q[1];
sx q[1];
rz(-0.3124803) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0346261) q[0];
sx q[0];
rz(-2.7101332) q[0];
sx q[0];
rz(1.2564684) q[0];
x q[1];
rz(0.58860345) q[2];
sx q[2];
rz(-0.75102511) q[2];
sx q[2];
rz(-2.5680755) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8487843) q[1];
sx q[1];
rz(-2.7393746) q[1];
sx q[1];
rz(2.4821698) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0279827) q[3];
sx q[3];
rz(-2.5590754) q[3];
sx q[3];
rz(1.5912671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3089402) q[2];
sx q[2];
rz(-2.4231484) q[2];
sx q[2];
rz(-0.18860513) q[2];
rz(-0.23379937) q[3];
sx q[3];
rz(-0.89582396) q[3];
sx q[3];
rz(2.5578267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
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
rz(0.69498649) q[0];
sx q[0];
rz(-1.8291031) q[0];
sx q[0];
rz(3.1332698) q[0];
rz(2.7024929) q[1];
sx q[1];
rz(-1.3689901) q[1];
sx q[1];
rz(1.2459374) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3942892) q[0];
sx q[0];
rz(-0.4957333) q[0];
sx q[0];
rz(-0.41252211) q[0];
rz(-pi) q[1];
rz(1.3338842) q[2];
sx q[2];
rz(-1.4258175) q[2];
sx q[2];
rz(-1.0558093) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.627927) q[1];
sx q[1];
rz(-2.6509319) q[1];
sx q[1];
rz(-0.8497683) q[1];
rz(-pi) q[2];
rz(0.28723081) q[3];
sx q[3];
rz(-0.93614791) q[3];
sx q[3];
rz(-1.3040445) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.7370854) q[2];
sx q[2];
rz(-0.046701996) q[2];
sx q[2];
rz(-1.5076293) q[2];
rz(-0.59219939) q[3];
sx q[3];
rz(-1.7630354) q[3];
sx q[3];
rz(-0.73970214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.020788766) q[0];
sx q[0];
rz(-1.4494267) q[0];
sx q[0];
rz(2.8756496) q[0];
rz(1.9256915) q[1];
sx q[1];
rz(-0.78548702) q[1];
sx q[1];
rz(1.9507834) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5191089) q[0];
sx q[0];
rz(-1.568123) q[0];
sx q[0];
rz(-0.83774211) q[0];
x q[1];
rz(1.0180151) q[2];
sx q[2];
rz(-1.5399122) q[2];
sx q[2];
rz(0.72252233) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.3490847) q[1];
sx q[1];
rz(-2.324632) q[1];
sx q[1];
rz(-1.5813827) q[1];
rz(-pi) q[2];
rz(-2.9407752) q[3];
sx q[3];
rz(-0.28017032) q[3];
sx q[3];
rz(-2.7602989) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.33209458) q[2];
sx q[2];
rz(-1.7569434) q[2];
sx q[2];
rz(0.14796999) q[2];
rz(-2.68908) q[3];
sx q[3];
rz(-1.0871525) q[3];
sx q[3];
rz(0.40542671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0026523503) q[0];
sx q[0];
rz(-1.2661221) q[0];
sx q[0];
rz(1.7953605) q[0];
rz(-3.1345308) q[1];
sx q[1];
rz(-0.66771737) q[1];
sx q[1];
rz(2.6341569) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1333727) q[0];
sx q[0];
rz(-2.3485561) q[0];
sx q[0];
rz(0.80312463) q[0];
rz(1.2260004) q[2];
sx q[2];
rz(-1.424768) q[2];
sx q[2];
rz(-0.70521077) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.7120613) q[1];
sx q[1];
rz(-1.1768394) q[1];
sx q[1];
rz(1.5809356) q[1];
x q[2];
rz(1.2984914) q[3];
sx q[3];
rz(-1.57988) q[3];
sx q[3];
rz(-0.42382012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.70600447) q[2];
sx q[2];
rz(-1.0554577) q[2];
sx q[2];
rz(2.9662507) q[2];
rz(0.38241479) q[3];
sx q[3];
rz(-1.9491111) q[3];
sx q[3];
rz(-0.46149883) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9984556) q[0];
sx q[0];
rz(-2.1303506) q[0];
sx q[0];
rz(1.829041) q[0];
rz(3.0328499) q[1];
sx q[1];
rz(-2.5332632) q[1];
sx q[1];
rz(2.463602) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0339542) q[0];
sx q[0];
rz(-1.0773412) q[0];
sx q[0];
rz(-0.02520589) q[0];
rz(-pi) q[1];
x q[1];
rz(0.021804734) q[2];
sx q[2];
rz(-0.74912375) q[2];
sx q[2];
rz(2.5072012) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.7398518) q[1];
sx q[1];
rz(-2.0900407) q[1];
sx q[1];
rz(-0.3335764) q[1];
x q[2];
rz(-2.4245913) q[3];
sx q[3];
rz(-1.1537408) q[3];
sx q[3];
rz(1.6959697) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.8678191) q[2];
sx q[2];
rz(-1.5673693) q[2];
sx q[2];
rz(-1.1567953) q[2];
rz(2.6371238) q[3];
sx q[3];
rz(-1.0320832) q[3];
sx q[3];
rz(0.57197905) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93541637) q[0];
sx q[0];
rz(-1.529006) q[0];
sx q[0];
rz(0.61844283) q[0];
rz(2.2930324) q[1];
sx q[1];
rz(-2.6237374) q[1];
sx q[1];
rz(-2.3458164) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1815565) q[0];
sx q[0];
rz(-2.040876) q[0];
sx q[0];
rz(-3.114469) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3168542) q[2];
sx q[2];
rz(-1.6055664) q[2];
sx q[2];
rz(2.8512213) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.3251212) q[1];
sx q[1];
rz(-0.42342788) q[1];
sx q[1];
rz(-0.85514833) q[1];
rz(-0.14988092) q[3];
sx q[3];
rz(-1.1544041) q[3];
sx q[3];
rz(-0.021377953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.985118) q[2];
sx q[2];
rz(-1.6343583) q[2];
sx q[2];
rz(-1.1953243) q[2];
rz(2.7965503) q[3];
sx q[3];
rz(-0.86193591) q[3];
sx q[3];
rz(-2.2752458) q[3];
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
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7095551) q[0];
sx q[0];
rz(-2.9032752) q[0];
sx q[0];
rz(-3.0639783) q[0];
rz(1.2331102) q[1];
sx q[1];
rz(-2.1014919) q[1];
sx q[1];
rz(0.79201039) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4302942) q[0];
sx q[0];
rz(-2.2195507) q[0];
sx q[0];
rz(2.4374967) q[0];
x q[1];
rz(-0.15969212) q[2];
sx q[2];
rz(-1.4174685) q[2];
sx q[2];
rz(-0.13736573) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.4153667) q[1];
sx q[1];
rz(-1.109352) q[1];
sx q[1];
rz(-2.9699294) q[1];
rz(-pi) q[2];
rz(0.80371334) q[3];
sx q[3];
rz(-1.6559542) q[3];
sx q[3];
rz(-2.5309895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.9218257) q[2];
sx q[2];
rz(-0.91138387) q[2];
sx q[2];
rz(-1.0731953) q[2];
rz(-0.24751599) q[3];
sx q[3];
rz(-1.211834) q[3];
sx q[3];
rz(0.016935067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3794004) q[0];
sx q[0];
rz(-1.301855) q[0];
sx q[0];
rz(-0.72574885) q[0];
rz(-2.2629867) q[1];
sx q[1];
rz(-2.286372) q[1];
sx q[1];
rz(1.1927037) q[1];
rz(-1.1996126) q[2];
sx q[2];
rz(-0.51358583) q[2];
sx q[2];
rz(-2.5805342) q[2];
rz(-1.6346651) q[3];
sx q[3];
rz(-2.4136834) q[3];
sx q[3];
rz(-0.53993213) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
