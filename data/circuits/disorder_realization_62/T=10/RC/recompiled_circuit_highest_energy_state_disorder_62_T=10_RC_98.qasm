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
rz(2.565413) q[0];
sx q[0];
rz(-1.6771069) q[0];
sx q[0];
rz(0.65281868) q[0];
rz(-0.44261143) q[1];
sx q[1];
rz(-2.0191329) q[1];
sx q[1];
rz(2.519156) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25678062) q[0];
sx q[0];
rz(-0.33397618) q[0];
sx q[0];
rz(-1.2080753) q[0];
x q[1];
rz(-0.60127778) q[2];
sx q[2];
rz(-1.8493422) q[2];
sx q[2];
rz(2.8097866) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.5099401) q[1];
sx q[1];
rz(-1.3871683) q[1];
sx q[1];
rz(-2.0765523) q[1];
rz(-pi) q[2];
rz(2.9838324) q[3];
sx q[3];
rz(-0.96786849) q[3];
sx q[3];
rz(-1.2186528) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6301253) q[2];
sx q[2];
rz(-1.7897391) q[2];
sx q[2];
rz(-0.54568616) q[2];
rz(2.4557579) q[3];
sx q[3];
rz(-0.67982173) q[3];
sx q[3];
rz(0.95664501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0551374) q[0];
sx q[0];
rz(-0.70239037) q[0];
sx q[0];
rz(1.0954683) q[0];
rz(-2.144004) q[1];
sx q[1];
rz(-1.2350524) q[1];
sx q[1];
rz(1.9445317) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8923949) q[0];
sx q[0];
rz(-2.4955242) q[0];
sx q[0];
rz(-1.1419673) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3648974) q[2];
sx q[2];
rz(-2.0414845) q[2];
sx q[2];
rz(1.350251) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.8833739) q[1];
sx q[1];
rz(-1.1158459) q[1];
sx q[1];
rz(2.8243498) q[1];
x q[2];
rz(-1.7577111) q[3];
sx q[3];
rz(-1.2469562) q[3];
sx q[3];
rz(1.8139386) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.7257488) q[2];
sx q[2];
rz(-1.7472605) q[2];
sx q[2];
rz(-1.7572629) q[2];
rz(1.7978801) q[3];
sx q[3];
rz(-1.5735156) q[3];
sx q[3];
rz(1.869092) q[3];
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
rz(pi/2) q[3];
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
rz(0.62190732) q[0];
sx q[0];
rz(-0.78751957) q[0];
sx q[0];
rz(0.4271048) q[0];
rz(-1.9097795) q[1];
sx q[1];
rz(-1.6424664) q[1];
sx q[1];
rz(-0.61418358) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1290951) q[0];
sx q[0];
rz(-1.3007727) q[0];
sx q[0];
rz(-2.8010023) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9236058) q[2];
sx q[2];
rz(-2.3828016) q[2];
sx q[2];
rz(1.2365149) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.28084785) q[1];
sx q[1];
rz(-0.13361803) q[1];
sx q[1];
rz(-2.0693151) q[1];
rz(-pi) q[2];
rz(-2.6355686) q[3];
sx q[3];
rz(-1.3576686) q[3];
sx q[3];
rz(2.0863078) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.4517639) q[2];
sx q[2];
rz(-0.58388766) q[2];
sx q[2];
rz(2.8705987) q[2];
rz(2.6594035) q[3];
sx q[3];
rz(-0.8420344) q[3];
sx q[3];
rz(-1.8577925) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86376205) q[0];
sx q[0];
rz(-1.8067124) q[0];
sx q[0];
rz(2.7226287) q[0];
rz(-0.22467443) q[1];
sx q[1];
rz(-1.9326262) q[1];
sx q[1];
rz(3.0640501) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52151742) q[0];
sx q[0];
rz(-1.8635611) q[0];
sx q[0];
rz(1.0241246) q[0];
rz(-0.77247844) q[2];
sx q[2];
rz(-1.315552) q[2];
sx q[2];
rz(-2.8359063) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.3540707) q[1];
sx q[1];
rz(-1.7131299) q[1];
sx q[1];
rz(0.44146367) q[1];
x q[2];
rz(0.61309149) q[3];
sx q[3];
rz(-2.7889851) q[3];
sx q[3];
rz(-2.4168454) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.055723995) q[2];
sx q[2];
rz(-0.33129498) q[2];
sx q[2];
rz(-1.9191939) q[2];
rz(-0.27082768) q[3];
sx q[3];
rz(-1.9041678) q[3];
sx q[3];
rz(2.6868611) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5772575) q[0];
sx q[0];
rz(-2.1463647) q[0];
sx q[0];
rz(1.3539535) q[0];
rz(2.1063781) q[1];
sx q[1];
rz(-1.926492) q[1];
sx q[1];
rz(-1.6114906) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3116118) q[0];
sx q[0];
rz(-0.84048827) q[0];
sx q[0];
rz(2.3118224) q[0];
x q[1];
rz(0.45299977) q[2];
sx q[2];
rz(-2.7794547) q[2];
sx q[2];
rz(2.1299794) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2097047) q[1];
sx q[1];
rz(-1.6242289) q[1];
sx q[1];
rz(-2.7393952) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5132929) q[3];
sx q[3];
rz(-1.981145) q[3];
sx q[3];
rz(-1.3259322) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.80663854) q[2];
sx q[2];
rz(-0.97794509) q[2];
sx q[2];
rz(3.0777001) q[2];
rz(0.70819267) q[3];
sx q[3];
rz(-1.4539366) q[3];
sx q[3];
rz(-1.8074869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42234364) q[0];
sx q[0];
rz(-3.1184734) q[0];
sx q[0];
rz(1.0759906) q[0];
rz(-0.05323449) q[1];
sx q[1];
rz(-2.2884171) q[1];
sx q[1];
rz(0.3784953) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.196825) q[0];
sx q[0];
rz(-1.440597) q[0];
sx q[0];
rz(0.2567592) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0113651) q[2];
sx q[2];
rz(-0.097245596) q[2];
sx q[2];
rz(-2.2878316) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6798181) q[1];
sx q[1];
rz(-1.8289806) q[1];
sx q[1];
rz(-1.3854909) q[1];
rz(-pi) q[2];
rz(-0.58121292) q[3];
sx q[3];
rz(-2.2848115) q[3];
sx q[3];
rz(-2.9670144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.1899015) q[2];
sx q[2];
rz(-1.4611763) q[2];
sx q[2];
rz(-1.0339197) q[2];
rz(2.2377491) q[3];
sx q[3];
rz(-0.90141064) q[3];
sx q[3];
rz(-0.044053642) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0251625) q[0];
sx q[0];
rz(-1.6950386) q[0];
sx q[0];
rz(2.548581) q[0];
rz(-3.016839) q[1];
sx q[1];
rz(-1.9826823) q[1];
sx q[1];
rz(-0.4932901) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.497555) q[0];
sx q[0];
rz(-1.8480267) q[0];
sx q[0];
rz(-1.9009186) q[0];
rz(-pi) q[1];
x q[1];
rz(0.22809314) q[2];
sx q[2];
rz(-1.5465294) q[2];
sx q[2];
rz(-2.9920141) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.4932781) q[1];
sx q[1];
rz(-1.9946792) q[1];
sx q[1];
rz(-2.5878169) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4629355) q[3];
sx q[3];
rz(-1.2965805) q[3];
sx q[3];
rz(-1.3317684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.3844246) q[2];
sx q[2];
rz(-1.1606471) q[2];
sx q[2];
rz(-1.9942795) q[2];
rz(0.82289639) q[3];
sx q[3];
rz(-1.1742679) q[3];
sx q[3];
rz(2.4790922) q[3];
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
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61854521) q[0];
sx q[0];
rz(-1.86919) q[0];
sx q[0];
rz(-0.4185032) q[0];
rz(0.11323994) q[1];
sx q[1];
rz(-1.6721882) q[1];
sx q[1];
rz(1.07771) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1215873) q[0];
sx q[0];
rz(-1.5888056) q[0];
sx q[0];
rz(-1.2822582) q[0];
rz(-1.244721) q[2];
sx q[2];
rz(-2.3387944) q[2];
sx q[2];
rz(1.7434415) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5716963) q[1];
sx q[1];
rz(-1.2087617) q[1];
sx q[1];
rz(0.052459474) q[1];
rz(-3.0642068) q[3];
sx q[3];
rz(-0.65801453) q[3];
sx q[3];
rz(-0.24694996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.8354127) q[2];
sx q[2];
rz(-0.70859185) q[2];
sx q[2];
rz(-1.4554679) q[2];
rz(-0.16737394) q[3];
sx q[3];
rz(-0.51515976) q[3];
sx q[3];
rz(2.0696056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48285943) q[0];
sx q[0];
rz(-1.3128244) q[0];
sx q[0];
rz(0.40153781) q[0];
rz(-2.7032779) q[1];
sx q[1];
rz(-2.3500748) q[1];
sx q[1];
rz(1.4567136) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0204937) q[0];
sx q[0];
rz(-1.5803845) q[0];
sx q[0];
rz(1.5367979) q[0];
rz(2.4859547) q[2];
sx q[2];
rz(-1.3961424) q[2];
sx q[2];
rz(1.8675592) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.7359384) q[1];
sx q[1];
rz(-0.88684139) q[1];
sx q[1];
rz(0.64534646) q[1];
rz(-pi) q[2];
rz(1.0965804) q[3];
sx q[3];
rz(-1.1062396) q[3];
sx q[3];
rz(1.43917) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.7598286) q[2];
sx q[2];
rz(-0.20918736) q[2];
sx q[2];
rz(-0.71167243) q[2];
rz(3.107374) q[3];
sx q[3];
rz(-2.4331369) q[3];
sx q[3];
rz(-1.9862407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
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
rz(1.3592247) q[0];
sx q[0];
rz(-1.490626) q[0];
sx q[0];
rz(-2.3427298) q[0];
rz(1.0460188) q[1];
sx q[1];
rz(-1.9452399) q[1];
sx q[1];
rz(2.9050713) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72457658) q[0];
sx q[0];
rz(-1.5931221) q[0];
sx q[0];
rz(-1.6315559) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3600618) q[2];
sx q[2];
rz(-2.3252914) q[2];
sx q[2];
rz(2.3280509) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.7812209) q[1];
sx q[1];
rz(-2.0016333) q[1];
sx q[1];
rz(-1.9442417) q[1];
rz(-pi) q[2];
rz(0.71098401) q[3];
sx q[3];
rz(-1.4812638) q[3];
sx q[3];
rz(2.3481365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.025042621) q[2];
sx q[2];
rz(-0.5566842) q[2];
sx q[2];
rz(-2.4994948) q[2];
rz(-2.5721278) q[3];
sx q[3];
rz(-0.48782188) q[3];
sx q[3];
rz(3.0552982) q[3];
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
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6713329) q[0];
sx q[0];
rz(-1.8076121) q[0];
sx q[0];
rz(-2.1847834) q[0];
rz(1.1915462) q[1];
sx q[1];
rz(-1.5793431) q[1];
sx q[1];
rz(1.6021077) q[1];
rz(-1.446522) q[2];
sx q[2];
rz(-2.3588603) q[2];
sx q[2];
rz(-1.0609577) q[2];
rz(2.2254734) q[3];
sx q[3];
rz(-2.6226433) q[3];
sx q[3];
rz(-0.15437689) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
