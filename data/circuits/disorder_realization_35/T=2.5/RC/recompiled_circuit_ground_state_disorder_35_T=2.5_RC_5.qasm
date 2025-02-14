OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.2090981) q[0];
sx q[0];
rz(-3.0734835) q[0];
sx q[0];
rz(2.3647302) q[0];
rz(-1.3070973) q[1];
sx q[1];
rz(-0.96944648) q[1];
sx q[1];
rz(1.3219272) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37047077) q[0];
sx q[0];
rz(-1.7092472) q[0];
sx q[0];
rz(-0.8531424) q[0];
rz(-pi) q[1];
rz(1.2250617) q[2];
sx q[2];
rz(-1.099473) q[2];
sx q[2];
rz(-2.2668348) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.894134) q[1];
sx q[1];
rz(-1.1620635) q[1];
sx q[1];
rz(-0.031886727) q[1];
rz(-pi) q[2];
rz(1.8613937) q[3];
sx q[3];
rz(-1.7589671) q[3];
sx q[3];
rz(-0.0083323697) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.063735828) q[2];
sx q[2];
rz(-1.3602164) q[2];
sx q[2];
rz(0.35749164) q[2];
rz(0.56719559) q[3];
sx q[3];
rz(-1.4128069) q[3];
sx q[3];
rz(2.7020057) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95328632) q[0];
sx q[0];
rz(-2.8585241) q[0];
sx q[0];
rz(1.333492) q[0];
rz(-0.4370583) q[1];
sx q[1];
rz(-0.60085618) q[1];
sx q[1];
rz(0.83998799) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49522049) q[0];
sx q[0];
rz(-0.53421181) q[0];
sx q[0];
rz(1.87508) q[0];
rz(-1.399038) q[2];
sx q[2];
rz(-2.1585625) q[2];
sx q[2];
rz(-1.5422271) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.31038302) q[1];
sx q[1];
rz(-1.0070756) q[1];
sx q[1];
rz(0.68117627) q[1];
rz(3.1026353) q[3];
sx q[3];
rz(-2.1202728) q[3];
sx q[3];
rz(0.4645068) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.4240894) q[2];
sx q[2];
rz(-1.6987897) q[2];
sx q[2];
rz(1.4060414) q[2];
rz(-2.9794335) q[3];
sx q[3];
rz(-1.5806961) q[3];
sx q[3];
rz(-2.6774008) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9185987) q[0];
sx q[0];
rz(-0.080568947) q[0];
sx q[0];
rz(-2.719847) q[0];
rz(-2.2987507) q[1];
sx q[1];
rz(-1.3858567) q[1];
sx q[1];
rz(2.8657894) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1037791) q[0];
sx q[0];
rz(-0.40291726) q[0];
sx q[0];
rz(0.62312868) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6709072) q[2];
sx q[2];
rz(-1.7690725) q[2];
sx q[2];
rz(2.7338365) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.9107863) q[1];
sx q[1];
rz(-1.2503997) q[1];
sx q[1];
rz(-1.5159392) q[1];
rz(-pi) q[2];
x q[2];
rz(0.86902501) q[3];
sx q[3];
rz(-1.6459673) q[3];
sx q[3];
rz(-1.6628671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.80321035) q[2];
sx q[2];
rz(-2.3440177) q[2];
sx q[2];
rz(0.69022834) q[2];
rz(-2.7817182) q[3];
sx q[3];
rz(-1.3709603) q[3];
sx q[3];
rz(-2.7580822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2279219) q[0];
sx q[0];
rz(-2.0991195) q[0];
sx q[0];
rz(-0.87164718) q[0];
rz(2.7323885) q[1];
sx q[1];
rz(-2.6458461) q[1];
sx q[1];
rz(0.31235487) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.384639) q[0];
sx q[0];
rz(-1.1273545) q[0];
sx q[0];
rz(1.4405946) q[0];
rz(-1.6705475) q[2];
sx q[2];
rz(-1.776223) q[2];
sx q[2];
rz(-2.2769986) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.67955454) q[1];
sx q[1];
rz(-1.7741388) q[1];
sx q[1];
rz(2.773079) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3757703) q[3];
sx q[3];
rz(-0.93149501) q[3];
sx q[3];
rz(0.89278683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.13545869) q[2];
sx q[2];
rz(-1.4486855) q[2];
sx q[2];
rz(0.82497605) q[2];
rz(-1.7168335) q[3];
sx q[3];
rz(-2.8202839) q[3];
sx q[3];
rz(0.43535522) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.41998) q[0];
sx q[0];
rz(-1.5897607) q[0];
sx q[0];
rz(0.87442526) q[0];
rz(-0.32886109) q[1];
sx q[1];
rz(-1.3855653) q[1];
sx q[1];
rz(1.0985451) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1060396) q[0];
sx q[0];
rz(-1.4604202) q[0];
sx q[0];
rz(-2.7644185) q[0];
rz(-0.036089049) q[2];
sx q[2];
rz(-2.3054625) q[2];
sx q[2];
rz(-0.058908894) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.23197996) q[1];
sx q[1];
rz(-1.8706338) q[1];
sx q[1];
rz(1.058387) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9727398) q[3];
sx q[3];
rz(-0.62905772) q[3];
sx q[3];
rz(-0.22543104) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.141779) q[2];
sx q[2];
rz(-1.1835316) q[2];
sx q[2];
rz(-2.4675274) q[2];
rz(-1.4874124) q[3];
sx q[3];
rz(-1.6964648) q[3];
sx q[3];
rz(-2.8580247) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11033002) q[0];
sx q[0];
rz(-2.1188348) q[0];
sx q[0];
rz(-1.3856101) q[0];
rz(1.1194057) q[1];
sx q[1];
rz(-1.6740572) q[1];
sx q[1];
rz(-1.3853692) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.114217) q[0];
sx q[0];
rz(-1.3189335) q[0];
sx q[0];
rz(-0.92433874) q[0];
x q[1];
rz(-1.4848861) q[2];
sx q[2];
rz(-1.0738147) q[2];
sx q[2];
rz(-0.53705207) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.0138356) q[1];
sx q[1];
rz(-2.7633939) q[1];
sx q[1];
rz(2.5950123) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6072261) q[3];
sx q[3];
rz(-1.7251937) q[3];
sx q[3];
rz(1.2702219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.063227) q[2];
sx q[2];
rz(-0.93457064) q[2];
sx q[2];
rz(2.7654977) q[2];
rz(-1.1345351) q[3];
sx q[3];
rz(-0.36342707) q[3];
sx q[3];
rz(-0.80662066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.019808708) q[0];
sx q[0];
rz(-0.60096318) q[0];
sx q[0];
rz(-3.0390749) q[0];
rz(-0.58492297) q[1];
sx q[1];
rz(-2.1939317) q[1];
sx q[1];
rz(1.9580511) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2823806) q[0];
sx q[0];
rz(-1.4420274) q[0];
sx q[0];
rz(-0.14315258) q[0];
rz(1.905422) q[2];
sx q[2];
rz(-0.93572223) q[2];
sx q[2];
rz(-2.2872567) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.592748) q[1];
sx q[1];
rz(-1.3164489) q[1];
sx q[1];
rz(2.2266822) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3041586) q[3];
sx q[3];
rz(-2.023306) q[3];
sx q[3];
rz(-0.072991144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.79503235) q[2];
sx q[2];
rz(-2.2117895) q[2];
sx q[2];
rz(2.9417876) q[2];
rz(-1.0605158) q[3];
sx q[3];
rz(-2.6001055) q[3];
sx q[3];
rz(0.04960355) q[3];
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
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26315966) q[0];
sx q[0];
rz(-1.6596154) q[0];
sx q[0];
rz(-0.31295452) q[0];
rz(-2.1921659) q[1];
sx q[1];
rz(-1.7890309) q[1];
sx q[1];
rz(1.688028) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67927406) q[0];
sx q[0];
rz(-1.7574213) q[0];
sx q[0];
rz(1.4786722) q[0];
x q[1];
rz(0.11067783) q[2];
sx q[2];
rz(-2.5191436) q[2];
sx q[2];
rz(-1.770293) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.18583567) q[1];
sx q[1];
rz(-1.794853) q[1];
sx q[1];
rz(0.70835374) q[1];
rz(-pi) q[2];
rz(-1.9096776) q[3];
sx q[3];
rz(-1.1383071) q[3];
sx q[3];
rz(0.56604715) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.35385418) q[2];
sx q[2];
rz(-2.165803) q[2];
sx q[2];
rz(-1.3346416) q[2];
rz(1.3245964) q[3];
sx q[3];
rz(-0.66671222) q[3];
sx q[3];
rz(1.9273531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.036309328) q[0];
sx q[0];
rz(-0.61573017) q[0];
sx q[0];
rz(-0.69806725) q[0];
rz(2.2659194) q[1];
sx q[1];
rz(-2.0798637) q[1];
sx q[1];
rz(-0.36044136) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1669641) q[0];
sx q[0];
rz(-2.1910843) q[0];
sx q[0];
rz(-1.0247158) q[0];
rz(-pi) q[1];
rz(-1.6351885) q[2];
sx q[2];
rz(-1.1887822) q[2];
sx q[2];
rz(-2.7329993) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.97353444) q[1];
sx q[1];
rz(-0.72859166) q[1];
sx q[1];
rz(-0.97961564) q[1];
rz(-pi) q[2];
rz(0.27020755) q[3];
sx q[3];
rz(-2.5761309) q[3];
sx q[3];
rz(-0.61494857) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2299819) q[2];
sx q[2];
rz(-1.4260099) q[2];
sx q[2];
rz(-0.81725517) q[2];
rz(2.3115555) q[3];
sx q[3];
rz(-0.54098141) q[3];
sx q[3];
rz(-2.9221007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63477409) q[0];
sx q[0];
rz(-0.84832484) q[0];
sx q[0];
rz(0.032055227) q[0];
rz(1.8219148) q[1];
sx q[1];
rz(-1.3751043) q[1];
sx q[1];
rz(-2.4874036) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3921627) q[0];
sx q[0];
rz(-1.6902958) q[0];
sx q[0];
rz(0.197535) q[0];
rz(-1.1172574) q[2];
sx q[2];
rz(-0.83570601) q[2];
sx q[2];
rz(-2.105956) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.54683751) q[1];
sx q[1];
rz(-2.0528626) q[1];
sx q[1];
rz(-2.099311) q[1];
rz(1.9496253) q[3];
sx q[3];
rz(-0.66756581) q[3];
sx q[3];
rz(-2.6258056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.0111982) q[2];
sx q[2];
rz(-2.0989213) q[2];
sx q[2];
rz(1.0065669) q[2];
rz(2.448163) q[3];
sx q[3];
rz(-1.3207685) q[3];
sx q[3];
rz(1.8600195) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
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
rz(1.4257767) q[0];
sx q[0];
rz(-2.4632813) q[0];
sx q[0];
rz(3.0671469) q[0];
rz(-0.88059942) q[1];
sx q[1];
rz(-1.1136628) q[1];
sx q[1];
rz(-2.640092) q[1];
rz(-2.1660317) q[2];
sx q[2];
rz(-0.79339334) q[2];
sx q[2];
rz(-1.7421772) q[2];
rz(-1.952259) q[3];
sx q[3];
rz(-0.91100024) q[3];
sx q[3];
rz(-0.80445214) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
