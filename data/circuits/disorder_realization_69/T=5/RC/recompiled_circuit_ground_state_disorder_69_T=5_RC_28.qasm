OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.8853814) q[0];
sx q[0];
rz(-2.4282832) q[0];
sx q[0];
rz(-1.927884) q[0];
rz(-1.3656536) q[1];
sx q[1];
rz(-2.8627099) q[1];
sx q[1];
rz(-0.20067781) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4532625) q[0];
sx q[0];
rz(-1.4879259) q[0];
sx q[0];
rz(2.0327264) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4083449) q[2];
sx q[2];
rz(-1.3562358) q[2];
sx q[2];
rz(-1.868737) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.37113443) q[1];
sx q[1];
rz(-1.7672897) q[1];
sx q[1];
rz(1.3875828) q[1];
x q[2];
rz(-0.49978017) q[3];
sx q[3];
rz(-2.793514) q[3];
sx q[3];
rz(1.6655371) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.063244907) q[2];
sx q[2];
rz(-1.1996317) q[2];
sx q[2];
rz(0.71301785) q[2];
rz(-1.2837422) q[3];
sx q[3];
rz(-2.4132437) q[3];
sx q[3];
rz(0.89948765) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9707608) q[0];
sx q[0];
rz(-1.66667) q[0];
sx q[0];
rz(-1.6677888) q[0];
rz(-2.5693192) q[1];
sx q[1];
rz(-1.0621366) q[1];
sx q[1];
rz(1.3776113) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.361918) q[0];
sx q[0];
rz(-2.5847844) q[0];
sx q[0];
rz(2.1614055) q[0];
x q[1];
rz(-2.4430635) q[2];
sx q[2];
rz(-2.4870092) q[2];
sx q[2];
rz(-1.6561001) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.0999801) q[1];
sx q[1];
rz(-2.6925107) q[1];
sx q[1];
rz(-0.84347165) q[1];
rz(-pi) q[2];
x q[2];
rz(0.99239852) q[3];
sx q[3];
rz(-2.7790894) q[3];
sx q[3];
rz(0.9622919) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.1231692) q[2];
sx q[2];
rz(-1.5095242) q[2];
sx q[2];
rz(-1.337576) q[2];
rz(-2.7959974) q[3];
sx q[3];
rz(-0.59185043) q[3];
sx q[3];
rz(2.2071297) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.086394101) q[0];
sx q[0];
rz(-0.18562695) q[0];
sx q[0];
rz(-2.5523972) q[0];
rz(-2.9335754) q[1];
sx q[1];
rz(-2.5841525) q[1];
sx q[1];
rz(0.20588188) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.170144) q[0];
sx q[0];
rz(-1.2177694) q[0];
sx q[0];
rz(0.63944177) q[0];
rz(-pi) q[1];
rz(-0.85742204) q[2];
sx q[2];
rz(-1.1829585) q[2];
sx q[2];
rz(-1.1042386) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.2370473) q[1];
sx q[1];
rz(-0.77213192) q[1];
sx q[1];
rz(0.44293483) q[1];
rz(1.9774164) q[3];
sx q[3];
rz(-1.6515035) q[3];
sx q[3];
rz(-0.77116283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.0516633) q[2];
sx q[2];
rz(-1.0445003) q[2];
sx q[2];
rz(0.54452407) q[2];
rz(-0.45474592) q[3];
sx q[3];
rz(-0.62362042) q[3];
sx q[3];
rz(-3.1375942) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.398448) q[0];
sx q[0];
rz(-2.5938617) q[0];
sx q[0];
rz(-2.4416583) q[0];
rz(-1.8327389) q[1];
sx q[1];
rz(-2.0715641) q[1];
sx q[1];
rz(-1.7226137) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5154079) q[0];
sx q[0];
rz(-1.1061449) q[0];
sx q[0];
rz(-2.9166469) q[0];
x q[1];
rz(-3.1221462) q[2];
sx q[2];
rz(-2.0692973) q[2];
sx q[2];
rz(-0.74048238) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.8992865) q[1];
sx q[1];
rz(-1.3315836) q[1];
sx q[1];
rz(-0.71372791) q[1];
rz(1.254564) q[3];
sx q[3];
rz(-1.4319311) q[3];
sx q[3];
rz(-0.75611243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.2122638) q[2];
sx q[2];
rz(-2.6876891) q[2];
sx q[2];
rz(-1.7960499) q[2];
rz(-2.4228607) q[3];
sx q[3];
rz(-0.74145442) q[3];
sx q[3];
rz(-2.45347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5162002) q[0];
sx q[0];
rz(-1.9560408) q[0];
sx q[0];
rz(-0.54310435) q[0];
rz(2.886046) q[1];
sx q[1];
rz(-0.4823904) q[1];
sx q[1];
rz(1.7593174) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7092752) q[0];
sx q[0];
rz(-1.0472282) q[0];
sx q[0];
rz(-2.0366336) q[0];
x q[1];
rz(-2.3747895) q[2];
sx q[2];
rz(-1.4652952) q[2];
sx q[2];
rz(0.86955429) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.2321736) q[1];
sx q[1];
rz(-2.5188753) q[1];
sx q[1];
rz(-2.0270585) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3997283) q[3];
sx q[3];
rz(-1.9661964) q[3];
sx q[3];
rz(0.65034894) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.9690659) q[2];
sx q[2];
rz(-1.669599) q[2];
sx q[2];
rz(-0.36953163) q[2];
rz(1.8116123) q[3];
sx q[3];
rz(-0.79052916) q[3];
sx q[3];
rz(-1.7867521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1177649) q[0];
sx q[0];
rz(-0.33482877) q[0];
sx q[0];
rz(-0.086300015) q[0];
rz(-1.6465126) q[1];
sx q[1];
rz(-1.1777271) q[1];
sx q[1];
rz(-0.42974791) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5983026) q[0];
sx q[0];
rz(-1.5006646) q[0];
sx q[0];
rz(1.381641) q[0];
rz(-pi) q[1];
rz(-0.65799539) q[2];
sx q[2];
rz(-1.3086196) q[2];
sx q[2];
rz(-0.10181759) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.4757929) q[1];
sx q[1];
rz(-1.1000888) q[1];
sx q[1];
rz(1.0757955) q[1];
rz(-2.4709701) q[3];
sx q[3];
rz(-1.792893) q[3];
sx q[3];
rz(-2.5012453) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.1090082) q[2];
sx q[2];
rz(-1.0823559) q[2];
sx q[2];
rz(-1.2004131) q[2];
rz(-0.16610185) q[3];
sx q[3];
rz(-2.3717272) q[3];
sx q[3];
rz(-2.2782245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2849041) q[0];
sx q[0];
rz(-2.3785474) q[0];
sx q[0];
rz(-0.82409182) q[0];
rz(-1.9169982) q[1];
sx q[1];
rz(-1.9173744) q[1];
sx q[1];
rz(-1.0138938) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4201502) q[0];
sx q[0];
rz(-0.80810302) q[0];
sx q[0];
rz(-2.0049014) q[0];
rz(-0.61182558) q[2];
sx q[2];
rz(-2.6227747) q[2];
sx q[2];
rz(0.045836115) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.8543267) q[1];
sx q[1];
rz(-1.1842898) q[1];
sx q[1];
rz(-0.58341656) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1070293) q[3];
sx q[3];
rz(-1.2295897) q[3];
sx q[3];
rz(-3.0337014) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.27807221) q[2];
sx q[2];
rz(-0.72276989) q[2];
sx q[2];
rz(1.2449167) q[2];
rz(-2.909929) q[3];
sx q[3];
rz(-1.04117) q[3];
sx q[3];
rz(-2.9932573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0949377) q[0];
sx q[0];
rz(-1.9355087) q[0];
sx q[0];
rz(2.1493602) q[0];
rz(-0.16920432) q[1];
sx q[1];
rz(-0.84183401) q[1];
sx q[1];
rz(1.7281035) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.468576) q[0];
sx q[0];
rz(-1.4961188) q[0];
sx q[0];
rz(-0.59665307) q[0];
rz(1.0186853) q[2];
sx q[2];
rz(-1.8183961) q[2];
sx q[2];
rz(1.1820861) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.38138776) q[1];
sx q[1];
rz(-1.6026114) q[1];
sx q[1];
rz(-0.29812584) q[1];
rz(-pi) q[2];
rz(1.0905209) q[3];
sx q[3];
rz(-1.7485707) q[3];
sx q[3];
rz(-2.4644574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.8168872) q[2];
sx q[2];
rz(-1.6605261) q[2];
sx q[2];
rz(1.4229577) q[2];
rz(-3.0286246) q[3];
sx q[3];
rz(-2.8739417) q[3];
sx q[3];
rz(2.512219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1896352) q[0];
sx q[0];
rz(-1.676214) q[0];
sx q[0];
rz(0.54963175) q[0];
rz(-0.59898392) q[1];
sx q[1];
rz(-0.87268972) q[1];
sx q[1];
rz(-1.6302861) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0534759) q[0];
sx q[0];
rz(-1.5742745) q[0];
sx q[0];
rz(-2.9649078) q[0];
x q[1];
rz(-0.51038607) q[2];
sx q[2];
rz(-2.5742219) q[2];
sx q[2];
rz(-2.9387143) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.0508742) q[1];
sx q[1];
rz(-1.1696739) q[1];
sx q[1];
rz(-1.0222597) q[1];
rz(-pi) q[2];
rz(-1.263932) q[3];
sx q[3];
rz(-1.0990541) q[3];
sx q[3];
rz(0.8620756) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1134243) q[2];
sx q[2];
rz(-1.3415965) q[2];
sx q[2];
rz(2.8070731) q[2];
rz(-2.4962375) q[3];
sx q[3];
rz(-1.0048451) q[3];
sx q[3];
rz(2.9042802) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.065011218) q[0];
sx q[0];
rz(-0.33877057) q[0];
sx q[0];
rz(0.6231128) q[0];
rz(0.30162853) q[1];
sx q[1];
rz(-2.5239065) q[1];
sx q[1];
rz(2.1004486) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.96466) q[0];
sx q[0];
rz(-2.3752691) q[0];
sx q[0];
rz(2.4936952) q[0];
rz(-pi) q[1];
rz(1.4782466) q[2];
sx q[2];
rz(-2.0299533) q[2];
sx q[2];
rz(2.8435304) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.99771254) q[1];
sx q[1];
rz(-1.6015983) q[1];
sx q[1];
rz(-0.7672337) q[1];
x q[2];
rz(-1.214056) q[3];
sx q[3];
rz(-1.3700273) q[3];
sx q[3];
rz(3.1266969) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4447896) q[2];
sx q[2];
rz(-0.6044693) q[2];
sx q[2];
rz(-2.1984072) q[2];
rz(-1.7275564) q[3];
sx q[3];
rz(-0.90641886) q[3];
sx q[3];
rz(-2.3448155) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4010451) q[0];
sx q[0];
rz(-2.7475806) q[0];
sx q[0];
rz(-0.077234118) q[0];
rz(-0.41863353) q[1];
sx q[1];
rz(-1.5593465) q[1];
sx q[1];
rz(-2.9895463) q[1];
rz(-3.0882193) q[2];
sx q[2];
rz(-0.57562258) q[2];
sx q[2];
rz(-2.876566) q[2];
rz(-2.7100415) q[3];
sx q[3];
rz(-1.8181556) q[3];
sx q[3];
rz(-2.8209742) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
