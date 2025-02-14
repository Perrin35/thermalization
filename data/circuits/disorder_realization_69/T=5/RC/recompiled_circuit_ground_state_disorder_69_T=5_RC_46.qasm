OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.2562113) q[0];
sx q[0];
rz(-0.71330944) q[0];
sx q[0];
rz(-1.2137086) q[0];
rz(-1.3656536) q[1];
sx q[1];
rz(-2.8627099) q[1];
sx q[1];
rz(-0.20067781) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4532625) q[0];
sx q[0];
rz(-1.4879259) q[0];
sx q[0];
rz(1.1088662) q[0];
rz(-2.9242582) q[2];
sx q[2];
rz(-1.4121018) q[2];
sx q[2];
rz(-2.8087698) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.37113443) q[1];
sx q[1];
rz(-1.7672897) q[1];
sx q[1];
rz(-1.7540098) q[1];
x q[2];
rz(0.49978017) q[3];
sx q[3];
rz(-2.793514) q[3];
sx q[3];
rz(1.4760555) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.063244907) q[2];
sx q[2];
rz(-1.941961) q[2];
sx q[2];
rz(-0.71301785) q[2];
rz(1.2837422) q[3];
sx q[3];
rz(-2.4132437) q[3];
sx q[3];
rz(-0.89948765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9707608) q[0];
sx q[0];
rz(-1.4749227) q[0];
sx q[0];
rz(1.4738039) q[0];
rz(2.5693192) q[1];
sx q[1];
rz(-2.0794561) q[1];
sx q[1];
rz(1.3776113) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77967465) q[0];
sx q[0];
rz(-0.55680823) q[0];
sx q[0];
rz(-2.1614055) q[0];
rz(-pi) q[1];
rz(-0.69852918) q[2];
sx q[2];
rz(-2.4870092) q[2];
sx q[2];
rz(1.6561001) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.2050537) q[1];
sx q[1];
rz(-1.2779826) q[1];
sx q[1];
rz(1.2252625) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2632971) q[3];
sx q[3];
rz(-1.3756985) q[3];
sx q[3];
rz(1.1565735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.1231692) q[2];
sx q[2];
rz(-1.6320684) q[2];
sx q[2];
rz(1.337576) q[2];
rz(2.7959974) q[3];
sx q[3];
rz(-0.59185043) q[3];
sx q[3];
rz(0.93446294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0551986) q[0];
sx q[0];
rz(-2.9559657) q[0];
sx q[0];
rz(2.5523972) q[0];
rz(2.9335754) q[1];
sx q[1];
rz(-2.5841525) q[1];
sx q[1];
rz(-0.20588188) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.170144) q[0];
sx q[0];
rz(-1.2177694) q[0];
sx q[0];
rz(0.63944177) q[0];
rz(-1.0127064) q[2];
sx q[2];
rz(-0.79540247) q[2];
sx q[2];
rz(-2.263042) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.66135397) q[1];
sx q[1];
rz(-1.8744555) q[1];
sx q[1];
rz(-2.420029) q[1];
rz(-pi) q[2];
x q[2];
rz(3.053756) q[3];
sx q[3];
rz(-1.9760152) q[3];
sx q[3];
rz(-2.3072568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0899293) q[2];
sx q[2];
rz(-1.0445003) q[2];
sx q[2];
rz(2.5970686) q[2];
rz(2.6868467) q[3];
sx q[3];
rz(-0.62362042) q[3];
sx q[3];
rz(-3.1375942) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74314463) q[0];
sx q[0];
rz(-2.5938617) q[0];
sx q[0];
rz(0.69993436) q[0];
rz(1.3088538) q[1];
sx q[1];
rz(-1.0700285) q[1];
sx q[1];
rz(1.7226137) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5154079) q[0];
sx q[0];
rz(-1.1061449) q[0];
sx q[0];
rz(2.9166469) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5350902) q[2];
sx q[2];
rz(-0.49884819) q[2];
sx q[2];
rz(0.6998261) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.0154889) q[1];
sx q[1];
rz(-2.2601193) q[1];
sx q[1];
rz(-1.8828805) q[1];
rz(-2.9955825) q[3];
sx q[3];
rz(-1.8838804) q[3];
sx q[3];
rz(2.2816471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.2122638) q[2];
sx q[2];
rz(-0.45390359) q[2];
sx q[2];
rz(-1.3455428) q[2];
rz(-0.71873194) q[3];
sx q[3];
rz(-0.74145442) q[3];
sx q[3];
rz(2.45347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6253925) q[0];
sx q[0];
rz(-1.9560408) q[0];
sx q[0];
rz(-0.54310435) q[0];
rz(-0.25554666) q[1];
sx q[1];
rz(-0.4823904) q[1];
sx q[1];
rz(1.7593174) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92111787) q[0];
sx q[0];
rz(-2.4555523) q[0];
sx q[0];
rz(-0.6612079) q[0];
rz(0.76680317) q[2];
sx q[2];
rz(-1.4652952) q[2];
sx q[2];
rz(-2.2720384) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.68866036) q[1];
sx q[1];
rz(-2.1218461) q[1];
sx q[1];
rz(2.8351985) q[1];
x q[2];
rz(0.38726728) q[3];
sx q[3];
rz(-0.42902374) q[3];
sx q[3];
rz(-0.22873345) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.9690659) q[2];
sx q[2];
rz(-1.4719937) q[2];
sx q[2];
rz(2.772061) q[2];
rz(1.3299804) q[3];
sx q[3];
rz(-2.3510635) q[3];
sx q[3];
rz(1.3548405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(-2.1177649) q[0];
sx q[0];
rz(-0.33482877) q[0];
sx q[0];
rz(3.0552926) q[0];
rz(-1.6465126) q[1];
sx q[1];
rz(-1.1777271) q[1];
sx q[1];
rz(-0.42974791) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5983026) q[0];
sx q[0];
rz(-1.6409281) q[0];
sx q[0];
rz(1.7599517) q[0];
rz(-pi) q[1];
rz(-2.7280732) q[2];
sx q[2];
rz(-2.4405688) q[2];
sx q[2];
rz(-1.9961693) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.60372231) q[1];
sx q[1];
rz(-0.66927823) q[1];
sx q[1];
rz(-2.3905588) q[1];
x q[2];
rz(1.8514351) q[3];
sx q[3];
rz(-2.2220816) q[3];
sx q[3];
rz(0.7574581) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1090082) q[2];
sx q[2];
rz(-1.0823559) q[2];
sx q[2];
rz(-1.2004131) q[2];
rz(-2.9754908) q[3];
sx q[3];
rz(-0.76986543) q[3];
sx q[3];
rz(-2.2782245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85668856) q[0];
sx q[0];
rz(-2.3785474) q[0];
sx q[0];
rz(-2.3175008) q[0];
rz(-1.2245945) q[1];
sx q[1];
rz(-1.2242182) q[1];
sx q[1];
rz(-1.0138938) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7214425) q[0];
sx q[0];
rz(-0.80810302) q[0];
sx q[0];
rz(2.0049014) q[0];
rz(0.61182558) q[2];
sx q[2];
rz(-2.6227747) q[2];
sx q[2];
rz(-0.045836115) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.9062324) q[1];
sx q[1];
rz(-2.4543833) q[1];
sx q[1];
rz(2.5053124) q[1];
rz(-1.0345633) q[3];
sx q[3];
rz(-1.9120029) q[3];
sx q[3];
rz(3.0337014) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8635204) q[2];
sx q[2];
rz(-2.4188228) q[2];
sx q[2];
rz(1.2449167) q[2];
rz(-0.23166367) q[3];
sx q[3];
rz(-2.1004227) q[3];
sx q[3];
rz(0.1483354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-2.046655) q[0];
sx q[0];
rz(-1.9355087) q[0];
sx q[0];
rz(-0.99223247) q[0];
rz(-0.16920432) q[1];
sx q[1];
rz(-0.84183401) q[1];
sx q[1];
rz(-1.4134891) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1344101) q[0];
sx q[0];
rz(-2.5408486) q[0];
sx q[0];
rz(-3.0092165) q[0];
rz(2.8529819) q[2];
sx q[2];
rz(-1.0373652) q[2];
sx q[2];
rz(2.6030428) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7602049) q[1];
sx q[1];
rz(-1.5389812) q[1];
sx q[1];
rz(0.29812584) q[1];
rz(-1.0905209) q[3];
sx q[3];
rz(-1.3930219) q[3];
sx q[3];
rz(-2.4644574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.32470545) q[2];
sx q[2];
rz(-1.6605261) q[2];
sx q[2];
rz(-1.7186349) q[2];
rz(-3.0286246) q[3];
sx q[3];
rz(-2.8739417) q[3];
sx q[3];
rz(-0.62937361) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1896352) q[0];
sx q[0];
rz(-1.4653787) q[0];
sx q[0];
rz(-2.5919609) q[0];
rz(-2.5426087) q[1];
sx q[1];
rz(-0.87268972) q[1];
sx q[1];
rz(-1.5113066) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6236512) q[0];
sx q[0];
rz(-1.7474801) q[0];
sx q[0];
rz(1.5743295) q[0];
rz(2.6312066) q[2];
sx q[2];
rz(-0.56737075) q[2];
sx q[2];
rz(-0.20287831) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.0907184) q[1];
sx q[1];
rz(-1.1696739) q[1];
sx q[1];
rz(-2.119333) q[1];
rz(-pi) q[2];
rz(2.6069943) q[3];
sx q[3];
rz(-2.5851997) q[3];
sx q[3];
rz(2.8884187) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.0281684) q[2];
sx q[2];
rz(-1.3415965) q[2];
sx q[2];
rz(-0.33451954) q[2];
rz(2.4962375) q[3];
sx q[3];
rz(-1.0048451) q[3];
sx q[3];
rz(0.2373124) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-0.065011218) q[0];
sx q[0];
rz(-0.33877057) q[0];
sx q[0];
rz(0.6231128) q[0];
rz(0.30162853) q[1];
sx q[1];
rz(-2.5239065) q[1];
sx q[1];
rz(-1.041144) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.96466) q[0];
sx q[0];
rz(-2.3752691) q[0];
sx q[0];
rz(-0.64789741) q[0];
rz(-pi) q[1];
rz(2.6807296) q[2];
sx q[2];
rz(-1.653737) q[2];
sx q[2];
rz(-1.2316224) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.99771254) q[1];
sx q[1];
rz(-1.6015983) q[1];
sx q[1];
rz(-2.374359) q[1];
rz(2.0984394) q[3];
sx q[3];
rz(-2.7343547) q[3];
sx q[3];
rz(1.094363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.4447896) q[2];
sx q[2];
rz(-2.5371234) q[2];
sx q[2];
rz(0.94318548) q[2];
rz(1.7275564) q[3];
sx q[3];
rz(-2.2351738) q[3];
sx q[3];
rz(-2.3448155) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74054756) q[0];
sx q[0];
rz(-2.7475806) q[0];
sx q[0];
rz(-0.077234118) q[0];
rz(0.41863353) q[1];
sx q[1];
rz(-1.5822462) q[1];
sx q[1];
rz(0.15204631) q[1];
rz(1.6054016) q[2];
sx q[2];
rz(-2.1454951) q[2];
sx q[2];
rz(-2.9401671) q[2];
rz(2.5984305) q[3];
sx q[3];
rz(-0.49351963) q[3];
sx q[3];
rz(-1.7388572) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
