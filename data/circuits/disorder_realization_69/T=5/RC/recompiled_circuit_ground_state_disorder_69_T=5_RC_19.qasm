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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6883301) q[0];
sx q[0];
rz(-1.6536667) q[0];
sx q[0];
rz(-2.0327264) q[0];
x q[1];
rz(0.21733445) q[2];
sx q[2];
rz(-1.7294908) q[2];
sx q[2];
rz(-0.33282285) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7704582) q[1];
sx q[1];
rz(-1.7672897) q[1];
sx q[1];
rz(-1.7540098) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8332769) q[3];
sx q[3];
rz(-1.7349958) q[3];
sx q[3];
rz(-2.7621244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.063244907) q[2];
sx q[2];
rz(-1.941961) q[2];
sx q[2];
rz(0.71301785) q[2];
rz(1.2837422) q[3];
sx q[3];
rz(-0.72834891) q[3];
sx q[3];
rz(-2.242105) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17083183) q[0];
sx q[0];
rz(-1.66667) q[0];
sx q[0];
rz(-1.6677888) q[0];
rz(2.5693192) q[1];
sx q[1];
rz(-2.0794561) q[1];
sx q[1];
rz(1.3776113) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77967465) q[0];
sx q[0];
rz(-0.55680823) q[0];
sx q[0];
rz(0.98018719) q[0];
rz(-pi) q[1];
x q[1];
rz(0.53133659) q[2];
sx q[2];
rz(-1.9730933) q[2];
sx q[2];
rz(0.50237331) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.26220825) q[1];
sx q[1];
rz(-1.9010494) q[1];
sx q[1];
rz(-0.31008215) q[1];
rz(-pi) q[2];
x q[2];
rz(0.20444025) q[3];
sx q[3];
rz(-1.8722765) q[3];
sx q[3];
rz(-2.7888576) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.01842347) q[2];
sx q[2];
rz(-1.6320684) q[2];
sx q[2];
rz(-1.8040166) q[2];
rz(0.34559524) q[3];
sx q[3];
rz(-2.5497422) q[3];
sx q[3];
rz(0.93446294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0551986) q[0];
sx q[0];
rz(-2.9559657) q[0];
sx q[0];
rz(2.5523972) q[0];
rz(-2.9335754) q[1];
sx q[1];
rz(-2.5841525) q[1];
sx q[1];
rz(-2.9357108) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1070646) q[0];
sx q[0];
rz(-2.42332) q[0];
sx q[0];
rz(0.55315259) q[0];
rz(-pi) q[1];
rz(0.85742204) q[2];
sx q[2];
rz(-1.1829585) q[2];
sx q[2];
rz(-2.037354) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2370473) q[1];
sx q[1];
rz(-0.77213192) q[1];
sx q[1];
rz(-0.44293483) q[1];
x q[2];
rz(1.1641762) q[3];
sx q[3];
rz(-1.6515035) q[3];
sx q[3];
rz(0.77116283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.0899293) q[2];
sx q[2];
rz(-2.0970924) q[2];
sx q[2];
rz(-2.5970686) q[2];
rz(-2.6868467) q[3];
sx q[3];
rz(-0.62362042) q[3];
sx q[3];
rz(-0.0039984306) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74314463) q[0];
sx q[0];
rz(-2.5938617) q[0];
sx q[0];
rz(-0.69993436) q[0];
rz(1.8327389) q[1];
sx q[1];
rz(-2.0715641) q[1];
sx q[1];
rz(-1.4189789) q[1];
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
rz(-1.6065025) q[2];
sx q[2];
rz(-2.6427445) q[2];
sx q[2];
rz(2.4417666) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.1261038) q[1];
sx q[1];
rz(-2.2601193) q[1];
sx q[1];
rz(-1.2587121) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1484234) q[3];
sx q[3];
rz(-2.7971533) q[3];
sx q[3];
rz(0.41448739) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2122638) q[2];
sx q[2];
rz(-2.6876891) q[2];
sx q[2];
rz(1.7960499) q[2];
rz(-2.4228607) q[3];
sx q[3];
rz(-2.4001382) q[3];
sx q[3];
rz(2.45347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5162002) q[0];
sx q[0];
rz(-1.9560408) q[0];
sx q[0];
rz(0.54310435) q[0];
rz(-2.886046) q[1];
sx q[1];
rz(-0.4823904) q[1];
sx q[1];
rz(1.3822752) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92111787) q[0];
sx q[0];
rz(-0.68604031) q[0];
sx q[0];
rz(-2.4803848) q[0];
rz(-pi) q[1];
rz(0.76680317) q[2];
sx q[2];
rz(-1.6762974) q[2];
sx q[2];
rz(2.2720384) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.4529323) q[1];
sx q[1];
rz(-1.0197465) q[1];
sx q[1];
rz(0.30639415) q[1];
rz(-0.38726728) q[3];
sx q[3];
rz(-2.7125689) q[3];
sx q[3];
rz(2.9128592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.9690659) q[2];
sx q[2];
rz(-1.4719937) q[2];
sx q[2];
rz(-0.36953163) q[2];
rz(1.3299804) q[3];
sx q[3];
rz(-2.3510635) q[3];
sx q[3];
rz(-1.7867521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1177649) q[0];
sx q[0];
rz(-2.8067639) q[0];
sx q[0];
rz(-3.0552926) q[0];
rz(-1.49508) q[1];
sx q[1];
rz(-1.1777271) q[1];
sx q[1];
rz(-2.7118447) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0140918) q[0];
sx q[0];
rz(-1.7594811) q[0];
sx q[0];
rz(-0.071401113) q[0];
x q[1];
rz(-2.4835973) q[2];
sx q[2];
rz(-1.832973) q[2];
sx q[2];
rz(-0.10181759) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.60372231) q[1];
sx q[1];
rz(-2.4723144) q[1];
sx q[1];
rz(-0.75103384) q[1];
rz(-pi) q[2];
rz(2.7930611) q[3];
sx q[3];
rz(-0.70100871) q[3];
sx q[3];
rz(1.9402869) q[3];
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
rz(2.9754908) q[3];
sx q[3];
rz(-2.3717272) q[3];
sx q[3];
rz(0.86336819) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2849041) q[0];
sx q[0];
rz(-0.76304522) q[0];
sx q[0];
rz(2.3175008) q[0];
rz(1.9169982) q[1];
sx q[1];
rz(-1.2242182) q[1];
sx q[1];
rz(2.1276988) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82912927) q[0];
sx q[0];
rz(-2.2861963) q[0];
sx q[0];
rz(2.7269643) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2538919) q[2];
sx q[2];
rz(-1.9887668) q[2];
sx q[2];
rz(0.63374139) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.287266) q[1];
sx q[1];
rz(-1.1842898) q[1];
sx q[1];
rz(2.5581761) q[1];
x q[2];
rz(2.7498662) q[3];
sx q[3];
rz(-2.0731032) q[3];
sx q[3];
rz(1.6592178) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.27807221) q[2];
sx q[2];
rz(-0.72276989) q[2];
sx q[2];
rz(1.2449167) q[2];
rz(2.909929) q[3];
sx q[3];
rz(-1.04117) q[3];
sx q[3];
rz(-0.1483354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0949377) q[0];
sx q[0];
rz(-1.2060839) q[0];
sx q[0];
rz(0.99223247) q[0];
rz(-0.16920432) q[1];
sx q[1];
rz(-0.84183401) q[1];
sx q[1];
rz(1.7281035) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0071825562) q[0];
sx q[0];
rz(-2.5408486) q[0];
sx q[0];
rz(3.0092165) q[0];
rz(2.1229073) q[2];
sx q[2];
rz(-1.8183961) q[2];
sx q[2];
rz(1.9595065) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.0553356) q[1];
sx q[1];
rz(-2.8418243) q[1];
sx q[1];
rz(0.10792984) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0510717) q[3];
sx q[3];
rz(-1.3930219) q[3];
sx q[3];
rz(2.4644574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.32470545) q[2];
sx q[2];
rz(-1.4810666) q[2];
sx q[2];
rz(-1.7186349) q[2];
rz(-0.11296806) q[3];
sx q[3];
rz(-0.26765099) q[3];
sx q[3];
rz(-0.62937361) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9519575) q[0];
sx q[0];
rz(-1.4653787) q[0];
sx q[0];
rz(-0.54963175) q[0];
rz(0.59898392) q[1];
sx q[1];
rz(-2.2689029) q[1];
sx q[1];
rz(1.5113066) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.64375) q[0];
sx q[0];
rz(-0.17671876) q[0];
sx q[0];
rz(0.019785893) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8725996) q[2];
sx q[2];
rz(-1.0827218) q[2];
sx q[2];
rz(0.38319626) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.0907184) q[1];
sx q[1];
rz(-1.9719187) q[1];
sx q[1];
rz(2.119333) q[1];
rz(-2.6069943) q[3];
sx q[3];
rz(-2.5851997) q[3];
sx q[3];
rz(0.25317395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.0281684) q[2];
sx q[2];
rz(-1.3415965) q[2];
sx q[2];
rz(0.33451954) q[2];
rz(-2.4962375) q[3];
sx q[3];
rz(-1.0048451) q[3];
sx q[3];
rz(2.9042802) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.065011218) q[0];
sx q[0];
rz(-2.8028221) q[0];
sx q[0];
rz(-0.6231128) q[0];
rz(-2.8399641) q[1];
sx q[1];
rz(-0.61768618) q[1];
sx q[1];
rz(-2.1004486) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2469676) q[0];
sx q[0];
rz(-1.1389705) q[0];
sx q[0];
rz(0.65460848) q[0];
rz(2.6807296) q[2];
sx q[2];
rz(-1.653737) q[2];
sx q[2];
rz(1.9099703) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.5981982) q[1];
sx q[1];
rz(-0.80401995) q[1];
sx q[1];
rz(1.5280185) q[1];
rz(-pi) q[2];
rz(-1.9275366) q[3];
sx q[3];
rz(-1.7715653) q[3];
sx q[3];
rz(-0.014895766) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.69680301) q[2];
sx q[2];
rz(-2.5371234) q[2];
sx q[2];
rz(2.1984072) q[2];
rz(1.4140363) q[3];
sx q[3];
rz(-2.2351738) q[3];
sx q[3];
rz(2.3448155) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-0.05337333) q[2];
sx q[2];
rz(-2.5659701) q[2];
sx q[2];
rz(0.26502668) q[2];
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
