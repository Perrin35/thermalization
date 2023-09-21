OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.3611276) q[0];
sx q[0];
rz(-0.68929231) q[0];
sx q[0];
rz(2.8110992) q[0];
rz(0.32980546) q[1];
sx q[1];
rz(5.4332241) q[1];
sx q[1];
rz(10.133893) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56624428) q[0];
sx q[0];
rz(-2.4465804) q[0];
sx q[0];
rz(-2.038875) q[0];
rz(-pi) q[1];
x q[1];
rz(0.032151392) q[2];
sx q[2];
rz(-1.2574147) q[2];
sx q[2];
rz(1.1653479) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.82422772) q[1];
sx q[1];
rz(-1.0365651) q[1];
sx q[1];
rz(-0.51170106) q[1];
rz(-0.37681864) q[3];
sx q[3];
rz(-2.0599277) q[3];
sx q[3];
rz(-0.34987846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.4101397) q[2];
sx q[2];
rz(-0.4099161) q[2];
sx q[2];
rz(1.5343792) q[2];
rz(0.93506995) q[3];
sx q[3];
rz(-1.3053223) q[3];
sx q[3];
rz(0.7888166) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7222897) q[0];
sx q[0];
rz(-3.0794444) q[0];
sx q[0];
rz(2.5193135) q[0];
rz(-0.17624804) q[1];
sx q[1];
rz(-1.9269678) q[1];
sx q[1];
rz(2.2252749) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57740649) q[0];
sx q[0];
rz(-0.86125492) q[0];
sx q[0];
rz(0.90155154) q[0];
rz(-3.0099478) q[2];
sx q[2];
rz(-2.3028214) q[2];
sx q[2];
rz(-1.4721578) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.037109) q[1];
sx q[1];
rz(-2.6625405) q[1];
sx q[1];
rz(-0.80161174) q[1];
x q[2];
rz(-1.6244435) q[3];
sx q[3];
rz(-2.7430153) q[3];
sx q[3];
rz(2.9159391) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.9006485) q[2];
sx q[2];
rz(-2.1449461) q[2];
sx q[2];
rz(2.3201578) q[2];
rz(0.017283043) q[3];
sx q[3];
rz(-2.343943) q[3];
sx q[3];
rz(0.78330529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18773742) q[0];
sx q[0];
rz(-1.5783577) q[0];
sx q[0];
rz(1.0082555) q[0];
rz(-3.1058274) q[1];
sx q[1];
rz(-1.6556031) q[1];
sx q[1];
rz(0.52454138) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4315955) q[0];
sx q[0];
rz(-1.734373) q[0];
sx q[0];
rz(3.1113935) q[0];
x q[1];
rz(0.67851615) q[2];
sx q[2];
rz(-1.3635474) q[2];
sx q[2];
rz(3.0845272) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.14568612) q[1];
sx q[1];
rz(-1.5090764) q[1];
sx q[1];
rz(-1.5476336) q[1];
rz(-1.2479765) q[3];
sx q[3];
rz(-0.4399235) q[3];
sx q[3];
rz(-0.80252121) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.3699469) q[2];
sx q[2];
rz(-1.6235855) q[2];
sx q[2];
rz(-1.4952205) q[2];
rz(-1.8203991) q[3];
sx q[3];
rz(-1.4097872) q[3];
sx q[3];
rz(-0.43740073) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8111073) q[0];
sx q[0];
rz(-0.91902584) q[0];
sx q[0];
rz(-0.088949732) q[0];
rz(0.51070172) q[1];
sx q[1];
rz(-0.82534868) q[1];
sx q[1];
rz(-0.68960062) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21259637) q[0];
sx q[0];
rz(-1.4082452) q[0];
sx q[0];
rz(0.29851229) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7933649) q[2];
sx q[2];
rz(-0.59154445) q[2];
sx q[2];
rz(-1.0500184) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.00011132414) q[1];
sx q[1];
rz(-2.3173216) q[1];
sx q[1];
rz(0.4517171) q[1];
x q[2];
rz(1.425399) q[3];
sx q[3];
rz(-2.0510011) q[3];
sx q[3];
rz(-2.8499545) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.4758063) q[2];
sx q[2];
rz(-0.67725724) q[2];
sx q[2];
rz(2.518667) q[2];
rz(1.1359435) q[3];
sx q[3];
rz(-0.1853075) q[3];
sx q[3];
rz(2.6749271) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85161197) q[0];
sx q[0];
rz(-1.2427793) q[0];
sx q[0];
rz(-2.5581397) q[0];
rz(1.1460229) q[1];
sx q[1];
rz(-1.6505516) q[1];
sx q[1];
rz(1.4978283) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4125815) q[0];
sx q[0];
rz(-0.81775613) q[0];
sx q[0];
rz(1.2110932) q[0];
x q[1];
rz(0.096502467) q[2];
sx q[2];
rz(-2.3415903) q[2];
sx q[2];
rz(-2.2068791) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.0651107) q[1];
sx q[1];
rz(-1.5981734) q[1];
sx q[1];
rz(-0.70365023) q[1];
x q[2];
rz(-1.4086401) q[3];
sx q[3];
rz(-2.080653) q[3];
sx q[3];
rz(1.4354524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.4862165) q[2];
sx q[2];
rz(-1.5465753) q[2];
sx q[2];
rz(1.5636469) q[2];
rz(2.2359713) q[3];
sx q[3];
rz(-2.8712397) q[3];
sx q[3];
rz(0.63846987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49801302) q[0];
sx q[0];
rz(-1.6828515) q[0];
sx q[0];
rz(-3.0946099) q[0];
rz(0.14818305) q[1];
sx q[1];
rz(-0.89887416) q[1];
sx q[1];
rz(1.4354338) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5934138) q[0];
sx q[0];
rz(-2.3313064) q[0];
sx q[0];
rz(-2.9368648) q[0];
rz(1.4749182) q[2];
sx q[2];
rz(-2.5898858) q[2];
sx q[2];
rz(0.72223896) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6378577) q[1];
sx q[1];
rz(-1.4459472) q[1];
sx q[1];
rz(-1.1605074) q[1];
rz(-pi) q[2];
rz(0.027408882) q[3];
sx q[3];
rz(-1.3331183) q[3];
sx q[3];
rz(-0.72565597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.0753714) q[2];
sx q[2];
rz(-2.1263289) q[2];
sx q[2];
rz(-3.0701239) q[2];
rz(1.6890769) q[3];
sx q[3];
rz(-1.444933) q[3];
sx q[3];
rz(0.92648363) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28618318) q[0];
sx q[0];
rz(-1.8137285) q[0];
sx q[0];
rz(-2.563971) q[0];
rz(-1.8619934) q[1];
sx q[1];
rz(-1.3459233) q[1];
sx q[1];
rz(2.1320027) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9855232) q[0];
sx q[0];
rz(-2.0445604) q[0];
sx q[0];
rz(-1.9457293) q[0];
rz(2.7570653) q[2];
sx q[2];
rz(-1.3587917) q[2];
sx q[2];
rz(0.42696135) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.13751444) q[1];
sx q[1];
rz(-0.61811781) q[1];
sx q[1];
rz(0.60933463) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.61543492) q[3];
sx q[3];
rz(-2.2781567) q[3];
sx q[3];
rz(0.59188852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.0499095) q[2];
sx q[2];
rz(-2.78806) q[2];
sx q[2];
rz(1.1996777) q[2];
rz(-0.47232929) q[3];
sx q[3];
rz(-1.4940558) q[3];
sx q[3];
rz(-1.002731) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6427479) q[0];
sx q[0];
rz(-1.5002748) q[0];
sx q[0];
rz(2.902466) q[0];
rz(0.7827951) q[1];
sx q[1];
rz(-2.6233964) q[1];
sx q[1];
rz(-2.696864) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1040161) q[0];
sx q[0];
rz(-0.74066478) q[0];
sx q[0];
rz(-2.2292577) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2112446) q[2];
sx q[2];
rz(-1.010251) q[2];
sx q[2];
rz(-2.408574) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.4152894) q[1];
sx q[1];
rz(-2.2518034) q[1];
sx q[1];
rz(-0.16420941) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4374251) q[3];
sx q[3];
rz(-1.9786454) q[3];
sx q[3];
rz(2.9123902) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.42177054) q[2];
sx q[2];
rz(-0.82227102) q[2];
sx q[2];
rz(-0.81531173) q[2];
rz(2.1067965) q[3];
sx q[3];
rz(-1.5765604) q[3];
sx q[3];
rz(1.8243272) q[3];
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
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3141044) q[0];
sx q[0];
rz(-2.1056471) q[0];
sx q[0];
rz(0.28717336) q[0];
rz(0.18889591) q[1];
sx q[1];
rz(-2.4177528) q[1];
sx q[1];
rz(2.8093991) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6824326) q[0];
sx q[0];
rz(-2.3258665) q[0];
sx q[0];
rz(-2.0072323) q[0];
rz(-pi) q[1];
rz(0.8719445) q[2];
sx q[2];
rz(-1.6064062) q[2];
sx q[2];
rz(-1.9277089) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.48391446) q[1];
sx q[1];
rz(-2.6473443) q[1];
sx q[1];
rz(2.5792522) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.513047) q[3];
sx q[3];
rz(-2.5124031) q[3];
sx q[3];
rz(2.7043846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2087848) q[2];
sx q[2];
rz(-2.1366426) q[2];
sx q[2];
rz(0.83958158) q[2];
rz(1.8509289) q[3];
sx q[3];
rz(-1.3495812) q[3];
sx q[3];
rz(0.65657842) q[3];
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
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.055450913) q[0];
sx q[0];
rz(-0.76776004) q[0];
sx q[0];
rz(-2.0822051) q[0];
rz(-2.1620031) q[1];
sx q[1];
rz(-1.9444119) q[1];
sx q[1];
rz(0.25451452) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8043038) q[0];
sx q[0];
rz(-0.24233195) q[0];
sx q[0];
rz(-1.2605577) q[0];
x q[1];
rz(-0.48760957) q[2];
sx q[2];
rz(-1.7874103) q[2];
sx q[2];
rz(0.85607869) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.49452457) q[1];
sx q[1];
rz(-1.8956603) q[1];
sx q[1];
rz(1.4855794) q[1];
x q[2];
rz(-1.47154) q[3];
sx q[3];
rz(-2.0072862) q[3];
sx q[3];
rz(2.1193159) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.09482) q[2];
sx q[2];
rz(-0.54868713) q[2];
sx q[2];
rz(-1.0478896) q[2];
rz(-1.7808328) q[3];
sx q[3];
rz(-2.4436617) q[3];
sx q[3];
rz(-0.60539436) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.027325252) q[0];
sx q[0];
rz(-2.0226759) q[0];
sx q[0];
rz(-0.080060536) q[0];
rz(-2.7813773) q[1];
sx q[1];
rz(-1.6811014) q[1];
sx q[1];
rz(-0.95492687) q[1];
rz(-1.9706456) q[2];
sx q[2];
rz(-2.5794537) q[2];
sx q[2];
rz(-0.07515547) q[2];
rz(-1.8923106) q[3];
sx q[3];
rz(-0.63890639) q[3];
sx q[3];
rz(-1.3756868) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
