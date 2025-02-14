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
rz(-1.7464632) q[0];
sx q[0];
rz(-0.75463086) q[0];
sx q[0];
rz(-2.057743) q[0];
rz(2.2639182) q[1];
sx q[1];
rz(4.35507) q[1];
sx q[1];
rz(9.9014643) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81945588) q[0];
sx q[0];
rz(-1.9662204) q[0];
sx q[0];
rz(-2.8614392) q[0];
rz(-2.0332912) q[2];
sx q[2];
rz(-1.1574189) q[2];
sx q[2];
rz(-0.55962901) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(3.1281761) q[1];
sx q[1];
rz(-2.8833296) q[1];
sx q[1];
rz(1.5034666) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.45248078) q[3];
sx q[3];
rz(-0.47305952) q[3];
sx q[3];
rz(2.5291131) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-3.1276663) q[2];
sx q[2];
rz(-1.4907336) q[2];
sx q[2];
rz(-2.9909383) q[2];
rz(1.539218) q[3];
sx q[3];
rz(-2.7645002) q[3];
sx q[3];
rz(-0.25130513) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(0.056331228) q[0];
sx q[0];
rz(-1.769861) q[0];
sx q[0];
rz(-2.7714609) q[0];
rz(-2.6689957) q[1];
sx q[1];
rz(-2.8234146) q[1];
sx q[1];
rz(-2.1841689) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9345624) q[0];
sx q[0];
rz(-1.0732722) q[0];
sx q[0];
rz(-0.091188641) q[0];
rz(2.967011) q[2];
sx q[2];
rz(-2.0160297) q[2];
sx q[2];
rz(-3.0246459) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1831187) q[1];
sx q[1];
rz(-1.4856824) q[1];
sx q[1];
rz(-2.4007011) q[1];
x q[2];
rz(2.4647275) q[3];
sx q[3];
rz(-1.3086623) q[3];
sx q[3];
rz(-1.8755975) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-3.1364253) q[2];
sx q[2];
rz(-0.74446669) q[2];
sx q[2];
rz(-2.1180507) q[2];
rz(1.3905904) q[3];
sx q[3];
rz(-1.7460881) q[3];
sx q[3];
rz(1.8243779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35109529) q[0];
sx q[0];
rz(-2.1320765) q[0];
sx q[0];
rz(2.5787831) q[0];
rz(-1.5142745) q[1];
sx q[1];
rz(-1.7104251) q[1];
sx q[1];
rz(3.0624342) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3252727) q[0];
sx q[0];
rz(-0.59945852) q[0];
sx q[0];
rz(-2.5599203) q[0];
x q[1];
rz(-3.0653238) q[2];
sx q[2];
rz(-2.6496716) q[2];
sx q[2];
rz(-2.9894376) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.9197977) q[1];
sx q[1];
rz(-0.88558965) q[1];
sx q[1];
rz(0.55660875) q[1];
rz(-pi) q[2];
rz(-0.8757365) q[3];
sx q[3];
rz(-2.1453224) q[3];
sx q[3];
rz(-0.72573001) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.2341653) q[2];
sx q[2];
rz(-1.558446) q[2];
sx q[2];
rz(1.3321715) q[2];
rz(0.35587707) q[3];
sx q[3];
rz(-1.9477113) q[3];
sx q[3];
rz(-0.98085421) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45977649) q[0];
sx q[0];
rz(-1.1695319) q[0];
sx q[0];
rz(2.0303149) q[0];
rz(-2.2214644) q[1];
sx q[1];
rz(-1.5891113) q[1];
sx q[1];
rz(2.8880602) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64719114) q[0];
sx q[0];
rz(-1.2526985) q[0];
sx q[0];
rz(0.53089106) q[0];
rz(0.16040921) q[2];
sx q[2];
rz(-1.4575851) q[2];
sx q[2];
rz(0.28217523) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.56139466) q[1];
sx q[1];
rz(-0.78402482) q[1];
sx q[1];
rz(-2.8785588) q[1];
rz(-pi) q[2];
x q[2];
rz(0.30999581) q[3];
sx q[3];
rz(-1.0862203) q[3];
sx q[3];
rz(0.44919168) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.8186875) q[2];
sx q[2];
rz(-0.56537586) q[2];
sx q[2];
rz(1.6820071) q[2];
rz(-2.5679576) q[3];
sx q[3];
rz(-1.2021659) q[3];
sx q[3];
rz(0.044376317) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6673073) q[0];
sx q[0];
rz(-1.0695589) q[0];
sx q[0];
rz(1.3473508) q[0];
rz(-0.59066331) q[1];
sx q[1];
rz(-1.9213516) q[1];
sx q[1];
rz(1.4264872) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0027699) q[0];
sx q[0];
rz(-1.213824) q[0];
sx q[0];
rz(2.7223552) q[0];
x q[1];
rz(-1.9587742) q[2];
sx q[2];
rz(-0.96599865) q[2];
sx q[2];
rz(-1.6590349) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.8914681) q[1];
sx q[1];
rz(-1.3169177) q[1];
sx q[1];
rz(1.6185544) q[1];
rz(1.5065881) q[3];
sx q[3];
rz(-1.2163463) q[3];
sx q[3];
rz(-1.7946769) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.9023989) q[2];
sx q[2];
rz(-1.0603443) q[2];
sx q[2];
rz(-1.9690465) q[2];
rz(-3.0555365) q[3];
sx q[3];
rz(-0.9525758) q[3];
sx q[3];
rz(0.15611592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0758783) q[0];
sx q[0];
rz(-1.5944163) q[0];
sx q[0];
rz(-2.2789047) q[0];
rz(0.24738303) q[1];
sx q[1];
rz(-1.3037953) q[1];
sx q[1];
rz(2.6695796) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9111689) q[0];
sx q[0];
rz(-0.8967451) q[0];
sx q[0];
rz(-1.1328902) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8389147) q[2];
sx q[2];
rz(-0.52137085) q[2];
sx q[2];
rz(2.9857295) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.5322383) q[1];
sx q[1];
rz(-1.2684457) q[1];
sx q[1];
rz(1.7560033) q[1];
rz(2.1930027) q[3];
sx q[3];
rz(-1.7441809) q[3];
sx q[3];
rz(-0.21188785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.90999675) q[2];
sx q[2];
rz(-0.96698499) q[2];
sx q[2];
rz(-1.7426573) q[2];
rz(1.4553962) q[3];
sx q[3];
rz(-2.3536436) q[3];
sx q[3];
rz(-0.56070352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48651925) q[0];
sx q[0];
rz(-1.9955248) q[0];
sx q[0];
rz(-2.6229677) q[0];
rz(2.4229166) q[1];
sx q[1];
rz(-2.6287754) q[1];
sx q[1];
rz(-1.3291043) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3603921) q[0];
sx q[0];
rz(-0.5628764) q[0];
sx q[0];
rz(0.84873523) q[0];
rz(-pi) q[1];
rz(1.2935712) q[2];
sx q[2];
rz(-1.6192064) q[2];
sx q[2];
rz(1.2836266) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.79256637) q[1];
sx q[1];
rz(-1.6036803) q[1];
sx q[1];
rz(0.83582857) q[1];
rz(-pi) q[2];
rz(0.98436004) q[3];
sx q[3];
rz(-1.5873261) q[3];
sx q[3];
rz(-0.23197385) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.79080498) q[2];
sx q[2];
rz(-2.4420276) q[2];
sx q[2];
rz(-3.1034071) q[2];
rz(-0.83032483) q[3];
sx q[3];
rz(-1.3874976) q[3];
sx q[3];
rz(0.047688095) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6787978) q[0];
sx q[0];
rz(-0.21634732) q[0];
sx q[0];
rz(1.7713254) q[0];
rz(2.7554152) q[1];
sx q[1];
rz(-1.4048856) q[1];
sx q[1];
rz(2.4716299) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6813107) q[0];
sx q[0];
rz(-0.71397266) q[0];
sx q[0];
rz(-0.80912368) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6514844) q[2];
sx q[2];
rz(-0.85570691) q[2];
sx q[2];
rz(-2.2632368) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.48310623) q[1];
sx q[1];
rz(-1.6838594) q[1];
sx q[1];
rz(2.1674564) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3588328) q[3];
sx q[3];
rz(-1.6073213) q[3];
sx q[3];
rz(-1.9741457) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.92114821) q[2];
sx q[2];
rz(-2.8577652) q[2];
sx q[2];
rz(2.0849126) q[2];
rz(1.4165261) q[3];
sx q[3];
rz(-1.8959911) q[3];
sx q[3];
rz(1.8511124) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0678299) q[0];
sx q[0];
rz(-2.1011598) q[0];
sx q[0];
rz(1.006806) q[0];
rz(2.6489068) q[1];
sx q[1];
rz(-1.7308116) q[1];
sx q[1];
rz(-0.83008343) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3842889) q[0];
sx q[0];
rz(-0.28097935) q[0];
sx q[0];
rz(3.021394) q[0];
rz(-2.5428204) q[2];
sx q[2];
rz(-1.3922924) q[2];
sx q[2];
rz(-1.1019966) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.54129564) q[1];
sx q[1];
rz(-1.4182555) q[1];
sx q[1];
rz(2.9593854) q[1];
x q[2];
rz(-1.1277783) q[3];
sx q[3];
rz(-1.5823871) q[3];
sx q[3];
rz(-3.1168337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.3229708) q[2];
sx q[2];
rz(-1.3016204) q[2];
sx q[2];
rz(1.0969561) q[2];
rz(-1.8638301) q[3];
sx q[3];
rz(-0.68456972) q[3];
sx q[3];
rz(0.6440312) q[3];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85035664) q[0];
sx q[0];
rz(-2.0391897) q[0];
sx q[0];
rz(2.7665603) q[0];
rz(-3.0229783) q[1];
sx q[1];
rz(-1.7053441) q[1];
sx q[1];
rz(-0.92432252) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7345006) q[0];
sx q[0];
rz(-1.0043) q[0];
sx q[0];
rz(-1.2248216) q[0];
rz(-pi) q[1];
rz(2.4690041) q[2];
sx q[2];
rz(-1.9832356) q[2];
sx q[2];
rz(2.6809566) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.754318) q[1];
sx q[1];
rz(-1.5994834) q[1];
sx q[1];
rz(2.7806675) q[1];
rz(-1.9162991) q[3];
sx q[3];
rz(-1.9959108) q[3];
sx q[3];
rz(0.88750741) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.2529926) q[2];
sx q[2];
rz(-1.000095) q[2];
sx q[2];
rz(-0.89317733) q[2];
rz(-1.1183974) q[3];
sx q[3];
rz(-2.1233386) q[3];
sx q[3];
rz(2.4848188) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3746344) q[0];
sx q[0];
rz(-1.081291) q[0];
sx q[0];
rz(-1.3670856) q[0];
rz(-2.959666) q[1];
sx q[1];
rz(-0.55768273) q[1];
sx q[1];
rz(2.2348977) q[1];
rz(-1.295455) q[2];
sx q[2];
rz(-0.32651797) q[2];
sx q[2];
rz(1.4985195) q[2];
rz(-0.652486) q[3];
sx q[3];
rz(-0.62118821) q[3];
sx q[3];
rz(-1.4493829) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
