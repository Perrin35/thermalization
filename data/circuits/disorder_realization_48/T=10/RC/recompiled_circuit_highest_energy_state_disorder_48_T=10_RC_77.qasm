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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9636757) q[0];
sx q[0];
rz(-2.6612894) q[0];
sx q[0];
rz(0.98573523) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.79390766) q[2];
sx q[2];
rz(-2.53144) q[2];
sx q[2];
rz(-2.8086503) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4922766) q[1];
sx q[1];
rz(-1.5879803) q[1];
sx q[1];
rz(-1.3130929) q[1];
x q[2];
rz(0.4313978) q[3];
sx q[3];
rz(-1.3702624) q[3];
sx q[3];
rz(-1.3667184) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.013926355) q[2];
sx q[2];
rz(-1.650859) q[2];
sx q[2];
rz(0.15065436) q[2];
rz(-1.539218) q[3];
sx q[3];
rz(-2.7645002) q[3];
sx q[3];
rz(-2.8902875) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0852614) q[0];
sx q[0];
rz(-1.3717317) q[0];
sx q[0];
rz(-0.37013176) q[0];
rz(2.6689957) q[1];
sx q[1];
rz(-2.8234146) q[1];
sx q[1];
rz(2.1841689) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8214398) q[0];
sx q[0];
rz(-1.4906881) q[0];
sx q[0];
rz(2.0700685) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9198869) q[2];
sx q[2];
rz(-2.6654976) q[2];
sx q[2];
rz(-2.8698336) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.53471662) q[1];
sx q[1];
rz(-0.83321111) q[1];
sx q[1];
rz(-1.4556769) q[1];
rz(-2.4647275) q[3];
sx q[3];
rz(-1.8329304) q[3];
sx q[3];
rz(1.2659951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.0051673278) q[2];
sx q[2];
rz(-0.74446669) q[2];
sx q[2];
rz(1.0235419) q[2];
rz(-1.7510022) q[3];
sx q[3];
rz(-1.7460881) q[3];
sx q[3];
rz(-1.3172147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7904974) q[0];
sx q[0];
rz(-2.1320765) q[0];
sx q[0];
rz(2.5787831) q[0];
rz(-1.6273181) q[1];
sx q[1];
rz(-1.7104251) q[1];
sx q[1];
rz(-3.0624342) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.889733) q[0];
sx q[0];
rz(-1.8859698) q[0];
sx q[0];
rz(-2.6227975) q[0];
rz(-pi) q[1];
rz(3.0653238) q[2];
sx q[2];
rz(-2.6496716) q[2];
sx q[2];
rz(-0.15215506) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.2217949) q[1];
sx q[1];
rz(-0.88558965) q[1];
sx q[1];
rz(2.5849839) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.78000809) q[3];
sx q[3];
rz(-0.87015188) q[3];
sx q[3];
rz(2.8740945) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.2341653) q[2];
sx q[2];
rz(-1.5831466) q[2];
sx q[2];
rz(1.3321715) q[2];
rz(2.7857156) q[3];
sx q[3];
rz(-1.9477113) q[3];
sx q[3];
rz(0.98085421) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6818162) q[0];
sx q[0];
rz(-1.1695319) q[0];
sx q[0];
rz(-2.0303149) q[0];
rz(2.2214644) q[1];
sx q[1];
rz(-1.5524813) q[1];
sx q[1];
rz(2.8880602) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0363844) q[0];
sx q[0];
rz(-2.0724793) q[0];
sx q[0];
rz(1.2060449) q[0];
x q[1];
rz(1.6854671) q[2];
sx q[2];
rz(-1.7301699) q[2];
sx q[2];
rz(-1.3068975) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.82099709) q[1];
sx q[1];
rz(-1.3861457) q[1];
sx q[1];
rz(2.3750633) q[1];
rz(-2.0757789) q[3];
sx q[3];
rz(-1.2974713) q[3];
sx q[3];
rz(1.8718639) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.8186875) q[2];
sx q[2];
rz(-2.5762168) q[2];
sx q[2];
rz(1.6820071) q[2];
rz(-2.5679576) q[3];
sx q[3];
rz(-1.9394268) q[3];
sx q[3];
rz(3.0972163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6673073) q[0];
sx q[0];
rz(-1.0695589) q[0];
sx q[0];
rz(1.7942418) q[0];
rz(2.5509293) q[1];
sx q[1];
rz(-1.2202411) q[1];
sx q[1];
rz(-1.4264872) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7280557) q[0];
sx q[0];
rz(-1.1794834) q[0];
sx q[0];
rz(1.9584459) q[0];
x q[1];
rz(1.1828184) q[2];
sx q[2];
rz(-2.175594) q[2];
sx q[2];
rz(-1.4825578) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.3326753) q[1];
sx q[1];
rz(-1.6170224) q[1];
sx q[1];
rz(-2.8874365) q[1];
x q[2];
rz(2.9699202) q[3];
sx q[3];
rz(-0.35997691) q[3];
sx q[3];
rz(1.6115007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.9023989) q[2];
sx q[2];
rz(-1.0603443) q[2];
sx q[2];
rz(1.9690465) q[2];
rz(3.0555365) q[3];
sx q[3];
rz(-0.9525758) q[3];
sx q[3];
rz(2.9854767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0657144) q[0];
sx q[0];
rz(-1.5471764) q[0];
sx q[0];
rz(-0.86268798) q[0];
rz(2.8942096) q[1];
sx q[1];
rz(-1.3037953) q[1];
sx q[1];
rz(-2.6695796) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87403546) q[0];
sx q[0];
rz(-2.3568601) q[0];
sx q[0];
rz(-2.6536056) q[0];
rz(-0.15100592) q[2];
sx q[2];
rz(-1.0698294) q[2];
sx q[2];
rz(-0.46268625) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.5322383) q[1];
sx q[1];
rz(-1.8731469) q[1];
sx q[1];
rz(1.3855893) q[1];
rz(-pi) q[2];
rz(-2.1930027) q[3];
sx q[3];
rz(-1.7441809) q[3];
sx q[3];
rz(-2.9297048) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.90999675) q[2];
sx q[2];
rz(-2.1746077) q[2];
sx q[2];
rz(1.3989353) q[2];
rz(-1.4553962) q[3];
sx q[3];
rz(-0.78794909) q[3];
sx q[3];
rz(-0.56070352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6550734) q[0];
sx q[0];
rz(-1.1460679) q[0];
sx q[0];
rz(0.51862496) q[0];
rz(0.71867603) q[1];
sx q[1];
rz(-2.6287754) q[1];
sx q[1];
rz(1.3291043) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78120056) q[0];
sx q[0];
rz(-0.5628764) q[0];
sx q[0];
rz(-0.84873523) q[0];
rz(1.8480214) q[2];
sx q[2];
rz(-1.5223862) q[2];
sx q[2];
rz(-1.8579661) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.74187255) q[1];
sx q[1];
rz(-2.406027) q[1];
sx q[1];
rz(-1.5217785) q[1];
x q[2];
rz(-1.5409327) q[3];
sx q[3];
rz(-0.58664188) q[3];
sx q[3];
rz(1.3636953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3507877) q[2];
sx q[2];
rz(-0.69956508) q[2];
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
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
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
rz(1.4627948) q[0];
sx q[0];
rz(-0.21634732) q[0];
sx q[0];
rz(1.3702673) q[0];
rz(-2.7554152) q[1];
sx q[1];
rz(-1.4048856) q[1];
sx q[1];
rz(0.6699627) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6813107) q[0];
sx q[0];
rz(-0.71397266) q[0];
sx q[0];
rz(2.332469) q[0];
rz(-1.6514844) q[2];
sx q[2];
rz(-0.85570691) q[2];
sx q[2];
rz(-0.87835588) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.6584864) q[1];
sx q[1];
rz(-1.4577333) q[1];
sx q[1];
rz(2.1674564) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3988204) q[3];
sx q[3];
rz(-2.9265518) q[3];
sx q[3];
rz(-2.5701534) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.92114821) q[2];
sx q[2];
rz(-0.28382742) q[2];
sx q[2];
rz(-2.0849126) q[2];
rz(1.7250666) q[3];
sx q[3];
rz(-1.8959911) q[3];
sx q[3];
rz(1.2904803) q[3];
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
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.073762745) q[0];
sx q[0];
rz(-2.1011598) q[0];
sx q[0];
rz(-1.006806) q[0];
rz(-2.6489068) q[1];
sx q[1];
rz(-1.7308116) q[1];
sx q[1];
rz(-2.3115092) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7573038) q[0];
sx q[0];
rz(-0.28097935) q[0];
sx q[0];
rz(-0.12019867) q[0];
x q[1];
rz(-1.3557498) q[2];
sx q[2];
rz(-2.1587662) q[2];
sx q[2];
rz(0.34823349) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.600297) q[1];
sx q[1];
rz(-1.4182555) q[1];
sx q[1];
rz(-0.18220724) q[1];
rz(-pi) q[2];
rz(2.0138143) q[3];
sx q[3];
rz(-1.5823871) q[3];
sx q[3];
rz(-3.1168337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.81862187) q[2];
sx q[2];
rz(-1.3016204) q[2];
sx q[2];
rz(-1.0969561) q[2];
rz(1.2777626) q[3];
sx q[3];
rz(-2.4570229) q[3];
sx q[3];
rz(2.4975615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.291236) q[0];
sx q[0];
rz(-2.0391897) q[0];
sx q[0];
rz(-2.7665603) q[0];
rz(-0.11861435) q[1];
sx q[1];
rz(-1.7053441) q[1];
sx q[1];
rz(-2.2172701) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7345006) q[0];
sx q[0];
rz(-2.1372927) q[0];
sx q[0];
rz(-1.9167711) q[0];
x q[1];
rz(0.61225981) q[2];
sx q[2];
rz(-2.3697402) q[2];
sx q[2];
rz(-1.5764232) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3872747) q[1];
sx q[1];
rz(-1.5994834) q[1];
sx q[1];
rz(-2.7806675) q[1];
rz(1.9162991) q[3];
sx q[3];
rz(-1.1456819) q[3];
sx q[3];
rz(-2.2540852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.2529926) q[2];
sx q[2];
rz(-2.1414976) q[2];
sx q[2];
rz(-0.89317733) q[2];
rz(2.0231953) q[3];
sx q[3];
rz(-1.018254) q[3];
sx q[3];
rz(0.65677381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7669582) q[0];
sx q[0];
rz(-2.0603016) q[0];
sx q[0];
rz(1.7745071) q[0];
rz(2.959666) q[1];
sx q[1];
rz(-2.5839099) q[1];
sx q[1];
rz(-0.90669496) q[1];
rz(1.8461377) q[2];
sx q[2];
rz(-0.32651797) q[2];
sx q[2];
rz(1.4985195) q[2];
rz(1.160865) q[3];
sx q[3];
rz(-2.0515473) q[3];
sx q[3];
rz(-0.69507364) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
