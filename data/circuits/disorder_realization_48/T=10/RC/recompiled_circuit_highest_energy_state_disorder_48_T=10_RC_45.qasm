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
rz(1.3951294) q[0];
sx q[0];
rz(-2.3869618) q[0];
sx q[0];
rz(-1.0838497) q[0];
rz(-0.87767449) q[1];
sx q[1];
rz(-1.2134774) q[1];
sx q[1];
rz(2.6649063) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81945588) q[0];
sx q[0];
rz(-1.1753723) q[0];
sx q[0];
rz(2.8614392) q[0];
rz(-pi) q[1];
rz(1.1083015) q[2];
sx q[2];
rz(-1.9841737) q[2];
sx q[2];
rz(0.55962901) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.4922766) q[1];
sx q[1];
rz(-1.5879803) q[1];
sx q[1];
rz(-1.8284998) q[1];
rz(1.3506557) q[3];
sx q[3];
rz(-1.992989) q[3];
sx q[3];
rz(3.0289502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-3.1276663) q[2];
sx q[2];
rz(-1.650859) q[2];
sx q[2];
rz(2.9909383) q[2];
rz(-1.6023747) q[3];
sx q[3];
rz(-2.7645002) q[3];
sx q[3];
rz(2.8902875) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0852614) q[0];
sx q[0];
rz(-1.3717317) q[0];
sx q[0];
rz(-0.37013176) q[0];
rz(0.47259694) q[1];
sx q[1];
rz(-0.31817803) q[1];
sx q[1];
rz(2.1841689) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8214398) q[0];
sx q[0];
rz(-1.4906881) q[0];
sx q[0];
rz(-1.0715241) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1195807) q[2];
sx q[2];
rz(-1.4133845) q[2];
sx q[2];
rz(-1.6119286) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.606876) q[1];
sx q[1];
rz(-0.83321111) q[1];
sx q[1];
rz(-1.4556769) q[1];
rz(-0.67686512) q[3];
sx q[3];
rz(-1.3086623) q[3];
sx q[3];
rz(1.2659951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.0051673278) q[2];
sx q[2];
rz(-2.397126) q[2];
sx q[2];
rz(-1.0235419) q[2];
rz(-1.7510022) q[3];
sx q[3];
rz(-1.7460881) q[3];
sx q[3];
rz(1.8243779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35109529) q[0];
sx q[0];
rz(-1.0095162) q[0];
sx q[0];
rz(0.56280953) q[0];
rz(1.6273181) q[1];
sx q[1];
rz(-1.4311675) q[1];
sx q[1];
rz(-3.0624342) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9978299) q[0];
sx q[0];
rz(-1.0799066) q[0];
sx q[0];
rz(-1.9299555) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6508826) q[2];
sx q[2];
rz(-1.5348002) q[2];
sx q[2];
rz(-1.7902059) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.1428433) q[1];
sx q[1];
rz(-2.2881737) q[1];
sx q[1];
rz(-0.99695506) q[1];
x q[2];
rz(-0.8757365) q[3];
sx q[3];
rz(-0.99627021) q[3];
sx q[3];
rz(0.72573001) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.2341653) q[2];
sx q[2];
rz(-1.5831466) q[2];
sx q[2];
rz(-1.8094212) q[2];
rz(2.7857156) q[3];
sx q[3];
rz(-1.9477113) q[3];
sx q[3];
rz(0.98085421) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(0.45977649) q[0];
sx q[0];
rz(-1.1695319) q[0];
sx q[0];
rz(-2.0303149) q[0];
rz(0.92012826) q[1];
sx q[1];
rz(-1.5524813) q[1];
sx q[1];
rz(-2.8880602) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43413517) q[0];
sx q[0];
rz(-0.61096707) q[0];
sx q[0];
rz(-2.5649628) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5229613) q[2];
sx q[2];
rz(-0.19605532) q[2];
sx q[2];
rz(-2.4624937) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2168) q[1];
sx q[1];
rz(-2.3210822) q[1];
sx q[1];
rz(-1.8245068) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0757789) q[3];
sx q[3];
rz(-1.2974713) q[3];
sx q[3];
rz(-1.2697288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3229052) q[2];
sx q[2];
rz(-2.5762168) q[2];
sx q[2];
rz(-1.4595855) q[2];
rz(-2.5679576) q[3];
sx q[3];
rz(-1.9394268) q[3];
sx q[3];
rz(-0.044376317) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6673073) q[0];
sx q[0];
rz(-2.0720338) q[0];
sx q[0];
rz(1.7942418) q[0];
rz(-0.59066331) q[1];
sx q[1];
rz(-1.9213516) q[1];
sx q[1];
rz(1.4264872) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41353696) q[0];
sx q[0];
rz(-1.9621092) q[0];
sx q[0];
rz(-1.1831468) q[0];
rz(-1.9587742) q[2];
sx q[2];
rz(-0.96599865) q[2];
sx q[2];
rz(1.4825578) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.8089174) q[1];
sx q[1];
rz(-1.5245702) q[1];
sx q[1];
rz(0.25415616) q[1];
rz(-pi) q[2];
rz(2.9699202) q[3];
sx q[3];
rz(-2.7816157) q[3];
sx q[3];
rz(1.5300919) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.23919375) q[2];
sx q[2];
rz(-1.0603443) q[2];
sx q[2];
rz(1.1725461) q[2];
rz(-3.0555365) q[3];
sx q[3];
rz(-0.9525758) q[3];
sx q[3];
rz(-2.9854767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0657144) q[0];
sx q[0];
rz(-1.5471764) q[0];
sx q[0];
rz(2.2789047) q[0];
rz(0.24738303) q[1];
sx q[1];
rz(-1.3037953) q[1];
sx q[1];
rz(2.6695796) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.056045) q[0];
sx q[0];
rz(-1.9084832) q[0];
sx q[0];
rz(-0.72283904) q[0];
x q[1];
rz(0.15100592) q[2];
sx q[2];
rz(-2.0717632) q[2];
sx q[2];
rz(-0.46268625) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6093543) q[1];
sx q[1];
rz(-1.2684457) q[1];
sx q[1];
rz(1.7560033) q[1];
rz(2.9293044) q[3];
sx q[3];
rz(-2.1822896) q[3];
sx q[3];
rz(1.4820198) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.90999675) q[2];
sx q[2];
rz(-0.96698499) q[2];
sx q[2];
rz(-1.7426573) q[2];
rz(-1.6861964) q[3];
sx q[3];
rz(-0.78794909) q[3];
sx q[3];
rz(0.56070352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48651925) q[0];
sx q[0];
rz(-1.9955248) q[0];
sx q[0];
rz(-2.6229677) q[0];
rz(-0.71867603) q[1];
sx q[1];
rz(-0.51281723) q[1];
sx q[1];
rz(-1.8124883) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14940093) q[0];
sx q[0];
rz(-1.9312381) q[0];
sx q[0];
rz(-2.0130231) q[0];
rz(-pi) q[1];
rz(1.2935712) q[2];
sx q[2];
rz(-1.6192064) q[2];
sx q[2];
rz(1.2836266) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.74187255) q[1];
sx q[1];
rz(-0.73556566) q[1];
sx q[1];
rz(-1.6198141) q[1];
rz(-pi) q[2];
rz(0.98436004) q[3];
sx q[3];
rz(-1.5542665) q[3];
sx q[3];
rz(0.23197385) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3507877) q[2];
sx q[2];
rz(-0.69956508) q[2];
sx q[2];
rz(0.038185509) q[2];
rz(-2.3112678) q[3];
sx q[3];
rz(-1.7540951) q[3];
sx q[3];
rz(0.047688095) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6787978) q[0];
sx q[0];
rz(-0.21634732) q[0];
sx q[0];
rz(-1.7713254) q[0];
rz(-2.7554152) q[1];
sx q[1];
rz(-1.4048856) q[1];
sx q[1];
rz(0.6699627) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7350175) q[0];
sx q[0];
rz(-1.1018714) q[0];
sx q[0];
rz(2.1308628) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4248895) q[2];
sx q[2];
rz(-1.5099025) q[2];
sx q[2];
rz(-2.3961807) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.6584864) q[1];
sx q[1];
rz(-1.6838594) q[1];
sx q[1];
rz(-0.97413625) q[1];
rz(-pi) q[2];
rz(0.037360351) q[3];
sx q[3];
rz(-1.3589763) q[3];
sx q[3];
rz(-0.39549124) q[3];
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
rz(-1.2456015) q[3];
sx q[3];
rz(-1.8511124) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0678299) q[0];
sx q[0];
rz(-2.1011598) q[0];
sx q[0];
rz(2.1347866) q[0];
rz(-0.49268588) q[1];
sx q[1];
rz(-1.410781) q[1];
sx q[1];
rz(-2.3115092) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0706129) q[0];
sx q[0];
rz(-1.5375397) q[0];
sx q[0];
rz(-0.279056) q[0];
rz(-pi) q[1];
rz(1.7858429) q[2];
sx q[2];
rz(-2.1587662) q[2];
sx q[2];
rz(-2.7933592) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.54129564) q[1];
sx q[1];
rz(-1.7233371) q[1];
sx q[1];
rz(2.9593854) q[1];
x q[2];
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
rz(2.3229708) q[2];
sx q[2];
rz(-1.8399723) q[2];
sx q[2];
rz(-1.0969561) q[2];
rz(-1.8638301) q[3];
sx q[3];
rz(-2.4570229) q[3];
sx q[3];
rz(-0.6440312) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85035664) q[0];
sx q[0];
rz(-2.0391897) q[0];
sx q[0];
rz(-2.7665603) q[0];
rz(-3.0229783) q[1];
sx q[1];
rz(-1.4362486) q[1];
sx q[1];
rz(0.92432252) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9572231) q[0];
sx q[0];
rz(-0.65378377) q[0];
sx q[0];
rz(2.6517646) q[0];
x q[1];
rz(-2.0807939) q[2];
sx q[2];
rz(-2.1782403) q[2];
sx q[2];
rz(-0.80112544) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.1076775) q[1];
sx q[1];
rz(-0.36201358) q[1];
sx q[1];
rz(-0.081078366) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.499324) q[3];
sx q[3];
rz(-0.54107252) q[3];
sx q[3];
rz(-0.16995811) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.88860005) q[2];
sx q[2];
rz(-1.000095) q[2];
sx q[2];
rz(-2.2484153) q[2];
rz(-2.0231953) q[3];
sx q[3];
rz(-1.018254) q[3];
sx q[3];
rz(2.4848188) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-1.7669582) q[0];
sx q[0];
rz(-2.0603016) q[0];
sx q[0];
rz(1.7745071) q[0];
rz(0.18192667) q[1];
sx q[1];
rz(-0.55768273) q[1];
sx q[1];
rz(2.2348977) q[1];
rz(0.091808783) q[2];
sx q[2];
rz(-1.8845857) q[2];
sx q[2];
rz(1.2086445) q[2];
rz(-1.160865) q[3];
sx q[3];
rz(-1.0900453) q[3];
sx q[3];
rz(2.446519) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
