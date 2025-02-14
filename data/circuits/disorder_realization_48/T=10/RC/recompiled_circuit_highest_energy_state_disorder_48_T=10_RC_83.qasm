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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9636757) q[0];
sx q[0];
rz(-0.48030329) q[0];
sx q[0];
rz(-0.98573523) q[0];
rz(-pi) q[1];
x q[1];
rz(0.79390766) q[2];
sx q[2];
rz(-0.61015266) q[2];
sx q[2];
rz(0.33294233) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(3.0585441) q[1];
sx q[1];
rz(-1.3131318) q[1];
sx q[1];
rz(-0.017770692) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3506557) q[3];
sx q[3];
rz(-1.992989) q[3];
sx q[3];
rz(3.0289502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-3.1276663) q[2];
sx q[2];
rz(-1.650859) q[2];
sx q[2];
rz(2.9909383) q[2];
rz(-1.6023747) q[3];
sx q[3];
rz(-2.7645002) q[3];
sx q[3];
rz(-0.25130513) q[3];
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
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0852614) q[0];
sx q[0];
rz(-1.769861) q[0];
sx q[0];
rz(2.7714609) q[0];
rz(2.6689957) q[1];
sx q[1];
rz(-0.31817803) q[1];
sx q[1];
rz(-2.1841689) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32015285) q[0];
sx q[0];
rz(-1.4906881) q[0];
sx q[0];
rz(-2.0700685) q[0];
rz(-pi) q[1];
rz(0.17458169) q[2];
sx q[2];
rz(-2.0160297) q[2];
sx q[2];
rz(-0.11694678) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.4366001) q[1];
sx q[1];
rz(-0.7448405) q[1];
sx q[1];
rz(-0.1257433) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4647275) q[3];
sx q[3];
rz(-1.3086623) q[3];
sx q[3];
rz(-1.8755975) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(3.1364253) q[2];
sx q[2];
rz(-2.397126) q[2];
sx q[2];
rz(1.0235419) q[2];
rz(-1.7510022) q[3];
sx q[3];
rz(-1.3955045) q[3];
sx q[3];
rz(1.3172147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9978299) q[0];
sx q[0];
rz(-1.0799066) q[0];
sx q[0];
rz(-1.9299555) q[0];
rz(-pi) q[1];
rz(3.0653238) q[2];
sx q[2];
rz(-0.4919211) q[2];
sx q[2];
rz(0.15215506) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.9197977) q[1];
sx q[1];
rz(-0.88558965) q[1];
sx q[1];
rz(-2.5849839) q[1];
x q[2];
rz(0.7003673) q[3];
sx q[3];
rz(-2.1383965) q[3];
sx q[3];
rz(-1.2705453) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2341653) q[2];
sx q[2];
rz(-1.5831466) q[2];
sx q[2];
rz(-1.3321715) q[2];
rz(2.7857156) q[3];
sx q[3];
rz(-1.9477113) q[3];
sx q[3];
rz(0.98085421) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(2.6818162) q[0];
sx q[0];
rz(-1.1695319) q[0];
sx q[0];
rz(-1.1112777) q[0];
rz(2.2214644) q[1];
sx q[1];
rz(-1.5891113) q[1];
sx q[1];
rz(-2.8880602) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7074575) q[0];
sx q[0];
rz(-0.61096707) q[0];
sx q[0];
rz(-0.57662983) q[0];
rz(-pi) q[1];
rz(-1.6854671) q[2];
sx q[2];
rz(-1.4114228) q[2];
sx q[2];
rz(-1.3068975) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.82099709) q[1];
sx q[1];
rz(-1.7554469) q[1];
sx q[1];
rz(0.76652938) q[1];
x q[2];
rz(0.30999581) q[3];
sx q[3];
rz(-1.0862203) q[3];
sx q[3];
rz(-2.692401) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.8186875) q[2];
sx q[2];
rz(-0.56537586) q[2];
sx q[2];
rz(1.6820071) q[2];
rz(0.5736351) q[3];
sx q[3];
rz(-1.2021659) q[3];
sx q[3];
rz(0.044376317) q[3];
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
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6673073) q[0];
sx q[0];
rz(-1.0695589) q[0];
sx q[0];
rz(1.3473508) q[0];
rz(-2.5509293) q[1];
sx q[1];
rz(-1.2202411) q[1];
sx q[1];
rz(1.4264872) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1388228) q[0];
sx q[0];
rz(-1.9277687) q[0];
sx q[0];
rz(0.41923747) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1828184) q[2];
sx q[2];
rz(-2.175594) q[2];
sx q[2];
rz(-1.6590349) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.8914681) q[1];
sx q[1];
rz(-1.3169177) q[1];
sx q[1];
rz(-1.5230383) q[1];
rz(2.9699202) q[3];
sx q[3];
rz(-0.35997691) q[3];
sx q[3];
rz(-1.5300919) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.9023989) q[2];
sx q[2];
rz(-1.0603443) q[2];
sx q[2];
rz(-1.1725461) q[2];
rz(-0.086056195) q[3];
sx q[3];
rz(-0.9525758) q[3];
sx q[3];
rz(2.9854767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0657144) q[0];
sx q[0];
rz(-1.5471764) q[0];
sx q[0];
rz(-0.86268798) q[0];
rz(0.24738303) q[1];
sx q[1];
rz(-1.3037953) q[1];
sx q[1];
rz(-0.47201306) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9111689) q[0];
sx q[0];
rz(-0.8967451) q[0];
sx q[0];
rz(-2.0087025) q[0];
x q[1];
rz(1.0649933) q[2];
sx q[2];
rz(-1.7031295) q[2];
sx q[2];
rz(1.1810609) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.6093543) q[1];
sx q[1];
rz(-1.2684457) q[1];
sx q[1];
rz(1.7560033) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2788762) q[3];
sx q[3];
rz(-0.64281598) q[3];
sx q[3];
rz(1.1228648) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.2315959) q[2];
sx q[2];
rz(-2.1746077) q[2];
sx q[2];
rz(-1.7426573) q[2];
rz(1.4553962) q[3];
sx q[3];
rz(-0.78794909) q[3];
sx q[3];
rz(0.56070352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6550734) q[0];
sx q[0];
rz(-1.1460679) q[0];
sx q[0];
rz(-2.6229677) q[0];
rz(2.4229166) q[1];
sx q[1];
rz(-2.6287754) q[1];
sx q[1];
rz(1.8124883) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9921917) q[0];
sx q[0];
rz(-1.2103545) q[0];
sx q[0];
rz(-2.0130231) q[0];
x q[1];
rz(-1.3955922) q[2];
sx q[2];
rz(-0.28131286) q[2];
sx q[2];
rz(3.0228721) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.3997201) q[1];
sx q[1];
rz(-2.406027) q[1];
sx q[1];
rz(-1.5217785) q[1];
x q[2];
rz(0.98436004) q[3];
sx q[3];
rz(-1.5542665) q[3];
sx q[3];
rz(0.23197385) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.79080498) q[2];
sx q[2];
rz(-0.69956508) q[2];
sx q[2];
rz(-0.038185509) q[2];
rz(0.83032483) q[3];
sx q[3];
rz(-1.3874976) q[3];
sx q[3];
rz(3.0939046) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6787978) q[0];
sx q[0];
rz(-2.9252453) q[0];
sx q[0];
rz(1.3702673) q[0];
rz(2.7554152) q[1];
sx q[1];
rz(-1.4048856) q[1];
sx q[1];
rz(-0.6699627) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7012335) q[0];
sx q[0];
rz(-2.0645077) q[0];
sx q[0];
rz(-0.53892737) q[0];
x q[1];
rz(-3.0490446) q[2];
sx q[2];
rz(-2.4227648) q[2];
sx q[2];
rz(0.75564849) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.252287) q[1];
sx q[1];
rz(-0.60599594) q[1];
sx q[1];
rz(1.3713981) q[1];
rz(-pi) q[2];
rz(-0.037360351) q[3];
sx q[3];
rz(-1.7826163) q[3];
sx q[3];
rz(-0.39549124) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.2204444) q[2];
sx q[2];
rz(-2.8577652) q[2];
sx q[2];
rz(1.05668) q[2];
rz(-1.7250666) q[3];
sx q[3];
rz(-1.2456015) q[3];
sx q[3];
rz(-1.8511124) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0678299) q[0];
sx q[0];
rz(-2.1011598) q[0];
sx q[0];
rz(-2.1347866) q[0];
rz(2.6489068) q[1];
sx q[1];
rz(-1.7308116) q[1];
sx q[1];
rz(2.3115092) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3842889) q[0];
sx q[0];
rz(-2.8606133) q[0];
sx q[0];
rz(-3.021394) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7858429) q[2];
sx q[2];
rz(-2.1587662) q[2];
sx q[2];
rz(-0.34823349) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.54129564) q[1];
sx q[1];
rz(-1.7233371) q[1];
sx q[1];
rz(0.18220724) q[1];
x q[2];
rz(2.0138143) q[3];
sx q[3];
rz(-1.5592056) q[3];
sx q[3];
rz(3.1168337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.3229708) q[2];
sx q[2];
rz(-1.3016204) q[2];
sx q[2];
rz(-2.0446365) q[2];
rz(1.8638301) q[3];
sx q[3];
rz(-2.4570229) q[3];
sx q[3];
rz(-2.4975615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85035664) q[0];
sx q[0];
rz(-1.1024029) q[0];
sx q[0];
rz(-2.7665603) q[0];
rz(-3.0229783) q[1];
sx q[1];
rz(-1.4362486) q[1];
sx q[1];
rz(-2.2172701) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7345006) q[0];
sx q[0];
rz(-1.0043) q[0];
sx q[0];
rz(1.9167711) q[0];
rz(-0.61225981) q[2];
sx q[2];
rz(-2.3697402) q[2];
sx q[2];
rz(1.5764232) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.1076775) q[1];
sx q[1];
rz(-0.36201358) q[1];
sx q[1];
rz(-0.081078366) q[1];
rz(2.499324) q[3];
sx q[3];
rz(-0.54107252) q[3];
sx q[3];
rz(-2.9716345) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.2529926) q[2];
sx q[2];
rz(-2.1414976) q[2];
sx q[2];
rz(0.89317733) q[2];
rz(-1.1183974) q[3];
sx q[3];
rz(-2.1233386) q[3];
sx q[3];
rz(-0.65677381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
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
rz(-3.0497839) q[2];
sx q[2];
rz(-1.8845857) q[2];
sx q[2];
rz(1.2086445) q[2];
rz(-0.51707324) q[3];
sx q[3];
rz(-1.2096249) q[3];
sx q[3];
rz(0.67740868) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
