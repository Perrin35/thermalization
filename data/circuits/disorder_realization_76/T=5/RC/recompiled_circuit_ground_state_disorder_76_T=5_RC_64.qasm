OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.4467093) q[0];
sx q[0];
rz(-2.0599685) q[0];
sx q[0];
rz(-2.4021436) q[0];
rz(2.3820355) q[1];
sx q[1];
rz(-1.3146725) q[1];
sx q[1];
rz(-1.7159599) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6524871) q[0];
sx q[0];
rz(-1.8579626) q[0];
sx q[0];
rz(-2.1225568) q[0];
x q[1];
rz(-0.47503586) q[2];
sx q[2];
rz(-0.54752195) q[2];
sx q[2];
rz(0.56625596) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.718218) q[1];
sx q[1];
rz(-1.6483232) q[1];
sx q[1];
rz(-1.7192057) q[1];
rz(-pi) q[2];
rz(-0.043685901) q[3];
sx q[3];
rz(-0.85245347) q[3];
sx q[3];
rz(1.9909315) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.2322959) q[2];
sx q[2];
rz(-0.84315073) q[2];
sx q[2];
rz(-0.24093957) q[2];
rz(-0.029189261) q[3];
sx q[3];
rz(-1.8029282) q[3];
sx q[3];
rz(-1.1791112) q[3];
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
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1643243) q[0];
sx q[0];
rz(-1.9376396) q[0];
sx q[0];
rz(-2.6569195) q[0];
rz(-2.4618705) q[1];
sx q[1];
rz(-1.8599963) q[1];
sx q[1];
rz(1.1601123) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3697333) q[0];
sx q[0];
rz(-1.5892649) q[0];
sx q[0];
rz(0.30402495) q[0];
rz(-pi) q[1];
rz(2.3789139) q[2];
sx q[2];
rz(-0.87717036) q[2];
sx q[2];
rz(0.1887624) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.94494941) q[1];
sx q[1];
rz(-1.3232051) q[1];
sx q[1];
rz(0.61243122) q[1];
rz(2.6187702) q[3];
sx q[3];
rz(-1.8588052) q[3];
sx q[3];
rz(3.0436181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.29740563) q[2];
sx q[2];
rz(-2.83941) q[2];
sx q[2];
rz(-2.6837132) q[2];
rz(2.0186021) q[3];
sx q[3];
rz(-2.1419958) q[3];
sx q[3];
rz(0.56639731) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26894012) q[0];
sx q[0];
rz(-2.432423) q[0];
sx q[0];
rz(-0.031524468) q[0];
rz(0.28745502) q[1];
sx q[1];
rz(-0.87702409) q[1];
sx q[1];
rz(-1.215975) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4327964) q[0];
sx q[0];
rz(-2.2026718) q[0];
sx q[0];
rz(-2.2461476) q[0];
x q[1];
rz(2.6570733) q[2];
sx q[2];
rz(-1.3860102) q[2];
sx q[2];
rz(-0.69714025) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.029473631) q[1];
sx q[1];
rz(-0.39327565) q[1];
sx q[1];
rz(2.359819) q[1];
x q[2];
rz(-0.76283703) q[3];
sx q[3];
rz(-0.82218542) q[3];
sx q[3];
rz(-1.959182) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-3.1027801) q[2];
sx q[2];
rz(-1.5422042) q[2];
sx q[2];
rz(-0.83596027) q[2];
rz(-1.3577667) q[3];
sx q[3];
rz(-2.416555) q[3];
sx q[3];
rz(0.74762216) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2826071) q[0];
sx q[0];
rz(-1.405412) q[0];
sx q[0];
rz(3.0614241) q[0];
rz(2.0591586) q[1];
sx q[1];
rz(-2.9083462) q[1];
sx q[1];
rz(1.3287883) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69705582) q[0];
sx q[0];
rz(-1.2460684) q[0];
sx q[0];
rz(2.6022807) q[0];
rz(-pi) q[1];
rz(-1.3983509) q[2];
sx q[2];
rz(-1.2524464) q[2];
sx q[2];
rz(-2.6686252) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.7419974) q[1];
sx q[1];
rz(-0.4542225) q[1];
sx q[1];
rz(1.9059559) q[1];
x q[2];
rz(1.7747709) q[3];
sx q[3];
rz(-1.2414724) q[3];
sx q[3];
rz(-2.9993111) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.3749915) q[2];
sx q[2];
rz(-2.1583755) q[2];
sx q[2];
rz(-0.90744606) q[2];
rz(-0.67029101) q[3];
sx q[3];
rz(-0.82796103) q[3];
sx q[3];
rz(1.7214187) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9012673) q[0];
sx q[0];
rz(-2.2690052) q[0];
sx q[0];
rz(2.7101044) q[0];
rz(1.1314499) q[1];
sx q[1];
rz(-1.1791469) q[1];
sx q[1];
rz(2.7010837) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1898855) q[0];
sx q[0];
rz(-1.110297) q[0];
sx q[0];
rz(0.13114625) q[0];
rz(2.9350403) q[2];
sx q[2];
rz(-1.9777386) q[2];
sx q[2];
rz(2.664704) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7311598) q[1];
sx q[1];
rz(-1.9010547) q[1];
sx q[1];
rz(-1.5291269) q[1];
rz(-2.2045361) q[3];
sx q[3];
rz(-2.0645294) q[3];
sx q[3];
rz(-0.99668324) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.12209192) q[2];
sx q[2];
rz(-2.2152405) q[2];
sx q[2];
rz(2.7276373) q[2];
rz(1.3881989) q[3];
sx q[3];
rz(-1.5475169) q[3];
sx q[3];
rz(3.0164914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51467657) q[0];
sx q[0];
rz(-1.7358945) q[0];
sx q[0];
rz(0.93389121) q[0];
rz(-0.67032188) q[1];
sx q[1];
rz(-1.2443845) q[1];
sx q[1];
rz(1.3353039) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0224755) q[0];
sx q[0];
rz(-0.60998218) q[0];
sx q[0];
rz(-1.7458292) q[0];
rz(2.983903) q[2];
sx q[2];
rz(-1.5957498) q[2];
sx q[2];
rz(1.4098997) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.090744734) q[1];
sx q[1];
rz(-3.0266592) q[1];
sx q[1];
rz(-1.7687709) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5794831) q[3];
sx q[3];
rz(-0.93378769) q[3];
sx q[3];
rz(-0.52906936) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.2495217) q[2];
sx q[2];
rz(-0.27270174) q[2];
sx q[2];
rz(1.2707233) q[2];
rz(-1.3696085) q[3];
sx q[3];
rz(-1.9000051) q[3];
sx q[3];
rz(2.6197267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(0.18913604) q[0];
sx q[0];
rz(-2.554775) q[0];
sx q[0];
rz(0.15368803) q[0];
rz(2.9391089) q[1];
sx q[1];
rz(-1.9053562) q[1];
sx q[1];
rz(-0.94857803) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21009203) q[0];
sx q[0];
rz(-2.4903959) q[0];
sx q[0];
rz(-0.33588846) q[0];
rz(-pi) q[1];
rz(-2.2099451) q[2];
sx q[2];
rz(-1.8338565) q[2];
sx q[2];
rz(-1.5342889) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.13742971) q[1];
sx q[1];
rz(-1.8297126) q[1];
sx q[1];
rz(0.24000411) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0074869) q[3];
sx q[3];
rz(-0.92536345) q[3];
sx q[3];
rz(-0.2303309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.004868) q[2];
sx q[2];
rz(-2.0437045) q[2];
sx q[2];
rz(0.0014121545) q[2];
rz(-0.062156113) q[3];
sx q[3];
rz(-1.4096189) q[3];
sx q[3];
rz(2.1803161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75981265) q[0];
sx q[0];
rz(-0.75730046) q[0];
sx q[0];
rz(1.9267474) q[0];
rz(-0.75792056) q[1];
sx q[1];
rz(-0.52752033) q[1];
sx q[1];
rz(2.5835999) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4014977) q[0];
sx q[0];
rz(-1.3245653) q[0];
sx q[0];
rz(-1.2241062) q[0];
x q[1];
rz(-0.55155374) q[2];
sx q[2];
rz(-1.6948912) q[2];
sx q[2];
rz(0.53260224) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.8415738) q[1];
sx q[1];
rz(-1.8487367) q[1];
sx q[1];
rz(-1.0701847) q[1];
x q[2];
rz(3.1147546) q[3];
sx q[3];
rz(-1.8478113) q[3];
sx q[3];
rz(-2.7153496) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.5672292) q[2];
sx q[2];
rz(-1.9476674) q[2];
sx q[2];
rz(-0.26926678) q[2];
rz(1.1632129) q[3];
sx q[3];
rz(-2.5389157) q[3];
sx q[3];
rz(-2.9009624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6729386) q[0];
sx q[0];
rz(-1.0373632) q[0];
sx q[0];
rz(1.0733676) q[0];
rz(-2.0138373) q[1];
sx q[1];
rz(-2.914371) q[1];
sx q[1];
rz(2.5849297) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7391602) q[0];
sx q[0];
rz(-0.35050979) q[0];
sx q[0];
rz(1.1465766) q[0];
x q[1];
rz(0.9999335) q[2];
sx q[2];
rz(-2.0682671) q[2];
sx q[2];
rz(-1.2116878) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.4924188) q[1];
sx q[1];
rz(-2.2606876) q[1];
sx q[1];
rz(2.0418732) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4707556) q[3];
sx q[3];
rz(-1.3373858) q[3];
sx q[3];
rz(-1.6897214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.90159455) q[2];
sx q[2];
rz(-1.0866714) q[2];
sx q[2];
rz(1.1617917) q[2];
rz(-2.6084172) q[3];
sx q[3];
rz(-0.52459255) q[3];
sx q[3];
rz(0.1575135) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6530782) q[0];
sx q[0];
rz(-2.1670659) q[0];
sx q[0];
rz(2.5467806) q[0];
rz(2.5626903) q[1];
sx q[1];
rz(-1.4756823) q[1];
sx q[1];
rz(2.4822809) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4730277) q[0];
sx q[0];
rz(-0.94946948) q[0];
sx q[0];
rz(-1.7354119) q[0];
rz(2.012678) q[2];
sx q[2];
rz(-2.127248) q[2];
sx q[2];
rz(-2.93612) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.54592268) q[1];
sx q[1];
rz(-2.8575142) q[1];
sx q[1];
rz(1.0417263) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3385184) q[3];
sx q[3];
rz(-0.37130203) q[3];
sx q[3];
rz(-2.7367531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.79002964) q[2];
sx q[2];
rz(-1.1976676) q[2];
sx q[2];
rz(3.0885546) q[2];
rz(1.3430345) q[3];
sx q[3];
rz(-1.8648632) q[3];
sx q[3];
rz(-0.0742577) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2863083) q[0];
sx q[0];
rz(-1.7779779) q[0];
sx q[0];
rz(1.5386982) q[0];
rz(-0.098943624) q[1];
sx q[1];
rz(-1.7715441) q[1];
sx q[1];
rz(-3.0861707) q[1];
rz(0.020053645) q[2];
sx q[2];
rz(-0.40553025) q[2];
sx q[2];
rz(2.5237906) q[2];
rz(0.74138882) q[3];
sx q[3];
rz(-1.6869873) q[3];
sx q[3];
rz(1.6968873) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
