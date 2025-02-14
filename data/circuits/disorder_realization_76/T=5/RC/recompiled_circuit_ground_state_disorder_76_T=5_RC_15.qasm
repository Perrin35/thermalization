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
rz(-0.75955716) q[1];
sx q[1];
rz(-1.8269202) q[1];
sx q[1];
rz(-1.4256328) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6524871) q[0];
sx q[0];
rz(-1.28363) q[0];
sx q[0];
rz(1.0190359) q[0];
rz(-pi) q[1];
rz(-0.47503586) q[2];
sx q[2];
rz(-0.54752195) q[2];
sx q[2];
rz(0.56625596) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5162075) q[1];
sx q[1];
rz(-2.974286) q[1];
sx q[1];
rz(-2.0545261) q[1];
x q[2];
rz(0.043685901) q[3];
sx q[3];
rz(-2.2891392) q[3];
sx q[3];
rz(1.9909315) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.90929675) q[2];
sx q[2];
rz(-2.2984419) q[2];
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
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1643243) q[0];
sx q[0];
rz(-1.9376396) q[0];
sx q[0];
rz(-2.6569195) q[0];
rz(0.67972216) q[1];
sx q[1];
rz(-1.8599963) q[1];
sx q[1];
rz(1.1601123) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80473111) q[0];
sx q[0];
rz(-1.2668249) q[0];
sx q[0];
rz(-1.5514403) q[0];
rz(-pi) q[1];
x q[1];
rz(0.76267879) q[2];
sx q[2];
rz(-0.87717036) q[2];
sx q[2];
rz(-0.1887624) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.1966432) q[1];
sx q[1];
rz(-1.8183876) q[1];
sx q[1];
rz(2.5291614) q[1];
x q[2];
rz(-2.6187702) q[3];
sx q[3];
rz(-1.2827875) q[3];
sx q[3];
rz(3.0436181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.29740563) q[2];
sx q[2];
rz(-0.30218267) q[2];
sx q[2];
rz(2.6837132) q[2];
rz(-1.1229905) q[3];
sx q[3];
rz(-0.99959683) q[3];
sx q[3];
rz(2.5751953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8726525) q[0];
sx q[0];
rz(-0.70916969) q[0];
sx q[0];
rz(-3.1100682) q[0];
rz(0.28745502) q[1];
sx q[1];
rz(-2.2645686) q[1];
sx q[1];
rz(-1.9256176) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7087962) q[0];
sx q[0];
rz(-2.2026718) q[0];
sx q[0];
rz(-2.2461476) q[0];
x q[1];
rz(-1.7789677) q[2];
sx q[2];
rz(-1.0952172) q[2];
sx q[2];
rz(2.1715234) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-3.112119) q[1];
sx q[1];
rz(-2.748317) q[1];
sx q[1];
rz(0.78177364) q[1];
rz(0.76283703) q[3];
sx q[3];
rz(-0.82218542) q[3];
sx q[3];
rz(1.959182) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.03881255) q[2];
sx q[2];
rz(-1.5422042) q[2];
sx q[2];
rz(0.83596027) q[2];
rz(-1.7838259) q[3];
sx q[3];
rz(-0.72503763) q[3];
sx q[3];
rz(-2.3939705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-1.8589856) q[0];
sx q[0];
rz(-1.7361807) q[0];
sx q[0];
rz(0.080168515) q[0];
rz(2.0591586) q[1];
sx q[1];
rz(-2.9083462) q[1];
sx q[1];
rz(1.3287883) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4565312) q[0];
sx q[0];
rz(-1.0624806) q[0];
sx q[0];
rz(1.1969181) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3983509) q[2];
sx q[2];
rz(-1.8891462) q[2];
sx q[2];
rz(-0.47296745) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.3722244) q[1];
sx q[1];
rz(-1.1435678) q[1];
sx q[1];
rz(-0.15924304) q[1];
rz(-pi) q[2];
x q[2];
rz(0.53502632) q[3];
sx q[3];
rz(-0.38540977) q[3];
sx q[3];
rz(0.42675323) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.7666011) q[2];
sx q[2];
rz(-2.1583755) q[2];
sx q[2];
rz(-0.90744606) q[2];
rz(0.67029101) q[3];
sx q[3];
rz(-2.3136316) q[3];
sx q[3];
rz(1.7214187) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2403253) q[0];
sx q[0];
rz(-2.2690052) q[0];
sx q[0];
rz(-2.7101044) q[0];
rz(-2.0101428) q[1];
sx q[1];
rz(-1.1791469) q[1];
sx q[1];
rz(-0.44050899) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2402396) q[0];
sx q[0];
rz(-0.47751891) q[0];
sx q[0];
rz(1.8285455) q[0];
rz(-1.155987) q[2];
sx q[2];
rz(-1.7602663) q[2];
sx q[2];
rz(-2.1304325) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.8590282) q[1];
sx q[1];
rz(-0.33278123) q[1];
sx q[1];
rz(-3.0206693) q[1];
rz(-pi) q[2];
rz(-2.308485) q[3];
sx q[3];
rz(-0.78189497) q[3];
sx q[3];
rz(1.9946919) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.0195007) q[2];
sx q[2];
rz(-2.2152405) q[2];
sx q[2];
rz(0.41395536) q[2];
rz(1.7533938) q[3];
sx q[3];
rz(-1.5940758) q[3];
sx q[3];
rz(3.0164914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51467657) q[0];
sx q[0];
rz(-1.4056982) q[0];
sx q[0];
rz(-2.2077014) q[0];
rz(-0.67032188) q[1];
sx q[1];
rz(-1.2443845) q[1];
sx q[1];
rz(-1.8062887) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7372197) q[0];
sx q[0];
rz(-1.4708733) q[0];
sx q[0];
rz(-2.1735682) q[0];
rz(-pi) q[1];
rz(2.9839758) q[2];
sx q[2];
rz(-0.15963563) q[2];
sx q[2];
rz(-3.1363413) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.090744734) q[1];
sx q[1];
rz(-0.11493348) q[1];
sx q[1];
rz(1.3728218) q[1];
rz(-1.5621095) q[3];
sx q[3];
rz(-2.207805) q[3];
sx q[3];
rz(0.52906936) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.2495217) q[2];
sx q[2];
rz(-0.27270174) q[2];
sx q[2];
rz(-1.2707233) q[2];
rz(1.3696085) q[3];
sx q[3];
rz(-1.9000051) q[3];
sx q[3];
rz(0.52186596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
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
rz(-0.18913604) q[0];
sx q[0];
rz(-2.554775) q[0];
sx q[0];
rz(-0.15368803) q[0];
rz(-0.20248374) q[1];
sx q[1];
rz(-1.9053562) q[1];
sx q[1];
rz(-0.94857803) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62371333) q[0];
sx q[0];
rz(-0.9615295) q[0];
sx q[0];
rz(-1.3246956) q[0];
x q[1];
rz(-0.3237299) q[2];
sx q[2];
rz(-2.1846131) q[2];
sx q[2];
rz(0.15440369) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.7708009) q[1];
sx q[1];
rz(-1.8026514) q[1];
sx q[1];
rz(1.3045909) q[1];
x q[2];
rz(-0.92102401) q[3];
sx q[3];
rz(-1.6778086) q[3];
sx q[3];
rz(1.2594852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.004868) q[2];
sx q[2];
rz(-2.0437045) q[2];
sx q[2];
rz(3.1401805) q[2];
rz(3.0794365) q[3];
sx q[3];
rz(-1.4096189) q[3];
sx q[3];
rz(-0.96127659) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.38178) q[0];
sx q[0];
rz(-2.3842922) q[0];
sx q[0];
rz(1.2148452) q[0];
rz(-0.75792056) q[1];
sx q[1];
rz(-0.52752033) q[1];
sx q[1];
rz(-0.55799276) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7400949) q[0];
sx q[0];
rz(-1.8170274) q[0];
sx q[0];
rz(1.9174865) q[0];
rz(-pi) q[1];
rz(-1.7162157) q[2];
sx q[2];
rz(-2.1176257) q[2];
sx q[2];
rz(1.1142004) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.3000189) q[1];
sx q[1];
rz(-1.8487367) q[1];
sx q[1];
rz(2.071408) q[1];
rz(-1.2936866) q[3];
sx q[3];
rz(-1.596611) q[3];
sx q[3];
rz(-1.1518948) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.5743635) q[2];
sx q[2];
rz(-1.1939253) q[2];
sx q[2];
rz(-2.8723259) q[2];
rz(-1.9783798) q[3];
sx q[3];
rz(-2.5389157) q[3];
sx q[3];
rz(-2.9009624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4686541) q[0];
sx q[0];
rz(-1.0373632) q[0];
sx q[0];
rz(-1.0733676) q[0];
rz(2.0138373) q[1];
sx q[1];
rz(-0.22722166) q[1];
sx q[1];
rz(2.5849297) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9087785) q[0];
sx q[0];
rz(-1.7126084) q[0];
sx q[0];
rz(1.2491666) q[0];
rz(-pi) q[1];
rz(2.3586541) q[2];
sx q[2];
rz(-2.4030444) q[2];
sx q[2];
rz(-0.99817456) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.60815114) q[1];
sx q[1];
rz(-1.2131629) q[1];
sx q[1];
rz(0.74700345) q[1];
x q[2];
rz(-1.4707556) q[3];
sx q[3];
rz(-1.8042068) q[3];
sx q[3];
rz(-1.6897214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.90159455) q[2];
sx q[2];
rz(-1.0866714) q[2];
sx q[2];
rz(-1.1617917) q[2];
rz(-2.6084172) q[3];
sx q[3];
rz(-0.52459255) q[3];
sx q[3];
rz(-2.9840792) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48851442) q[0];
sx q[0];
rz(-2.1670659) q[0];
sx q[0];
rz(-0.59481204) q[0];
rz(-2.5626903) q[1];
sx q[1];
rz(-1.4756823) q[1];
sx q[1];
rz(0.65931177) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1950451) q[0];
sx q[0];
rz(-2.5016245) q[0];
sx q[0];
rz(-2.9165687) q[0];
rz(-2.5393007) q[2];
sx q[2];
rz(-2.4458439) q[2];
sx q[2];
rz(2.2056923) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.54592268) q[1];
sx q[1];
rz(-0.28407846) q[1];
sx q[1];
rz(-1.0417263) q[1];
x q[2];
rz(1.8439383) q[3];
sx q[3];
rz(-1.3160664) q[3];
sx q[3];
rz(-1.8984853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.79002964) q[2];
sx q[2];
rz(-1.943925) q[2];
sx q[2];
rz(-0.053038049) q[2];
rz(1.3430345) q[3];
sx q[3];
rz(-1.2767295) q[3];
sx q[3];
rz(-3.067335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85528436) q[0];
sx q[0];
rz(-1.3636148) q[0];
sx q[0];
rz(-1.6028945) q[0];
rz(-0.098943624) q[1];
sx q[1];
rz(-1.7715441) q[1];
sx q[1];
rz(-3.0861707) q[1];
rz(3.121539) q[2];
sx q[2];
rz(-2.7360624) q[2];
sx q[2];
rz(-0.61780203) q[2];
rz(0.17114279) q[3];
sx q[3];
rz(-2.3928693) q[3];
sx q[3];
rz(0.00015043845) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
