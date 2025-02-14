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
rz(1.3069557) q[0];
sx q[0];
rz(4.0210273) q[0];
sx q[0];
rz(9.9183912) q[0];
rz(-1.2797132) q[1];
sx q[1];
rz(3.7682025) q[1];
sx q[1];
rz(11.087853) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0644306) q[0];
sx q[0];
rz(-2.0666984) q[0];
sx q[0];
rz(2.2999963) q[0];
rz(-pi) q[1];
rz(-0.25716146) q[2];
sx q[2];
rz(-0.92338054) q[2];
sx q[2];
rz(-2.0055298) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.6329816) q[1];
sx q[1];
rz(-2.2618083) q[1];
sx q[1];
rz(1.2862842) q[1];
rz(-pi) q[2];
rz(0.36840393) q[3];
sx q[3];
rz(-2.1737636) q[3];
sx q[3];
rz(-1.6507698) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.0066068) q[2];
sx q[2];
rz(-2.0015494) q[2];
sx q[2];
rz(2.0331649) q[2];
rz(1.3125575) q[3];
sx q[3];
rz(-2.8969942) q[3];
sx q[3];
rz(-0.31497064) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7600064) q[0];
sx q[0];
rz(-0.8257603) q[0];
sx q[0];
rz(-3.1257358) q[0];
rz(-0.99041692) q[1];
sx q[1];
rz(-0.76737338) q[1];
sx q[1];
rz(-2.2784065) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8548633) q[0];
sx q[0];
rz(-0.83896381) q[0];
sx q[0];
rz(-1.4409723) q[0];
x q[1];
rz(0.77570373) q[2];
sx q[2];
rz(-0.45368567) q[2];
sx q[2];
rz(-1.1248447) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.9941123) q[1];
sx q[1];
rz(-1.0801472) q[1];
sx q[1];
rz(-1.7685686) q[1];
rz(1.7402418) q[3];
sx q[3];
rz(-0.72899216) q[3];
sx q[3];
rz(0.21641009) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.70011675) q[2];
sx q[2];
rz(-2.8642544) q[2];
sx q[2];
rz(2.6727943) q[2];
rz(0.68764728) q[3];
sx q[3];
rz(-1.629849) q[3];
sx q[3];
rz(-2.8414753) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1402635) q[0];
sx q[0];
rz(-2.2183473) q[0];
sx q[0];
rz(0.73626751) q[0];
rz(-2.7983792) q[1];
sx q[1];
rz(-2.2540269) q[1];
sx q[1];
rz(0.18377486) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.60584468) q[0];
sx q[0];
rz(-2.3130406) q[0];
sx q[0];
rz(1.4808307) q[0];
rz(-0.69664069) q[2];
sx q[2];
rz(-1.9844779) q[2];
sx q[2];
rz(0.51332849) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.0597093) q[1];
sx q[1];
rz(-2.3569111) q[1];
sx q[1];
rz(-2.7126524) q[1];
rz(0.36005693) q[3];
sx q[3];
rz(-2.0814133) q[3];
sx q[3];
rz(1.0086446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.78202128) q[2];
sx q[2];
rz(-2.6247793) q[2];
sx q[2];
rz(2.6256631) q[2];
rz(-1.1082209) q[3];
sx q[3];
rz(-2.5907232) q[3];
sx q[3];
rz(1.9903323) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16294031) q[0];
sx q[0];
rz(-1.9300224) q[0];
sx q[0];
rz(-0.51280713) q[0];
rz(2.6817952) q[1];
sx q[1];
rz(-1.5492946) q[1];
sx q[1];
rz(2.3777681) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99256809) q[0];
sx q[0];
rz(-1.6852753) q[0];
sx q[0];
rz(-3.0106198) q[0];
rz(-pi) q[1];
rz(-0.41641367) q[2];
sx q[2];
rz(-1.8685307) q[2];
sx q[2];
rz(-2.1418051) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.3556174) q[1];
sx q[1];
rz(-1.2316868) q[1];
sx q[1];
rz(0.29735844) q[1];
x q[2];
rz(-0.34448907) q[3];
sx q[3];
rz(-1.7055344) q[3];
sx q[3];
rz(1.8014993) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.1226471) q[2];
sx q[2];
rz(-2.985869) q[2];
sx q[2];
rz(-2.8141008) q[2];
rz(-2.859595) q[3];
sx q[3];
rz(-1.3982541) q[3];
sx q[3];
rz(1.5871083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9649488) q[0];
sx q[0];
rz(-2.7556941) q[0];
sx q[0];
rz(-2.1642245) q[0];
rz(-2.7727959) q[1];
sx q[1];
rz(-0.86124033) q[1];
sx q[1];
rz(-1.4126973) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57713014) q[0];
sx q[0];
rz(-2.6106195) q[0];
sx q[0];
rz(-1.1891548) q[0];
rz(-3.1277282) q[2];
sx q[2];
rz(-2.0875476) q[2];
sx q[2];
rz(2.2114829) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1407847) q[1];
sx q[1];
rz(-2.85726) q[1];
sx q[1];
rz(-1.6823014) q[1];
rz(-1.2431954) q[3];
sx q[3];
rz(-1.0169815) q[3];
sx q[3];
rz(-2.6462951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.29350027) q[2];
sx q[2];
rz(-0.81587452) q[2];
sx q[2];
rz(0.14981848) q[2];
rz(0.39854974) q[3];
sx q[3];
rz(-1.5236676) q[3];
sx q[3];
rz(0.15628763) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77001101) q[0];
sx q[0];
rz(-3.0122029) q[0];
sx q[0];
rz(2.5883801) q[0];
rz(-2.5999056) q[1];
sx q[1];
rz(-2.0952416) q[1];
sx q[1];
rz(-2.6936626) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49419379) q[0];
sx q[0];
rz(-1.1851536) q[0];
sx q[0];
rz(-1.4146039) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7121353) q[2];
sx q[2];
rz(-1.9782559) q[2];
sx q[2];
rz(0.40312672) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.74400157) q[1];
sx q[1];
rz(-2.7251456) q[1];
sx q[1];
rz(-1.2214425) q[1];
rz(-pi) q[2];
rz(1.5138835) q[3];
sx q[3];
rz(-1.2384733) q[3];
sx q[3];
rz(0.11385465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.64122671) q[2];
sx q[2];
rz(-1.1643103) q[2];
sx q[2];
rz(-0.79916239) q[2];
rz(0.3847807) q[3];
sx q[3];
rz(-0.32618263) q[3];
sx q[3];
rz(1.0001812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7810829) q[0];
sx q[0];
rz(-1.0436844) q[0];
sx q[0];
rz(-0.59180301) q[0];
rz(2.7599755) q[1];
sx q[1];
rz(-0.38002574) q[1];
sx q[1];
rz(2.9891678) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29983175) q[0];
sx q[0];
rz(-0.61474568) q[0];
sx q[0];
rz(-2.1257945) q[0];
x q[1];
rz(-0.47758684) q[2];
sx q[2];
rz(-1.7268333) q[2];
sx q[2];
rz(-2.0618771) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.85411863) q[1];
sx q[1];
rz(-1.2240613) q[1];
sx q[1];
rz(1.5902014) q[1];
rz(-pi) q[2];
x q[2];
rz(0.82181886) q[3];
sx q[3];
rz(-0.83335175) q[3];
sx q[3];
rz(2.5488473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.3064208) q[2];
sx q[2];
rz(-0.99205899) q[2];
sx q[2];
rz(-1.4527808) q[2];
rz(-2.6460904) q[3];
sx q[3];
rz(-0.82363868) q[3];
sx q[3];
rz(-0.56408322) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38967663) q[0];
sx q[0];
rz(-0.037411995) q[0];
sx q[0];
rz(-2.9801242) q[0];
rz(3.1096733) q[1];
sx q[1];
rz(-2.4999764) q[1];
sx q[1];
rz(1.8643103) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91457462) q[0];
sx q[0];
rz(-0.76138568) q[0];
sx q[0];
rz(-2.7622591) q[0];
rz(-1.0890245) q[2];
sx q[2];
rz(-2.6725997) q[2];
sx q[2];
rz(0.8587786) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.1254221) q[1];
sx q[1];
rz(-1.6787663) q[1];
sx q[1];
rz(2.8362464) q[1];
x q[2];
rz(1.444449) q[3];
sx q[3];
rz(-0.29623162) q[3];
sx q[3];
rz(2.031523) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.6892467) q[2];
sx q[2];
rz(-2.2301058) q[2];
sx q[2];
rz(2.7455043) q[2];
rz(-2.7416157) q[3];
sx q[3];
rz(-0.59498274) q[3];
sx q[3];
rz(0.8766492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9006627) q[0];
sx q[0];
rz(-0.018095896) q[0];
sx q[0];
rz(-2.990429) q[0];
rz(-0.97686544) q[1];
sx q[1];
rz(-2.7604389) q[1];
sx q[1];
rz(2.3968598) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3522986) q[0];
sx q[0];
rz(-1.8478824) q[0];
sx q[0];
rz(-0.21488551) q[0];
x q[1];
rz(1.9836203) q[2];
sx q[2];
rz(-1.7016439) q[2];
sx q[2];
rz(1.7826338) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.7661644) q[1];
sx q[1];
rz(-0.64830983) q[1];
sx q[1];
rz(0.0023771087) q[1];
rz(-pi) q[2];
rz(1.8119393) q[3];
sx q[3];
rz(-1.2281115) q[3];
sx q[3];
rz(1.8543275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.62290827) q[2];
sx q[2];
rz(-1.9025977) q[2];
sx q[2];
rz(-2.1935479) q[2];
rz(-2.0840123) q[3];
sx q[3];
rz(-1.6653929) q[3];
sx q[3];
rz(1.074056) q[3];
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
x q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0483765) q[0];
sx q[0];
rz(-0.26234782) q[0];
sx q[0];
rz(-0.23396215) q[0];
rz(-1.9562862) q[1];
sx q[1];
rz(-2.2133841) q[1];
sx q[1];
rz(-1.8918461) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.043086655) q[0];
sx q[0];
rz(-1.0191425) q[0];
sx q[0];
rz(-1.2292525) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7777257) q[2];
sx q[2];
rz(-0.96038681) q[2];
sx q[2];
rz(2.0774942) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.1703012) q[1];
sx q[1];
rz(-2.2066322) q[1];
sx q[1];
rz(2.2595992) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3951282) q[3];
sx q[3];
rz(-2.2445956) q[3];
sx q[3];
rz(1.2487981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.11448161) q[2];
sx q[2];
rz(-1.9805084) q[2];
sx q[2];
rz(1.4410045) q[2];
rz(3.0103185) q[3];
sx q[3];
rz(-0.461853) q[3];
sx q[3];
rz(-2.89768) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.092939) q[0];
sx q[0];
rz(-0.96554148) q[0];
sx q[0];
rz(1.1337793) q[0];
rz(-1.0406021) q[1];
sx q[1];
rz(-2.2941209) q[1];
sx q[1];
rz(2.6170731) q[1];
rz(2.8362989) q[2];
sx q[2];
rz(-0.96571453) q[2];
sx q[2];
rz(-2.1681186) q[2];
rz(-0.63325044) q[3];
sx q[3];
rz(-1.4489104) q[3];
sx q[3];
rz(-0.62319402) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
