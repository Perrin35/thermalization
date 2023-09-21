OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.52842927) q[0];
sx q[0];
rz(-1.0597205) q[0];
sx q[0];
rz(-2.4106195) q[0];
rz(1.641474) q[1];
sx q[1];
rz(5.2483622) q[1];
sx q[1];
rz(11.622826) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4817754) q[0];
sx q[0];
rz(-3.0660015) q[0];
sx q[0];
rz(-0.48150058) q[0];
rz(-pi) q[1];
x q[1];
rz(1.655683) q[2];
sx q[2];
rz(-0.92157084) q[2];
sx q[2];
rz(0.50253403) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.700625) q[1];
sx q[1];
rz(-1.96083) q[1];
sx q[1];
rz(0.36867152) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5991873) q[3];
sx q[3];
rz(-0.35202682) q[3];
sx q[3];
rz(2.5792518) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.8865108) q[2];
sx q[2];
rz(-1.3811029) q[2];
sx q[2];
rz(1.8908267) q[2];
rz(-1.4261774) q[3];
sx q[3];
rz(-0.91606796) q[3];
sx q[3];
rz(2.1616518) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.132906) q[0];
sx q[0];
rz(-2.0693021) q[0];
sx q[0];
rz(2.5426478) q[0];
rz(1.8006181) q[1];
sx q[1];
rz(-2.1913765) q[1];
sx q[1];
rz(-0.96639955) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9430267) q[0];
sx q[0];
rz(-2.1417924) q[0];
sx q[0];
rz(2.1023554) q[0];
x q[1];
rz(2.5194089) q[2];
sx q[2];
rz(-1.5303648) q[2];
sx q[2];
rz(-1.7379023) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.4233154) q[1];
sx q[1];
rz(-1.6498955) q[1];
sx q[1];
rz(-3.0273816) q[1];
rz(0.74921272) q[3];
sx q[3];
rz(-1.2470761) q[3];
sx q[3];
rz(2.2183228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.0559343) q[2];
sx q[2];
rz(-2.3264383) q[2];
sx q[2];
rz(2.8404964) q[2];
rz(1.1931233) q[3];
sx q[3];
rz(-1.5914702) q[3];
sx q[3];
rz(-1.955207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10467228) q[0];
sx q[0];
rz(-1.5043229) q[0];
sx q[0];
rz(-1.954129) q[0];
rz(-1.9056412) q[1];
sx q[1];
rz(-2.1042447) q[1];
sx q[1];
rz(1.8240066) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31908195) q[0];
sx q[0];
rz(-1.3586449) q[0];
sx q[0];
rz(2.4427419) q[0];
rz(-1.7682398) q[2];
sx q[2];
rz(-1.3609481) q[2];
sx q[2];
rz(-0.91453493) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9356404) q[1];
sx q[1];
rz(-1.9899568) q[1];
sx q[1];
rz(1.6816891) q[1];
rz(1.3665974) q[3];
sx q[3];
rz(-2.6150828) q[3];
sx q[3];
rz(-0.029475676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.1304156) q[2];
sx q[2];
rz(-1.7480363) q[2];
sx q[2];
rz(2.9349566) q[2];
rz(0.7080428) q[3];
sx q[3];
rz(-2.9338624) q[3];
sx q[3];
rz(2.2487683) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3669423) q[0];
sx q[0];
rz(-0.21454021) q[0];
sx q[0];
rz(-0.88622093) q[0];
rz(2.1318502) q[1];
sx q[1];
rz(-2.2354398) q[1];
sx q[1];
rz(1.2264235) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0013989) q[0];
sx q[0];
rz(-1.5045907) q[0];
sx q[0];
rz(-2.2982236) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2573651) q[2];
sx q[2];
rz(-1.9366169) q[2];
sx q[2];
rz(-2.5522752) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7734079) q[1];
sx q[1];
rz(-0.28973026) q[1];
sx q[1];
rz(-0.75399953) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6552116) q[3];
sx q[3];
rz(-1.2380621) q[3];
sx q[3];
rz(-2.2330724) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.6440789) q[2];
sx q[2];
rz(-1.3138694) q[2];
sx q[2];
rz(0.99299661) q[2];
rz(-1.3126866) q[3];
sx q[3];
rz(-2.1195181) q[3];
sx q[3];
rz(2.657857) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.118367) q[0];
sx q[0];
rz(-2.826773) q[0];
sx q[0];
rz(-1.9556048) q[0];
rz(-1.4233308) q[1];
sx q[1];
rz(-1.6512197) q[1];
sx q[1];
rz(0.58247724) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3956446) q[0];
sx q[0];
rz(-2.1714604) q[0];
sx q[0];
rz(3.0546741) q[0];
x q[1];
rz(0.8930348) q[2];
sx q[2];
rz(-1.4544011) q[2];
sx q[2];
rz(1.900577) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.40618784) q[1];
sx q[1];
rz(-1.5899854) q[1];
sx q[1];
rz(0.71449844) q[1];
x q[2];
rz(-2.7353103) q[3];
sx q[3];
rz(-1.632382) q[3];
sx q[3];
rz(0.53698925) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.4867268) q[2];
sx q[2];
rz(-1.5348237) q[2];
sx q[2];
rz(3.0598818) q[2];
rz(-2.667526) q[3];
sx q[3];
rz(-1.2636377) q[3];
sx q[3];
rz(-1.4985532) q[3];
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
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.896647) q[0];
sx q[0];
rz(-1.1279673) q[0];
sx q[0];
rz(0.0078049302) q[0];
rz(-1.7410949) q[1];
sx q[1];
rz(-0.85406071) q[1];
sx q[1];
rz(2.0369464) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1305599) q[0];
sx q[0];
rz(-1.4078119) q[0];
sx q[0];
rz(0.66793229) q[0];
rz(-0.64013021) q[2];
sx q[2];
rz(-1.374561) q[2];
sx q[2];
rz(3.047608) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.066597477) q[1];
sx q[1];
rz(-1.0269594) q[1];
sx q[1];
rz(0.72504136) q[1];
rz(1.2475345) q[3];
sx q[3];
rz(-0.48719104) q[3];
sx q[3];
rz(-2.4647453) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.2509987) q[2];
sx q[2];
rz(-2.4119174) q[2];
sx q[2];
rz(-0.55523038) q[2];
rz(-0.17272078) q[3];
sx q[3];
rz(-1.3135066) q[3];
sx q[3];
rz(1.593332) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5531439) q[0];
sx q[0];
rz(-0.48184904) q[0];
sx q[0];
rz(0.061766457) q[0];
rz(-0.24208367) q[1];
sx q[1];
rz(-0.37529072) q[1];
sx q[1];
rz(-1.1118836) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0430849) q[0];
sx q[0];
rz(-2.3102009) q[0];
sx q[0];
rz(0.99494536) q[0];
rz(-1.3948453) q[2];
sx q[2];
rz(-0.07882747) q[2];
sx q[2];
rz(1.4836756) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.8916787) q[1];
sx q[1];
rz(-2.2644682) q[1];
sx q[1];
rz(1.2077043) q[1];
rz(-pi) q[2];
rz(0.32825177) q[3];
sx q[3];
rz(-1.8040856) q[3];
sx q[3];
rz(1.0321898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3322488) q[2];
sx q[2];
rz(-2.377254) q[2];
sx q[2];
rz(-0.81364441) q[2];
rz(1.404473) q[3];
sx q[3];
rz(-2.9128894) q[3];
sx q[3];
rz(-0.61541921) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5161045) q[0];
sx q[0];
rz(-1.7913211) q[0];
sx q[0];
rz(-1.425449) q[0];
rz(-1.6199934) q[1];
sx q[1];
rz(-0.75192538) q[1];
sx q[1];
rz(2.5040748) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2190327) q[0];
sx q[0];
rz(-2.1171283) q[0];
sx q[0];
rz(2.4541897) q[0];
rz(1.2136739) q[2];
sx q[2];
rz(-2.2631858) q[2];
sx q[2];
rz(-1.1536319) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9962822) q[1];
sx q[1];
rz(-1.1728371) q[1];
sx q[1];
rz(-0.63509649) q[1];
rz(-pi) q[2];
x q[2];
rz(0.57070891) q[3];
sx q[3];
rz(-1.4548886) q[3];
sx q[3];
rz(1.4012208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1311538) q[2];
sx q[2];
rz(-0.70019478) q[2];
sx q[2];
rz(2.0054224) q[2];
rz(1.4853959) q[3];
sx q[3];
rz(-2.6172726) q[3];
sx q[3];
rz(1.9320528) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
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
rz(1.5159601) q[0];
sx q[0];
rz(-0.78173286) q[0];
sx q[0];
rz(-2.4654454) q[0];
rz(0.82538429) q[1];
sx q[1];
rz(-2.717658) q[1];
sx q[1];
rz(2.1264145) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7069106) q[0];
sx q[0];
rz(-0.43723956) q[0];
sx q[0];
rz(2.8291563) q[0];
rz(-pi) q[1];
x q[1];
rz(0.19998156) q[2];
sx q[2];
rz(-2.4790384) q[2];
sx q[2];
rz(-2.8560864) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.2296914) q[1];
sx q[1];
rz(-1.9209314) q[1];
sx q[1];
rz(-1.9446816) q[1];
rz(2.3750651) q[3];
sx q[3];
rz(-1.818728) q[3];
sx q[3];
rz(-0.31162308) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4201346) q[2];
sx q[2];
rz(-2.9253503) q[2];
sx q[2];
rz(0.96735111) q[2];
rz(1.5445276) q[3];
sx q[3];
rz(-1.8287851) q[3];
sx q[3];
rz(2.7887204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6185146) q[0];
sx q[0];
rz(-1.5627562) q[0];
sx q[0];
rz(0.05649795) q[0];
rz(2.0691195) q[1];
sx q[1];
rz(-1.2697376) q[1];
sx q[1];
rz(-1.4046232) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0427986) q[0];
sx q[0];
rz(-2.7569175) q[0];
sx q[0];
rz(-1.2205475) q[0];
rz(2.8374412) q[2];
sx q[2];
rz(-1.4545822) q[2];
sx q[2];
rz(-2.846037) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.23552588) q[1];
sx q[1];
rz(-2.0350254) q[1];
sx q[1];
rz(-2.3266351) q[1];
rz(-pi) q[2];
rz(0.15053648) q[3];
sx q[3];
rz(-1.5621462) q[3];
sx q[3];
rz(2.0899525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.29356062) q[2];
sx q[2];
rz(-0.75389391) q[2];
sx q[2];
rz(0.30612293) q[2];
rz(-0.39811578) q[3];
sx q[3];
rz(-1.476215) q[3];
sx q[3];
rz(1.0242296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8067779) q[0];
sx q[0];
rz(-1.0659185) q[0];
sx q[0];
rz(0.48068) q[0];
rz(-1.1595935) q[1];
sx q[1];
rz(-1.0212785) q[1];
sx q[1];
rz(-2.8550128) q[1];
rz(2.2904916) q[2];
sx q[2];
rz(-2.7534178) q[2];
sx q[2];
rz(-0.22321246) q[2];
rz(-1.0767827) q[3];
sx q[3];
rz(-1.2772588) q[3];
sx q[3];
rz(0.96989934) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];