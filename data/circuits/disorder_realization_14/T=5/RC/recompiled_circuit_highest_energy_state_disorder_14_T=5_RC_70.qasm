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
rz(1.5481663) q[0];
sx q[0];
rz(3.610008) q[0];
sx q[0];
rz(9.2863291) q[0];
rz(-1.2568714) q[1];
sx q[1];
rz(-1.0705907) q[1];
sx q[1];
rz(0.021477403) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81651743) q[0];
sx q[0];
rz(-1.3680352) q[0];
sx q[0];
rz(1.9295503) q[0];
rz(-pi) q[1];
x q[1];
rz(0.98928605) q[2];
sx q[2];
rz(-2.044319) q[2];
sx q[2];
rz(0.038528942) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.8285458) q[1];
sx q[1];
rz(-1.4171825) q[1];
sx q[1];
rz(1.7474354) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6329166) q[3];
sx q[3];
rz(-2.4812316) q[3];
sx q[3];
rz(2.3213399) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.59637493) q[2];
sx q[2];
rz(-0.45520982) q[2];
sx q[2];
rz(-1.5111766) q[2];
rz(-3.0916072) q[3];
sx q[3];
rz(-1.2286011) q[3];
sx q[3];
rz(2.889192) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2216457) q[0];
sx q[0];
rz(-1.9939461) q[0];
sx q[0];
rz(-3.1044712) q[0];
rz(0.014017398) q[1];
sx q[1];
rz(-2.5633096) q[1];
sx q[1];
rz(0.23606539) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0229313) q[0];
sx q[0];
rz(-0.076181024) q[0];
sx q[0];
rz(-0.19019048) q[0];
rz(-0.71546258) q[2];
sx q[2];
rz(-1.9505378) q[2];
sx q[2];
rz(-1.9850784) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.46690698) q[1];
sx q[1];
rz(-1.0471736) q[1];
sx q[1];
rz(-0.44498131) q[1];
x q[2];
rz(-1.2108324) q[3];
sx q[3];
rz(-1.3356555) q[3];
sx q[3];
rz(2.1699303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.0125121) q[2];
sx q[2];
rz(-1.7725638) q[2];
sx q[2];
rz(3.1079666) q[2];
rz(-2.8513837) q[3];
sx q[3];
rz(-0.84607327) q[3];
sx q[3];
rz(0.39053759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.520312) q[0];
sx q[0];
rz(-0.97977591) q[0];
sx q[0];
rz(2.8424971) q[0];
rz(-1.5215123) q[1];
sx q[1];
rz(-3.1144996) q[1];
sx q[1];
rz(-3.1153968) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8854302) q[0];
sx q[0];
rz(-3.1195398) q[0];
sx q[0];
rz(-1.2718448) q[0];
rz(-1.8511591) q[2];
sx q[2];
rz(-1.4329264) q[2];
sx q[2];
rz(3.1119135) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.44134162) q[1];
sx q[1];
rz(-2.1192351) q[1];
sx q[1];
rz(2.3943108) q[1];
x q[2];
rz(-0.45460017) q[3];
sx q[3];
rz(-0.48991436) q[3];
sx q[3];
rz(-0.59441209) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.4574778) q[2];
sx q[2];
rz(-0.44197765) q[2];
sx q[2];
rz(2.5341865) q[2];
rz(2.7720747) q[3];
sx q[3];
rz(-1.642546) q[3];
sx q[3];
rz(-3.1336237) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58005106) q[0];
sx q[0];
rz(-0.20195584) q[0];
sx q[0];
rz(-1.1308905) q[0];
rz(2.789403) q[1];
sx q[1];
rz(-0.49030855) q[1];
sx q[1];
rz(2.9255548) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0982042) q[0];
sx q[0];
rz(-1.3281315) q[0];
sx q[0];
rz(-2.674874) q[0];
rz(-pi) q[1];
rz(-0.44682002) q[2];
sx q[2];
rz(-1.4654612) q[2];
sx q[2];
rz(-0.0086431816) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.46680784) q[1];
sx q[1];
rz(-2.4744316) q[1];
sx q[1];
rz(2.2391367) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0918255) q[3];
sx q[3];
rz(-0.19102016) q[3];
sx q[3];
rz(2.614792) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.617368) q[2];
sx q[2];
rz(-0.91297954) q[2];
sx q[2];
rz(0.53140223) q[2];
rz(1.4517387) q[3];
sx q[3];
rz(-2.6302591) q[3];
sx q[3];
rz(-1.0303729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3568929) q[0];
sx q[0];
rz(-2.0578616) q[0];
sx q[0];
rz(-2.8392131) q[0];
rz(1.1131635) q[1];
sx q[1];
rz(-1.0013564) q[1];
sx q[1];
rz(-2.0804292) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36398023) q[0];
sx q[0];
rz(-1.5060695) q[0];
sx q[0];
rz(1.510468) q[0];
x q[1];
rz(0.2246173) q[2];
sx q[2];
rz(-1.6664522) q[2];
sx q[2];
rz(0.10570733) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.1984242) q[1];
sx q[1];
rz(-2.3454002) q[1];
sx q[1];
rz(-1.7502341) q[1];
rz(-pi) q[2];
rz(2.4441884) q[3];
sx q[3];
rz(-0.38129574) q[3];
sx q[3];
rz(0.21496102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.314996) q[2];
sx q[2];
rz(-1.3773409) q[2];
sx q[2];
rz(-2.946741) q[2];
rz(2.3991614) q[3];
sx q[3];
rz(-2.4104379) q[3];
sx q[3];
rz(-2.8661695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9575397) q[0];
sx q[0];
rz(-2.4803211) q[0];
sx q[0];
rz(1.1141962) q[0];
rz(2.8320352) q[1];
sx q[1];
rz(-0.93300262) q[1];
sx q[1];
rz(-1.385744) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4969869) q[0];
sx q[0];
rz(-2.806555) q[0];
sx q[0];
rz(-2.1315639) q[0];
x q[1];
rz(-1.0631752) q[2];
sx q[2];
rz(-1.4709899) q[2];
sx q[2];
rz(-0.55528477) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.09090611) q[1];
sx q[1];
rz(-2.202817) q[1];
sx q[1];
rz(1.8952096) q[1];
rz(-pi) q[2];
x q[2];
rz(0.65151229) q[3];
sx q[3];
rz(-0.88884547) q[3];
sx q[3];
rz(-1.3523952) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.4625357) q[2];
sx q[2];
rz(-0.44687301) q[2];
sx q[2];
rz(0.40979579) q[2];
rz(-3.0637686) q[3];
sx q[3];
rz(-1.2160559) q[3];
sx q[3];
rz(-0.76243824) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78793144) q[0];
sx q[0];
rz(-0.15058148) q[0];
sx q[0];
rz(1.6702363) q[0];
rz(0.81368601) q[1];
sx q[1];
rz(-2.4206471) q[1];
sx q[1];
rz(-2.241316) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6532458) q[0];
sx q[0];
rz(-1.6890959) q[0];
sx q[0];
rz(1.5100046) q[0];
x q[1];
rz(1.288432) q[2];
sx q[2];
rz(-0.29855628) q[2];
sx q[2];
rz(2.3207842) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.011411) q[1];
sx q[1];
rz(-2.8754042) q[1];
sx q[1];
rz(-1.5397416) q[1];
rz(-0.26746349) q[3];
sx q[3];
rz(-0.24175343) q[3];
sx q[3];
rz(3.0694244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.6206996) q[2];
sx q[2];
rz(-0.40105477) q[2];
sx q[2];
rz(-0.14191423) q[2];
rz(2.9698931) q[3];
sx q[3];
rz(-1.2407691) q[3];
sx q[3];
rz(-2.5920674) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5080344) q[0];
sx q[0];
rz(-1.1548076) q[0];
sx q[0];
rz(1.5800193) q[0];
rz(2.436807) q[1];
sx q[1];
rz(-0.90122688) q[1];
sx q[1];
rz(0.27552342) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6632623) q[0];
sx q[0];
rz(-1.5092634) q[0];
sx q[0];
rz(-1.4429773) q[0];
rz(0.24532206) q[2];
sx q[2];
rz(-1.6082014) q[2];
sx q[2];
rz(2.268301) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.49373473) q[1];
sx q[1];
rz(-2.0701412) q[1];
sx q[1];
rz(2.8060421) q[1];
rz(-pi) q[2];
rz(-0.51361689) q[3];
sx q[3];
rz(-0.73459638) q[3];
sx q[3];
rz(-2.9111322) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.0602818) q[2];
sx q[2];
rz(-1.5014481) q[2];
sx q[2];
rz(0.29067972) q[2];
rz(-2.3447013) q[3];
sx q[3];
rz(-2.6698038) q[3];
sx q[3];
rz(-2.1577788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31442916) q[0];
sx q[0];
rz(-1.4845347) q[0];
sx q[0];
rz(2.2737801) q[0];
rz(2.9293291) q[1];
sx q[1];
rz(-2.2607195) q[1];
sx q[1];
rz(0.27264047) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72948706) q[0];
sx q[0];
rz(-0.22815591) q[0];
sx q[0];
rz(1.9081997) q[0];
rz(-pi) q[1];
rz(-0.10137239) q[2];
sx q[2];
rz(-2.3194312) q[2];
sx q[2];
rz(-3.0790975) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.3301047) q[1];
sx q[1];
rz(-1.5327785) q[1];
sx q[1];
rz(2.2056715) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8707859) q[3];
sx q[3];
rz(-1.9686896) q[3];
sx q[3];
rz(0.34664819) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(3.0923882) q[2];
sx q[2];
rz(-0.35085756) q[2];
sx q[2];
rz(1.997088) q[2];
rz(-0.62659621) q[3];
sx q[3];
rz(-2.383039) q[3];
sx q[3];
rz(0.315061) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
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
rz(0.7476244) q[0];
sx q[0];
rz(-2.4691041) q[0];
sx q[0];
rz(-2.5482063) q[0];
rz(2.4001832) q[1];
sx q[1];
rz(-0.64677042) q[1];
sx q[1];
rz(1.5945565) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55047148) q[0];
sx q[0];
rz(-2.897399) q[0];
sx q[0];
rz(2.0376192) q[0];
rz(-0.41489652) q[2];
sx q[2];
rz(-1.7977568) q[2];
sx q[2];
rz(1.0912053) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.98982544) q[1];
sx q[1];
rz(-1.5011468) q[1];
sx q[1];
rz(-2.4098957) q[1];
rz(-0.029742777) q[3];
sx q[3];
rz(-1.3689201) q[3];
sx q[3];
rz(-1.0591398) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.048162248) q[2];
sx q[2];
rz(-2.5967345) q[2];
sx q[2];
rz(2.1989934) q[2];
rz(-1.7132828) q[3];
sx q[3];
rz(-1.2497679) q[3];
sx q[3];
rz(1.1552756) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1930595) q[0];
sx q[0];
rz(-1.5567224) q[0];
sx q[0];
rz(-1.3514883) q[0];
rz(1.9652741) q[1];
sx q[1];
rz(-0.90570025) q[1];
sx q[1];
rz(-0.66014231) q[1];
rz(-1.0084441) q[2];
sx q[2];
rz(-2.8165419) q[2];
sx q[2];
rz(-1.4680924) q[2];
rz(0.42849937) q[3];
sx q[3];
rz(-0.57572031) q[3];
sx q[3];
rz(-2.4888103) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
