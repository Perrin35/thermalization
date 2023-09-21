OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.4361753) q[0];
sx q[0];
rz(-0.56646148) q[0];
sx q[0];
rz(0.17106549) q[0];
rz(-1.6879727) q[1];
sx q[1];
rz(3.4847335) q[1];
sx q[1];
rz(7.6141678) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91141191) q[0];
sx q[0];
rz(-2.4898306) q[0];
sx q[0];
rz(-0.08641152) q[0];
rz(-0.51585977) q[2];
sx q[2];
rz(-1.5896279) q[2];
sx q[2];
rz(-3.0992103) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.060974412) q[1];
sx q[1];
rz(-1.8218826) q[1];
sx q[1];
rz(0.017647839) q[1];
x q[2];
rz(-1.347723) q[3];
sx q[3];
rz(-2.2017456) q[3];
sx q[3];
rz(2.0032721) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.3657637) q[2];
sx q[2];
rz(-0.77314955) q[2];
sx q[2];
rz(1.8120871) q[2];
rz(-1.7261516) q[3];
sx q[3];
rz(-1.5751782) q[3];
sx q[3];
rz(-2.9392488) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99130327) q[0];
sx q[0];
rz(-2.5973899) q[0];
sx q[0];
rz(1.4527028) q[0];
rz(-2.9406722) q[1];
sx q[1];
rz(-1.0849489) q[1];
sx q[1];
rz(-0.30028775) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18314221) q[0];
sx q[0];
rz(-2.169325) q[0];
sx q[0];
rz(1.4248225) q[0];
rz(-pi) q[1];
rz(0.57057256) q[2];
sx q[2];
rz(-1.3010446) q[2];
sx q[2];
rz(2.3938993) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-3.049388) q[1];
sx q[1];
rz(-2.7687679) q[1];
sx q[1];
rz(2.2224109) q[1];
rz(-pi) q[2];
rz(-1.3817915) q[3];
sx q[3];
rz(-1.4100037) q[3];
sx q[3];
rz(-0.27634987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.369027) q[2];
sx q[2];
rz(-0.3521266) q[2];
sx q[2];
rz(1.0380113) q[2];
rz(2.299262) q[3];
sx q[3];
rz(-1.4146283) q[3];
sx q[3];
rz(0.37030181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
rz(-1.5791941) q[0];
sx q[0];
rz(-0.47214046) q[0];
sx q[0];
rz(2.753479) q[0];
rz(-0.072470486) q[1];
sx q[1];
rz(-1.4265172) q[1];
sx q[1];
rz(-0.31633502) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64992031) q[0];
sx q[0];
rz(-2.8965817) q[0];
sx q[0];
rz(-1.578376) q[0];
x q[1];
rz(-2.5306273) q[2];
sx q[2];
rz(-0.66662153) q[2];
sx q[2];
rz(1.1952458) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.12831941) q[1];
sx q[1];
rz(-2.0187906) q[1];
sx q[1];
rz(2.0112579) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2315138) q[3];
sx q[3];
rz(-1.8456036) q[3];
sx q[3];
rz(2.2846082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.0324273) q[2];
sx q[2];
rz(-0.18523231) q[2];
sx q[2];
rz(-0.16080984) q[2];
rz(-3.0155449) q[3];
sx q[3];
rz(-1.7536609) q[3];
sx q[3];
rz(-2.1179312) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9008824) q[0];
sx q[0];
rz(-0.70403376) q[0];
sx q[0];
rz(0.8316935) q[0];
rz(-1.3446993) q[1];
sx q[1];
rz(-0.99935499) q[1];
sx q[1];
rz(1.5136738) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8320223) q[0];
sx q[0];
rz(-2.4674468) q[0];
sx q[0];
rz(-0.073622965) q[0];
rz(-pi) q[1];
rz(2.1349147) q[2];
sx q[2];
rz(-1.5022105) q[2];
sx q[2];
rz(2.2401631) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.0537783) q[1];
sx q[1];
rz(-0.44856854) q[1];
sx q[1];
rz(1.5320369) q[1];
x q[2];
rz(-1.1553331) q[3];
sx q[3];
rz(-0.88298015) q[3];
sx q[3];
rz(0.26879877) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.7164798) q[2];
sx q[2];
rz(-1.0093062) q[2];
sx q[2];
rz(1.9173737) q[2];
rz(-2.7246357) q[3];
sx q[3];
rz(-1.0427534) q[3];
sx q[3];
rz(-2.6127889) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.058218) q[0];
sx q[0];
rz(-2.871802) q[0];
sx q[0];
rz(-1.8378687) q[0];
rz(2.6858792) q[1];
sx q[1];
rz(-0.2625176) q[1];
sx q[1];
rz(3.0900893) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18729678) q[0];
sx q[0];
rz(-1.4547537) q[0];
sx q[0];
rz(3.0412578) q[0];
rz(0.15437834) q[2];
sx q[2];
rz(-1.2882243) q[2];
sx q[2];
rz(0.64235657) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0806549) q[1];
sx q[1];
rz(-0.47532156) q[1];
sx q[1];
rz(1.2259543) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0052161) q[3];
sx q[3];
rz(-0.73835056) q[3];
sx q[3];
rz(0.34463681) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.45433989) q[2];
sx q[2];
rz(-1.3663102) q[2];
sx q[2];
rz(-2.5435737) q[2];
rz(2.5915742) q[3];
sx q[3];
rz(-2.9019182) q[3];
sx q[3];
rz(-2.0050744) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0710058) q[0];
sx q[0];
rz(-0.07659176) q[0];
sx q[0];
rz(2.9396074) q[0];
rz(-0.96549353) q[1];
sx q[1];
rz(-1.0243203) q[1];
sx q[1];
rz(-0.083267033) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9986388) q[0];
sx q[0];
rz(-2.8790701) q[0];
sx q[0];
rz(-0.72857626) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8667738) q[2];
sx q[2];
rz(-1.8289369) q[2];
sx q[2];
rz(-1.9354265) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.49230089) q[1];
sx q[1];
rz(-1.2052844) q[1];
sx q[1];
rz(-3.027012) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5842651) q[3];
sx q[3];
rz(-1.3152272) q[3];
sx q[3];
rz(-0.42423466) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.10963708) q[2];
sx q[2];
rz(-1.457931) q[2];
sx q[2];
rz(-1.3827682) q[2];
rz(-0.2494732) q[3];
sx q[3];
rz(-1.8981257) q[3];
sx q[3];
rz(-1.680254) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4270585) q[0];
sx q[0];
rz(-2.336851) q[0];
sx q[0];
rz(-0.89685857) q[0];
rz(-0.25009051) q[1];
sx q[1];
rz(-0.18083328) q[1];
sx q[1];
rz(1.1175964) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16545907) q[0];
sx q[0];
rz(-1.4384067) q[0];
sx q[0];
rz(0.73029851) q[0];
rz(-pi) q[1];
x q[1];
rz(2.932981) q[2];
sx q[2];
rz(-2.2715855) q[2];
sx q[2];
rz(2.029062) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.84278216) q[1];
sx q[1];
rz(-0.59814765) q[1];
sx q[1];
rz(-0.36069718) q[1];
rz(-pi) q[2];
rz(-2.696633) q[3];
sx q[3];
rz(-0.58044725) q[3];
sx q[3];
rz(0.85948932) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.0573132) q[2];
sx q[2];
rz(-1.3998312) q[2];
sx q[2];
rz(2.9220667) q[2];
rz(-2.7548742) q[3];
sx q[3];
rz(-0.76824776) q[3];
sx q[3];
rz(2.1462671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0893843) q[0];
sx q[0];
rz(-1.5445671) q[0];
sx q[0];
rz(2.9242933) q[0];
rz(-0.030933881) q[1];
sx q[1];
rz(-0.63806454) q[1];
sx q[1];
rz(-0.6589748) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17651672) q[0];
sx q[0];
rz(-0.79916164) q[0];
sx q[0];
rz(1.6060711) q[0];
rz(-pi) q[1];
rz(2.4202004) q[2];
sx q[2];
rz(-2.112769) q[2];
sx q[2];
rz(1.1162356) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.98098552) q[1];
sx q[1];
rz(-1.3541344) q[1];
sx q[1];
rz(-0.76088455) q[1];
x q[2];
rz(-0.068816618) q[3];
sx q[3];
rz(-2.2127082) q[3];
sx q[3];
rz(1.2115692) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.4387681) q[2];
sx q[2];
rz(-1.7610901) q[2];
sx q[2];
rz(-0.32021114) q[2];
rz(-1.5252339) q[3];
sx q[3];
rz(-1.9459008) q[3];
sx q[3];
rz(1.8922136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3808688) q[0];
sx q[0];
rz(-0.18146935) q[0];
sx q[0];
rz(0.6828126) q[0];
rz(-2.071351) q[1];
sx q[1];
rz(-1.2425334) q[1];
sx q[1];
rz(-1.210093) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8378728) q[0];
sx q[0];
rz(-1.923773) q[0];
sx q[0];
rz(-2.9205802) q[0];
x q[1];
rz(-2.1754873) q[2];
sx q[2];
rz(-1.7793057) q[2];
sx q[2];
rz(-1.0143806) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.51325071) q[1];
sx q[1];
rz(-1.2770137) q[1];
sx q[1];
rz(-2.252584) q[1];
rz(-1.1990511) q[3];
sx q[3];
rz(-0.98402714) q[3];
sx q[3];
rz(1.1128807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.5658297) q[2];
sx q[2];
rz(-1.016022) q[2];
sx q[2];
rz(-0.56662095) q[2];
rz(-2.2120655) q[3];
sx q[3];
rz(-0.98171392) q[3];
sx q[3];
rz(1.7738336) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.728445) q[0];
sx q[0];
rz(-0.25756535) q[0];
sx q[0];
rz(1.43191) q[0];
rz(-2.6152949) q[1];
sx q[1];
rz(-2.6142575) q[1];
sx q[1];
rz(0.79968232) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4096217) q[0];
sx q[0];
rz(-2.7831166) q[0];
sx q[0];
rz(1.8147857) q[0];
x q[1];
rz(-1.6522371) q[2];
sx q[2];
rz(-1.7310206) q[2];
sx q[2];
rz(2.3022431) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.77196808) q[1];
sx q[1];
rz(-2.5110285) q[1];
sx q[1];
rz(-2.1715013) q[1];
rz(-pi) q[2];
rz(-1.4950072) q[3];
sx q[3];
rz(-1.7160551) q[3];
sx q[3];
rz(-2.4491572) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.8563103) q[2];
sx q[2];
rz(-0.44390634) q[2];
sx q[2];
rz(2.4463859) q[2];
rz(0.55784145) q[3];
sx q[3];
rz(-1.3643684) q[3];
sx q[3];
rz(2.4485596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0859062) q[0];
sx q[0];
rz(-1.9491371) q[0];
sx q[0];
rz(0.51167713) q[0];
rz(-1.2790537) q[1];
sx q[1];
rz(-2.3976354) q[1];
sx q[1];
rz(2.4639113) q[1];
rz(0.33591253) q[2];
sx q[2];
rz(-1.5699785) q[2];
sx q[2];
rz(1.7287398) q[2];
rz(-2.9813319) q[3];
sx q[3];
rz(-0.80453034) q[3];
sx q[3];
rz(2.9321032) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];