OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.6736421) q[0];
sx q[0];
rz(-1.8704432) q[0];
sx q[0];
rz(-0.76572642) q[0];
rz(5.8869047) q[1];
sx q[1];
rz(6.1904391) q[1];
sx q[1];
rz(9.9775597) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7460913) q[0];
sx q[0];
rz(-1.1679018) q[0];
sx q[0];
rz(-0.42568107) q[0];
rz(0.80443126) q[2];
sx q[2];
rz(-1.0519769) q[2];
sx q[2];
rz(-3.0592953) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8914767) q[1];
sx q[1];
rz(-1.0179903) q[1];
sx q[1];
rz(2.7673036) q[1];
rz(-pi) q[2];
rz(-2.2329198) q[3];
sx q[3];
rz(-1.931802) q[3];
sx q[3];
rz(-0.26546016) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.3679686) q[2];
sx q[2];
rz(-1.5017193) q[2];
sx q[2];
rz(-0.088851301) q[2];
rz(0.39696524) q[3];
sx q[3];
rz(-1.0395972) q[3];
sx q[3];
rz(1.6840434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6615768) q[0];
sx q[0];
rz(-2.6341697) q[0];
sx q[0];
rz(0.85451025) q[0];
rz(-0.62141934) q[1];
sx q[1];
rz(-2.8050551) q[1];
sx q[1];
rz(-1.052676) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1922469) q[0];
sx q[0];
rz(-1.7070012) q[0];
sx q[0];
rz(-1.6736567) q[0];
rz(-pi) q[1];
rz(-1.1243049) q[2];
sx q[2];
rz(-1.1928176) q[2];
sx q[2];
rz(1.93731) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.8206975) q[1];
sx q[1];
rz(-1.7716494) q[1];
sx q[1];
rz(-0.83495514) q[1];
x q[2];
rz(1.5731811) q[3];
sx q[3];
rz(-2.5376476) q[3];
sx q[3];
rz(1.2177474) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.4738327) q[2];
sx q[2];
rz(-0.86060539) q[2];
sx q[2];
rz(-2.1916126) q[2];
rz(-2.1330323) q[3];
sx q[3];
rz(-0.63298321) q[3];
sx q[3];
rz(2.8620201) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94473332) q[0];
sx q[0];
rz(-1.9254528) q[0];
sx q[0];
rz(-1.7945633) q[0];
rz(-1.0579717) q[1];
sx q[1];
rz(-2.8949013) q[1];
sx q[1];
rz(-0.11014858) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63907097) q[0];
sx q[0];
rz(-0.57146954) q[0];
sx q[0];
rz(-0.95614627) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4227432) q[2];
sx q[2];
rz(-1.6534936) q[2];
sx q[2];
rz(0.18872866) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.2825605) q[1];
sx q[1];
rz(-0.60884005) q[1];
sx q[1];
rz(0.18688272) q[1];
rz(-pi) q[2];
rz(2.6950257) q[3];
sx q[3];
rz(-1.1249229) q[3];
sx q[3];
rz(-1.2620516) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.18016732) q[2];
sx q[2];
rz(-1.5632997) q[2];
sx q[2];
rz(0.041672826) q[2];
rz(-2.9336119) q[3];
sx q[3];
rz(-2.9198923) q[3];
sx q[3];
rz(0.25377932) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5623077) q[0];
sx q[0];
rz(-1.8203745) q[0];
sx q[0];
rz(2.1606309) q[0];
rz(2.7085069) q[1];
sx q[1];
rz(-1.667645) q[1];
sx q[1];
rz(1.513419) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0312322) q[0];
sx q[0];
rz(-1.1735532) q[0];
sx q[0];
rz(1.375694) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.01362205) q[2];
sx q[2];
rz(-2.2421466) q[2];
sx q[2];
rz(-1.7126132) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.7119496) q[1];
sx q[1];
rz(-0.50772655) q[1];
sx q[1];
rz(-1.7832548) q[1];
x q[2];
rz(0.28339213) q[3];
sx q[3];
rz(-0.95855721) q[3];
sx q[3];
rz(1.3536842) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.5341586) q[2];
sx q[2];
rz(-2.6223493) q[2];
sx q[2];
rz(-0.81238166) q[2];
rz(1.2477929) q[3];
sx q[3];
rz(-1.2500074) q[3];
sx q[3];
rz(1.8947424) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3355376) q[0];
sx q[0];
rz(-1.0220818) q[0];
sx q[0];
rz(-1.442765) q[0];
rz(-0.45817786) q[1];
sx q[1];
rz(-1.5135601) q[1];
sx q[1];
rz(0.75622574) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9423507) q[0];
sx q[0];
rz(-2.1926542) q[0];
sx q[0];
rz(-2.950475) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.64018546) q[2];
sx q[2];
rz(-2.019503) q[2];
sx q[2];
rz(1.9404836) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.56179433) q[1];
sx q[1];
rz(-1.7057912) q[1];
sx q[1];
rz(1.6588103) q[1];
rz(1.4706572) q[3];
sx q[3];
rz(-2.7960145) q[3];
sx q[3];
rz(-2.341193) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.4592287) q[2];
sx q[2];
rz(-1.3882779) q[2];
sx q[2];
rz(2.1759822) q[2];
rz(0.56973488) q[3];
sx q[3];
rz(-2.0163048) q[3];
sx q[3];
rz(-1.4558571) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4023034) q[0];
sx q[0];
rz(-1.9434513) q[0];
sx q[0];
rz(-1.7560316) q[0];
rz(2.3784474) q[1];
sx q[1];
rz(-1.4475854) q[1];
sx q[1];
rz(-3.0398583) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3107233) q[0];
sx q[0];
rz(-2.1889157) q[0];
sx q[0];
rz(0.11316664) q[0];
rz(-pi) q[1];
rz(-2.7587682) q[2];
sx q[2];
rz(-1.6144132) q[2];
sx q[2];
rz(0.62784615) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.49470529) q[1];
sx q[1];
rz(-0.27783074) q[1];
sx q[1];
rz(1.3459413) q[1];
rz(-pi) q[2];
rz(2.3699573) q[3];
sx q[3];
rz(-2.2766678) q[3];
sx q[3];
rz(0.40654686) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(3.0330641) q[2];
sx q[2];
rz(-1.6988924) q[2];
sx q[2];
rz(0.70890439) q[2];
rz(-2.3809643) q[3];
sx q[3];
rz(-1.9899188) q[3];
sx q[3];
rz(2.2061548) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0660504) q[0];
sx q[0];
rz(-1.9882555) q[0];
sx q[0];
rz(-0.72108889) q[0];
rz(-2.1636294) q[1];
sx q[1];
rz(-2.5826192) q[1];
sx q[1];
rz(1.5616547) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6128591) q[0];
sx q[0];
rz(-1.4943143) q[0];
sx q[0];
rz(-1.4974817) q[0];
rz(-pi) q[1];
rz(-0.97550895) q[2];
sx q[2];
rz(-1.7726608) q[2];
sx q[2];
rz(-1.5306115) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.4750907) q[1];
sx q[1];
rz(-1.9618841) q[1];
sx q[1];
rz(-0.60152454) q[1];
rz(-1.1317113) q[3];
sx q[3];
rz(-2.4222322) q[3];
sx q[3];
rz(-2.2124825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.28479031) q[2];
sx q[2];
rz(-2.6140116) q[2];
sx q[2];
rz(-1.6214726) q[2];
rz(2.704845) q[3];
sx q[3];
rz(-2.2468061) q[3];
sx q[3];
rz(-0.53957087) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0590416) q[0];
sx q[0];
rz(-2.3024547) q[0];
sx q[0];
rz(0.82373291) q[0];
rz(2.0027022) q[1];
sx q[1];
rz(-1.7887807) q[1];
sx q[1];
rz(-1.1762071) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7438854) q[0];
sx q[0];
rz(-1.526014) q[0];
sx q[0];
rz(-1.0683879) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5912639) q[2];
sx q[2];
rz(-2.4541353) q[2];
sx q[2];
rz(2.2000809) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.215428) q[1];
sx q[1];
rz(-1.3326982) q[1];
sx q[1];
rz(1.5808616) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5309991) q[3];
sx q[3];
rz(-2.1615513) q[3];
sx q[3];
rz(-1.9692957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.99522432) q[2];
sx q[2];
rz(-2.5547042) q[2];
sx q[2];
rz(-0.5874908) q[2];
rz(2.7557709) q[3];
sx q[3];
rz(-0.84857517) q[3];
sx q[3];
rz(-2.2390656) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1189482) q[0];
sx q[0];
rz(-2.1653439) q[0];
sx q[0];
rz(-2.5318085) q[0];
rz(1.7976286) q[1];
sx q[1];
rz(-1.1577497) q[1];
sx q[1];
rz(-2.9885898) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6082108) q[0];
sx q[0];
rz(-2.675288) q[0];
sx q[0];
rz(-0.043688579) q[0];
x q[1];
rz(-2.346368) q[2];
sx q[2];
rz(-1.5780026) q[2];
sx q[2];
rz(-2.0237271) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.24125946) q[1];
sx q[1];
rz(-1.5521084) q[1];
sx q[1];
rz(2.2128723) q[1];
rz(-2.2629635) q[3];
sx q[3];
rz(-1.9464916) q[3];
sx q[3];
rz(-0.13334783) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.0387705) q[2];
sx q[2];
rz(-1.5800579) q[2];
sx q[2];
rz(-2.5193396) q[2];
rz(1.9410939) q[3];
sx q[3];
rz(-1.4193204) q[3];
sx q[3];
rz(2.5672369) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1822561) q[0];
sx q[0];
rz(-3.0756364) q[0];
sx q[0];
rz(0.35520735) q[0];
rz(-1.1025053) q[1];
sx q[1];
rz(-3.1246154) q[1];
sx q[1];
rz(0.65974081) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3400354) q[0];
sx q[0];
rz(-1.7424955) q[0];
sx q[0];
rz(-0.6562018) q[0];
rz(1.9842568) q[2];
sx q[2];
rz(-1.454118) q[2];
sx q[2];
rz(1.9422187) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7463508) q[1];
sx q[1];
rz(-1.244427) q[1];
sx q[1];
rz(2.8119948) q[1];
rz(-1.5054613) q[3];
sx q[3];
rz(-2.4628277) q[3];
sx q[3];
rz(-0.15798727) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.3416662) q[2];
sx q[2];
rz(-0.081192668) q[2];
sx q[2];
rz(2.9426835) q[2];
rz(-0.79814664) q[3];
sx q[3];
rz(-1.9559559) q[3];
sx q[3];
rz(-1.2822436) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2220919) q[0];
sx q[0];
rz(-1.5380479) q[0];
sx q[0];
rz(2.0684239) q[0];
rz(2.3566698) q[1];
sx q[1];
rz(-2.2833318) q[1];
sx q[1];
rz(-2.2699184) q[1];
rz(-1.9368108) q[2];
sx q[2];
rz(-0.13199619) q[2];
sx q[2];
rz(0.97313626) q[2];
rz(-1.296923) q[3];
sx q[3];
rz(-0.92007888) q[3];
sx q[3];
rz(-2.32213) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
