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
rz(-2.3209651) q[0];
sx q[0];
rz(-2.4790915) q[0];
sx q[0];
rz(-2.878046) q[0];
rz(1.1072371) q[1];
sx q[1];
rz(4.4237408) q[1];
sx q[1];
rz(10.10034) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0209262) q[0];
sx q[0];
rz(-1.5772737) q[0];
sx q[0];
rz(-1.9064039) q[0];
x q[1];
rz(2.3604806) q[2];
sx q[2];
rz(-0.75133577) q[2];
sx q[2];
rz(1.9066115) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.88982595) q[1];
sx q[1];
rz(-1.7210362) q[1];
sx q[1];
rz(1.9637693) q[1];
rz(2.1249902) q[3];
sx q[3];
rz(-0.57376353) q[3];
sx q[3];
rz(0.91780182) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.88966173) q[2];
sx q[2];
rz(-1.2222922) q[2];
sx q[2];
rz(2.7737854) q[2];
rz(-2.8273072) q[3];
sx q[3];
rz(-0.29044423) q[3];
sx q[3];
rz(0.17816003) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80938584) q[0];
sx q[0];
rz(-0.089501373) q[0];
sx q[0];
rz(-1.7820763) q[0];
rz(-1.2844194) q[1];
sx q[1];
rz(-1.9764683) q[1];
sx q[1];
rz(3.0899835) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29362677) q[0];
sx q[0];
rz(-0.32362263) q[0];
sx q[0];
rz(-0.96903649) q[0];
rz(-pi) q[1];
rz(1.2003635) q[2];
sx q[2];
rz(-1.3784172) q[2];
sx q[2];
rz(-0.34483257) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7588531) q[1];
sx q[1];
rz(-1.570228) q[1];
sx q[1];
rz(1.9447536) q[1];
rz(-pi) q[2];
x q[2];
rz(0.69689023) q[3];
sx q[3];
rz(-2.5514388) q[3];
sx q[3];
rz(1.3499638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-3.0525558) q[2];
sx q[2];
rz(-2.8812228) q[2];
sx q[2];
rz(-0.54711771) q[2];
rz(-1.8735006) q[3];
sx q[3];
rz(-1.3601466) q[3];
sx q[3];
rz(0.33856302) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7608305) q[0];
sx q[0];
rz(-2.1050504) q[0];
sx q[0];
rz(1.8007675) q[0];
rz(-2.2876168) q[1];
sx q[1];
rz(-1.2869336) q[1];
sx q[1];
rz(2.8391848) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5216999) q[0];
sx q[0];
rz(-2.1850153) q[0];
sx q[0];
rz(2.8015561) q[0];
rz(-pi) q[1];
rz(-1.339389) q[2];
sx q[2];
rz(-1.1303899) q[2];
sx q[2];
rz(-2.6803859) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.9871618) q[1];
sx q[1];
rz(-0.87244528) q[1];
sx q[1];
rz(-2.4157546) q[1];
rz(-pi) q[2];
x q[2];
rz(0.27389483) q[3];
sx q[3];
rz(-2.0467826) q[3];
sx q[3];
rz(2.7416503) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.2053908) q[2];
sx q[2];
rz(-1.256029) q[2];
sx q[2];
rz(2.5993247) q[2];
rz(-0.99644709) q[3];
sx q[3];
rz(-2.9139329) q[3];
sx q[3];
rz(-0.12797932) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9484613) q[0];
sx q[0];
rz(-1.4840935) q[0];
sx q[0];
rz(1.2633854) q[0];
rz(2.7073233) q[1];
sx q[1];
rz(-2.3691005) q[1];
sx q[1];
rz(1.451452) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93218373) q[0];
sx q[0];
rz(-1.1157633) q[0];
sx q[0];
rz(-1.8952151) q[0];
rz(-pi) q[1];
rz(-2.5927438) q[2];
sx q[2];
rz(-1.9611231) q[2];
sx q[2];
rz(0.091920741) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.8122708) q[1];
sx q[1];
rz(-1.656161) q[1];
sx q[1];
rz(1.4567767) q[1];
rz(-pi) q[2];
rz(-2.583458) q[3];
sx q[3];
rz(-1.0252748) q[3];
sx q[3];
rz(1.5212052) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0749391) q[2];
sx q[2];
rz(-2.1982919) q[2];
sx q[2];
rz(2.2554876) q[2];
rz(-0.0063627176) q[3];
sx q[3];
rz(-1.8013835) q[3];
sx q[3];
rz(1.1191012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56115443) q[0];
sx q[0];
rz(-2.6851324) q[0];
sx q[0];
rz(-1.891267) q[0];
rz(-0.23288947) q[1];
sx q[1];
rz(-2.3857375) q[1];
sx q[1];
rz(-0.09045352) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88921493) q[0];
sx q[0];
rz(-1.4486635) q[0];
sx q[0];
rz(1.8218763) q[0];
rz(1.8402507) q[2];
sx q[2];
rz(-1.1134992) q[2];
sx q[2];
rz(-2.2895165) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.4793746) q[1];
sx q[1];
rz(-1.5081128) q[1];
sx q[1];
rz(0.65356837) q[1];
rz(1.8985604) q[3];
sx q[3];
rz(-1.6521427) q[3];
sx q[3];
rz(-1.3535172) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.49454409) q[2];
sx q[2];
rz(-0.19379751) q[2];
sx q[2];
rz(1.9336644) q[2];
rz(-0.9211933) q[3];
sx q[3];
rz(-2.2974206) q[3];
sx q[3];
rz(-0.45879656) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10814609) q[0];
sx q[0];
rz(-1.7195846) q[0];
sx q[0];
rz(2.6500927) q[0];
rz(2.4132701) q[1];
sx q[1];
rz(-2.0305384) q[1];
sx q[1];
rz(-0.19457766) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.174343) q[0];
sx q[0];
rz(-1.3493269) q[0];
sx q[0];
rz(0.99507777) q[0];
x q[1];
rz(-0.80912239) q[2];
sx q[2];
rz(-0.97130191) q[2];
sx q[2];
rz(-0.15669151) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.3272596) q[1];
sx q[1];
rz(-2.2047298) q[1];
sx q[1];
rz(-0.64822207) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4846912) q[3];
sx q[3];
rz(-1.3852775) q[3];
sx q[3];
rz(-0.82478648) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3538094) q[2];
sx q[2];
rz(-2.1259191) q[2];
sx q[2];
rz(0.0018399012) q[2];
rz(-1.6483824) q[3];
sx q[3];
rz(-2.4108678) q[3];
sx q[3];
rz(-2.8572237) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88347411) q[0];
sx q[0];
rz(-0.26171568) q[0];
sx q[0];
rz(-1.0650241) q[0];
rz(0.26271543) q[1];
sx q[1];
rz(-2.5698667) q[1];
sx q[1];
rz(-0.46331847) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.036788851) q[0];
sx q[0];
rz(-0.98829174) q[0];
sx q[0];
rz(0.1881442) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.58775522) q[2];
sx q[2];
rz(-2.3238733) q[2];
sx q[2];
rz(0.14434926) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.1541248) q[1];
sx q[1];
rz(-0.5041669) q[1];
sx q[1];
rz(-0.28480633) q[1];
rz(-0.024302146) q[3];
sx q[3];
rz(-1.6858471) q[3];
sx q[3];
rz(-0.76434154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-3.1411983) q[2];
sx q[2];
rz(-1.980915) q[2];
sx q[2];
rz(2.311643) q[2];
rz(-1.0478033) q[3];
sx q[3];
rz(-2.6959097) q[3];
sx q[3];
rz(-2.9629663) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8038427) q[0];
sx q[0];
rz(-0.99030817) q[0];
sx q[0];
rz(-2.3867699) q[0];
rz(2.7524475) q[1];
sx q[1];
rz(-1.6837589) q[1];
sx q[1];
rz(0.39774242) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5311814) q[0];
sx q[0];
rz(-1.2875778) q[0];
sx q[0];
rz(-2.8736658) q[0];
rz(-pi) q[1];
x q[1];
rz(0.55126528) q[2];
sx q[2];
rz(-1.3756075) q[2];
sx q[2];
rz(2.8697687) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.3263742) q[1];
sx q[1];
rz(-2.7644016) q[1];
sx q[1];
rz(0.81881028) q[1];
x q[2];
rz(0.74716178) q[3];
sx q[3];
rz(-1.0667757) q[3];
sx q[3];
rz(0.51591831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.8336746) q[2];
sx q[2];
rz(-2.8454744) q[2];
sx q[2];
rz(2.8274242) q[2];
rz(1.2416154) q[3];
sx q[3];
rz(-1.5528468) q[3];
sx q[3];
rz(-1.9023021) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32387963) q[0];
sx q[0];
rz(-0.68809026) q[0];
sx q[0];
rz(-0.24539465) q[0];
rz(-1.2303526) q[1];
sx q[1];
rz(-3.0407258) q[1];
sx q[1];
rz(0.51602498) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9629899) q[0];
sx q[0];
rz(-1.7699452) q[0];
sx q[0];
rz(0.38802035) q[0];
rz(-0.86271) q[2];
sx q[2];
rz(-2.1893775) q[2];
sx q[2];
rz(-2.1052306) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0210683) q[1];
sx q[1];
rz(-1.85441) q[1];
sx q[1];
rz(0.3801769) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0050659) q[3];
sx q[3];
rz(-1.6322005) q[3];
sx q[3];
rz(1.9207973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.5838991) q[2];
sx q[2];
rz(-1.9198753) q[2];
sx q[2];
rz(-2.058513) q[2];
rz(0.22100581) q[3];
sx q[3];
rz(-2.998304) q[3];
sx q[3];
rz(-2.3451282) q[3];
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
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1204656) q[0];
sx q[0];
rz(-0.080862232) q[0];
sx q[0];
rz(-0.13588151) q[0];
rz(1.7120301) q[1];
sx q[1];
rz(-2.0202961) q[1];
sx q[1];
rz(-0.28724614) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0025072) q[0];
sx q[0];
rz(-1.814743) q[0];
sx q[0];
rz(2.9480145) q[0];
rz(1.2340698) q[2];
sx q[2];
rz(-2.8904466) q[2];
sx q[2];
rz(0.8539356) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.1767856) q[1];
sx q[1];
rz(-2.5075956) q[1];
sx q[1];
rz(0.37668677) q[1];
x q[2];
rz(-0.98669085) q[3];
sx q[3];
rz(-1.3249448) q[3];
sx q[3];
rz(2.33146) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2222432) q[2];
sx q[2];
rz(-1.9025849) q[2];
sx q[2];
rz(2.3075721) q[2];
rz(-2.833994) q[3];
sx q[3];
rz(-0.37720507) q[3];
sx q[3];
rz(-0.26267499) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1766227) q[0];
sx q[0];
rz(-1.3055834) q[0];
sx q[0];
rz(2.1156043) q[0];
rz(2.7550244) q[1];
sx q[1];
rz(-1.3810806) q[1];
sx q[1];
rz(2.5234533) q[1];
rz(2.7170277) q[2];
sx q[2];
rz(-1.4506315) q[2];
sx q[2];
rz(-1.8898024) q[2];
rz(0.97053075) q[3];
sx q[3];
rz(-2.8875792) q[3];
sx q[3];
rz(-2.447788) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
