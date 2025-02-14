OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.9417579) q[0];
sx q[0];
rz(-2.4500442) q[0];
sx q[0];
rz(0.77342311) q[0];
rz(3.050488) q[1];
sx q[1];
rz(2.3478822) q[1];
sx q[1];
rz(4.4749727) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9851935) q[0];
sx q[0];
rz(-0.65843136) q[0];
sx q[0];
rz(-0.014617217) q[0];
x q[1];
rz(1.0716419) q[2];
sx q[2];
rz(-1.7737651) q[2];
sx q[2];
rz(-0.90107337) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.65705196) q[1];
sx q[1];
rz(-0.96772497) q[1];
sx q[1];
rz(-2.6588427) q[1];
rz(1.5035787) q[3];
sx q[3];
rz(-1.5537694) q[3];
sx q[3];
rz(0.57651455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.29204631) q[2];
sx q[2];
rz(-1.8201733) q[2];
sx q[2];
rz(-0.41479659) q[2];
rz(-1.8997806) q[3];
sx q[3];
rz(-1.1153699) q[3];
sx q[3];
rz(2.950086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5135797) q[0];
sx q[0];
rz(-0.11592557) q[0];
sx q[0];
rz(-2.7017748) q[0];
rz(-0.10305931) q[1];
sx q[1];
rz(-0.28918806) q[1];
sx q[1];
rz(-1.9724847) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4415574) q[0];
sx q[0];
rz(-0.9481144) q[0];
sx q[0];
rz(0.17130987) q[0];
x q[1];
rz(-2.1610519) q[2];
sx q[2];
rz(-1.5276747) q[2];
sx q[2];
rz(2.0415963) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.00094571908) q[1];
sx q[1];
rz(-0.92921153) q[1];
sx q[1];
rz(-2.916275) q[1];
x q[2];
rz(-1.7984932) q[3];
sx q[3];
rz(-2.4191354) q[3];
sx q[3];
rz(1.8060773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.49587265) q[2];
sx q[2];
rz(-1.3912667) q[2];
sx q[2];
rz(-2.8953841) q[2];
rz(-1.5496893) q[3];
sx q[3];
rz(-0.61384765) q[3];
sx q[3];
rz(-1.7483819) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98591268) q[0];
sx q[0];
rz(-1.14013) q[0];
sx q[0];
rz(-0.21173665) q[0];
rz(0.011292975) q[1];
sx q[1];
rz(-0.79835049) q[1];
sx q[1];
rz(-1.6976154) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31682184) q[0];
sx q[0];
rz(-1.5708692) q[0];
sx q[0];
rz(0.035158872) q[0];
rz(3.090834) q[2];
sx q[2];
rz(-1.3302667) q[2];
sx q[2];
rz(1.5367791) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.1423751) q[1];
sx q[1];
rz(-0.96227598) q[1];
sx q[1];
rz(-2.0318248) q[1];
rz(-pi) q[2];
rz(-2.6908247) q[3];
sx q[3];
rz(-1.6524466) q[3];
sx q[3];
rz(-1.9729561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.4766562) q[2];
sx q[2];
rz(-1.0756476) q[2];
sx q[2];
rz(-0.54534379) q[2];
rz(-0.90965811) q[3];
sx q[3];
rz(-2.6792512) q[3];
sx q[3];
rz(-1.9285412) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.720541) q[0];
sx q[0];
rz(-0.67890972) q[0];
sx q[0];
rz(-2.3140267) q[0];
rz(-1.958581) q[1];
sx q[1];
rz(-1.8667826) q[1];
sx q[1];
rz(1.2659198) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0166548) q[0];
sx q[0];
rz(-1.6724574) q[0];
sx q[0];
rz(-3.0502425) q[0];
x q[1];
rz(-3.1304743) q[2];
sx q[2];
rz(-1.3633565) q[2];
sx q[2];
rz(0.34726663) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.3490679) q[1];
sx q[1];
rz(-1.5390087) q[1];
sx q[1];
rz(-2.1963901) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2692415) q[3];
sx q[3];
rz(-2.8954801) q[3];
sx q[3];
rz(1.2239561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.8852692) q[2];
sx q[2];
rz(-1.4726535) q[2];
sx q[2];
rz(0.081776865) q[2];
rz(2.8335588) q[3];
sx q[3];
rz(-2.6576198) q[3];
sx q[3];
rz(-1.5380194) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1643739) q[0];
sx q[0];
rz(-0.27537167) q[0];
sx q[0];
rz(0.070505738) q[0];
rz(-2.8071857) q[1];
sx q[1];
rz(-2.6102378) q[1];
sx q[1];
rz(-0.45323429) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8716259) q[0];
sx q[0];
rz(-2.0701417) q[0];
sx q[0];
rz(0.64887394) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9137827) q[2];
sx q[2];
rz(-1.1880298) q[2];
sx q[2];
rz(-1.5441976) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.6649225) q[1];
sx q[1];
rz(-1.3241265) q[1];
sx q[1];
rz(1.4566521) q[1];
rz(-pi) q[2];
rz(2.5284834) q[3];
sx q[3];
rz(-0.96115784) q[3];
sx q[3];
rz(-0.10559374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.8411023) q[2];
sx q[2];
rz(-1.6359685) q[2];
sx q[2];
rz(0.21031586) q[2];
rz(2.7503843) q[3];
sx q[3];
rz(-2.936383) q[3];
sx q[3];
rz(2.6989663) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5000358) q[0];
sx q[0];
rz(-2.0099202) q[0];
sx q[0];
rz(-2.9697707) q[0];
rz(1.3459282) q[1];
sx q[1];
rz(-0.95564061) q[1];
sx q[1];
rz(1.9926434) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35691914) q[0];
sx q[0];
rz(-1.5789512) q[0];
sx q[0];
rz(-0.11155931) q[0];
x q[1];
rz(0.8900155) q[2];
sx q[2];
rz(-2.7820754) q[2];
sx q[2];
rz(2.7744164) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.51275245) q[1];
sx q[1];
rz(-1.0003371) q[1];
sx q[1];
rz(-2.9243584) q[1];
x q[2];
rz(-1.2158209) q[3];
sx q[3];
rz(-2.6067511) q[3];
sx q[3];
rz(-2.8876784) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.12080869) q[2];
sx q[2];
rz(-0.9641996) q[2];
sx q[2];
rz(2.7210893) q[2];
rz(-1.6546107) q[3];
sx q[3];
rz(-1.1543115) q[3];
sx q[3];
rz(-1.1544863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3625951) q[0];
sx q[0];
rz(-1.2263466) q[0];
sx q[0];
rz(0.37780365) q[0];
rz(-1.3899577) q[1];
sx q[1];
rz(-1.4327587) q[1];
sx q[1];
rz(2.7094944) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88277757) q[0];
sx q[0];
rz(-1.3963789) q[0];
sx q[0];
rz(2.9091878) q[0];
rz(-pi) q[1];
rz(2.7977706) q[2];
sx q[2];
rz(-3.0431586) q[2];
sx q[2];
rz(1.7119031) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.3258265) q[1];
sx q[1];
rz(-1.6066505) q[1];
sx q[1];
rz(-2.6795866) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7338846) q[3];
sx q[3];
rz(-2.2033436) q[3];
sx q[3];
rz(3.0303101) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.3720588) q[2];
sx q[2];
rz(-2.5199315) q[2];
sx q[2];
rz(-2.9080234) q[2];
rz(-1.4522878) q[3];
sx q[3];
rz(-1.7100311) q[3];
sx q[3];
rz(-1.5805894) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5643144) q[0];
sx q[0];
rz(-2.1121341) q[0];
sx q[0];
rz(2.8218063) q[0];
rz(-2.2794967) q[1];
sx q[1];
rz(-1.2513221) q[1];
sx q[1];
rz(1.3592985) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0687953) q[0];
sx q[0];
rz(-2.2900343) q[0];
sx q[0];
rz(2.3510975) q[0];
rz(-pi) q[1];
rz(2.1807387) q[2];
sx q[2];
rz(-2.8981588) q[2];
sx q[2];
rz(-3.0185901) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.7118468) q[1];
sx q[1];
rz(-2.2825438) q[1];
sx q[1];
rz(-1.3414308) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8639761) q[3];
sx q[3];
rz(-0.95562387) q[3];
sx q[3];
rz(2.6108116) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.2916145) q[2];
sx q[2];
rz(-2.7225967) q[2];
sx q[2];
rz(-0.46763793) q[2];
rz(-0.48401287) q[3];
sx q[3];
rz(-1.368112) q[3];
sx q[3];
rz(1.8638994) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17396946) q[0];
sx q[0];
rz(-1.4894217) q[0];
sx q[0];
rz(-2.0066579) q[0];
rz(0.060215503) q[1];
sx q[1];
rz(-2.4048012) q[1];
sx q[1];
rz(1.6103475) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0400684) q[0];
sx q[0];
rz(-2.1520242) q[0];
sx q[0];
rz(-2.5826497) q[0];
x q[1];
rz(-0.68708276) q[2];
sx q[2];
rz(-1.1456881) q[2];
sx q[2];
rz(-1.8532422) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.3432777) q[1];
sx q[1];
rz(-0.28099842) q[1];
sx q[1];
rz(-0.35759371) q[1];
x q[2];
rz(-0.78218145) q[3];
sx q[3];
rz(-1.5986658) q[3];
sx q[3];
rz(-3.0620344) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.327329) q[2];
sx q[2];
rz(-1.9244104) q[2];
sx q[2];
rz(-2.4129756) q[2];
rz(-0.21555756) q[3];
sx q[3];
rz(-1.7451127) q[3];
sx q[3];
rz(-1.0936776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2586867) q[0];
sx q[0];
rz(-3.1104493) q[0];
sx q[0];
rz(0.31797847) q[0];
rz(1.7440354) q[1];
sx q[1];
rz(-1.3293068) q[1];
sx q[1];
rz(0.2831645) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8737333) q[0];
sx q[0];
rz(-1.443131) q[0];
sx q[0];
rz(-3.0588849) q[0];
rz(-pi) q[1];
rz(-1.4941169) q[2];
sx q[2];
rz(-0.81816219) q[2];
sx q[2];
rz(-0.61329182) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.4589112) q[1];
sx q[1];
rz(-2.030917) q[1];
sx q[1];
rz(2.3691872) q[1];
rz(-pi) q[2];
rz(1.1569144) q[3];
sx q[3];
rz(-0.22651895) q[3];
sx q[3];
rz(3.1030797) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.91750034) q[2];
sx q[2];
rz(-1.516569) q[2];
sx q[2];
rz(3.0523114) q[2];
rz(-1.8381522) q[3];
sx q[3];
rz(-2.3035514) q[3];
sx q[3];
rz(-0.97344056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3598809) q[0];
sx q[0];
rz(-0.64890535) q[0];
sx q[0];
rz(-2.1606408) q[0];
rz(2.5557062) q[1];
sx q[1];
rz(-1.8460907) q[1];
sx q[1];
rz(1.5150217) q[1];
rz(2.9286251) q[2];
sx q[2];
rz(-1.3942547) q[2];
sx q[2];
rz(-0.19993776) q[2];
rz(2.8845434) q[3];
sx q[3];
rz(-1.9405196) q[3];
sx q[3];
rz(-1.6911239) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
