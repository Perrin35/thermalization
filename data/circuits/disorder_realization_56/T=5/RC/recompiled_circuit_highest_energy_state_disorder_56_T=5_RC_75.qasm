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
rz(1.4752969) q[0];
sx q[0];
rz(-1.2694321) q[0];
sx q[0];
rz(0.67212927) q[0];
rz(-2.693306) q[1];
sx q[1];
rz(-1.6192133) q[1];
sx q[1];
rz(-0.7575922) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2494932) q[0];
sx q[0];
rz(-1.6803015) q[0];
sx q[0];
rz(0.076105781) q[0];
x q[1];
rz(2.5301928) q[2];
sx q[2];
rz(-0.58908236) q[2];
sx q[2];
rz(0.13797274) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.5930773) q[1];
sx q[1];
rz(-2.6679949) q[1];
sx q[1];
rz(0.51215902) q[1];
rz(2.1610307) q[3];
sx q[3];
rz(-0.96348982) q[3];
sx q[3];
rz(3.0890298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.0483094) q[2];
sx q[2];
rz(-1.8123241) q[2];
sx q[2];
rz(1.7993571) q[2];
rz(-1.7736769) q[3];
sx q[3];
rz(-1.0932837) q[3];
sx q[3];
rz(-1.5526519) q[3];
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
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9369478) q[0];
sx q[0];
rz(-2.1277675) q[0];
sx q[0];
rz(-2.192705) q[0];
rz(0.10143796) q[1];
sx q[1];
rz(-1.0600435) q[1];
sx q[1];
rz(0.92996517) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0596493) q[0];
sx q[0];
rz(-0.23632061) q[0];
sx q[0];
rz(-0.44775072) q[0];
rz(-pi) q[1];
rz(-3.1279966) q[2];
sx q[2];
rz(-2.6363723) q[2];
sx q[2];
rz(1.9590953) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.2163712) q[1];
sx q[1];
rz(-2.9489046) q[1];
sx q[1];
rz(0.57740258) q[1];
rz(-2.3306777) q[3];
sx q[3];
rz(-1.0958015) q[3];
sx q[3];
rz(0.37372933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.0019504) q[2];
sx q[2];
rz(-1.6763687) q[2];
sx q[2];
rz(-0.742221) q[2];
rz(-0.46418515) q[3];
sx q[3];
rz(-1.3451385) q[3];
sx q[3];
rz(-1.5531042) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4699698) q[0];
sx q[0];
rz(-2.9556584) q[0];
sx q[0];
rz(-1.6999014) q[0];
rz(-0.16547671) q[1];
sx q[1];
rz(-2.4022357) q[1];
sx q[1];
rz(-0.86732078) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21737032) q[0];
sx q[0];
rz(-1.5708367) q[0];
sx q[0];
rz(1.5702308) q[0];
x q[1];
rz(1.1050842) q[2];
sx q[2];
rz(-2.5869479) q[2];
sx q[2];
rz(1.1797734) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.2677541) q[1];
sx q[1];
rz(-1.8513362) q[1];
sx q[1];
rz(-2.8024017) q[1];
rz(-0.8115143) q[3];
sx q[3];
rz(-1.1574739) q[3];
sx q[3];
rz(-2.1275525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1778339) q[2];
sx q[2];
rz(-1.9671665) q[2];
sx q[2];
rz(0.7589232) q[2];
rz(0.80398503) q[3];
sx q[3];
rz(-2.1920429) q[3];
sx q[3];
rz(2.4132531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3114965) q[0];
sx q[0];
rz(-1.8599956) q[0];
sx q[0];
rz(0.48402825) q[0];
rz(1.6054035) q[1];
sx q[1];
rz(-1.7264629) q[1];
sx q[1];
rz(-1.7162292) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7361476) q[0];
sx q[0];
rz(-2.510294) q[0];
sx q[0];
rz(2.5456356) q[0];
x q[1];
rz(-2.8292848) q[2];
sx q[2];
rz(-1.3538401) q[2];
sx q[2];
rz(-1.7157451) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.0850544) q[1];
sx q[1];
rz(-0.94496545) q[1];
sx q[1];
rz(-0.27295785) q[1];
rz(-pi) q[2];
rz(0.45141545) q[3];
sx q[3];
rz(-1.2299996) q[3];
sx q[3];
rz(1.4489074) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.75590762) q[2];
sx q[2];
rz(-1.0681095) q[2];
sx q[2];
rz(2.8483086) q[2];
rz(3.0934603) q[3];
sx q[3];
rz(-1.3061413) q[3];
sx q[3];
rz(0.3046681) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3047979) q[0];
sx q[0];
rz(-1.7339107) q[0];
sx q[0];
rz(2.0378713) q[0];
rz(-2.2301105) q[1];
sx q[1];
rz(-1.5244923) q[1];
sx q[1];
rz(0.10890659) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86065642) q[0];
sx q[0];
rz(-2.0514279) q[0];
sx q[0];
rz(-2.3970766) q[0];
rz(-pi) q[1];
x q[1];
rz(0.53498603) q[2];
sx q[2];
rz(-1.4145383) q[2];
sx q[2];
rz(-1.6339169) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.8629294) q[1];
sx q[1];
rz(-0.68611523) q[1];
sx q[1];
rz(1.8300232) q[1];
x q[2];
rz(1.2690684) q[3];
sx q[3];
rz(-1.0142039) q[3];
sx q[3];
rz(-2.3315786) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.4917422) q[2];
sx q[2];
rz(-1.3567341) q[2];
sx q[2];
rz(-0.25807992) q[2];
rz(0.38170013) q[3];
sx q[3];
rz(-0.93770599) q[3];
sx q[3];
rz(2.5842353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3968762) q[0];
sx q[0];
rz(-0.98137403) q[0];
sx q[0];
rz(1.9816403) q[0];
rz(-0.57811919) q[1];
sx q[1];
rz(-1.4719529) q[1];
sx q[1];
rz(-2.8181308) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1581381) q[0];
sx q[0];
rz(-1.1295415) q[0];
sx q[0];
rz(0.47636124) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6450364) q[2];
sx q[2];
rz(-1.7944031) q[2];
sx q[2];
rz(-1.3204492) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.0486794) q[1];
sx q[1];
rz(-1.7968751) q[1];
sx q[1];
rz(2.5552555) q[1];
x q[2];
rz(2.6810535) q[3];
sx q[3];
rz(-1.8611188) q[3];
sx q[3];
rz(3.0825305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.20299992) q[2];
sx q[2];
rz(-2.3679569) q[2];
sx q[2];
rz(2.2933551) q[2];
rz(1.8741459) q[3];
sx q[3];
rz(-1.7402612) q[3];
sx q[3];
rz(-0.69376865) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11442014) q[0];
sx q[0];
rz(-2.4878451) q[0];
sx q[0];
rz(-0.87345901) q[0];
rz(1.9050441) q[1];
sx q[1];
rz(-2.1658587) q[1];
sx q[1];
rz(3.0580318) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.75514166) q[0];
sx q[0];
rz(-1.6928732) q[0];
sx q[0];
rz(2.8194619) q[0];
rz(-pi) q[1];
rz(-1.0413076) q[2];
sx q[2];
rz(-1.7779967) q[2];
sx q[2];
rz(2.4251314) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8712743) q[1];
sx q[1];
rz(-1.9650998) q[1];
sx q[1];
rz(-0.81812268) q[1];
rz(-pi) q[2];
rz(-0.6912937) q[3];
sx q[3];
rz(-0.45540998) q[3];
sx q[3];
rz(0.56870715) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.7053335) q[2];
sx q[2];
rz(-1.6062364) q[2];
sx q[2];
rz(-1.7748888) q[2];
rz(-0.27240917) q[3];
sx q[3];
rz(-2.1209769) q[3];
sx q[3];
rz(2.3830856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(0.90683872) q[0];
sx q[0];
rz(-0.82643569) q[0];
sx q[0];
rz(-0.93836623) q[0];
rz(-0.80288184) q[1];
sx q[1];
rz(-2.1736841) q[1];
sx q[1];
rz(-2.3209007) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9296449) q[0];
sx q[0];
rz(-1.4370455) q[0];
sx q[0];
rz(1.9424214) q[0];
x q[1];
rz(1.9782045) q[2];
sx q[2];
rz(-2.738224) q[2];
sx q[2];
rz(-0.58248108) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.31466111) q[1];
sx q[1];
rz(-0.95035997) q[1];
sx q[1];
rz(-0.13038306) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9441767) q[3];
sx q[3];
rz(-2.7146517) q[3];
sx q[3];
rz(-2.3641242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.20251033) q[2];
sx q[2];
rz(-0.15190092) q[2];
sx q[2];
rz(-2.0738156) q[2];
rz(-1.1311401) q[3];
sx q[3];
rz(-0.83429566) q[3];
sx q[3];
rz(-0.080032674) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.985567) q[0];
sx q[0];
rz(-1.9859059) q[0];
sx q[0];
rz(2.735403) q[0];
rz(1.4093026) q[1];
sx q[1];
rz(-2.1235762) q[1];
sx q[1];
rz(-1.0850151) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70124088) q[0];
sx q[0];
rz(-0.95537649) q[0];
sx q[0];
rz(-1.3568272) q[0];
x q[1];
rz(-3.135097) q[2];
sx q[2];
rz(-0.61031872) q[2];
sx q[2];
rz(1.6808866) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.4915732) q[1];
sx q[1];
rz(-2.2597092) q[1];
sx q[1];
rz(0.42999173) q[1];
rz(-pi) q[2];
rz(-2.3906443) q[3];
sx q[3];
rz(-2.4118858) q[3];
sx q[3];
rz(-2.2580604) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.5558418) q[2];
sx q[2];
rz(-2.5297574) q[2];
sx q[2];
rz(0.93070585) q[2];
rz(1.7454923) q[3];
sx q[3];
rz(-1.7788818) q[3];
sx q[3];
rz(-3.0220368) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68468204) q[0];
sx q[0];
rz(-0.97701183) q[0];
sx q[0];
rz(-1.8344301) q[0];
rz(0.3859418) q[1];
sx q[1];
rz(-1.394505) q[1];
sx q[1];
rz(1.0345667) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32337727) q[0];
sx q[0];
rz(-1.3971097) q[0];
sx q[0];
rz(-3.0370569) q[0];
rz(-pi) q[1];
rz(-2.7427086) q[2];
sx q[2];
rz(-0.49408696) q[2];
sx q[2];
rz(1.2436155) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9683796) q[1];
sx q[1];
rz(-0.82260859) q[1];
sx q[1];
rz(2.2065225) q[1];
rz(-pi) q[2];
x q[2];
rz(0.8108658) q[3];
sx q[3];
rz(-1.4812143) q[3];
sx q[3];
rz(0.60096622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.6941541) q[2];
sx q[2];
rz(-2.4093781) q[2];
sx q[2];
rz(0.36699692) q[2];
rz(2.0224109) q[3];
sx q[3];
rz(-0.39803353) q[3];
sx q[3];
rz(1.6163577) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.781453) q[0];
sx q[0];
rz(-1.9522788) q[0];
sx q[0];
rz(2.0576394) q[0];
rz(-3.0837334) q[1];
sx q[1];
rz(-1.2670988) q[1];
sx q[1];
rz(1.7115464) q[1];
rz(-2.62769) q[2];
sx q[2];
rz(-2.9436417) q[2];
sx q[2];
rz(-1.8463989) q[2];
rz(-0.41856159) q[3];
sx q[3];
rz(-0.7444612) q[3];
sx q[3];
rz(0.85519467) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
