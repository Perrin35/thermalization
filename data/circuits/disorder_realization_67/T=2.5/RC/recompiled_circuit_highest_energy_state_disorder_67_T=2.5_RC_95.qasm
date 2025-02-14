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
rz(1.2713852) q[0];
sx q[0];
rz(-0.013590824) q[0];
sx q[0];
rz(3.115227) q[0];
rz(0.68323505) q[1];
sx q[1];
rz(4.9024138) q[1];
sx q[1];
rz(9.496357) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9893104) q[0];
sx q[0];
rz(-1.5610862) q[0];
sx q[0];
rz(-0.025061314) q[0];
rz(-pi) q[1];
rz(-2.6144892) q[2];
sx q[2];
rz(-0.44473916) q[2];
sx q[2];
rz(1.9183733) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.6759807) q[1];
sx q[1];
rz(-0.019400228) q[1];
sx q[1];
rz(1.3127203) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2504225) q[3];
sx q[3];
rz(-1.8719058) q[3];
sx q[3];
rz(1.6973349) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.8611531) q[2];
sx q[2];
rz(-2.4344378) q[2];
sx q[2];
rz(0.46671483) q[2];
rz(0.47433445) q[3];
sx q[3];
rz(-0.021183906) q[3];
sx q[3];
rz(0.02136136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83653432) q[0];
sx q[0];
rz(-2.6439522) q[0];
sx q[0];
rz(-0.019813892) q[0];
rz(1.5927947) q[1];
sx q[1];
rz(-2.9217547) q[1];
sx q[1];
rz(1.6682495) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9813741) q[0];
sx q[0];
rz(-0.88079534) q[0];
sx q[0];
rz(-2.2771858) q[0];
rz(-0.25924087) q[2];
sx q[2];
rz(-1.8564285) q[2];
sx q[2];
rz(-2.223857) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.32343601) q[1];
sx q[1];
rz(-1.4940573) q[1];
sx q[1];
rz(-1.3763001) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.73592379) q[3];
sx q[3];
rz(-0.97565996) q[3];
sx q[3];
rz(1.5335961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.3026498) q[2];
sx q[2];
rz(-0.5984211) q[2];
sx q[2];
rz(1.2897276) q[2];
rz(-1.2368115) q[3];
sx q[3];
rz(-2.808282) q[3];
sx q[3];
rz(-0.70241565) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1692093) q[0];
sx q[0];
rz(-1.9820259) q[0];
sx q[0];
rz(-1.5706536) q[0];
rz(-1.4717357) q[1];
sx q[1];
rz(-1.5262628) q[1];
sx q[1];
rz(-2.6872046) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7189302) q[0];
sx q[0];
rz(-2.5276428) q[0];
sx q[0];
rz(1.3387247) q[0];
rz(-pi) q[1];
rz(1.7128031) q[2];
sx q[2];
rz(-1.5446086) q[2];
sx q[2];
rz(-2.533503) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.4207698) q[1];
sx q[1];
rz(-1.5481564) q[1];
sx q[1];
rz(-1.4525502) q[1];
rz(-pi) q[2];
rz(-0.79409365) q[3];
sx q[3];
rz(-2.3183868) q[3];
sx q[3];
rz(2.2633207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.74139524) q[2];
sx q[2];
rz(-3.1252842) q[2];
sx q[2];
rz(-0.34747094) q[2];
rz(-0.45626429) q[3];
sx q[3];
rz(-3.1268692) q[3];
sx q[3];
rz(-0.99304503) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72661138) q[0];
sx q[0];
rz(-1.9009637) q[0];
sx q[0];
rz(1.4062784) q[0];
rz(-0.44334385) q[1];
sx q[1];
rz(-2.1163546) q[1];
sx q[1];
rz(-1.571507) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9260028) q[0];
sx q[0];
rz(-2.2099751) q[0];
sx q[0];
rz(1.502468) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.072567) q[2];
sx q[2];
rz(-1.6352425) q[2];
sx q[2];
rz(0.2914337) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.0896596) q[1];
sx q[1];
rz(-1.5610236) q[1];
sx q[1];
rz(-1.2925413) q[1];
rz(-2.1277694) q[3];
sx q[3];
rz(-2.644745) q[3];
sx q[3];
rz(0.58409494) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7683679) q[2];
sx q[2];
rz(-0.39373213) q[2];
sx q[2];
rz(3.0623398) q[2];
rz(1.1894038) q[3];
sx q[3];
rz(-1.3711843) q[3];
sx q[3];
rz(-1.6325379) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4382512) q[0];
sx q[0];
rz(-0.51333135) q[0];
sx q[0];
rz(-2.3355423) q[0];
rz(0.84000677) q[1];
sx q[1];
rz(-0.012931074) q[1];
sx q[1];
rz(-0.77617019) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0426725) q[0];
sx q[0];
rz(-1.4399043) q[0];
sx q[0];
rz(-0.58986536) q[0];
rz(-pi) q[1];
rz(3.1311536) q[2];
sx q[2];
rz(-1.5722646) q[2];
sx q[2];
rz(1.2900342) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.35368109) q[1];
sx q[1];
rz(-3.0056142) q[1];
sx q[1];
rz(0.087894364) q[1];
rz(-pi) q[2];
rz(2.8819894) q[3];
sx q[3];
rz(-1.4624551) q[3];
sx q[3];
rz(-1.3110127) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.025803056) q[2];
sx q[2];
rz(-1.5610521) q[2];
sx q[2];
rz(2.4157794) q[2];
rz(-0.18802655) q[3];
sx q[3];
rz(-0.060421061) q[3];
sx q[3];
rz(2.4316725) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1793154) q[0];
sx q[0];
rz(-0.55914068) q[0];
sx q[0];
rz(-2.6002) q[0];
rz(0.18556449) q[1];
sx q[1];
rz(-1.5927529) q[1];
sx q[1];
rz(-3.013179) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.550116) q[0];
sx q[0];
rz(-1.2696365) q[0];
sx q[0];
rz(-0.37806435) q[0];
x q[1];
rz(-1.6096079) q[2];
sx q[2];
rz(-0.12141849) q[2];
sx q[2];
rz(1.5281637) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.7306108) q[1];
sx q[1];
rz(-1.1514947) q[1];
sx q[1];
rz(1.4813408) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5153993) q[3];
sx q[3];
rz(-2.9528816) q[3];
sx q[3];
rz(-0.95079899) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.7554756) q[2];
sx q[2];
rz(-0.057689276) q[2];
sx q[2];
rz(-2.3003787) q[2];
rz(2.9403213) q[3];
sx q[3];
rz(-1.5320675) q[3];
sx q[3];
rz(-2.8819528) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5614618) q[0];
sx q[0];
rz(-0.78829563) q[0];
sx q[0];
rz(1.5808251) q[0];
rz(-2.4687817) q[1];
sx q[1];
rz(-1.4038059) q[1];
sx q[1];
rz(3.1127081) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88850856) q[0];
sx q[0];
rz(-2.7612503) q[0];
sx q[0];
rz(2.4230291) q[0];
x q[1];
rz(-2.7682226) q[2];
sx q[2];
rz(-0.50818102) q[2];
sx q[2];
rz(2.1477107) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.4116482) q[1];
sx q[1];
rz(-1.333451) q[1];
sx q[1];
rz(-1.9945595) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0402803) q[3];
sx q[3];
rz(-2.2757752) q[3];
sx q[3];
rz(-1.353273) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.9081356) q[2];
sx q[2];
rz(-2.5448749) q[2];
sx q[2];
rz(0.92823589) q[2];
rz(0.2975896) q[3];
sx q[3];
rz(-2.987515) q[3];
sx q[3];
rz(-0.84460622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0723648) q[0];
sx q[0];
rz(-0.2437676) q[0];
sx q[0];
rz(3.0425161) q[0];
rz(2.15436) q[1];
sx q[1];
rz(-1.8376553) q[1];
sx q[1];
rz(2.6027021) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19463704) q[0];
sx q[0];
rz(-1.4879003) q[0];
sx q[0];
rz(0.006587365) q[0];
rz(-pi) q[1];
x q[1];
rz(1.107223) q[2];
sx q[2];
rz(-3.0863783) q[2];
sx q[2];
rz(2.5718073) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.0090333825) q[1];
sx q[1];
rz(-1.0485639) q[1];
sx q[1];
rz(-1.9886677) q[1];
rz(-pi) q[2];
rz(1.5032926) q[3];
sx q[3];
rz(-1.5482404) q[3];
sx q[3];
rz(1.8543138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.146356) q[2];
sx q[2];
rz(-0.015492798) q[2];
sx q[2];
rz(-2.8323979) q[2];
rz(-2.501798) q[3];
sx q[3];
rz(-3.1412509) q[3];
sx q[3];
rz(3.0777847) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8362506) q[0];
sx q[0];
rz(-0.59665614) q[0];
sx q[0];
rz(-3.0869361) q[0];
rz(-2.0089741) q[1];
sx q[1];
rz(-1.9839958) q[1];
sx q[1];
rz(-1.8554912) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0564467) q[0];
sx q[0];
rz(-2.9582199) q[0];
sx q[0];
rz(1.3297179) q[0];
x q[1];
rz(-0.80100061) q[2];
sx q[2];
rz(-3.082976) q[2];
sx q[2];
rz(-0.80551565) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.9735582) q[1];
sx q[1];
rz(-1.4806627) q[1];
sx q[1];
rz(0.22613495) q[1];
x q[2];
rz(-1.083582) q[3];
sx q[3];
rz(-1.5069252) q[3];
sx q[3];
rz(-2.4364249) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.22141156) q[2];
sx q[2];
rz(-2.5765918) q[2];
sx q[2];
rz(-2.0378225) q[2];
rz(1.6086027) q[3];
sx q[3];
rz(-0.04403232) q[3];
sx q[3];
rz(0.65582961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.00103818) q[0];
sx q[0];
rz(-2.9619205) q[0];
sx q[0];
rz(3.1371064) q[0];
rz(1.5890315) q[1];
sx q[1];
rz(-1.6922502) q[1];
sx q[1];
rz(-3.0844614) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.321314) q[0];
sx q[0];
rz(-1.6666659) q[0];
sx q[0];
rz(1.3906327) q[0];
rz(-pi) q[1];
x q[1];
rz(0.56978858) q[2];
sx q[2];
rz(-2.8789911) q[2];
sx q[2];
rz(0.71209967) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.88567121) q[1];
sx q[1];
rz(-2.2388865) q[1];
sx q[1];
rz(2.9873965) q[1];
rz(-0.1520429) q[3];
sx q[3];
rz(-0.83629823) q[3];
sx q[3];
rz(-1.1610069) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.39157465) q[2];
sx q[2];
rz(-3.1142758) q[2];
sx q[2];
rz(-2.2726783) q[2];
rz(0.9817552) q[3];
sx q[3];
rz(-3.1120286) q[3];
sx q[3];
rz(2.6386007) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0960196) q[0];
sx q[0];
rz(-1.6572784) q[0];
sx q[0];
rz(-1.4838765) q[0];
rz(-0.45687301) q[1];
sx q[1];
rz(-0.15468205) q[1];
sx q[1];
rz(-0.044943132) q[1];
rz(-0.76404608) q[2];
sx q[2];
rz(-1.7052787) q[2];
sx q[2];
rz(1.6324001) q[2];
rz(-0.66309912) q[3];
sx q[3];
rz(-1.6840006) q[3];
sx q[3];
rz(-1.4406289) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
