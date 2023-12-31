OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.47473946) q[0];
sx q[0];
rz(-0.82959509) q[0];
sx q[0];
rz(-2.9876246) q[0];
rz(0.83377588) q[1];
sx q[1];
rz(-2.1492465) q[1];
sx q[1];
rz(-0.33831236) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96903893) q[0];
sx q[0];
rz(-1.2533029) q[0];
sx q[0];
rz(0.25704268) q[0];
rz(-1.1233653) q[2];
sx q[2];
rz(-0.72927232) q[2];
sx q[2];
rz(-0.73993081) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.073804341) q[1];
sx q[1];
rz(-1.6238302) q[1];
sx q[1];
rz(0.28633134) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.1258718) q[3];
sx q[3];
rz(-2.0831046) q[3];
sx q[3];
rz(-0.44954625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.14264318) q[2];
sx q[2];
rz(-0.34037408) q[2];
sx q[2];
rz(-1.1738698) q[2];
rz(3.0657892) q[3];
sx q[3];
rz(-1.9971763) q[3];
sx q[3];
rz(0.092806667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8006111) q[0];
sx q[0];
rz(-1.0656463) q[0];
sx q[0];
rz(0.064963438) q[0];
rz(-0.57463542) q[1];
sx q[1];
rz(-0.42962933) q[1];
sx q[1];
rz(1.8992791) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.500538) q[0];
sx q[0];
rz(-0.28028742) q[0];
sx q[0];
rz(2.6602402) q[0];
rz(-2.20349) q[2];
sx q[2];
rz(-1.6993076) q[2];
sx q[2];
rz(-3.1046257) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.6537135) q[1];
sx q[1];
rz(-2.5839845) q[1];
sx q[1];
rz(2.8370268) q[1];
rz(-2.4317125) q[3];
sx q[3];
rz(-1.1871561) q[3];
sx q[3];
rz(0.82304728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.80766455) q[2];
sx q[2];
rz(-2.0662722) q[2];
sx q[2];
rz(-2.5578965) q[2];
rz(-0.57404533) q[3];
sx q[3];
rz(-1.1255001) q[3];
sx q[3];
rz(-3.0103502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7211001) q[0];
sx q[0];
rz(-0.89389602) q[0];
sx q[0];
rz(2.4131391) q[0];
rz(1.6473673) q[1];
sx q[1];
rz(-0.39847001) q[1];
sx q[1];
rz(-1.0167936) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8249614) q[0];
sx q[0];
rz(-1.6065856) q[0];
sx q[0];
rz(0.081598452) q[0];
x q[1];
rz(1.0810034) q[2];
sx q[2];
rz(-1.0989597) q[2];
sx q[2];
rz(-2.8420198) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6779855) q[1];
sx q[1];
rz(-1.4336168) q[1];
sx q[1];
rz(-2.3198818) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4025027) q[3];
sx q[3];
rz(-2.8561391) q[3];
sx q[3];
rz(-1.9608378) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.8016522) q[2];
sx q[2];
rz(-1.6093971) q[2];
sx q[2];
rz(-1.0495079) q[2];
rz(-2.5028051) q[3];
sx q[3];
rz(-2.5103266) q[3];
sx q[3];
rz(-1.9558186) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3291572) q[0];
sx q[0];
rz(-1.8742467) q[0];
sx q[0];
rz(1.6695492) q[0];
rz(-0.73515785) q[1];
sx q[1];
rz(-0.77886326) q[1];
sx q[1];
rz(-2.8947815) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26537672) q[0];
sx q[0];
rz(-0.88337684) q[0];
sx q[0];
rz(2.9486297) q[0];
x q[1];
rz(-1.7814126) q[2];
sx q[2];
rz(-1.1014551) q[2];
sx q[2];
rz(2.5460555) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.88627316) q[1];
sx q[1];
rz(-1.5181932) q[1];
sx q[1];
rz(-2.7707151) q[1];
x q[2];
rz(-2.187192) q[3];
sx q[3];
rz(-1.3351001) q[3];
sx q[3];
rz(2.1637722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.4776769) q[2];
sx q[2];
rz(-2.0229979) q[2];
sx q[2];
rz(1.5412615) q[2];
rz(-0.70704308) q[3];
sx q[3];
rz(-1.0975081) q[3];
sx q[3];
rz(-2.2533806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33070579) q[0];
sx q[0];
rz(-2.4208477) q[0];
sx q[0];
rz(1.8141618) q[0];
rz(1.56303) q[1];
sx q[1];
rz(-2.6674318) q[1];
sx q[1];
rz(-0.24838233) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7029593) q[0];
sx q[0];
rz(-2.6002433) q[0];
sx q[0];
rz(-0.066141733) q[0];
rz(-pi) q[1];
rz(-2.4400473) q[2];
sx q[2];
rz(-0.24214673) q[2];
sx q[2];
rz(-2.0602351) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.69265428) q[1];
sx q[1];
rz(-0.9874953) q[1];
sx q[1];
rz(1.3650236) q[1];
rz(-pi) q[2];
rz(1.3316657) q[3];
sx q[3];
rz(-2.2056747) q[3];
sx q[3];
rz(-2.2374416) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.43626943) q[2];
sx q[2];
rz(-2.1537809) q[2];
sx q[2];
rz(-0.79745897) q[2];
rz(-2.752839) q[3];
sx q[3];
rz(-2.5377486) q[3];
sx q[3];
rz(-2.6388772) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1335063) q[0];
sx q[0];
rz(-3.0651423) q[0];
sx q[0];
rz(-1.7957934) q[0];
rz(-1.0812409) q[1];
sx q[1];
rz(-1.2370279) q[1];
sx q[1];
rz(-0.1246917) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9799177) q[0];
sx q[0];
rz(-1.7470164) q[0];
sx q[0];
rz(-1.6582703) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5325534) q[2];
sx q[2];
rz(-1.8134724) q[2];
sx q[2];
rz(-0.86953029) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.60285073) q[1];
sx q[1];
rz(-2.3067143) q[1];
sx q[1];
rz(-1.49453) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.74495875) q[3];
sx q[3];
rz(-0.31674851) q[3];
sx q[3];
rz(-1.484364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.51320118) q[2];
sx q[2];
rz(-1.1867563) q[2];
sx q[2];
rz(-2.52264) q[2];
rz(2.0882873) q[3];
sx q[3];
rz(-2.9635933) q[3];
sx q[3];
rz(-0.9296023) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58182794) q[0];
sx q[0];
rz(-1.7511837) q[0];
sx q[0];
rz(1.0429617) q[0];
rz(2.6783121) q[1];
sx q[1];
rz(-2.0279341) q[1];
sx q[1];
rz(-1.0707062) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2830414) q[0];
sx q[0];
rz(-1.0867449) q[0];
sx q[0];
rz(-1.5714684) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7037017) q[2];
sx q[2];
rz(-1.9153567) q[2];
sx q[2];
rz(0.29495707) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(3.0106525) q[1];
sx q[1];
rz(-0.7796692) q[1];
sx q[1];
rz(-0.14203771) q[1];
rz(3.0827818) q[3];
sx q[3];
rz(-2.8088514) q[3];
sx q[3];
rz(0.22418338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.36859194) q[2];
sx q[2];
rz(-0.7395145) q[2];
sx q[2];
rz(-2.8179742) q[2];
rz(2.1598024) q[3];
sx q[3];
rz(-0.86172813) q[3];
sx q[3];
rz(1.0872844) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6535646) q[0];
sx q[0];
rz(-1.7249148) q[0];
sx q[0];
rz(1.9352242) q[0];
rz(1.2127097) q[1];
sx q[1];
rz(-0.85314631) q[1];
sx q[1];
rz(-2.1941197) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8038717) q[0];
sx q[0];
rz(-2.5944355) q[0];
sx q[0];
rz(-1.7659811) q[0];
rz(-1.4119554) q[2];
sx q[2];
rz(-1.7454595) q[2];
sx q[2];
rz(0.58089248) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.265043) q[1];
sx q[1];
rz(-1.8076841) q[1];
sx q[1];
rz(-2.151728) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2097589) q[3];
sx q[3];
rz(-2.3449538) q[3];
sx q[3];
rz(0.79209298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.56132135) q[2];
sx q[2];
rz(-1.7611971) q[2];
sx q[2];
rz(-0.46978152) q[2];
rz(1.3011159) q[3];
sx q[3];
rz(-1.4353292) q[3];
sx q[3];
rz(0.27967134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25205055) q[0];
sx q[0];
rz(-2.7520576) q[0];
sx q[0];
rz(-1.8126194) q[0];
rz(2.3503616) q[1];
sx q[1];
rz(-0.33114094) q[1];
sx q[1];
rz(-2.9387617) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3644667) q[0];
sx q[0];
rz(-1.5317894) q[0];
sx q[0];
rz(1.5968621) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9183396) q[2];
sx q[2];
rz(-1.3666144) q[2];
sx q[2];
rz(0.60165652) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.087346615) q[1];
sx q[1];
rz(-0.11212238) q[1];
sx q[1];
rz(-2.0263158) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4167452) q[3];
sx q[3];
rz(-1.1136347) q[3];
sx q[3];
rz(-0.16897133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7245076) q[2];
sx q[2];
rz(-0.26792002) q[2];
sx q[2];
rz(0.99651304) q[2];
rz(0.35342446) q[3];
sx q[3];
rz(-0.74520183) q[3];
sx q[3];
rz(-0.80741185) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.080169454) q[0];
sx q[0];
rz(-2.3286979) q[0];
sx q[0];
rz(0.18173519) q[0];
rz(3.0985447) q[1];
sx q[1];
rz(-2.4964066) q[1];
sx q[1];
rz(2.8607686) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40598665) q[0];
sx q[0];
rz(-1.9970903) q[0];
sx q[0];
rz(2.6398525) q[0];
rz(-0.73239399) q[2];
sx q[2];
rz(-2.8576982) q[2];
sx q[2];
rz(0.64901272) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.2320497) q[1];
sx q[1];
rz(-3.0735525) q[1];
sx q[1];
rz(-1.0104936) q[1];
x q[2];
rz(-2.9452768) q[3];
sx q[3];
rz(-1.2442949) q[3];
sx q[3];
rz(1.452009) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3165555) q[2];
sx q[2];
rz(-1.8871769) q[2];
sx q[2];
rz(0.62310702) q[2];
rz(-1.0021707) q[3];
sx q[3];
rz(-1.3391756) q[3];
sx q[3];
rz(0.56308693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6476718) q[0];
sx q[0];
rz(-1.5734084) q[0];
sx q[0];
rz(-1.5403803) q[0];
rz(-2.2676246) q[1];
sx q[1];
rz(-1.0653492) q[1];
sx q[1];
rz(-3.031562) q[1];
rz(2.3360844) q[2];
sx q[2];
rz(-1.1136354) q[2];
sx q[2];
rz(-1.8352933) q[2];
rz(-1.1602041) q[3];
sx q[3];
rz(-1.2413597) q[3];
sx q[3];
rz(1.6551457) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
