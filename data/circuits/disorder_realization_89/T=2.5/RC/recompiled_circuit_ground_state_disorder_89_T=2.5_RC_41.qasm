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
rz(-2.3681695) q[0];
rz(3.050488) q[1];
sx q[1];
rz(2.3478822) q[1];
sx q[1];
rz(4.4749727) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42595902) q[0];
sx q[0];
rz(-1.5618526) q[0];
sx q[0];
rz(-0.65837966) q[0];
rz(2.9113468) q[2];
sx q[2];
rz(-2.0587927) q[2];
sx q[2];
rz(-0.56026087) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.7385834) q[1];
sx q[1];
rz(-2.3883551) q[1];
sx q[1];
rz(-0.97767104) q[1];
rz(-pi) q[2];
rz(-3.1245272) q[3];
sx q[3];
rz(-1.6380042) q[3];
sx q[3];
rz(-2.1484571) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.8495463) q[2];
sx q[2];
rz(-1.3214194) q[2];
sx q[2];
rz(0.41479659) q[2];
rz(1.8997806) q[3];
sx q[3];
rz(-2.0262227) q[3];
sx q[3];
rz(2.950086) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5135797) q[0];
sx q[0];
rz(-0.11592557) q[0];
sx q[0];
rz(0.43981788) q[0];
rz(-3.0385333) q[1];
sx q[1];
rz(-0.28918806) q[1];
sx q[1];
rz(-1.1691079) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77020184) q[0];
sx q[0];
rz(-1.4318716) q[0];
sx q[0];
rz(-0.94111218) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1610519) q[2];
sx q[2];
rz(-1.613918) q[2];
sx q[2];
rz(1.0999964) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4354187) q[1];
sx q[1];
rz(-1.7507554) q[1];
sx q[1];
rz(2.2246996) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.86128791) q[3];
sx q[3];
rz(-1.7206186) q[3];
sx q[3];
rz(-0.063171841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.64572) q[2];
sx q[2];
rz(-1.7503259) q[2];
sx q[2];
rz(-0.24620852) q[2];
rz(1.5496893) q[3];
sx q[3];
rz(-0.61384765) q[3];
sx q[3];
rz(1.7483819) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.15568) q[0];
sx q[0];
rz(-2.0014626) q[0];
sx q[0];
rz(-0.21173665) q[0];
rz(-0.011292975) q[1];
sx q[1];
rz(-2.3432422) q[1];
sx q[1];
rz(1.4439772) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2560476) q[0];
sx q[0];
rz(-0.035158947) q[0];
sx q[0];
rz(-0.0020744046) q[0];
rz(-pi) q[1];
x q[1];
rz(0.050758661) q[2];
sx q[2];
rz(-1.3302667) q[2];
sx q[2];
rz(-1.5367791) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.1423751) q[1];
sx q[1];
rz(-0.96227598) q[1];
sx q[1];
rz(1.1097679) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6908247) q[3];
sx q[3];
rz(-1.6524466) q[3];
sx q[3];
rz(-1.9729561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.4766562) q[2];
sx q[2];
rz(-2.0659451) q[2];
sx q[2];
rz(0.54534379) q[2];
rz(2.2319345) q[3];
sx q[3];
rz(-2.6792512) q[3];
sx q[3];
rz(1.2130515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42105168) q[0];
sx q[0];
rz(-2.4626829) q[0];
sx q[0];
rz(-2.3140267) q[0];
rz(-1.1830117) q[1];
sx q[1];
rz(-1.8667826) q[1];
sx q[1];
rz(-1.2659198) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2823845) q[0];
sx q[0];
rz(-0.13656884) q[0];
sx q[0];
rz(0.84117667) q[0];
rz(1.363344) q[2];
sx q[2];
rz(-1.5599164) q[2];
sx q[2];
rz(-1.2258197) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.3862761) q[1];
sx q[1];
rz(-2.1960253) q[1];
sx q[1];
rz(-0.03920725) q[1];
x q[2];
rz(1.7608579) q[3];
sx q[3];
rz(-1.4134839) q[3];
sx q[3];
rz(2.8050131) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.8852692) q[2];
sx q[2];
rz(-1.4726535) q[2];
sx q[2];
rz(-0.081776865) q[2];
rz(-2.8335588) q[3];
sx q[3];
rz(-2.6576198) q[3];
sx q[3];
rz(1.5380194) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9772188) q[0];
sx q[0];
rz(-2.866221) q[0];
sx q[0];
rz(-0.070505738) q[0];
rz(2.8071857) q[1];
sx q[1];
rz(-0.53135482) q[1];
sx q[1];
rz(2.6883584) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8639979) q[0];
sx q[0];
rz(-2.3454614) q[0];
sx q[0];
rz(-0.73426883) q[0];
rz(-2.9137827) q[2];
sx q[2];
rz(-1.9535629) q[2];
sx q[2];
rz(1.5441976) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.2259648) q[1];
sx q[1];
rz(-0.27131042) q[1];
sx q[1];
rz(0.42479272) q[1];
rz(0.86403697) q[3];
sx q[3];
rz(-2.0620769) q[3];
sx q[3];
rz(-1.2934409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.3004904) q[2];
sx q[2];
rz(-1.6359685) q[2];
sx q[2];
rz(2.9312768) q[2];
rz(2.7503843) q[3];
sx q[3];
rz(-2.936383) q[3];
sx q[3];
rz(2.6989663) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6415569) q[0];
sx q[0];
rz(-2.0099202) q[0];
sx q[0];
rz(0.17182194) q[0];
rz(-1.3459282) q[1];
sx q[1];
rz(-2.185952) q[1];
sx q[1];
rz(-1.1489493) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1412107) q[0];
sx q[0];
rz(-3.0297369) q[0];
sx q[0];
rz(3.0684708) q[0];
rz(2.2515772) q[2];
sx q[2];
rz(-2.7820754) q[2];
sx q[2];
rz(0.36717626) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.51275245) q[1];
sx q[1];
rz(-2.1412555) q[1];
sx q[1];
rz(2.9243584) q[1];
x q[2];
rz(2.9385185) q[3];
sx q[3];
rz(-2.0690898) q[3];
sx q[3];
rz(-2.4808675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.12080869) q[2];
sx q[2];
rz(-0.9641996) q[2];
sx q[2];
rz(-2.7210893) q[2];
rz(-1.6546107) q[3];
sx q[3];
rz(-1.1543115) q[3];
sx q[3];
rz(1.9871064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3625951) q[0];
sx q[0];
rz(-1.2263466) q[0];
sx q[0];
rz(2.763789) q[0];
rz(-1.3899577) q[1];
sx q[1];
rz(-1.4327587) q[1];
sx q[1];
rz(-0.43209824) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8209131) q[0];
sx q[0];
rz(-0.28962505) q[0];
sx q[0];
rz(-2.4885213) q[0];
x q[1];
rz(3.0488857) q[2];
sx q[2];
rz(-1.5376629) q[2];
sx q[2];
rz(-2.9404145) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.3687146) q[1];
sx q[1];
rz(-2.0324824) q[1];
sx q[1];
rz(1.5307472) q[1];
x q[2];
rz(-0.89698741) q[3];
sx q[3];
rz(-1.8963061) q[3];
sx q[3];
rz(1.432076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3720588) q[2];
sx q[2];
rz(-2.5199315) q[2];
sx q[2];
rz(0.2335693) q[2];
rz(-1.6893049) q[3];
sx q[3];
rz(-1.4315616) q[3];
sx q[3];
rz(1.5610032) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5643144) q[0];
sx q[0];
rz(-1.0294585) q[0];
sx q[0];
rz(-0.3197864) q[0];
rz(0.86209595) q[1];
sx q[1];
rz(-1.8902706) q[1];
sx q[1];
rz(1.7822942) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0523543) q[0];
sx q[0];
rz(-2.1349094) q[0];
sx q[0];
rz(-2.4648239) q[0];
rz(-0.96085397) q[2];
sx q[2];
rz(-2.8981588) q[2];
sx q[2];
rz(-3.0185901) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.055147) q[1];
sx q[1];
rz(-2.3999955) q[1];
sx q[1];
rz(0.25772175) q[1];
x q[2];
rz(-1.9407746) q[3];
sx q[3];
rz(-2.4741459) q[3];
sx q[3];
rz(0.98946179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.2916145) q[2];
sx q[2];
rz(-2.7225967) q[2];
sx q[2];
rz(0.46763793) q[2];
rz(-2.6575798) q[3];
sx q[3];
rz(-1.7734807) q[3];
sx q[3];
rz(1.8638994) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9676232) q[0];
sx q[0];
rz(-1.4894217) q[0];
sx q[0];
rz(2.0066579) q[0];
rz(-3.0813772) q[1];
sx q[1];
rz(-2.4048012) q[1];
sx q[1];
rz(-1.5312451) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13847199) q[0];
sx q[0];
rz(-2.0299737) q[0];
sx q[0];
rz(-0.91158406) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5217053) q[2];
sx q[2];
rz(-0.7893749) q[2];
sx q[2];
rz(-0.7482341) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.57219124) q[1];
sx q[1];
rz(-1.6680155) q[1];
sx q[1];
rz(-2.8775294) q[1];
rz(-3.1020622) q[3];
sx q[3];
rz(-0.78257221) q[3];
sx q[3];
rz(-1.4631997) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.81426364) q[2];
sx q[2];
rz(-1.9244104) q[2];
sx q[2];
rz(2.4129756) q[2];
rz(-0.21555756) q[3];
sx q[3];
rz(-1.39648) q[3];
sx q[3];
rz(-2.047915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88290596) q[0];
sx q[0];
rz(-3.1104493) q[0];
sx q[0];
rz(0.31797847) q[0];
rz(-1.3975573) q[1];
sx q[1];
rz(-1.8122858) q[1];
sx q[1];
rz(2.8584282) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8737333) q[0];
sx q[0];
rz(-1.6984617) q[0];
sx q[0];
rz(-3.0588849) q[0];
rz(-pi) q[1];
rz(-2.3874902) q[2];
sx q[2];
rz(-1.5148544) q[2];
sx q[2];
rz(-0.90503446) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.6826815) q[1];
sx q[1];
rz(-2.030917) q[1];
sx q[1];
rz(2.3691872) q[1];
x q[2];
rz(-0.092425032) q[3];
sx q[3];
rz(-1.3636949) q[3];
sx q[3];
rz(2.7566119) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.2240923) q[2];
sx q[2];
rz(-1.516569) q[2];
sx q[2];
rz(-3.0523114) q[2];
rz(1.8381522) q[3];
sx q[3];
rz(-0.83804122) q[3];
sx q[3];
rz(2.1681521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3598809) q[0];
sx q[0];
rz(-2.4926873) q[0];
sx q[0];
rz(0.98095184) q[0];
rz(-2.5557062) q[1];
sx q[1];
rz(-1.2955019) q[1];
sx q[1];
rz(-1.6265709) q[1];
rz(-1.3902612) q[2];
sx q[2];
rz(-1.3611887) q[2];
sx q[2];
rz(1.4088189) q[2];
rz(-1.1896776) q[3];
sx q[3];
rz(-1.8101235) q[3];
sx q[3];
rz(2.9265612) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
