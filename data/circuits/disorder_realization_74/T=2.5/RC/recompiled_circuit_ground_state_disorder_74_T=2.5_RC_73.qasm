OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.7862406) q[0];
sx q[0];
rz(-2.8347637) q[0];
sx q[0];
rz(2.8428349) q[0];
rz(-0.6660676) q[1];
sx q[1];
rz(2.5112285) q[1];
sx q[1];
rz(11.093333) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.436384) q[0];
sx q[0];
rz(-1.602631) q[0];
sx q[0];
rz(-2.3408606) q[0];
x q[1];
rz(0.017919964) q[2];
sx q[2];
rz(-1.2700873) q[2];
sx q[2];
rz(1.199388) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.7466429) q[1];
sx q[1];
rz(-2.267845) q[1];
sx q[1];
rz(-0.66956981) q[1];
rz(-2.8979635) q[3];
sx q[3];
rz(-0.28377658) q[3];
sx q[3];
rz(0.42967969) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.1851958) q[2];
sx q[2];
rz(-0.38072017) q[2];
sx q[2];
rz(1.3107276) q[2];
rz(-0.091536097) q[3];
sx q[3];
rz(-0.59713489) q[3];
sx q[3];
rz(-2.6105647) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
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
rz(-0.071534261) q[0];
sx q[0];
rz(-2.0508843) q[0];
sx q[0];
rz(2.6614905) q[0];
rz(1.9942888) q[1];
sx q[1];
rz(-2.5105748) q[1];
sx q[1];
rz(2.4400585) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6433367) q[0];
sx q[0];
rz(-1.0717046) q[0];
sx q[0];
rz(-0.82478158) q[0];
rz(-pi) q[1];
rz(0.84757324) q[2];
sx q[2];
rz(-2.5511191) q[2];
sx q[2];
rz(-2.5440885) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.5614723) q[1];
sx q[1];
rz(-1.1271203) q[1];
sx q[1];
rz(2.2444112) q[1];
x q[2];
rz(1.2626889) q[3];
sx q[3];
rz(-2.1854221) q[3];
sx q[3];
rz(-2.7648787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.39701617) q[2];
sx q[2];
rz(-1.959356) q[2];
sx q[2];
rz(-1.5783295) q[2];
rz(-2.8619316) q[3];
sx q[3];
rz(-0.71664387) q[3];
sx q[3];
rz(-0.40951148) q[3];
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
rz(2.3438943) q[0];
sx q[0];
rz(-2.7293623) q[0];
sx q[0];
rz(-1.4233587) q[0];
rz(3.0176945) q[1];
sx q[1];
rz(-2.2975477) q[1];
sx q[1];
rz(-1.5843676) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7182962) q[0];
sx q[0];
rz(-2.3580821) q[0];
sx q[0];
rz(-1.9968724) q[0];
x q[1];
rz(0.84752797) q[2];
sx q[2];
rz(-1.28994) q[2];
sx q[2];
rz(0.47570634) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.8516006) q[1];
sx q[1];
rz(-2.8826635) q[1];
sx q[1];
rz(-1.2570791) q[1];
rz(-pi) q[2];
x q[2];
rz(0.73488124) q[3];
sx q[3];
rz(-0.85414825) q[3];
sx q[3];
rz(1.7295966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.0763756) q[2];
sx q[2];
rz(-2.4880444) q[2];
sx q[2];
rz(2.2067113) q[2];
rz(-0.30385083) q[3];
sx q[3];
rz(-0.6128208) q[3];
sx q[3];
rz(-2.2936308) q[3];
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
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76871753) q[0];
sx q[0];
rz(-2.3621552) q[0];
sx q[0];
rz(-0.26707643) q[0];
rz(-3.0067387) q[1];
sx q[1];
rz(-2.5649773) q[1];
sx q[1];
rz(0.83736247) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1601279) q[0];
sx q[0];
rz(-2.5257887) q[0];
sx q[0];
rz(-2.0439022) q[0];
x q[1];
rz(-2.1979273) q[2];
sx q[2];
rz(-0.84463464) q[2];
sx q[2];
rz(-0.86792714) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.76772753) q[1];
sx q[1];
rz(-0.30555913) q[1];
sx q[1];
rz(0.69626804) q[1];
x q[2];
rz(-2.2211406) q[3];
sx q[3];
rz(-3.0643769) q[3];
sx q[3];
rz(-1.6473351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.5460633) q[2];
sx q[2];
rz(-2.9661621) q[2];
sx q[2];
rz(2.9626633) q[2];
rz(1.1692125) q[3];
sx q[3];
rz(-1.365265) q[3];
sx q[3];
rz(0.21807142) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1201852) q[0];
sx q[0];
rz(-2.6299801) q[0];
sx q[0];
rz(-2.2688493) q[0];
rz(-1.9841638) q[1];
sx q[1];
rz(-2.2667784) q[1];
sx q[1];
rz(-3.1246368) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2568396) q[0];
sx q[0];
rz(-1.5829979) q[0];
sx q[0];
rz(0.73464616) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.59973195) q[2];
sx q[2];
rz(-1.5790734) q[2];
sx q[2];
rz(3.099583) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.3185721) q[1];
sx q[1];
rz(-1.8649532) q[1];
sx q[1];
rz(-0.78359722) q[1];
rz(-1.8618552) q[3];
sx q[3];
rz(-1.3123056) q[3];
sx q[3];
rz(0.84421009) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0852647) q[2];
sx q[2];
rz(-1.4339829) q[2];
sx q[2];
rz(-2.8223574) q[2];
rz(-2.4625835) q[3];
sx q[3];
rz(-1.7947936) q[3];
sx q[3];
rz(0.80176789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7769258) q[0];
sx q[0];
rz(-0.32061446) q[0];
sx q[0];
rz(2.4989682) q[0];
rz(-1.369426) q[1];
sx q[1];
rz(-0.6232999) q[1];
sx q[1];
rz(-2.1646037) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1851121) q[0];
sx q[0];
rz(-1.6688884) q[0];
sx q[0];
rz(3.0177659) q[0];
rz(-pi) q[1];
rz(-0.50708167) q[2];
sx q[2];
rz(-2.0424105) q[2];
sx q[2];
rz(0.44632775) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.7599267) q[1];
sx q[1];
rz(-1.5846415) q[1];
sx q[1];
rz(2.2626876) q[1];
rz(-pi) q[2];
rz(3.0877572) q[3];
sx q[3];
rz(-0.98568688) q[3];
sx q[3];
rz(-1.1424292) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.52648181) q[2];
sx q[2];
rz(-1.7340163) q[2];
sx q[2];
rz(1.2283481) q[2];
rz(-2.2743716) q[3];
sx q[3];
rz(-2.7286178) q[3];
sx q[3];
rz(0.76798463) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41384554) q[0];
sx q[0];
rz(-0.21738805) q[0];
sx q[0];
rz(-2.9687498) q[0];
rz(2.3364283) q[1];
sx q[1];
rz(-2.5925437) q[1];
sx q[1];
rz(3.0056312) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7522536) q[0];
sx q[0];
rz(-2.2704008) q[0];
sx q[0];
rz(-2.9851578) q[0];
rz(-2.3195761) q[2];
sx q[2];
rz(-1.4944585) q[2];
sx q[2];
rz(2.1427296) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.6579305) q[1];
sx q[1];
rz(-1.4560946) q[1];
sx q[1];
rz(1.5516267) q[1];
rz(2.7649859) q[3];
sx q[3];
rz(-0.35264243) q[3];
sx q[3];
rz(-0.4761179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2148296) q[2];
sx q[2];
rz(-1.5462993) q[2];
sx q[2];
rz(0.1845486) q[2];
rz(2.9522225) q[3];
sx q[3];
rz(-2.6034077) q[3];
sx q[3];
rz(-2.2034933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1592584) q[0];
sx q[0];
rz(-2.8026447) q[0];
sx q[0];
rz(-0.92639297) q[0];
rz(-3.1023846) q[1];
sx q[1];
rz(-2.4640633) q[1];
sx q[1];
rz(-2.1047986) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77039546) q[0];
sx q[0];
rz(-1.7465034) q[0];
sx q[0];
rz(-1.8699339) q[0];
x q[1];
rz(-2.762297) q[2];
sx q[2];
rz(-2.6808028) q[2];
sx q[2];
rz(-0.57578218) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.8421887) q[1];
sx q[1];
rz(-1.8129947) q[1];
sx q[1];
rz(-0.41917015) q[1];
rz(-2.751334) q[3];
sx q[3];
rz(-1.7446127) q[3];
sx q[3];
rz(2.9040608) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.2490273) q[2];
sx q[2];
rz(-1.1250863) q[2];
sx q[2];
rz(-0.79891515) q[2];
rz(-0.51698452) q[3];
sx q[3];
rz(-2.7345782) q[3];
sx q[3];
rz(-0.5504722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35140458) q[0];
sx q[0];
rz(-3.1365972) q[0];
sx q[0];
rz(-1.6222401) q[0];
rz(1.5478569) q[1];
sx q[1];
rz(-2.3212815) q[1];
sx q[1];
rz(2.4511852) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0102745) q[0];
sx q[0];
rz(-1.2206435) q[0];
sx q[0];
rz(1.8084333) q[0];
x q[1];
rz(1.4682653) q[2];
sx q[2];
rz(-0.77418113) q[2];
sx q[2];
rz(0.098086327) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.71296299) q[1];
sx q[1];
rz(-0.74364122) q[1];
sx q[1];
rz(2.3548467) q[1];
rz(-pi) q[2];
x q[2];
rz(1.810174) q[3];
sx q[3];
rz(-2.5898511) q[3];
sx q[3];
rz(0.61659471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.4385628) q[2];
sx q[2];
rz(-2.5327693) q[2];
sx q[2];
rz(-2.1717066) q[2];
rz(-0.25740933) q[3];
sx q[3];
rz(-2.9195547) q[3];
sx q[3];
rz(-1.359587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
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
rz(2.6064706) q[0];
sx q[0];
rz(-0.73000014) q[0];
sx q[0];
rz(-3.0560793) q[0];
rz(2.8657148) q[1];
sx q[1];
rz(-0.62092263) q[1];
sx q[1];
rz(-1.9524908) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6333503) q[0];
sx q[0];
rz(-1.9353284) q[0];
sx q[0];
rz(2.7056498) q[0];
rz(-pi) q[1];
rz(-3.1136572) q[2];
sx q[2];
rz(-0.31012529) q[2];
sx q[2];
rz(-1.5903697) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.738742) q[1];
sx q[1];
rz(-1.1711118) q[1];
sx q[1];
rz(0.12895361) q[1];
x q[2];
rz(-0.25190763) q[3];
sx q[3];
rz(-2.7828476) q[3];
sx q[3];
rz(0.22049604) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.4756061) q[2];
sx q[2];
rz(-2.2624272) q[2];
sx q[2];
rz(2.6574668) q[2];
rz(3.0020946) q[3];
sx q[3];
rz(-0.77553427) q[3];
sx q[3];
rz(0.16105306) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.465268) q[0];
sx q[0];
rz(-1.9857255) q[0];
sx q[0];
rz(2.3406512) q[0];
rz(-1.4576661) q[1];
sx q[1];
rz(-1.4977581) q[1];
sx q[1];
rz(-1.3118634) q[1];
rz(-1.7745849) q[2];
sx q[2];
rz(-2.2098354) q[2];
sx q[2];
rz(-0.60571203) q[2];
rz(-2.3845354) q[3];
sx q[3];
rz(-1.5499877) q[3];
sx q[3];
rz(-3.1201759) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
