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
rz(0.45646271) q[0];
sx q[0];
rz(-2.2678092) q[0];
sx q[0];
rz(1.1377347) q[0];
rz(-3.0955834) q[1];
sx q[1];
rz(-2.6061821) q[1];
sx q[1];
rz(1.5387662) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3483602) q[0];
sx q[0];
rz(-1.5377511) q[0];
sx q[0];
rz(-1.6152302) q[0];
x q[1];
rz(-1.4307664) q[2];
sx q[2];
rz(-0.7225572) q[2];
sx q[2];
rz(0.089499105) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.0592717) q[1];
sx q[1];
rz(-2.6387625) q[1];
sx q[1];
rz(2.7053653) q[1];
x q[2];
rz(-1.7102431) q[3];
sx q[3];
rz(-1.5155503) q[3];
sx q[3];
rz(-1.6193438) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.5753182) q[2];
sx q[2];
rz(-2.7988269) q[2];
sx q[2];
rz(0.23635593) q[2];
rz(-2.6929839) q[3];
sx q[3];
rz(-1.6581422) q[3];
sx q[3];
rz(-1.0033222) q[3];
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
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3818632) q[0];
sx q[0];
rz(-0.48668447) q[0];
sx q[0];
rz(-0.78616649) q[0];
rz(-2.3568514) q[1];
sx q[1];
rz(-1.4737543) q[1];
sx q[1];
rz(-3.0395708) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52962063) q[0];
sx q[0];
rz(-2.1054322) q[0];
sx q[0];
rz(-1.8055339) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8672075) q[2];
sx q[2];
rz(-0.80889946) q[2];
sx q[2];
rz(-2.1975225) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0151012) q[1];
sx q[1];
rz(-0.41200629) q[1];
sx q[1];
rz(0.73275779) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.24436538) q[3];
sx q[3];
rz(-2.3373342) q[3];
sx q[3];
rz(-1.3028631) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1138136) q[2];
sx q[2];
rz(-1.6543829) q[2];
sx q[2];
rz(-0.2629183) q[2];
rz(-1.3024088) q[3];
sx q[3];
rz(-1.115256) q[3];
sx q[3];
rz(-0.7836248) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1394434) q[0];
sx q[0];
rz(-0.0037010598) q[0];
sx q[0];
rz(-1.8078467) q[0];
rz(-1.3847146) q[1];
sx q[1];
rz(-0.92888558) q[1];
sx q[1];
rz(0.76470107) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.156835) q[0];
sx q[0];
rz(-1.7647224) q[0];
sx q[0];
rz(2.8976999) q[0];
x q[1];
rz(2.5909293) q[2];
sx q[2];
rz(-0.43573365) q[2];
sx q[2];
rz(-0.41415641) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.050753442) q[1];
sx q[1];
rz(-1.3982441) q[1];
sx q[1];
rz(-1.8295294) q[1];
rz(-2.8370492) q[3];
sx q[3];
rz(-2.2650654) q[3];
sx q[3];
rz(-2.3643431) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.1769522) q[2];
sx q[2];
rz(-1.2531589) q[2];
sx q[2];
rz(0.18356744) q[2];
rz(0.14487264) q[3];
sx q[3];
rz(-1.8795857) q[3];
sx q[3];
rz(2.9174793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95530987) q[0];
sx q[0];
rz(-1.1306385) q[0];
sx q[0];
rz(-0.282222) q[0];
rz(1.1497633) q[1];
sx q[1];
rz(-1.3590004) q[1];
sx q[1];
rz(-1.6414292) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.951913) q[0];
sx q[0];
rz(-1.743058) q[0];
sx q[0];
rz(-2.7776615) q[0];
rz(-pi) q[1];
rz(1.9740536) q[2];
sx q[2];
rz(-2.5563498) q[2];
sx q[2];
rz(-2.5257021) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.4373926) q[1];
sx q[1];
rz(-1.3199894) q[1];
sx q[1];
rz(1.9168684) q[1];
rz(-pi) q[2];
rz(2.8372466) q[3];
sx q[3];
rz(-1.4831717) q[3];
sx q[3];
rz(2.9612142) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.1993316) q[2];
sx q[2];
rz(-1.6959689) q[2];
sx q[2];
rz(1.4873827) q[2];
rz(-0.88996327) q[3];
sx q[3];
rz(-1.7421937) q[3];
sx q[3];
rz(0.14931211) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1509961) q[0];
sx q[0];
rz(-0.52495933) q[0];
sx q[0];
rz(-1.4599482) q[0];
rz(1.8563942) q[1];
sx q[1];
rz(-1.3672914) q[1];
sx q[1];
rz(2.6885414) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59555028) q[0];
sx q[0];
rz(-0.30883967) q[0];
sx q[0];
rz(1.7504391) q[0];
rz(-pi) q[1];
rz(0.052930367) q[2];
sx q[2];
rz(-0.9890511) q[2];
sx q[2];
rz(0.60110053) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.275251) q[1];
sx q[1];
rz(-1.3853933) q[1];
sx q[1];
rz(-0.37948186) q[1];
x q[2];
rz(-0.83638453) q[3];
sx q[3];
rz(-2.5646113) q[3];
sx q[3];
rz(-2.0744689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.65752658) q[2];
sx q[2];
rz(-0.47407293) q[2];
sx q[2];
rz(-1.5566114) q[2];
rz(-1.5786242) q[3];
sx q[3];
rz(-1.2789187) q[3];
sx q[3];
rz(2.1141619) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9180561) q[0];
sx q[0];
rz(-2.1503088) q[0];
sx q[0];
rz(-1.9687442) q[0];
rz(-0.46074834) q[1];
sx q[1];
rz(-1.2604424) q[1];
sx q[1];
rz(1.6430829) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39517477) q[0];
sx q[0];
rz(-0.9580637) q[0];
sx q[0];
rz(-0.00092555028) q[0];
rz(-pi) q[1];
x q[1];
rz(0.33289624) q[2];
sx q[2];
rz(-1.3182827) q[2];
sx q[2];
rz(1.1292924) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-3.0705432) q[1];
sx q[1];
rz(-2.8724246) q[1];
sx q[1];
rz(1.3530003) q[1];
rz(0.86885683) q[3];
sx q[3];
rz(-1.6770937) q[3];
sx q[3];
rz(-2.2095263) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.0074761) q[2];
sx q[2];
rz(-1.4223998) q[2];
sx q[2];
rz(2.9134992) q[2];
rz(2.5988233) q[3];
sx q[3];
rz(-2.1398862) q[3];
sx q[3];
rz(-0.91226474) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0246564) q[0];
sx q[0];
rz(-2.0967364) q[0];
sx q[0];
rz(-2.5352617) q[0];
rz(1.8527276) q[1];
sx q[1];
rz(-1.5325129) q[1];
sx q[1];
rz(2.2859763) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0887849) q[0];
sx q[0];
rz(-1.1644378) q[0];
sx q[0];
rz(-2.4929422) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3683041) q[2];
sx q[2];
rz(-1.6318562) q[2];
sx q[2];
rz(1.9715015) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2825477) q[1];
sx q[1];
rz(-2.5063516) q[1];
sx q[1];
rz(-2.7905491) q[1];
rz(1.110027) q[3];
sx q[3];
rz(-0.52692014) q[3];
sx q[3];
rz(0.44749081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.3768846) q[2];
sx q[2];
rz(-0.43262425) q[2];
sx q[2];
rz(3.063859) q[2];
rz(3.0149095) q[3];
sx q[3];
rz(-2.099791) q[3];
sx q[3];
rz(-0.83109394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7938101) q[0];
sx q[0];
rz(-0.3599444) q[0];
sx q[0];
rz(0.98156324) q[0];
rz(1.3579824) q[1];
sx q[1];
rz(-0.49109083) q[1];
sx q[1];
rz(-0.29409274) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2002073) q[0];
sx q[0];
rz(-1.6305171) q[0];
sx q[0];
rz(-0.30903791) q[0];
rz(-pi) q[1];
rz(0.63801672) q[2];
sx q[2];
rz(-0.45244452) q[2];
sx q[2];
rz(-2.5663301) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.8311685) q[1];
sx q[1];
rz(-0.71485315) q[1];
sx q[1];
rz(1.1352886) q[1];
rz(-2.1721187) q[3];
sx q[3];
rz(-0.088299835) q[3];
sx q[3];
rz(-1.2754319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.10245094) q[2];
sx q[2];
rz(-2.4852018) q[2];
sx q[2];
rz(-1.0605109) q[2];
rz(-2.5833526) q[3];
sx q[3];
rz(-2.5679936) q[3];
sx q[3];
rz(3.1225045) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7790262) q[0];
sx q[0];
rz(-2.9988852) q[0];
sx q[0];
rz(-0.78417626) q[0];
rz(-1.342429) q[1];
sx q[1];
rz(-1.289184) q[1];
sx q[1];
rz(-2.4450891) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7292228) q[0];
sx q[0];
rz(-1.3431088) q[0];
sx q[0];
rz(2.9231637) q[0];
rz(0.67217153) q[2];
sx q[2];
rz(-1.9709338) q[2];
sx q[2];
rz(-0.53022417) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.81791524) q[1];
sx q[1];
rz(-1.9888048) q[1];
sx q[1];
rz(-2.6004385) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0075865) q[3];
sx q[3];
rz(-1.7184765) q[3];
sx q[3];
rz(2.0509843) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.8037618) q[2];
sx q[2];
rz(-2.3904114) q[2];
sx q[2];
rz(1.8776228) q[2];
rz(-1.7654644) q[3];
sx q[3];
rz(-0.91457808) q[3];
sx q[3];
rz(-1.5553364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92397583) q[0];
sx q[0];
rz(-0.042984977) q[0];
sx q[0];
rz(-1.6643583) q[0];
rz(-0.90786511) q[1];
sx q[1];
rz(-1.6198747) q[1];
sx q[1];
rz(-2.1715865) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9912787) q[0];
sx q[0];
rz(-1.3320001) q[0];
sx q[0];
rz(-0.97183174) q[0];
rz(1.9978421) q[2];
sx q[2];
rz(-2.5055024) q[2];
sx q[2];
rz(1.9502806) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.86783389) q[1];
sx q[1];
rz(-2.1911439) q[1];
sx q[1];
rz(-1.378744) q[1];
rz(0.72204263) q[3];
sx q[3];
rz(-2.0835428) q[3];
sx q[3];
rz(1.3861314) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.9475391) q[2];
sx q[2];
rz(-2.2668362) q[2];
sx q[2];
rz(-2.8694895) q[2];
rz(2.2943606) q[3];
sx q[3];
rz(-0.50744414) q[3];
sx q[3];
rz(1.0890755) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67615164) q[0];
sx q[0];
rz(-1.1347102) q[0];
sx q[0];
rz(1.4665428) q[0];
rz(-0.33879406) q[1];
sx q[1];
rz(-1.6962961) q[1];
sx q[1];
rz(1.2512339) q[1];
rz(-1.9235545) q[2];
sx q[2];
rz(-2.259476) q[2];
sx q[2];
rz(-1.8720575) q[2];
rz(-1.1477039) q[3];
sx q[3];
rz(-1.6467027) q[3];
sx q[3];
rz(2.3859181) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
