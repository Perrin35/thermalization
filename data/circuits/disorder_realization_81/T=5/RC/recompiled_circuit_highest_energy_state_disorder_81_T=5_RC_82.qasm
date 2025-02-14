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
rz(4.0153761) q[0];
sx q[0];
rz(10.562513) q[0];
rz(-3.0955834) q[1];
sx q[1];
rz(0.53541056) q[1];
sx q[1];
rz(14.169197) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22390511) q[0];
sx q[0];
rz(-1.5263867) q[0];
sx q[0];
rz(3.1085148) q[0];
rz(-pi) q[1];
rz(0.85311546) q[2];
sx q[2];
rz(-1.6632281) q[2];
sx q[2];
rz(1.554956) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0178609) q[1];
sx q[1];
rz(-1.7758472) q[1];
sx q[1];
rz(-0.46242949) q[1];
x q[2];
rz(-3.0858062) q[3];
sx q[3];
rz(-1.4315637) q[3];
sx q[3];
rz(3.1007953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.5753182) q[2];
sx q[2];
rz(-2.7988269) q[2];
sx q[2];
rz(0.23635593) q[2];
rz(0.44860873) q[3];
sx q[3];
rz(-1.6581422) q[3];
sx q[3];
rz(2.1382704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7597294) q[0];
sx q[0];
rz(-2.6549082) q[0];
sx q[0];
rz(2.3554262) q[0];
rz(0.78474125) q[1];
sx q[1];
rz(-1.6678383) q[1];
sx q[1];
rz(3.0395708) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.611972) q[0];
sx q[0];
rz(-2.1054322) q[0];
sx q[0];
rz(1.8055339) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8445016) q[2];
sx q[2];
rz(-0.80657437) q[2];
sx q[2];
rz(-1.3606124) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.1338623) q[1];
sx q[1];
rz(-1.8419767) q[1];
sx q[1];
rz(-0.31409632) q[1];
rz(-2.8972273) q[3];
sx q[3];
rz(-0.8042585) q[3];
sx q[3];
rz(-1.3028631) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.027779) q[2];
sx q[2];
rz(-1.4872097) q[2];
sx q[2];
rz(0.2629183) q[2];
rz(1.8391838) q[3];
sx q[3];
rz(-1.115256) q[3];
sx q[3];
rz(2.3579679) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(-0.0021492783) q[0];
sx q[0];
rz(-3.1378916) q[0];
sx q[0];
rz(-1.8078467) q[0];
rz(-1.7568781) q[1];
sx q[1];
rz(-2.2127071) q[1];
sx q[1];
rz(0.76470107) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9847577) q[0];
sx q[0];
rz(-1.7647224) q[0];
sx q[0];
rz(-0.24389275) q[0];
rz(-pi) q[1];
rz(-2.5909293) q[2];
sx q[2];
rz(-2.705859) q[2];
sx q[2];
rz(2.7274362) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4746318) q[1];
sx q[1];
rz(-1.8256011) q[1];
sx q[1];
rz(-2.9632225) q[1];
x q[2];
rz(-0.85327236) q[3];
sx q[3];
rz(-1.8033334) q[3];
sx q[3];
rz(0.59508376) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.1769522) q[2];
sx q[2];
rz(-1.8884337) q[2];
sx q[2];
rz(-0.18356744) q[2];
rz(-0.14487264) q[3];
sx q[3];
rz(-1.262007) q[3];
sx q[3];
rz(-0.22411331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95530987) q[0];
sx q[0];
rz(-1.1306385) q[0];
sx q[0];
rz(2.8593707) q[0];
rz(-1.1497633) q[1];
sx q[1];
rz(-1.3590004) q[1];
sx q[1];
rz(-1.5001635) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.041752664) q[0];
sx q[0];
rz(-2.7405993) q[0];
sx q[0];
rz(2.6869511) q[0];
rz(-pi) q[1];
rz(1.9740536) q[2];
sx q[2];
rz(-2.5563498) q[2];
sx q[2];
rz(0.61589059) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.4373926) q[1];
sx q[1];
rz(-1.3199894) q[1];
sx q[1];
rz(1.2247242) q[1];
x q[2];
rz(2.8372466) q[3];
sx q[3];
rz(-1.658421) q[3];
sx q[3];
rz(-2.9612142) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.942261) q[2];
sx q[2];
rz(-1.6959689) q[2];
sx q[2];
rz(-1.65421) q[2];
rz(-2.2516294) q[3];
sx q[3];
rz(-1.3993989) q[3];
sx q[3];
rz(-2.9922805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99059659) q[0];
sx q[0];
rz(-2.6166333) q[0];
sx q[0];
rz(1.4599482) q[0];
rz(1.8563942) q[1];
sx q[1];
rz(-1.7743013) q[1];
sx q[1];
rz(0.45305124) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.734402) q[0];
sx q[0];
rz(-1.2670867) q[0];
sx q[0];
rz(3.0846473) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6510681) q[2];
sx q[2];
rz(-0.58387127) q[2];
sx q[2];
rz(-0.50498) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.275251) q[1];
sx q[1];
rz(-1.3853933) q[1];
sx q[1];
rz(-2.7621108) q[1];
rz(-pi) q[2];
rz(-0.41129859) q[3];
sx q[3];
rz(-1.1539475) q[3];
sx q[3];
rz(1.251879) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.65752658) q[2];
sx q[2];
rz(-2.6675197) q[2];
sx q[2];
rz(-1.5566114) q[2];
rz(-1.5629684) q[3];
sx q[3];
rz(-1.862674) q[3];
sx q[3];
rz(-1.0274308) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-2.9180561) q[0];
sx q[0];
rz(-2.1503088) q[0];
sx q[0];
rz(-1.1728485) q[0];
rz(2.6808443) q[1];
sx q[1];
rz(-1.2604424) q[1];
sx q[1];
rz(1.6430829) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1750893) q[0];
sx q[0];
rz(-1.5700392) q[0];
sx q[0];
rz(0.9580635) q[0];
rz(0.33289624) q[2];
sx q[2];
rz(-1.82331) q[2];
sx q[2];
rz(2.0123002) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.9869796) q[1];
sx q[1];
rz(-1.8334532) q[1];
sx q[1];
rz(0.059537454) q[1];
rz(-pi) q[2];
rz(-1.7345627) q[3];
sx q[3];
rz(-2.4330045) q[3];
sx q[3];
rz(2.6276789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.0074761) q[2];
sx q[2];
rz(-1.4223998) q[2];
sx q[2];
rz(0.22809347) q[2];
rz(-0.54276931) q[3];
sx q[3];
rz(-1.0017064) q[3];
sx q[3];
rz(-2.2293279) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11693624) q[0];
sx q[0];
rz(-2.0967364) q[0];
sx q[0];
rz(-2.5352617) q[0];
rz(-1.2888651) q[1];
sx q[1];
rz(-1.5325129) q[1];
sx q[1];
rz(-0.85561633) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3324749) q[0];
sx q[0];
rz(-0.98253378) q[0];
sx q[0];
rz(1.0756798) q[0];
x q[1];
rz(-2.3683041) q[2];
sx q[2];
rz(-1.6318562) q[2];
sx q[2];
rz(-1.9715015) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.2825477) q[1];
sx q[1];
rz(-0.63524109) q[1];
sx q[1];
rz(-2.7905491) q[1];
x q[2];
rz(2.0511914) q[3];
sx q[3];
rz(-1.3452936) q[3];
sx q[3];
rz(0.71798872) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.3768846) q[2];
sx q[2];
rz(-0.43262425) q[2];
sx q[2];
rz(-3.063859) q[2];
rz(-3.0149095) q[3];
sx q[3];
rz(-2.099791) q[3];
sx q[3];
rz(0.83109394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34778255) q[0];
sx q[0];
rz(-0.3599444) q[0];
sx q[0];
rz(-2.1600294) q[0];
rz(-1.7836102) q[1];
sx q[1];
rz(-2.6505018) q[1];
sx q[1];
rz(0.29409274) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9413853) q[0];
sx q[0];
rz(-1.5110755) q[0];
sx q[0];
rz(-0.30903791) q[0];
rz(-0.37224877) q[2];
sx q[2];
rz(-1.3073834) q[2];
sx q[2];
rz(-1.5579223) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.5433054) q[1];
sx q[1];
rz(-1.2906055) q[1];
sx q[1];
rz(-2.2375475) q[1];
x q[2];
rz(-3.0915501) q[3];
sx q[3];
rz(-1.4980157) q[3];
sx q[3];
rz(2.4693054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.10245094) q[2];
sx q[2];
rz(-0.65639085) q[2];
sx q[2];
rz(-2.0810818) q[2];
rz(-2.5833526) q[3];
sx q[3];
rz(-0.57359901) q[3];
sx q[3];
rz(0.019088117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3625665) q[0];
sx q[0];
rz(-0.14270742) q[0];
sx q[0];
rz(-0.78417626) q[0];
rz(1.7991637) q[1];
sx q[1];
rz(-1.289184) q[1];
sx q[1];
rz(-2.4450891) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7769515) q[0];
sx q[0];
rz(-0.31420194) q[0];
sx q[0];
rz(2.3227341) q[0];
x q[1];
rz(-1.0752468) q[2];
sx q[2];
rz(-0.96002561) q[2];
sx q[2];
rz(1.3411759) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.99217691) q[1];
sx q[1];
rz(-2.0609629) q[1];
sx q[1];
rz(2.0489245) q[1];
x q[2];
rz(0.83902766) q[3];
sx q[3];
rz(-0.19908842) q[3];
sx q[3];
rz(-1.8323048) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.33783087) q[2];
sx q[2];
rz(-2.3904114) q[2];
sx q[2];
rz(1.8776228) q[2];
rz(1.3761282) q[3];
sx q[3];
rz(-2.2270146) q[3];
sx q[3];
rz(1.5553364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92397583) q[0];
sx q[0];
rz(-3.0986077) q[0];
sx q[0];
rz(1.4772344) q[0];
rz(0.90786511) q[1];
sx q[1];
rz(-1.6198747) q[1];
sx q[1];
rz(-0.97000617) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0546717) q[0];
sx q[0];
rz(-0.63935125) q[0];
sx q[0];
rz(1.1631835) q[0];
rz(-2.8447609) q[2];
sx q[2];
rz(-2.1420711) q[2];
sx q[2];
rz(-2.4650857) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.2737588) q[1];
sx q[1];
rz(-0.9504488) q[1];
sx q[1];
rz(1.378744) q[1];
rz(2.41955) q[3];
sx q[3];
rz(-1.0580499) q[3];
sx q[3];
rz(-1.7554612) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.19405356) q[2];
sx q[2];
rz(-2.2668362) q[2];
sx q[2];
rz(-0.27210316) q[2];
rz(2.2943606) q[3];
sx q[3];
rz(-0.50744414) q[3];
sx q[3];
rz(1.0890755) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
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
rz(0.67615164) q[0];
sx q[0];
rz(-2.0068824) q[0];
sx q[0];
rz(-1.6750499) q[0];
rz(0.33879406) q[1];
sx q[1];
rz(-1.4452965) q[1];
sx q[1];
rz(-1.8903587) q[1];
rz(-1.2180381) q[2];
sx q[2];
rz(-0.88211664) q[2];
sx q[2];
rz(1.2695352) q[2];
rz(-1.9938888) q[3];
sx q[3];
rz(-1.4948899) q[3];
sx q[3];
rz(-0.75567452) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
