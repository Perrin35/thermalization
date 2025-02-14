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
rz(-2.6851299) q[0];
sx q[0];
rz(-0.87378341) q[0];
sx q[0];
rz(2.003858) q[0];
rz(-3.0955834) q[1];
sx q[1];
rz(-2.6061821) q[1];
sx q[1];
rz(1.5387662) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3483602) q[0];
sx q[0];
rz(-1.5377511) q[0];
sx q[0];
rz(-1.5263625) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0191604) q[2];
sx q[2];
rz(-0.85683595) q[2];
sx q[2];
rz(3.0453504) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.5482169) q[1];
sx q[1];
rz(-1.1187828) q[1];
sx q[1];
rz(1.3424681) q[1];
rz(-3.0858062) q[3];
sx q[3];
rz(-1.4315637) q[3];
sx q[3];
rz(-0.040797357) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.5753182) q[2];
sx q[2];
rz(-0.34276572) q[2];
sx q[2];
rz(0.23635593) q[2];
rz(-2.6929839) q[3];
sx q[3];
rz(-1.6581422) q[3];
sx q[3];
rz(2.1382704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7597294) q[0];
sx q[0];
rz(-0.48668447) q[0];
sx q[0];
rz(2.3554262) q[0];
rz(2.3568514) q[1];
sx q[1];
rz(-1.6678383) q[1];
sx q[1];
rz(0.10202185) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9791631) q[0];
sx q[0];
rz(-1.7722881) q[0];
sx q[0];
rz(-2.5946846) q[0];
rz(-pi) q[1];
rz(-2.3574102) q[2];
sx q[2];
rz(-1.3578556) q[2];
sx q[2];
rz(2.7226457) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.7914572) q[1];
sx q[1];
rz(-1.8730436) q[1];
sx q[1];
rz(1.2863897) q[1];
rz(-1.8169448) q[3];
sx q[3];
rz(-0.79691468) q[3];
sx q[3];
rz(-1.6479542) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.1138136) q[2];
sx q[2];
rz(-1.6543829) q[2];
sx q[2];
rz(-2.8786744) q[2];
rz(1.8391838) q[3];
sx q[3];
rz(-1.115256) q[3];
sx q[3];
rz(2.3579679) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0021492783) q[0];
sx q[0];
rz(-3.1378916) q[0];
sx q[0];
rz(-1.333746) q[0];
rz(-1.3847146) q[1];
sx q[1];
rz(-0.92888558) q[1];
sx q[1];
rz(0.76470107) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36604213) q[0];
sx q[0];
rz(-1.3315655) q[0];
sx q[0];
rz(-1.3711098) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.5506634) q[2];
sx q[2];
rz(-2.705859) q[2];
sx q[2];
rz(0.41415641) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(3.0908392) q[1];
sx q[1];
rz(-1.3982441) q[1];
sx q[1];
rz(-1.3120632) q[1];
rz(-pi) q[2];
rz(-2.2883203) q[3];
sx q[3];
rz(-1.8033334) q[3];
sx q[3];
rz(2.5465089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.9646405) q[2];
sx q[2];
rz(-1.8884337) q[2];
sx q[2];
rz(0.18356744) q[2];
rz(0.14487264) q[3];
sx q[3];
rz(-1.8795857) q[3];
sx q[3];
rz(-0.22411331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95530987) q[0];
sx q[0];
rz(-2.0109542) q[0];
sx q[0];
rz(-2.8593707) q[0];
rz(1.1497633) q[1];
sx q[1];
rz(-1.3590004) q[1];
sx q[1];
rz(-1.6414292) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.951913) q[0];
sx q[0];
rz(-1.743058) q[0];
sx q[0];
rz(-0.36393117) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1181951) q[2];
sx q[2];
rz(-1.3522902) q[2];
sx q[2];
rz(-0.61321875) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7042) q[1];
sx q[1];
rz(-1.8216032) q[1];
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
rz(-pi) q[1];
x q[1];
rz(1.942261) q[2];
sx q[2];
rz(-1.4456238) q[2];
sx q[2];
rz(-1.4873827) q[2];
rz(2.2516294) q[3];
sx q[3];
rz(-1.3993989) q[3];
sx q[3];
rz(2.9922805) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1509961) q[0];
sx q[0];
rz(-2.6166333) q[0];
sx q[0];
rz(-1.6816444) q[0];
rz(1.2851985) q[1];
sx q[1];
rz(-1.7743013) q[1];
sx q[1];
rz(-0.45305124) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40719068) q[0];
sx q[0];
rz(-1.8745059) q[0];
sx q[0];
rz(3.0846473) q[0];
rz(-pi) q[1];
rz(-3.0886623) q[2];
sx q[2];
rz(-0.9890511) q[2];
sx q[2];
rz(-2.5404921) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.0041607) q[1];
sx q[1];
rz(-2.7212226) q[1];
sx q[1];
rz(0.46868268) q[1];
rz(-pi) q[2];
rz(1.1207709) q[3];
sx q[3];
rz(-1.1965568) q[3];
sx q[3];
rz(-2.9974724) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.4840661) q[2];
sx q[2];
rz(-2.6675197) q[2];
sx q[2];
rz(1.5849812) q[2];
rz(1.5629684) q[3];
sx q[3];
rz(-1.2789187) q[3];
sx q[3];
rz(2.1141619) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9180561) q[0];
sx q[0];
rz(-2.1503088) q[0];
sx q[0];
rz(1.9687442) q[0];
rz(2.6808443) q[1];
sx q[1];
rz(-1.2604424) q[1];
sx q[1];
rz(-1.4985098) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7448085) q[0];
sx q[0];
rz(-0.61273324) q[0];
sx q[0];
rz(1.5694797) q[0];
rz(-2.8086964) q[2];
sx q[2];
rz(-1.82331) q[2];
sx q[2];
rz(-1.1292924) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.431659) q[1];
sx q[1];
rz(-1.5133031) q[1];
sx q[1];
rz(1.8338982) q[1];
rz(-pi) q[2];
x q[2];
rz(1.40703) q[3];
sx q[3];
rz(-2.4330045) q[3];
sx q[3];
rz(-0.51391376) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.1341165) q[2];
sx q[2];
rz(-1.7191929) q[2];
sx q[2];
rz(0.22809347) q[2];
rz(2.5988233) q[3];
sx q[3];
rz(-2.1398862) q[3];
sx q[3];
rz(-0.91226474) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11693624) q[0];
sx q[0];
rz(-1.0448562) q[0];
sx q[0];
rz(0.60633099) q[0];
rz(-1.8527276) q[1];
sx q[1];
rz(-1.6090798) q[1];
sx q[1];
rz(2.2859763) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3324749) q[0];
sx q[0];
rz(-0.98253378) q[0];
sx q[0];
rz(1.0756798) q[0];
rz(-0.77328859) q[2];
sx q[2];
rz(-1.5097364) q[2];
sx q[2];
rz(-1.9715015) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.8555755) q[1];
sx q[1];
rz(-0.97980503) q[1];
sx q[1];
rz(-1.3225287) q[1];
rz(-pi) q[2];
rz(-1.110027) q[3];
sx q[3];
rz(-2.6146725) q[3];
sx q[3];
rz(-2.6941018) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.764708) q[2];
sx q[2];
rz(-2.7089684) q[2];
sx q[2];
rz(3.063859) q[2];
rz(-0.12668315) q[3];
sx q[3];
rz(-2.099791) q[3];
sx q[3];
rz(-0.83109394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34778255) q[0];
sx q[0];
rz(-2.7816483) q[0];
sx q[0];
rz(-0.98156324) q[0];
rz(1.7836102) q[1];
sx q[1];
rz(-0.49109083) q[1];
sx q[1];
rz(0.29409274) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9413853) q[0];
sx q[0];
rz(-1.5110755) q[0];
sx q[0];
rz(-2.8325547) q[0];
rz(2.5035759) q[2];
sx q[2];
rz(-2.6891481) q[2];
sx q[2];
rz(-2.5663301) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.24178019) q[1];
sx q[1];
rz(-2.2072148) q[1];
sx q[1];
rz(0.35102014) q[1];
x q[2];
rz(-2.1721187) q[3];
sx q[3];
rz(-3.0532928) q[3];
sx q[3];
rz(1.2754319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.10245094) q[2];
sx q[2];
rz(-0.65639085) q[2];
sx q[2];
rz(1.0605109) q[2];
rz(-0.55824009) q[3];
sx q[3];
rz(-0.57359901) q[3];
sx q[3];
rz(3.1225045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-1.3625665) q[0];
sx q[0];
rz(-0.14270742) q[0];
sx q[0];
rz(-2.3574164) q[0];
rz(-1.7991637) q[1];
sx q[1];
rz(-1.8524086) q[1];
sx q[1];
rz(0.69650355) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2084889) q[0];
sx q[0];
rz(-1.7835) q[0];
sx q[0];
rz(1.8038294) q[0];
rz(-pi) q[1];
rz(2.5449341) q[2];
sx q[2];
rz(-0.76603973) q[2];
sx q[2];
rz(-2.5562364) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.9827658) q[1];
sx q[1];
rz(-2.4707795) q[1];
sx q[1];
rz(0.71158408) q[1];
rz(-0.83902766) q[3];
sx q[3];
rz(-0.19908842) q[3];
sx q[3];
rz(-1.3092878) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.33783087) q[2];
sx q[2];
rz(-0.75118128) q[2];
sx q[2];
rz(-1.8776228) q[2];
rz(-1.3761282) q[3];
sx q[3];
rz(-2.2270146) q[3];
sx q[3];
rz(-1.5553364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92397583) q[0];
sx q[0];
rz(-3.0986077) q[0];
sx q[0];
rz(1.4772344) q[0];
rz(-2.2337275) q[1];
sx q[1];
rz(-1.5217179) q[1];
sx q[1];
rz(0.97000617) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9912787) q[0];
sx q[0];
rz(-1.3320001) q[0];
sx q[0];
rz(2.1697609) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9978421) q[2];
sx q[2];
rz(-2.5055024) q[2];
sx q[2];
rz(1.191312) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.86783389) q[1];
sx q[1];
rz(-0.9504488) q[1];
sx q[1];
rz(-1.7628487) q[1];
x q[2];
rz(-0.92718683) q[3];
sx q[3];
rz(-2.1845316) q[3];
sx q[3];
rz(-0.2231547) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.19405356) q[2];
sx q[2];
rz(-0.87475646) q[2];
sx q[2];
rz(-2.8694895) q[2];
rz(-2.2943606) q[3];
sx q[3];
rz(-0.50744414) q[3];
sx q[3];
rz(2.0525172) q[3];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.465441) q[0];
sx q[0];
rz(-1.1347102) q[0];
sx q[0];
rz(1.4665428) q[0];
rz(2.8027986) q[1];
sx q[1];
rz(-1.6962961) q[1];
sx q[1];
rz(1.2512339) q[1];
rz(-2.4215578) q[2];
sx q[2];
rz(-1.8408114) q[2];
sx q[2];
rz(2.6104952) q[2];
rz(1.1477039) q[3];
sx q[3];
rz(-1.4948899) q[3];
sx q[3];
rz(-0.75567452) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
