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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22390511) q[0];
sx q[0];
rz(-1.5263867) q[0];
sx q[0];
rz(-3.1085148) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4307664) q[2];
sx q[2];
rz(-2.4190355) q[2];
sx q[2];
rz(0.089499105) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5482169) q[1];
sx q[1];
rz(-2.0228099) q[1];
sx q[1];
rz(1.7991245) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9494667) q[3];
sx q[3];
rz(-0.14992564) q[3];
sx q[3];
rz(0.42335864) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.5662745) q[2];
sx q[2];
rz(-0.34276572) q[2];
sx q[2];
rz(-2.9052367) q[2];
rz(2.6929839) q[3];
sx q[3];
rz(-1.6581422) q[3];
sx q[3];
rz(1.0033222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7597294) q[0];
sx q[0];
rz(-2.6549082) q[0];
sx q[0];
rz(-2.3554262) q[0];
rz(-2.3568514) q[1];
sx q[1];
rz(-1.6678383) q[1];
sx q[1];
rz(-0.10202185) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.090793153) q[0];
sx q[0];
rz(-2.5623119) q[0];
sx q[0];
rz(-2.7673278) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3574102) q[2];
sx q[2];
rz(-1.783737) q[2];
sx q[2];
rz(-0.41894693) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.35013546) q[1];
sx q[1];
rz(-1.8730436) q[1];
sx q[1];
rz(1.2863897) q[1];
x q[2];
rz(-0.24436538) q[3];
sx q[3];
rz(-2.3373342) q[3];
sx q[3];
rz(1.8387295) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.1138136) q[2];
sx q[2];
rz(-1.6543829) q[2];
sx q[2];
rz(2.8786744) q[2];
rz(-1.3024088) q[3];
sx q[3];
rz(-2.0263367) q[3];
sx q[3];
rz(0.7836248) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0021492783) q[0];
sx q[0];
rz(-0.0037010598) q[0];
sx q[0];
rz(1.333746) q[0];
rz(1.3847146) q[1];
sx q[1];
rz(-0.92888558) q[1];
sx q[1];
rz(-0.76470107) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9847577) q[0];
sx q[0];
rz(-1.3768702) q[0];
sx q[0];
rz(-0.24389275) q[0];
rz(-pi) q[1];
rz(2.7638859) q[2];
sx q[2];
rz(-1.3481082) q[2];
sx q[2];
rz(1.6645704) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6669608) q[1];
sx q[1];
rz(-1.8256011) q[1];
sx q[1];
rz(2.9632225) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9165048) q[3];
sx q[3];
rz(-0.74784333) q[3];
sx q[3];
rz(-1.9077099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.1769522) q[2];
sx q[2];
rz(-1.8884337) q[2];
sx q[2];
rz(-2.9580252) q[2];
rz(-2.99672) q[3];
sx q[3];
rz(-1.262007) q[3];
sx q[3];
rz(0.22411331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95530987) q[0];
sx q[0];
rz(-2.0109542) q[0];
sx q[0];
rz(2.8593707) q[0];
rz(1.1497633) q[1];
sx q[1];
rz(-1.7825922) q[1];
sx q[1];
rz(1.6414292) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44631413) q[0];
sx q[0];
rz(-1.9290961) q[0];
sx q[0];
rz(-1.3867239) q[0];
x q[1];
rz(1.167539) q[2];
sx q[2];
rz(-0.58524281) q[2];
sx q[2];
rz(-2.5257021) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4053661) q[1];
sx q[1];
rz(-0.42441503) q[1];
sx q[1];
rz(2.2176803) q[1];
rz(-0.28516523) q[3];
sx q[3];
rz(-2.8252606) q[3];
sx q[3];
rz(-1.118699) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.1993316) q[2];
sx q[2];
rz(-1.6959689) q[2];
sx q[2];
rz(1.65421) q[2];
rz(-2.2516294) q[3];
sx q[3];
rz(-1.3993989) q[3];
sx q[3];
rz(0.14931211) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99059659) q[0];
sx q[0];
rz(-2.6166333) q[0];
sx q[0];
rz(-1.4599482) q[0];
rz(-1.2851985) q[1];
sx q[1];
rz(-1.3672914) q[1];
sx q[1];
rz(-0.45305124) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40719068) q[0];
sx q[0];
rz(-1.2670867) q[0];
sx q[0];
rz(3.0846473) q[0];
rz(1.4905246) q[2];
sx q[2];
rz(-0.58387127) q[2];
sx q[2];
rz(0.50498) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0041607) q[1];
sx q[1];
rz(-2.7212226) q[1];
sx q[1];
rz(-2.67291) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0208218) q[3];
sx q[3];
rz(-1.9450359) q[3];
sx q[3];
rz(-2.9974724) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.65752658) q[2];
sx q[2];
rz(-2.6675197) q[2];
sx q[2];
rz(1.5566114) q[2];
rz(1.5629684) q[3];
sx q[3];
rz(-1.862674) q[3];
sx q[3];
rz(-2.1141619) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
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
rz(0.46074834) q[1];
sx q[1];
rz(-1.2604424) q[1];
sx q[1];
rz(1.4985098) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1750893) q[0];
sx q[0];
rz(-1.5700392) q[0];
sx q[0];
rz(2.1835292) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.837311) q[2];
sx q[2];
rz(-1.2488447) q[2];
sx q[2];
rz(2.6139174) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.431659) q[1];
sx q[1];
rz(-1.6282896) q[1];
sx q[1];
rz(1.8338982) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2727358) q[3];
sx q[3];
rz(-1.6770937) q[3];
sx q[3];
rz(-2.2095263) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.1341165) q[2];
sx q[2];
rz(-1.7191929) q[2];
sx q[2];
rz(-0.22809347) q[2];
rz(2.5988233) q[3];
sx q[3];
rz(-1.0017064) q[3];
sx q[3];
rz(0.91226474) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11693624) q[0];
sx q[0];
rz(-2.0967364) q[0];
sx q[0];
rz(0.60633099) q[0];
rz(1.2888651) q[1];
sx q[1];
rz(-1.6090798) q[1];
sx q[1];
rz(2.2859763) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3324749) q[0];
sx q[0];
rz(-2.1590589) q[0];
sx q[0];
rz(-1.0756798) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4855723) q[2];
sx q[2];
rz(-2.3422675) q[2];
sx q[2];
rz(-2.8003789) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.42509584) q[1];
sx q[1];
rz(-1.7762868) q[1];
sx q[1];
rz(-2.53611) q[1];
x q[2];
rz(2.8884573) q[3];
sx q[3];
rz(-2.0380424) q[3];
sx q[3];
rz(0.96880355) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7938101) q[0];
sx q[0];
rz(-0.3599444) q[0];
sx q[0];
rz(2.1600294) q[0];
rz(1.7836102) q[1];
sx q[1];
rz(-2.6505018) q[1];
sx q[1];
rz(2.8474999) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9413853) q[0];
sx q[0];
rz(-1.6305171) q[0];
sx q[0];
rz(-2.8325547) q[0];
x q[1];
rz(-0.37224877) q[2];
sx q[2];
rz(-1.8342092) q[2];
sx q[2];
rz(-1.5836704) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.8311685) q[1];
sx q[1];
rz(-0.71485315) q[1];
sx q[1];
rz(2.0063041) q[1];
x q[2];
rz(-0.050042583) q[3];
sx q[3];
rz(-1.4980157) q[3];
sx q[3];
rz(-2.4693054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0391417) q[2];
sx q[2];
rz(-2.4852018) q[2];
sx q[2];
rz(-1.0605109) q[2];
rz(2.5833526) q[3];
sx q[3];
rz(-0.57359901) q[3];
sx q[3];
rz(3.1225045) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3625665) q[0];
sx q[0];
rz(-2.9988852) q[0];
sx q[0];
rz(0.78417626) q[0];
rz(-1.342429) q[1];
sx q[1];
rz(-1.289184) q[1];
sx q[1];
rz(0.69650355) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7292228) q[0];
sx q[0];
rz(-1.3431088) q[0];
sx q[0];
rz(0.21842893) q[0];
rz(-2.0663459) q[2];
sx q[2];
rz(-2.181567) q[2];
sx q[2];
rz(1.3411759) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1494157) q[1];
sx q[1];
rz(-1.0806298) q[1];
sx q[1];
rz(-1.0926682) q[1];
rz(-pi) q[2];
x q[2];
rz(2.302565) q[3];
sx q[3];
rz(-0.19908842) q[3];
sx q[3];
rz(-1.3092878) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.8037618) q[2];
sx q[2];
rz(-0.75118128) q[2];
sx q[2];
rz(1.8776228) q[2];
rz(1.7654644) q[3];
sx q[3];
rz(-0.91457808) q[3];
sx q[3];
rz(-1.5862563) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2176168) q[0];
sx q[0];
rz(-0.042984977) q[0];
sx q[0];
rz(-1.4772344) q[0];
rz(2.2337275) q[1];
sx q[1];
rz(-1.5217179) q[1];
sx q[1];
rz(-0.97000617) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.086921) q[0];
sx q[0];
rz(-0.63935125) q[0];
sx q[0];
rz(-1.1631835) q[0];
rz(-2.1625948) q[2];
sx q[2];
rz(-1.819397) q[2];
sx q[2];
rz(-0.73038855) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.86783389) q[1];
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
rz(-pi/2) q[1];
sx q[1];
rz(2.9475391) q[2];
sx q[2];
rz(-2.2668362) q[2];
sx q[2];
rz(-0.27210316) q[2];
rz(2.2943606) q[3];
sx q[3];
rz(-2.6341485) q[3];
sx q[3];
rz(-1.0890755) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67615164) q[0];
sx q[0];
rz(-1.1347102) q[0];
sx q[0];
rz(1.4665428) q[0];
rz(-2.8027986) q[1];
sx q[1];
rz(-1.4452965) q[1];
sx q[1];
rz(-1.8903587) q[1];
rz(-2.7441944) q[2];
sx q[2];
rz(-2.3811679) q[2];
sx q[2];
rz(0.74447348) q[2];
rz(-1.7539514) q[3];
sx q[3];
rz(-0.42944254) q[3];
sx q[3];
rz(-2.4933168) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
