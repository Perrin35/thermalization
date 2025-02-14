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
rz(-2.7954571) q[0];
sx q[0];
rz(-1.6383642) q[0];
sx q[0];
rz(2.4099953) q[0];
rz(-2.6919964) q[1];
sx q[1];
rz(-0.1875339) q[1];
sx q[1];
rz(-0.43391689) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2303378) q[0];
sx q[0];
rz(-0.83090913) q[0];
sx q[0];
rz(3.0940542) q[0];
x q[1];
rz(0.5753168) q[2];
sx q[2];
rz(-2.0849209) q[2];
sx q[2];
rz(2.7935087) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.190769) q[1];
sx q[1];
rz(-0.94591037) q[1];
sx q[1];
rz(-0.48887078) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4106969) q[3];
sx q[3];
rz(-2.948048) q[3];
sx q[3];
rz(1.5805472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.4291541) q[2];
sx q[2];
rz(-1.6753847) q[2];
sx q[2];
rz(-2.1377371) q[2];
rz(-2.6058274) q[3];
sx q[3];
rz(-2.2771213) q[3];
sx q[3];
rz(-0.0011477688) q[3];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67343229) q[0];
sx q[0];
rz(-2.4591481) q[0];
sx q[0];
rz(-2.4144507) q[0];
rz(-0.06761059) q[1];
sx q[1];
rz(-1.7865684) q[1];
sx q[1];
rz(-2.0179857) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1557187) q[0];
sx q[0];
rz(-2.576568) q[0];
sx q[0];
rz(-0.54604097) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4482165) q[2];
sx q[2];
rz(-1.5059587) q[2];
sx q[2];
rz(-1.4139869) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.1084845) q[1];
sx q[1];
rz(-0.29769167) q[1];
sx q[1];
rz(-1.9161403) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3109834) q[3];
sx q[3];
rz(-2.5334483) q[3];
sx q[3];
rz(0.45765314) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.2850538) q[2];
sx q[2];
rz(-1.5736138) q[2];
sx q[2];
rz(-0.48173586) q[2];
rz(2.5887865) q[3];
sx q[3];
rz(-2.063664) q[3];
sx q[3];
rz(-0.089574561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8168617) q[0];
sx q[0];
rz(-0.4751927) q[0];
sx q[0];
rz(3.0480296) q[0];
rz(-1.3803253) q[1];
sx q[1];
rz(-2.1880136) q[1];
sx q[1];
rz(-0.63724744) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4555511) q[0];
sx q[0];
rz(-1.0555869) q[0];
sx q[0];
rz(-1.0710708) q[0];
rz(-pi) q[1];
rz(2.5334444) q[2];
sx q[2];
rz(-1.7575193) q[2];
sx q[2];
rz(-3.1352459) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.67575021) q[1];
sx q[1];
rz(-1.9681962) q[1];
sx q[1];
rz(-0.028156412) q[1];
rz(-pi) q[2];
rz(-0.1430238) q[3];
sx q[3];
rz(-0.29353729) q[3];
sx q[3];
rz(1.0931726) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.20643413) q[2];
sx q[2];
rz(-0.55525246) q[2];
sx q[2];
rz(-0.68378249) q[2];
rz(-0.23009662) q[3];
sx q[3];
rz(-1.4068539) q[3];
sx q[3];
rz(-2.7195215) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.049659599) q[0];
sx q[0];
rz(-0.50685087) q[0];
sx q[0];
rz(1.4403213) q[0];
rz(-0.47538844) q[1];
sx q[1];
rz(-2.5071867) q[1];
sx q[1];
rz(2.5604274) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0345349) q[0];
sx q[0];
rz(-2.2904582) q[0];
sx q[0];
rz(-0.11373539) q[0];
rz(-pi) q[1];
rz(-1.1085132) q[2];
sx q[2];
rz(-1.7812742) q[2];
sx q[2];
rz(-1.1365979) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.5161203) q[1];
sx q[1];
rz(-0.56159084) q[1];
sx q[1];
rz(0.59188868) q[1];
rz(-2.851296) q[3];
sx q[3];
rz(-0.50945849) q[3];
sx q[3];
rz(-2.6313033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0873969) q[2];
sx q[2];
rz(-0.18614686) q[2];
sx q[2];
rz(2.6206214) q[2];
rz(-1.6446796) q[3];
sx q[3];
rz(-1.4119586) q[3];
sx q[3];
rz(-1.514667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5533376) q[0];
sx q[0];
rz(-2.8659358) q[0];
sx q[0];
rz(2.0113373) q[0];
rz(2.0626119) q[1];
sx q[1];
rz(-1.1993473) q[1];
sx q[1];
rz(-2.030453) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2030226) q[0];
sx q[0];
rz(-2.6513854) q[0];
sx q[0];
rz(1.9287648) q[0];
rz(-pi) q[1];
rz(0.07569282) q[2];
sx q[2];
rz(-0.53097979) q[2];
sx q[2];
rz(-2.7493283) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.4823522) q[1];
sx q[1];
rz(-0.9352881) q[1];
sx q[1];
rz(-0.19584943) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2419315) q[3];
sx q[3];
rz(-2.3658345) q[3];
sx q[3];
rz(-2.69913) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.40372103) q[2];
sx q[2];
rz(-0.51509866) q[2];
sx q[2];
rz(2.1007288) q[2];
rz(-0.8199842) q[3];
sx q[3];
rz(-1.9085725) q[3];
sx q[3];
rz(-0.6238873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2417004) q[0];
sx q[0];
rz(-0.59881678) q[0];
sx q[0];
rz(-0.65648055) q[0];
rz(1.3735324) q[1];
sx q[1];
rz(-0.7936002) q[1];
sx q[1];
rz(-2.0097282) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29378018) q[0];
sx q[0];
rz(-2.4180857) q[0];
sx q[0];
rz(2.840994) q[0];
x q[1];
rz(-0.64986367) q[2];
sx q[2];
rz(-2.3520326) q[2];
sx q[2];
rz(0.87433483) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.060185) q[1];
sx q[1];
rz(-0.8601195) q[1];
sx q[1];
rz(-3.1050472) q[1];
rz(-pi) q[2];
rz(0.15040654) q[3];
sx q[3];
rz(-1.5664706) q[3];
sx q[3];
rz(-1.9459917) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.446283) q[2];
sx q[2];
rz(-0.60574564) q[2];
sx q[2];
rz(0.31965762) q[2];
rz(0.22932912) q[3];
sx q[3];
rz(-1.0234443) q[3];
sx q[3];
rz(2.2314609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0685773) q[0];
sx q[0];
rz(-1.9972766) q[0];
sx q[0];
rz(-1.9101494) q[0];
rz(3.0629509) q[1];
sx q[1];
rz(-1.4631203) q[1];
sx q[1];
rz(0.15422779) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32145661) q[0];
sx q[0];
rz(-0.38561441) q[0];
sx q[0];
rz(-0.48954757) q[0];
x q[1];
rz(-0.94701887) q[2];
sx q[2];
rz(-0.49206054) q[2];
sx q[2];
rz(0.90384441) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.63111607) q[1];
sx q[1];
rz(-2.0492024) q[1];
sx q[1];
rz(1.9892742) q[1];
x q[2];
rz(2.7307512) q[3];
sx q[3];
rz(-1.7663398) q[3];
sx q[3];
rz(-1.1943749) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.8057033) q[2];
sx q[2];
rz(-1.5259537) q[2];
sx q[2];
rz(1.4914782) q[2];
rz(-1.7283745) q[3];
sx q[3];
rz(-1.7468942) q[3];
sx q[3];
rz(-1.570805) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0989646) q[0];
sx q[0];
rz(-1.093981) q[0];
sx q[0];
rz(-1.4549103) q[0];
rz(-2.1658354) q[1];
sx q[1];
rz(-1.059633) q[1];
sx q[1];
rz(1.8029433) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.094269116) q[0];
sx q[0];
rz(-1.9587095) q[0];
sx q[0];
rz(0.83705018) q[0];
rz(-pi) q[1];
rz(-0.33390121) q[2];
sx q[2];
rz(-1.0907764) q[2];
sx q[2];
rz(1.0261818) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.2690791) q[1];
sx q[1];
rz(-0.95035205) q[1];
sx q[1];
rz(0.2618813) q[1];
rz(-pi) q[2];
rz(-2.4782031) q[3];
sx q[3];
rz(-1.8691392) q[3];
sx q[3];
rz(1.364352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.51155382) q[2];
sx q[2];
rz(-1.0147107) q[2];
sx q[2];
rz(2.3021728) q[2];
rz(-1.9073585) q[3];
sx q[3];
rz(-2.0525565) q[3];
sx q[3];
rz(-0.031410005) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79427528) q[0];
sx q[0];
rz(-1.3823771) q[0];
sx q[0];
rz(-2.7224702) q[0];
rz(1.6436815) q[1];
sx q[1];
rz(-1.7828015) q[1];
sx q[1];
rz(0.43325123) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99892761) q[0];
sx q[0];
rz(-1.6892528) q[0];
sx q[0];
rz(-1.1112369) q[0];
rz(-pi) q[1];
x q[1];
rz(1.584132) q[2];
sx q[2];
rz(-2.2277701) q[2];
sx q[2];
rz(1.8767534) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.309714) q[1];
sx q[1];
rz(-1.0597178) q[1];
sx q[1];
rz(-0.32541497) q[1];
x q[2];
rz(-2.8178704) q[3];
sx q[3];
rz(-2.5383213) q[3];
sx q[3];
rz(-2.3761185) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.7822632) q[2];
sx q[2];
rz(-1.9908345) q[2];
sx q[2];
rz(-0.19598728) q[2];
rz(-1.1228784) q[3];
sx q[3];
rz(-2.4298318) q[3];
sx q[3];
rz(-1.6897197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
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
rz(2.6561683) q[0];
sx q[0];
rz(-0.66101414) q[0];
sx q[0];
rz(1.0957023) q[0];
rz(2.6307259) q[1];
sx q[1];
rz(-0.26900649) q[1];
sx q[1];
rz(2.9313472) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73186648) q[0];
sx q[0];
rz(-1.7667701) q[0];
sx q[0];
rz(2.1541697) q[0];
x q[1];
rz(-0.11504193) q[2];
sx q[2];
rz(-0.91918531) q[2];
sx q[2];
rz(0.63331214) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(3.0970891) q[1];
sx q[1];
rz(-2.0538446) q[1];
sx q[1];
rz(1.0853069) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5737567) q[3];
sx q[3];
rz(-0.68151268) q[3];
sx q[3];
rz(-0.62076724) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.11253396) q[2];
sx q[2];
rz(-1.6363279) q[2];
sx q[2];
rz(-1.0851592) q[2];
rz(-2.8339913) q[3];
sx q[3];
rz(-2.8234973) q[3];
sx q[3];
rz(-0.65348452) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45162421) q[0];
sx q[0];
rz(-1.1545447) q[0];
sx q[0];
rz(1.92323) q[0];
rz(-2.4901509) q[1];
sx q[1];
rz(-2.1919498) q[1];
sx q[1];
rz(0.776074) q[1];
rz(2.2617359) q[2];
sx q[2];
rz(-1.5590038) q[2];
sx q[2];
rz(-1.8270983) q[2];
rz(0.93919803) q[3];
sx q[3];
rz(-1.9892577) q[3];
sx q[3];
rz(1.6792959) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
