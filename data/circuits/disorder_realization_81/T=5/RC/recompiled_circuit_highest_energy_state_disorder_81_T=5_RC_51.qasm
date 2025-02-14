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
rz(-1.1377347) q[0];
rz(-3.0955834) q[1];
sx q[1];
rz(-2.6061821) q[1];
sx q[1];
rz(-1.6028264) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3483602) q[0];
sx q[0];
rz(-1.6038415) q[0];
sx q[0];
rz(-1.6152302) q[0];
rz(-pi) q[1];
rz(0.85311546) q[2];
sx q[2];
rz(-1.4783646) q[2];
sx q[2];
rz(1.5866367) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.5482169) q[1];
sx q[1];
rz(-1.1187828) q[1];
sx q[1];
rz(-1.7991245) q[1];
rz(-pi) q[2];
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
rz(-2.7988269) q[2];
sx q[2];
rz(2.9052367) q[2];
rz(2.6929839) q[3];
sx q[3];
rz(-1.4834504) q[3];
sx q[3];
rz(-1.0033222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
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
rz(-2.3568514) q[1];
sx q[1];
rz(-1.4737543) q[1];
sx q[1];
rz(0.10202185) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52962063) q[0];
sx q[0];
rz(-2.1054322) q[0];
sx q[0];
rz(-1.8055339) q[0];
x q[1];
rz(1.2743852) q[2];
sx q[2];
rz(-0.80889946) q[2];
sx q[2];
rz(-2.1975225) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.35013546) q[1];
sx q[1];
rz(-1.268549) q[1];
sx q[1];
rz(1.2863897) q[1];
rz(-pi) q[2];
rz(2.3524093) q[3];
sx q[3];
rz(-1.7459646) q[3];
sx q[3];
rz(0.096668124) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.1138136) q[2];
sx q[2];
rz(-1.4872097) q[2];
sx q[2];
rz(-0.2629183) q[2];
rz(1.8391838) q[3];
sx q[3];
rz(-1.115256) q[3];
sx q[3];
rz(-0.7836248) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1394434) q[0];
sx q[0];
rz(-0.0037010598) q[0];
sx q[0];
rz(1.333746) q[0];
rz(1.3847146) q[1];
sx q[1];
rz(-2.2127071) q[1];
sx q[1];
rz(0.76470107) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.156835) q[0];
sx q[0];
rz(-1.3768702) q[0];
sx q[0];
rz(0.24389275) q[0];
x q[1];
rz(2.7638859) q[2];
sx q[2];
rz(-1.3481082) q[2];
sx q[2];
rz(-1.4770222) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.0955268) q[1];
sx q[1];
rz(-0.30991677) q[1];
sx q[1];
rz(0.97280963) q[1];
rz(-2.8370492) q[3];
sx q[3];
rz(-2.2650654) q[3];
sx q[3];
rz(0.77724953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.9646405) q[2];
sx q[2];
rz(-1.2531589) q[2];
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
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95530987) q[0];
sx q[0];
rz(-2.0109542) q[0];
sx q[0];
rz(2.8593707) q[0];
rz(1.9918293) q[1];
sx q[1];
rz(-1.3590004) q[1];
sx q[1];
rz(1.6414292) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.041752664) q[0];
sx q[0];
rz(-0.40099335) q[0];
sx q[0];
rz(-2.6869511) q[0];
rz(-0.25441443) q[2];
sx q[2];
rz(-2.1037648) q[2];
sx q[2];
rz(1.0889458) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7042) q[1];
sx q[1];
rz(-1.3199894) q[1];
sx q[1];
rz(1.2247242) q[1];
rz(-pi) q[2];
rz(0.30434609) q[3];
sx q[3];
rz(-1.658421) q[3];
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
rz(2.2516294) q[3];
sx q[3];
rz(-1.7421937) q[3];
sx q[3];
rz(-2.9922805) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99059659) q[0];
sx q[0];
rz(-0.52495933) q[0];
sx q[0];
rz(1.4599482) q[0];
rz(-1.2851985) q[1];
sx q[1];
rz(-1.3672914) q[1];
sx q[1];
rz(-0.45305124) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1465587) q[0];
sx q[0];
rz(-1.5164598) q[0];
sx q[0];
rz(-1.8749692) q[0];
rz(-pi) q[1];
x q[1];
rz(0.98840752) q[2];
sx q[2];
rz(-1.526579) q[2];
sx q[2];
rz(2.1427936) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.0041607) q[1];
sx q[1];
rz(-2.7212226) q[1];
sx q[1];
rz(-0.46868268) q[1];
rz(2.0208218) q[3];
sx q[3];
rz(-1.9450359) q[3];
sx q[3];
rz(0.14412021) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.65752658) q[2];
sx q[2];
rz(-0.47407293) q[2];
sx q[2];
rz(-1.5849812) q[2];
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
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22353657) q[0];
sx q[0];
rz(-0.99128381) q[0];
sx q[0];
rz(-1.9687442) q[0];
rz(2.6808443) q[1];
sx q[1];
rz(-1.2604424) q[1];
sx q[1];
rz(1.6430829) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7464179) q[0];
sx q[0];
rz(-0.9580637) q[0];
sx q[0];
rz(3.1406671) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8086964) q[2];
sx q[2];
rz(-1.82331) q[2];
sx q[2];
rz(2.0123002) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.15461309) q[1];
sx q[1];
rz(-1.8334532) q[1];
sx q[1];
rz(-0.059537454) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.40703) q[3];
sx q[3];
rz(-0.70858817) q[3];
sx q[3];
rz(2.6276789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.1341165) q[2];
sx q[2];
rz(-1.4223998) q[2];
sx q[2];
rz(2.9134992) q[2];
rz(-0.54276931) q[3];
sx q[3];
rz(-2.1398862) q[3];
sx q[3];
rz(-0.91226474) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11693624) q[0];
sx q[0];
rz(-2.0967364) q[0];
sx q[0];
rz(2.5352617) q[0];
rz(1.8527276) q[1];
sx q[1];
rz(-1.6090798) q[1];
sx q[1];
rz(0.85561633) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0528078) q[0];
sx q[0];
rz(-1.1644378) q[0];
sx q[0];
rz(0.64865048) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4855723) q[2];
sx q[2];
rz(-2.3422675) q[2];
sx q[2];
rz(0.34121379) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.7164968) q[1];
sx q[1];
rz(-1.7762868) q[1];
sx q[1];
rz(2.53611) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8884573) q[3];
sx q[3];
rz(-2.0380424) q[3];
sx q[3];
rz(2.1727891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.764708) q[2];
sx q[2];
rz(-0.43262425) q[2];
sx q[2];
rz(3.063859) q[2];
rz(-3.0149095) q[3];
sx q[3];
rz(-2.099791) q[3];
sx q[3];
rz(0.83109394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34778255) q[0];
sx q[0];
rz(-0.3599444) q[0];
sx q[0];
rz(-0.98156324) q[0];
rz(-1.3579824) q[1];
sx q[1];
rz(-0.49109083) q[1];
sx q[1];
rz(0.29409274) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61035778) q[0];
sx q[0];
rz(-1.2623275) q[0];
sx q[0];
rz(1.5081132) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7693439) q[2];
sx q[2];
rz(-1.3073834) q[2];
sx q[2];
rz(-1.5579223) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.31042415) q[1];
sx q[1];
rz(-0.71485315) q[1];
sx q[1];
rz(-1.1352886) q[1];
x q[2];
rz(-0.96947396) q[3];
sx q[3];
rz(-3.0532928) q[3];
sx q[3];
rz(1.8661608) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(3.0391417) q[2];
sx q[2];
rz(-2.4852018) q[2];
sx q[2];
rz(2.0810818) q[2];
rz(-2.5833526) q[3];
sx q[3];
rz(-0.57359901) q[3];
sx q[3];
rz(0.019088117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3625665) q[0];
sx q[0];
rz(-2.9988852) q[0];
sx q[0];
rz(2.3574164) q[0];
rz(-1.7991637) q[1];
sx q[1];
rz(-1.289184) q[1];
sx q[1];
rz(-0.69650355) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2084889) q[0];
sx q[0];
rz(-1.7835) q[0];
sx q[0];
rz(1.8038294) q[0];
rz(-pi) q[1];
rz(-2.0663459) q[2];
sx q[2];
rz(-2.181567) q[2];
sx q[2];
rz(1.3411759) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.99217691) q[1];
sx q[1];
rz(-1.0806298) q[1];
sx q[1];
rz(1.0926682) q[1];
rz(3.0075865) q[3];
sx q[3];
rz(-1.7184765) q[3];
sx q[3];
rz(2.0509843) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.8037618) q[2];
sx q[2];
rz(-2.3904114) q[2];
sx q[2];
rz(1.2639698) q[2];
rz(-1.7654644) q[3];
sx q[3];
rz(-2.2270146) q[3];
sx q[3];
rz(1.5553364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92397583) q[0];
sx q[0];
rz(-3.0986077) q[0];
sx q[0];
rz(1.6643583) q[0];
rz(0.90786511) q[1];
sx q[1];
rz(-1.5217179) q[1];
sx q[1];
rz(0.97000617) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15031397) q[0];
sx q[0];
rz(-1.3320001) q[0];
sx q[0];
rz(-2.1697609) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1625948) q[2];
sx q[2];
rz(-1.3221957) q[2];
sx q[2];
rz(-0.73038855) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.59040356) q[1];
sx q[1];
rz(-1.7267372) q[1];
sx q[1];
rz(2.5124418) q[1];
x q[2];
rz(0.72204263) q[3];
sx q[3];
rz(-2.0835428) q[3];
sx q[3];
rz(-1.7554612) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.9475391) q[2];
sx q[2];
rz(-2.2668362) q[2];
sx q[2];
rz(-2.8694895) q[2];
rz(2.2943606) q[3];
sx q[3];
rz(-2.6341485) q[3];
sx q[3];
rz(-1.0890755) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(2.8027986) q[1];
sx q[1];
rz(-1.6962961) q[1];
sx q[1];
rz(1.2512339) q[1];
rz(1.9235545) q[2];
sx q[2];
rz(-0.88211664) q[2];
sx q[2];
rz(1.2695352) q[2];
rz(1.7539514) q[3];
sx q[3];
rz(-2.7121501) q[3];
sx q[3];
rz(0.6482758) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
