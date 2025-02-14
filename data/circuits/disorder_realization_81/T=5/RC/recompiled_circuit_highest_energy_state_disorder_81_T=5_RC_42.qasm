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
rz(-2.6061821) q[1];
sx q[1];
rz(-1.6028264) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7249776) q[0];
sx q[0];
rz(-3.0862245) q[0];
sx q[0];
rz(-0.93101089) q[0];
x q[1];
rz(-2.2884772) q[2];
sx q[2];
rz(-1.4783646) q[2];
sx q[2];
rz(-1.554956) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5482169) q[1];
sx q[1];
rz(-2.0228099) q[1];
sx q[1];
rz(1.3424681) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7102431) q[3];
sx q[3];
rz(-1.5155503) q[3];
sx q[3];
rz(1.6193438) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.5753182) q[2];
sx q[2];
rz(-0.34276572) q[2];
sx q[2];
rz(-0.23635593) q[2];
rz(-0.44860873) q[3];
sx q[3];
rz(-1.6581422) q[3];
sx q[3];
rz(-2.1382704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7597294) q[0];
sx q[0];
rz(-0.48668447) q[0];
sx q[0];
rz(2.3554262) q[0];
rz(0.78474125) q[1];
sx q[1];
rz(-1.4737543) q[1];
sx q[1];
rz(-3.0395708) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.090793153) q[0];
sx q[0];
rz(-2.5623119) q[0];
sx q[0];
rz(-2.7673278) q[0];
rz(-0.29709105) q[2];
sx q[2];
rz(-0.80657437) q[2];
sx q[2];
rz(1.7809803) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.7914572) q[1];
sx q[1];
rz(-1.8730436) q[1];
sx q[1];
rz(1.8552029) q[1];
rz(0.78918334) q[3];
sx q[3];
rz(-1.7459646) q[3];
sx q[3];
rz(-0.096668124) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.027779) q[2];
sx q[2];
rz(-1.4872097) q[2];
sx q[2];
rz(-0.2629183) q[2];
rz(-1.8391838) q[3];
sx q[3];
rz(-2.0263367) q[3];
sx q[3];
rz(-0.7836248) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1394434) q[0];
sx q[0];
rz(-0.0037010598) q[0];
sx q[0];
rz(1.333746) q[0];
rz(-1.7568781) q[1];
sx q[1];
rz(-2.2127071) q[1];
sx q[1];
rz(-2.3768916) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7755505) q[0];
sx q[0];
rz(-1.3315655) q[0];
sx q[0];
rz(-1.7704829) q[0];
rz(-pi) q[1];
x q[1];
rz(0.5506634) q[2];
sx q[2];
rz(-0.43573365) q[2];
sx q[2];
rz(0.41415641) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.0908392) q[1];
sx q[1];
rz(-1.3982441) q[1];
sx q[1];
rz(-1.3120632) q[1];
rz(-pi) q[2];
rz(2.2883203) q[3];
sx q[3];
rz(-1.3382592) q[3];
sx q[3];
rz(-0.59508376) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.9646405) q[2];
sx q[2];
rz(-1.2531589) q[2];
sx q[2];
rz(2.9580252) q[2];
rz(0.14487264) q[3];
sx q[3];
rz(-1.262007) q[3];
sx q[3];
rz(0.22411331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1862828) q[0];
sx q[0];
rz(-1.1306385) q[0];
sx q[0];
rz(-2.8593707) q[0];
rz(-1.1497633) q[1];
sx q[1];
rz(-1.7825922) q[1];
sx q[1];
rz(1.5001635) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44631413) q[0];
sx q[0];
rz(-1.9290961) q[0];
sx q[0];
rz(-1.3867239) q[0];
rz(-1.167539) q[2];
sx q[2];
rz(-0.58524281) q[2];
sx q[2];
rz(2.5257021) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.044149787) q[1];
sx q[1];
rz(-1.9056093) q[1];
sx q[1];
rz(2.8756932) q[1];
rz(-pi) q[2];
rz(-0.28516523) q[3];
sx q[3];
rz(-2.8252606) q[3];
sx q[3];
rz(2.0228937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.1993316) q[2];
sx q[2];
rz(-1.4456238) q[2];
sx q[2];
rz(-1.4873827) q[2];
rz(0.88996327) q[3];
sx q[3];
rz(-1.3993989) q[3];
sx q[3];
rz(0.14931211) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1509961) q[0];
sx q[0];
rz(-2.6166333) q[0];
sx q[0];
rz(1.4599482) q[0];
rz(-1.8563942) q[1];
sx q[1];
rz(-1.3672914) q[1];
sx q[1];
rz(0.45305124) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59555028) q[0];
sx q[0];
rz(-0.30883967) q[0];
sx q[0];
rz(-1.3911535) q[0];
rz(-0.052930367) q[2];
sx q[2];
rz(-0.9890511) q[2];
sx q[2];
rz(-0.60110053) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.86634163) q[1];
sx q[1];
rz(-1.3853933) q[1];
sx q[1];
rz(2.7621108) q[1];
x q[2];
rz(-2.3052081) q[3];
sx q[3];
rz(-0.57698133) q[3];
sx q[3];
rz(1.0671237) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.65752658) q[2];
sx q[2];
rz(-2.6675197) q[2];
sx q[2];
rz(-1.5566114) q[2];
rz(1.5629684) q[3];
sx q[3];
rz(-1.862674) q[3];
sx q[3];
rz(1.0274308) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22353657) q[0];
sx q[0];
rz(-0.99128381) q[0];
sx q[0];
rz(1.1728485) q[0];
rz(2.6808443) q[1];
sx q[1];
rz(-1.2604424) q[1];
sx q[1];
rz(-1.4985098) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1750893) q[0];
sx q[0];
rz(-1.5715535) q[0];
sx q[0];
rz(-2.1835292) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.837311) q[2];
sx q[2];
rz(-1.892748) q[2];
sx q[2];
rz(0.52767524) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.071049404) q[1];
sx q[1];
rz(-0.26916801) q[1];
sx q[1];
rz(-1.7885923) q[1];
rz(-3.0027578) q[3];
sx q[3];
rz(-0.87362008) q[3];
sx q[3];
rz(-2.4133854) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.1341165) q[2];
sx q[2];
rz(-1.4223998) q[2];
sx q[2];
rz(-0.22809347) q[2];
rz(-2.5988233) q[3];
sx q[3];
rz(-2.1398862) q[3];
sx q[3];
rz(-2.2293279) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0246564) q[0];
sx q[0];
rz(-2.0967364) q[0];
sx q[0];
rz(2.5352617) q[0];
rz(1.2888651) q[1];
sx q[1];
rz(-1.6090798) q[1];
sx q[1];
rz(2.2859763) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0528078) q[0];
sx q[0];
rz(-1.9771549) q[0];
sx q[0];
rz(0.64865048) q[0];
rz(2.3683041) q[2];
sx q[2];
rz(-1.6318562) q[2];
sx q[2];
rz(-1.1700912) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.2825477) q[1];
sx q[1];
rz(-2.5063516) q[1];
sx q[1];
rz(2.7905491) q[1];
x q[2];
rz(1.0904013) q[3];
sx q[3];
rz(-1.7962991) q[3];
sx q[3];
rz(0.71798872) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3768846) q[2];
sx q[2];
rz(-0.43262425) q[2];
sx q[2];
rz(-0.077733668) q[2];
rz(-3.0149095) q[3];
sx q[3];
rz(-2.099791) q[3];
sx q[3];
rz(0.83109394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7938101) q[0];
sx q[0];
rz(-2.7816483) q[0];
sx q[0];
rz(-2.1600294) q[0];
rz(1.7836102) q[1];
sx q[1];
rz(-2.6505018) q[1];
sx q[1];
rz(-0.29409274) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61035778) q[0];
sx q[0];
rz(-1.8792651) q[0];
sx q[0];
rz(-1.5081132) q[0];
rz(1.852599) q[2];
sx q[2];
rz(-1.929612) q[2];
sx q[2];
rz(-0.11419088) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.24178019) q[1];
sx q[1];
rz(-0.93437785) q[1];
sx q[1];
rz(-2.7905725) q[1];
x q[2];
rz(1.4979248) q[3];
sx q[3];
rz(-1.5208863) q[3];
sx q[3];
rz(-2.2467256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.10245094) q[2];
sx q[2];
rz(-2.4852018) q[2];
sx q[2];
rz(-2.0810818) q[2];
rz(-2.5833526) q[3];
sx q[3];
rz(-2.5679936) q[3];
sx q[3];
rz(-0.019088117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7790262) q[0];
sx q[0];
rz(-0.14270742) q[0];
sx q[0];
rz(2.3574164) q[0];
rz(-1.7991637) q[1];
sx q[1];
rz(-1.289184) q[1];
sx q[1];
rz(-0.69650355) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41236988) q[0];
sx q[0];
rz(-1.3431088) q[0];
sx q[0];
rz(-2.9231637) q[0];
rz(-pi) q[1];
rz(2.5449341) q[2];
sx q[2];
rz(-2.3755529) q[2];
sx q[2];
rz(-0.5853563) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3236774) q[1];
sx q[1];
rz(-1.9888048) q[1];
sx q[1];
rz(0.54115414) q[1];
rz(-pi) q[2];
x q[2];
rz(2.302565) q[3];
sx q[3];
rz(-2.9425042) q[3];
sx q[3];
rz(1.3092878) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.33783087) q[2];
sx q[2];
rz(-0.75118128) q[2];
sx q[2];
rz(-1.2639698) q[2];
rz(-1.7654644) q[3];
sx q[3];
rz(-0.91457808) q[3];
sx q[3];
rz(-1.5553364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2176168) q[0];
sx q[0];
rz(-0.042984977) q[0];
sx q[0];
rz(-1.6643583) q[0];
rz(-2.2337275) q[1];
sx q[1];
rz(-1.5217179) q[1];
sx q[1];
rz(0.97000617) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5805626) q[0];
sx q[0];
rz(-2.1505158) q[0];
sx q[0];
rz(-2.8549578) q[0];
rz(-pi) q[1];
rz(2.8447609) q[2];
sx q[2];
rz(-2.1420711) q[2];
sx q[2];
rz(-0.67650696) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.1906311) q[1];
sx q[1];
rz(-2.4959557) q[1];
sx q[1];
rz(-0.26107045) q[1];
x q[2];
rz(-2.41955) q[3];
sx q[3];
rz(-1.0580499) q[3];
sx q[3];
rz(-1.3861314) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.19405356) q[2];
sx q[2];
rz(-2.2668362) q[2];
sx q[2];
rz(0.27210316) q[2];
rz(-2.2943606) q[3];
sx q[3];
rz(-0.50744414) q[3];
sx q[3];
rz(-1.0890755) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67615164) q[0];
sx q[0];
rz(-2.0068824) q[0];
sx q[0];
rz(-1.6750499) q[0];
rz(-2.8027986) q[1];
sx q[1];
rz(-1.4452965) q[1];
sx q[1];
rz(-1.8903587) q[1];
rz(-0.39739824) q[2];
sx q[2];
rz(-0.76042475) q[2];
sx q[2];
rz(-2.3971192) q[2];
rz(0.083214464) q[3];
sx q[3];
rz(-1.1490001) q[3];
sx q[3];
rz(0.84925539) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
