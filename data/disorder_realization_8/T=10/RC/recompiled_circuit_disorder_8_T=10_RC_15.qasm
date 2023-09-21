OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.8060057) q[0];
sx q[0];
rz(-0.94526362) q[0];
sx q[0];
rz(2.6160016) q[0];
rz(0.2431915) q[1];
sx q[1];
rz(-1.9089729) q[1];
sx q[1];
rz(0.90484172) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4558976) q[0];
sx q[0];
rz(-2.064961) q[0];
sx q[0];
rz(-2.6463638) q[0];
rz(-pi) q[1];
x q[1];
rz(0.22613871) q[2];
sx q[2];
rz(-1.7625426) q[2];
sx q[2];
rz(0.5118256) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.2658087) q[1];
sx q[1];
rz(-1.335653) q[1];
sx q[1];
rz(1.1633412) q[1];
x q[2];
rz(-1.8815133) q[3];
sx q[3];
rz(-1.7508535) q[3];
sx q[3];
rz(2.3678399) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.20377542) q[2];
sx q[2];
rz(-1.4062466) q[2];
sx q[2];
rz(3.0453483) q[2];
rz(1.0359267) q[3];
sx q[3];
rz(-0.38714287) q[3];
sx q[3];
rz(-0.15371418) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44089833) q[0];
sx q[0];
rz(-2.7504524) q[0];
sx q[0];
rz(-0.76517117) q[0];
rz(-1.8493429) q[1];
sx q[1];
rz(-0.48520979) q[1];
sx q[1];
rz(0.66295019) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.457263) q[0];
sx q[0];
rz(-1.508679) q[0];
sx q[0];
rz(-1.6371884) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2234369) q[2];
sx q[2];
rz(-1.895972) q[2];
sx q[2];
rz(2.6174389) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.029164974) q[1];
sx q[1];
rz(-1.6532835) q[1];
sx q[1];
rz(-0.44177456) q[1];
x q[2];
rz(-0.48682537) q[3];
sx q[3];
rz(-1.7142222) q[3];
sx q[3];
rz(1.6008582) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.048916653) q[2];
sx q[2];
rz(-2.0930223) q[2];
sx q[2];
rz(0.32901397) q[2];
rz(-0.66550231) q[3];
sx q[3];
rz(-0.21829675) q[3];
sx q[3];
rz(-1.3177419) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5927785) q[0];
sx q[0];
rz(-2.9187262) q[0];
sx q[0];
rz(2.9192525) q[0];
rz(-1.0173343) q[1];
sx q[1];
rz(-0.72128123) q[1];
sx q[1];
rz(-2.6229048) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.737239) q[0];
sx q[0];
rz(-2.3271932) q[0];
sx q[0];
rz(-2.0703719) q[0];
rz(-pi) q[1];
rz(2.8754183) q[2];
sx q[2];
rz(-1.6632348) q[2];
sx q[2];
rz(-3.0150974) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.4244713) q[1];
sx q[1];
rz(-0.79454225) q[1];
sx q[1];
rz(0.24309991) q[1];
x q[2];
rz(-3.0869811) q[3];
sx q[3];
rz(-1.5699937) q[3];
sx q[3];
rz(-2.4459248) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.8609994) q[2];
sx q[2];
rz(-0.36182797) q[2];
sx q[2];
rz(-0.068543531) q[2];
rz(2.5391501) q[3];
sx q[3];
rz(-0.76255637) q[3];
sx q[3];
rz(-3.0025735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.65072) q[0];
sx q[0];
rz(-2.3644709) q[0];
sx q[0];
rz(-2.9673476) q[0];
rz(2.6113367) q[1];
sx q[1];
rz(-1.4825876) q[1];
sx q[1];
rz(-2.6285016) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78631567) q[0];
sx q[0];
rz(-2.516541) q[0];
sx q[0];
rz(3.0062208) q[0];
x q[1];
rz(2.1999173) q[2];
sx q[2];
rz(-2.2721014) q[2];
sx q[2];
rz(-0.47052449) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.2345703) q[1];
sx q[1];
rz(-1.7090194) q[1];
sx q[1];
rz(-0.45781086) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4245093) q[3];
sx q[3];
rz(-2.9330367) q[3];
sx q[3];
rz(1.2775161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.47485581) q[2];
sx q[2];
rz(-0.36617827) q[2];
sx q[2];
rz(-0.22988698) q[2];
rz(0.41904467) q[3];
sx q[3];
rz(-1.34904) q[3];
sx q[3];
rz(0.45927799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4693562) q[0];
sx q[0];
rz(-0.75119632) q[0];
sx q[0];
rz(-0.67681926) q[0];
rz(0.49304402) q[1];
sx q[1];
rz(-0.9489916) q[1];
sx q[1];
rz(-0.61606032) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6989243) q[0];
sx q[0];
rz(-1.5183503) q[0];
sx q[0];
rz(-3.1057182) q[0];
rz(0.36523833) q[2];
sx q[2];
rz(-1.4589981) q[2];
sx q[2];
rz(1.0527843) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.7580326) q[1];
sx q[1];
rz(-2.5639503) q[1];
sx q[1];
rz(-0.010096117) q[1];
rz(-pi) q[2];
x q[2];
rz(0.89415278) q[3];
sx q[3];
rz(-1.812462) q[3];
sx q[3];
rz(1.5622996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4074576) q[2];
sx q[2];
rz(-2.5871758) q[2];
sx q[2];
rz(2.8862254) q[2];
rz(-1.5363961) q[3];
sx q[3];
rz(-2.1891749) q[3];
sx q[3];
rz(0.77409625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7191294) q[0];
sx q[0];
rz(-0.96619773) q[0];
sx q[0];
rz(2.6690924) q[0];
rz(2.6155112) q[1];
sx q[1];
rz(-0.20985797) q[1];
sx q[1];
rz(2.2568259) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9677744) q[0];
sx q[0];
rz(-2.7533555) q[0];
sx q[0];
rz(-1.7675179) q[0];
rz(-pi) q[1];
rz(2.9315345) q[2];
sx q[2];
rz(-1.5569599) q[2];
sx q[2];
rz(-2.7163598) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.92237597) q[1];
sx q[1];
rz(-1.1616542) q[1];
sx q[1];
rz(-0.31422305) q[1];
rz(0.24598908) q[3];
sx q[3];
rz(-2.0912598) q[3];
sx q[3];
rz(-2.1465079) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.3970967) q[2];
sx q[2];
rz(-2.3366191) q[2];
sx q[2];
rz(-0.3113783) q[2];
rz(-1.3686251) q[3];
sx q[3];
rz(-2.6840648) q[3];
sx q[3];
rz(0.5293203) q[3];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41806528) q[0];
sx q[0];
rz(-0.27856809) q[0];
sx q[0];
rz(3.080522) q[0];
rz(-0.04018499) q[1];
sx q[1];
rz(-1.1611074) q[1];
sx q[1];
rz(0.73289245) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31156763) q[0];
sx q[0];
rz(-2.4550779) q[0];
sx q[0];
rz(0.65450432) q[0];
rz(-2.6635025) q[2];
sx q[2];
rz(-2.2149137) q[2];
sx q[2];
rz(2.9113876) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.7072308) q[1];
sx q[1];
rz(-1.5773298) q[1];
sx q[1];
rz(-2.1965501) q[1];
rz(-pi) q[2];
x q[2];
rz(1.93768) q[3];
sx q[3];
rz(-3.0266579) q[3];
sx q[3];
rz(0.68655187) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.0733033) q[2];
sx q[2];
rz(-1.9446334) q[2];
sx q[2];
rz(2.627009) q[2];
rz(1.2375281) q[3];
sx q[3];
rz(-2.5585744) q[3];
sx q[3];
rz(2.5966743) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.085389) q[0];
sx q[0];
rz(-1.4852925) q[0];
sx q[0];
rz(0.7094267) q[0];
rz(1.5052694) q[1];
sx q[1];
rz(-1.0737597) q[1];
sx q[1];
rz(2.8628796) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80517171) q[0];
sx q[0];
rz(-1.780691) q[0];
sx q[0];
rz(1.4324485) q[0];
rz(0.3785554) q[2];
sx q[2];
rz(-2.5494529) q[2];
sx q[2];
rz(-2.2373667) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.41234327) q[1];
sx q[1];
rz(-1.7885498) q[1];
sx q[1];
rz(-1.2530243) q[1];
rz(-1.9023499) q[3];
sx q[3];
rz(-2.9058876) q[3];
sx q[3];
rz(1.9086259) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8401106) q[2];
sx q[2];
rz(-1.297941) q[2];
sx q[2];
rz(0.056079496) q[2];
rz(-0.85514832) q[3];
sx q[3];
rz(-0.44848281) q[3];
sx q[3];
rz(-2.7364031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99496019) q[0];
sx q[0];
rz(-2.9652847) q[0];
sx q[0];
rz(-0.96889281) q[0];
rz(2.6682207) q[1];
sx q[1];
rz(-2.3469766) q[1];
sx q[1];
rz(1.1425346) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56209598) q[0];
sx q[0];
rz(-1.8651433) q[0];
sx q[0];
rz(-3.0466945) q[0];
rz(-1.5683453) q[2];
sx q[2];
rz(-2.2741389) q[2];
sx q[2];
rz(-3.0424812) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.6856319) q[1];
sx q[1];
rz(-2.5870442) q[1];
sx q[1];
rz(-0.23994069) q[1];
x q[2];
rz(1.8248796) q[3];
sx q[3];
rz(-0.68613201) q[3];
sx q[3];
rz(2.7568267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.2450976) q[2];
sx q[2];
rz(-1.2197887) q[2];
sx q[2];
rz(-2.1208105) q[2];
rz(-0.3237237) q[3];
sx q[3];
rz(-2.3886069) q[3];
sx q[3];
rz(3.1304205) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97994119) q[0];
sx q[0];
rz(-0.027898235) q[0];
sx q[0];
rz(-0.7014057) q[0];
rz(0.91570634) q[1];
sx q[1];
rz(-2.1332824) q[1];
sx q[1];
rz(1.2385626) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1322051) q[0];
sx q[0];
rz(-1.3552249) q[0];
sx q[0];
rz(-2.0995887) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6221223) q[2];
sx q[2];
rz(-0.81846279) q[2];
sx q[2];
rz(-2.3527956) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-3.0760127) q[1];
sx q[1];
rz(-0.91846839) q[1];
sx q[1];
rz(-0.14821649) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7391316) q[3];
sx q[3];
rz(-0.59909648) q[3];
sx q[3];
rz(1.7408016) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.68676585) q[2];
sx q[2];
rz(-1.0964311) q[2];
sx q[2];
rz(-3.0977541) q[2];
rz(-1.2010126) q[3];
sx q[3];
rz(-2.4062556) q[3];
sx q[3];
rz(2.1380077) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4713521) q[0];
sx q[0];
rz(-0.72605194) q[0];
sx q[0];
rz(-1.3656021) q[0];
rz(0.025370601) q[1];
sx q[1];
rz(-1.8352958) q[1];
sx q[1];
rz(-1.8713554) q[1];
rz(0.92924835) q[2];
sx q[2];
rz(-2.6570448) q[2];
sx q[2];
rz(-1.6397283) q[2];
rz(2.9712215) q[3];
sx q[3];
rz(-0.80633612) q[3];
sx q[3];
rz(2.7663305) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];