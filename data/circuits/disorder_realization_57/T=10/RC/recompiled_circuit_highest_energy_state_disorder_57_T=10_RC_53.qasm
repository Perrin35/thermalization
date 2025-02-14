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
rz(1.1433831) q[0];
sx q[0];
rz(-1.5902061) q[0];
sx q[0];
rz(-2.8144612) q[0];
rz(1.7786572) q[1];
sx q[1];
rz(-2.5442446) q[1];
sx q[1];
rz(-0.30244952) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10714794) q[0];
sx q[0];
rz(-1.7361138) q[0];
sx q[0];
rz(1.3470696) q[0];
rz(-pi) q[1];
rz(-2.2277432) q[2];
sx q[2];
rz(-2.3656332) q[2];
sx q[2];
rz(2.1405107) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.55284269) q[1];
sx q[1];
rz(-1.8859914) q[1];
sx q[1];
rz(-2.9957459) q[1];
rz(0.47306319) q[3];
sx q[3];
rz(-1.6636208) q[3];
sx q[3];
rz(1.0582093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.3737619) q[2];
sx q[2];
rz(-3.0475898) q[2];
sx q[2];
rz(2.8685699) q[2];
rz(2.7772969) q[3];
sx q[3];
rz(-1.7184869) q[3];
sx q[3];
rz(0.020615904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.453422) q[0];
sx q[0];
rz(-2.0437129) q[0];
sx q[0];
rz(-1.706644) q[0];
rz(-0.20225987) q[1];
sx q[1];
rz(-0.63019284) q[1];
sx q[1];
rz(0.30466255) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.019783171) q[0];
sx q[0];
rz(-2.2643548) q[0];
sx q[0];
rz(2.8163363) q[0];
rz(2.8424047) q[2];
sx q[2];
rz(-1.169489) q[2];
sx q[2];
rz(-0.68986675) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.4604291) q[1];
sx q[1];
rz(-2.7507456) q[1];
sx q[1];
rz(-0.50182287) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6913139) q[3];
sx q[3];
rz(-1.4637865) q[3];
sx q[3];
rz(-1.7086687) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.8649586) q[2];
sx q[2];
rz(-0.43143299) q[2];
sx q[2];
rz(-1.6032093) q[2];
rz(-2.3254584) q[3];
sx q[3];
rz(-0.82681257) q[3];
sx q[3];
rz(-0.89474595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64033878) q[0];
sx q[0];
rz(-2.8479939) q[0];
sx q[0];
rz(-1.5119934) q[0];
rz(-1.8005796) q[1];
sx q[1];
rz(-0.87313849) q[1];
sx q[1];
rz(0.27892932) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19353774) q[0];
sx q[0];
rz(-2.9729261) q[0];
sx q[0];
rz(0.70341603) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.52742676) q[2];
sx q[2];
rz(-2.1373539) q[2];
sx q[2];
rz(-1.1184831) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3027569) q[1];
sx q[1];
rz(-2.347337) q[1];
sx q[1];
rz(1.1357186) q[1];
rz(-pi) q[2];
rz(-1.5666714) q[3];
sx q[3];
rz(-0.71625159) q[3];
sx q[3];
rz(-0.28820693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.5035847) q[2];
sx q[2];
rz(-1.4853442) q[2];
sx q[2];
rz(0.33835641) q[2];
rz(-1.407297) q[3];
sx q[3];
rz(-1.1275848) q[3];
sx q[3];
rz(0.97051632) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9426743) q[0];
sx q[0];
rz(-3.119454) q[0];
sx q[0];
rz(0.66377798) q[0];
rz(2.5860419) q[1];
sx q[1];
rz(-1.5063565) q[1];
sx q[1];
rz(-1.3551855) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5252939) q[0];
sx q[0];
rz(-0.64726603) q[0];
sx q[0];
rz(-3.0891524) q[0];
x q[1];
rz(-1.6333601) q[2];
sx q[2];
rz(-0.35919562) q[2];
sx q[2];
rz(-2.0421093) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.04482567) q[1];
sx q[1];
rz(-1.191865) q[1];
sx q[1];
rz(2.1842537) q[1];
x q[2];
rz(-1.438645) q[3];
sx q[3];
rz(-2.6717917) q[3];
sx q[3];
rz(-0.93549226) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.3163471) q[2];
sx q[2];
rz(-2.1630042) q[2];
sx q[2];
rz(-0.52581954) q[2];
rz(0.60872269) q[3];
sx q[3];
rz(-0.63568297) q[3];
sx q[3];
rz(-0.70613247) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71259251) q[0];
sx q[0];
rz(-1.1845931) q[0];
sx q[0];
rz(2.1767966) q[0];
rz(-0.4110128) q[1];
sx q[1];
rz(-2.1844468) q[1];
sx q[1];
rz(-2.029665) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0135213) q[0];
sx q[0];
rz(-3.1134634) q[0];
sx q[0];
rz(-0.97814409) q[0];
rz(-pi) q[1];
rz(2.0093727) q[2];
sx q[2];
rz(-2.2242332) q[2];
sx q[2];
rz(-0.48023047) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.7605684) q[1];
sx q[1];
rz(-0.79729092) q[1];
sx q[1];
rz(-0.83177318) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6451177) q[3];
sx q[3];
rz(-2.8395445) q[3];
sx q[3];
rz(-0.99629096) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.3045706) q[2];
sx q[2];
rz(-1.3304173) q[2];
sx q[2];
rz(-2.0988317) q[2];
rz(1.1834831) q[3];
sx q[3];
rz(-0.8684929) q[3];
sx q[3];
rz(-2.1709501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17827621) q[0];
sx q[0];
rz(-1.9177508) q[0];
sx q[0];
rz(0.29603145) q[0];
rz(3.1405247) q[1];
sx q[1];
rz(-1.3384534) q[1];
sx q[1];
rz(-2.1086878) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4460239) q[0];
sx q[0];
rz(-3.1349413) q[0];
sx q[0];
rz(-2.5435996) q[0];
x q[1];
rz(-1.270711) q[2];
sx q[2];
rz(-0.3870766) q[2];
sx q[2];
rz(-2.84969) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.7958837) q[1];
sx q[1];
rz(-2.8879037) q[1];
sx q[1];
rz(-0.2233875) q[1];
rz(2.7240678) q[3];
sx q[3];
rz(-2.333765) q[3];
sx q[3];
rz(-2.1958627) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.380015) q[2];
sx q[2];
rz(-1.4662611) q[2];
sx q[2];
rz(-1.638691) q[2];
rz(2.5401529) q[3];
sx q[3];
rz(-0.99965874) q[3];
sx q[3];
rz(-2.7433266) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10694207) q[0];
sx q[0];
rz(-1.5926188) q[0];
sx q[0];
rz(-2.4079127) q[0];
rz(-0.94379464) q[1];
sx q[1];
rz(-1.5422041) q[1];
sx q[1];
rz(-2.5234047) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3933936) q[0];
sx q[0];
rz(-1.3915147) q[0];
sx q[0];
rz(-0.050672942) q[0];
x q[1];
rz(-2.3011977) q[2];
sx q[2];
rz(-1.5638509) q[2];
sx q[2];
rz(-0.23887979) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.427278) q[1];
sx q[1];
rz(-0.66525092) q[1];
sx q[1];
rz(-1.6131496) q[1];
x q[2];
rz(0.47241601) q[3];
sx q[3];
rz(-0.69184408) q[3];
sx q[3];
rz(1.7418246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.9435297) q[2];
sx q[2];
rz(-1.1897831) q[2];
sx q[2];
rz(0.61723906) q[2];
rz(-1.9700358) q[3];
sx q[3];
rz(-0.37552437) q[3];
sx q[3];
rz(0.61081162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1307369) q[0];
sx q[0];
rz(-2.3541088) q[0];
sx q[0];
rz(-2.6686344) q[0];
rz(-1.4594151) q[1];
sx q[1];
rz(-1.727203) q[1];
sx q[1];
rz(0.032729538) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6115295) q[0];
sx q[0];
rz(-0.034688799) q[0];
sx q[0];
rz(-2.7682224) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1104537) q[2];
sx q[2];
rz(-2.2810069) q[2];
sx q[2];
rz(-1.4905765) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.53712691) q[1];
sx q[1];
rz(-2.6422738) q[1];
sx q[1];
rz(2.0313655) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2975912) q[3];
sx q[3];
rz(-3.0298067) q[3];
sx q[3];
rz(2.8120086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.1590283) q[2];
sx q[2];
rz(-0.50243598) q[2];
sx q[2];
rz(2.0459797) q[2];
rz(-0.76121965) q[3];
sx q[3];
rz(-1.8571564) q[3];
sx q[3];
rz(0.42263862) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.465076) q[0];
sx q[0];
rz(-3.0618771) q[0];
sx q[0];
rz(-0.7088784) q[0];
rz(2.830016) q[1];
sx q[1];
rz(-2.0496924) q[1];
sx q[1];
rz(-2.8332205) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.209359) q[0];
sx q[0];
rz(-1.4054125) q[0];
sx q[0];
rz(-0.39914058) q[0];
rz(-0.73268415) q[2];
sx q[2];
rz(-2.6044629) q[2];
sx q[2];
rz(2.2718475) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.2110097) q[1];
sx q[1];
rz(-0.66291729) q[1];
sx q[1];
rz(2.7410313) q[1];
rz(-pi) q[2];
x q[2];
rz(0.25777581) q[3];
sx q[3];
rz(-2.3376102) q[3];
sx q[3];
rz(-2.3641158) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.10820216) q[2];
sx q[2];
rz(-1.4806925) q[2];
sx q[2];
rz(-2.9208276) q[2];
rz(1.0734142) q[3];
sx q[3];
rz(-0.99075166) q[3];
sx q[3];
rz(-2.3555135) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49755001) q[0];
sx q[0];
rz(-2.2798517) q[0];
sx q[0];
rz(1.0368689) q[0];
rz(1.6715624) q[1];
sx q[1];
rz(-1.0187047) q[1];
sx q[1];
rz(1.9482313) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65914175) q[0];
sx q[0];
rz(-1.5406784) q[0];
sx q[0];
rz(2.6994266) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.5035053) q[2];
sx q[2];
rz(-0.76447884) q[2];
sx q[2];
rz(-1.6703005) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.5037704) q[1];
sx q[1];
rz(-0.71412266) q[1];
sx q[1];
rz(-2.8526866) q[1];
rz(-1.3383629) q[3];
sx q[3];
rz(-1.4475736) q[3];
sx q[3];
rz(-0.49613813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.6679473) q[2];
sx q[2];
rz(-1.2845984) q[2];
sx q[2];
rz(-3.0038707) q[2];
rz(-1.5406476) q[3];
sx q[3];
rz(-2.9100304) q[3];
sx q[3];
rz(-2.2102977) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66844983) q[0];
sx q[0];
rz(-1.4746329) q[0];
sx q[0];
rz(2.0869577) q[0];
rz(2.6558381) q[1];
sx q[1];
rz(-1.2243441) q[1];
sx q[1];
rz(-1.4996554) q[1];
rz(-0.68436868) q[2];
sx q[2];
rz(-1.3791313) q[2];
sx q[2];
rz(-2.287938) q[2];
rz(1.0578591) q[3];
sx q[3];
rz(-0.87659642) q[3];
sx q[3];
rz(2.07758) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
