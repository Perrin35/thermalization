OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.7282038) q[0];
sx q[0];
rz(-2.0079186) q[0];
sx q[0];
rz(1.5490305) q[0];
rz(1.6917317) q[1];
sx q[1];
rz(-0.65728846) q[1];
sx q[1];
rz(-2.5974098) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.777594) q[0];
sx q[0];
rz(-1.6382123) q[0];
sx q[0];
rz(-1.6075587) q[0];
rz(-pi) q[1];
rz(-1.5266225) q[2];
sx q[2];
rz(-1.6720811) q[2];
sx q[2];
rz(-0.12461187) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.7218329) q[1];
sx q[1];
rz(-0.83336035) q[1];
sx q[1];
rz(-2.1268197) q[1];
rz(-0.55239001) q[3];
sx q[3];
rz(-0.033621764) q[3];
sx q[3];
rz(0.53510964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.071775285) q[2];
sx q[2];
rz(-1.8775619) q[2];
sx q[2];
rz(-1.7791746) q[2];
rz(-3.1133364) q[3];
sx q[3];
rz(-1.3794206) q[3];
sx q[3];
rz(-0.79022592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4966999) q[0];
sx q[0];
rz(-2.4364478) q[0];
sx q[0];
rz(-1.9702966) q[0];
rz(0.21121875) q[1];
sx q[1];
rz(-0.44208458) q[1];
sx q[1];
rz(1.404095) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1275741) q[0];
sx q[0];
rz(-1.9919792) q[0];
sx q[0];
rz(0.76571) q[0];
rz(2.6589083) q[2];
sx q[2];
rz(-1.4233372) q[2];
sx q[2];
rz(-0.54373103) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.5215056) q[1];
sx q[1];
rz(-1.5937935) q[1];
sx q[1];
rz(-0.28316811) q[1];
rz(-2.8849765) q[3];
sx q[3];
rz(-0.49391541) q[3];
sx q[3];
rz(-1.4790505) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.30963787) q[2];
sx q[2];
rz(-1.6654623) q[2];
sx q[2];
rz(0.94397604) q[2];
rz(-0.55654636) q[3];
sx q[3];
rz(-2.6440933) q[3];
sx q[3];
rz(-1.336162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8495162) q[0];
sx q[0];
rz(-2.3777666) q[0];
sx q[0];
rz(1.7180432) q[0];
rz(0.81958333) q[1];
sx q[1];
rz(-2.3312566) q[1];
sx q[1];
rz(-0.56366411) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99479988) q[0];
sx q[0];
rz(-2.1437862) q[0];
sx q[0];
rz(1.6550199) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.59994772) q[2];
sx q[2];
rz(-1.7581345) q[2];
sx q[2];
rz(1.7768605) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.015805294) q[1];
sx q[1];
rz(-1.8965221) q[1];
sx q[1];
rz(-2.159352) q[1];
rz(-0.54791252) q[3];
sx q[3];
rz(-0.69763819) q[3];
sx q[3];
rz(1.1988977) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.6391969) q[2];
sx q[2];
rz(-1.1221308) q[2];
sx q[2];
rz(-1.0602661) q[2];
rz(-0.075573102) q[3];
sx q[3];
rz(-2.2836756) q[3];
sx q[3];
rz(-3.0269572) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.32325) q[0];
sx q[0];
rz(-0.30775726) q[0];
sx q[0];
rz(0.19317214) q[0];
rz(-1.5974143) q[1];
sx q[1];
rz(-1.9924106) q[1];
sx q[1];
rz(-0.65778041) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13947978) q[0];
sx q[0];
rz(-1.6008899) q[0];
sx q[0];
rz(0.17480236) q[0];
rz(-1.1409608) q[2];
sx q[2];
rz(-1.8550711) q[2];
sx q[2];
rz(0.62361275) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(3.0187877) q[1];
sx q[1];
rz(-2.0969166) q[1];
sx q[1];
rz(2.8469574) q[1];
rz(-1.2781906) q[3];
sx q[3];
rz(-1.8183823) q[3];
sx q[3];
rz(-2.2729371) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.0956991) q[2];
sx q[2];
rz(-2.7313488) q[2];
sx q[2];
rz(2.0945385) q[2];
rz(-2.0783157) q[3];
sx q[3];
rz(-1.9789109) q[3];
sx q[3];
rz(1.8803546) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7970153) q[0];
sx q[0];
rz(-0.319096) q[0];
sx q[0];
rz(2.1176594) q[0];
rz(-0.78760415) q[1];
sx q[1];
rz(-1.5658295) q[1];
sx q[1];
rz(-0.62686282) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.244827) q[0];
sx q[0];
rz(-1.5875495) q[0];
sx q[0];
rz(-1.6003952) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3847694) q[2];
sx q[2];
rz(-0.42381091) q[2];
sx q[2];
rz(-0.28048453) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.53351952) q[1];
sx q[1];
rz(-0.43577172) q[1];
sx q[1];
rz(2.3418952) q[1];
rz(-pi) q[2];
rz(1.2792148) q[3];
sx q[3];
rz(-0.40816187) q[3];
sx q[3];
rz(2.6274031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.91288599) q[2];
sx q[2];
rz(-1.9084946) q[2];
sx q[2];
rz(-2.1234925) q[2];
rz(0.61156887) q[3];
sx q[3];
rz(-1.8194149) q[3];
sx q[3];
rz(-0.95156804) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6927032) q[0];
sx q[0];
rz(-1.3494116) q[0];
sx q[0];
rz(1.1992136) q[0];
rz(1.9723643) q[1];
sx q[1];
rz(-1.0511845) q[1];
sx q[1];
rz(0.32454023) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31774662) q[0];
sx q[0];
rz(-2.4577603) q[0];
sx q[0];
rz(-0.86156396) q[0];
x q[1];
rz(-2.3338823) q[2];
sx q[2];
rz(-1.0683904) q[2];
sx q[2];
rz(-3.0903357) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3970916) q[1];
sx q[1];
rz(-2.607064) q[1];
sx q[1];
rz(-1.6194653) q[1];
x q[2];
rz(-0.59431521) q[3];
sx q[3];
rz(-2.6404877) q[3];
sx q[3];
rz(1.2831812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.67363182) q[2];
sx q[2];
rz(-1.8362074) q[2];
sx q[2];
rz(1.8661873) q[2];
rz(2.0868789) q[3];
sx q[3];
rz(-0.79189363) q[3];
sx q[3];
rz(-1.5462497) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4642898) q[0];
sx q[0];
rz(-1.3339366) q[0];
sx q[0];
rz(-1.7328847) q[0];
rz(-2.6858221) q[1];
sx q[1];
rz(-2.9401638) q[1];
sx q[1];
rz(1.2021525) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5802655) q[0];
sx q[0];
rz(-2.1939477) q[0];
sx q[0];
rz(-2.6448963) q[0];
rz(-pi) q[1];
rz(-1.9583148) q[2];
sx q[2];
rz(-0.85867184) q[2];
sx q[2];
rz(-0.14856635) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.051604328) q[1];
sx q[1];
rz(-0.12433908) q[1];
sx q[1];
rz(-2.405637) q[1];
rz(2.8771411) q[3];
sx q[3];
rz(-1.7489986) q[3];
sx q[3];
rz(-1.7314535) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.000164) q[2];
sx q[2];
rz(-1.8464073) q[2];
sx q[2];
rz(-2.2953575) q[2];
rz(-2.6464461) q[3];
sx q[3];
rz(-0.64619243) q[3];
sx q[3];
rz(-1.600986) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.071035944) q[0];
sx q[0];
rz(-1.2905916) q[0];
sx q[0];
rz(-2.0284247) q[0];
rz(-2.44599) q[1];
sx q[1];
rz(-0.38989392) q[1];
sx q[1];
rz(-1.6092469) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5672011) q[0];
sx q[0];
rz(-2.6041457) q[0];
sx q[0];
rz(2.6525932) q[0];
rz(-2.6534383) q[2];
sx q[2];
rz(-1.9794954) q[2];
sx q[2];
rz(0.28446928) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.53837) q[1];
sx q[1];
rz(-0.45062989) q[1];
sx q[1];
rz(1.5249114) q[1];
x q[2];
rz(-0.60204695) q[3];
sx q[3];
rz(-1.6685408) q[3];
sx q[3];
rz(-2.5412113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.9939076) q[2];
sx q[2];
rz(-0.2299749) q[2];
sx q[2];
rz(-2.2873986) q[2];
rz(1.2285852) q[3];
sx q[3];
rz(-1.5766141) q[3];
sx q[3];
rz(1.4860229) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5478741) q[0];
sx q[0];
rz(-0.49867189) q[0];
sx q[0];
rz(2.4771931) q[0];
rz(1.9539072) q[1];
sx q[1];
rz(-0.70078754) q[1];
sx q[1];
rz(1.857035) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8878471) q[0];
sx q[0];
rz(-2.3652024) q[0];
sx q[0];
rz(1.1573769) q[0];
x q[1];
rz(-2.8753488) q[2];
sx q[2];
rz(-1.9387445) q[2];
sx q[2];
rz(-0.14012303) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.43209546) q[1];
sx q[1];
rz(-2.2120683) q[1];
sx q[1];
rz(-2.1367367) q[1];
rz(-pi) q[2];
rz(-3.0751238) q[3];
sx q[3];
rz(-2.7507938) q[3];
sx q[3];
rz(1.9399411) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.61218843) q[2];
sx q[2];
rz(-2.2470784) q[2];
sx q[2];
rz(0.55076304) q[2];
rz(2.4380056) q[3];
sx q[3];
rz(-2.1313322) q[3];
sx q[3];
rz(-1.5065058) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4187014) q[0];
sx q[0];
rz(-0.35690618) q[0];
sx q[0];
rz(-3.0974467) q[0];
rz(-1.528953) q[1];
sx q[1];
rz(-1.2395369) q[1];
sx q[1];
rz(2.3619161) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44256193) q[0];
sx q[0];
rz(-2.583722) q[0];
sx q[0];
rz(-2.872422) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0585971) q[2];
sx q[2];
rz(-2.3239845) q[2];
sx q[2];
rz(1.6542733) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.4317707) q[1];
sx q[1];
rz(-0.55574544) q[1];
sx q[1];
rz(-2.7906448) q[1];
rz(-1.1053106) q[3];
sx q[3];
rz(-2.2959384) q[3];
sx q[3];
rz(0.5474962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.49446517) q[2];
sx q[2];
rz(-2.4202042) q[2];
sx q[2];
rz(-1.5987827) q[2];
rz(0.70458448) q[3];
sx q[3];
rz(-1.6274118) q[3];
sx q[3];
rz(3.1295479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71173944) q[0];
sx q[0];
rz(-1.585351) q[0];
sx q[0];
rz(2.4899695) q[0];
rz(1.7977057) q[1];
sx q[1];
rz(-1.4812891) q[1];
sx q[1];
rz(-0.67566009) q[1];
rz(-2.8316108) q[2];
sx q[2];
rz(-1.1520755) q[2];
sx q[2];
rz(0.68785695) q[2];
rz(0.6380973) q[3];
sx q[3];
rz(-1.4641855) q[3];
sx q[3];
rz(-2.2104213) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
