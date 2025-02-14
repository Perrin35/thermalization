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
rz(1.4242564) q[0];
sx q[0];
rz(-1.2196701) q[0];
sx q[0];
rz(-0.57504672) q[0];
rz(-2.9991034) q[1];
sx q[1];
rz(-1.2187076) q[1];
sx q[1];
rz(0.86427468) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48102114) q[0];
sx q[0];
rz(-0.29566524) q[0];
sx q[0];
rz(-0.85938248) q[0];
rz(-1.3623385) q[2];
sx q[2];
rz(-2.499806) q[2];
sx q[2];
rz(-1.2520977) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.083347224) q[1];
sx q[1];
rz(-1.5687404) q[1];
sx q[1];
rz(2.6115692) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4349917) q[3];
sx q[3];
rz(-2.4445783) q[3];
sx q[3];
rz(1.5829338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.064726) q[2];
sx q[2];
rz(-1.4996424) q[2];
sx q[2];
rz(-2.1323252) q[2];
rz(2.3276954) q[3];
sx q[3];
rz(-2.8129357) q[3];
sx q[3];
rz(0.3391268) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.086804248) q[0];
sx q[0];
rz(-1.4084933) q[0];
sx q[0];
rz(-2.8993697) q[0];
rz(-1.9107266) q[1];
sx q[1];
rz(-0.63175285) q[1];
sx q[1];
rz(-2.3435074) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4674007) q[0];
sx q[0];
rz(-0.62819203) q[0];
sx q[0];
rz(-2.9487361) q[0];
x q[1];
rz(1.1801486) q[2];
sx q[2];
rz(-2.7728348) q[2];
sx q[2];
rz(-0.40290305) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.4646513) q[1];
sx q[1];
rz(-1.92917) q[1];
sx q[1];
rz(0.48848571) q[1];
rz(0.40753461) q[3];
sx q[3];
rz(-2.7544602) q[3];
sx q[3];
rz(-0.77104309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.65608394) q[2];
sx q[2];
rz(-1.069243) q[2];
sx q[2];
rz(0.83267027) q[2];
rz(0.63140702) q[3];
sx q[3];
rz(-0.74149817) q[3];
sx q[3];
rz(-3.1413063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1239419) q[0];
sx q[0];
rz(-0.90921062) q[0];
sx q[0];
rz(-0.18774524) q[0];
rz(-0.84367696) q[1];
sx q[1];
rz(-1.2742821) q[1];
sx q[1];
rz(-0.63293308) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7596066) q[0];
sx q[0];
rz(-1.4999903) q[0];
sx q[0];
rz(1.9340865) q[0];
rz(-0.68667163) q[2];
sx q[2];
rz(-1.8455659) q[2];
sx q[2];
rz(-2.6426154) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9687815) q[1];
sx q[1];
rz(-2.079614) q[1];
sx q[1];
rz(-1.225586) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2237596) q[3];
sx q[3];
rz(-1.6926024) q[3];
sx q[3];
rz(-0.95637874) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.8741499) q[2];
sx q[2];
rz(-0.94620693) q[2];
sx q[2];
rz(1.2713185) q[2];
rz(1.4960131) q[3];
sx q[3];
rz(-1.6451903) q[3];
sx q[3];
rz(2.1700844) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8895759) q[0];
sx q[0];
rz(-0.75279623) q[0];
sx q[0];
rz(2.4821607) q[0];
rz(-1.8239498) q[1];
sx q[1];
rz(-1.3392071) q[1];
sx q[1];
rz(0.79016322) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1135318) q[0];
sx q[0];
rz(-1.694593) q[0];
sx q[0];
rz(-0.33756279) q[0];
x q[1];
rz(-1.6035346) q[2];
sx q[2];
rz(-1.5509836) q[2];
sx q[2];
rz(-2.0231501) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.15406628) q[1];
sx q[1];
rz(-0.48168698) q[1];
sx q[1];
rz(-2.9996526) q[1];
rz(-pi) q[2];
rz(-0.42415027) q[3];
sx q[3];
rz(-0.59643918) q[3];
sx q[3];
rz(3.0354478) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.059375199) q[2];
sx q[2];
rz(-1.5267812) q[2];
sx q[2];
rz(-2.8455632) q[2];
rz(-2.5791903) q[3];
sx q[3];
rz(-2.2735169) q[3];
sx q[3];
rz(-3.0276022) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.247308) q[0];
sx q[0];
rz(-1.1881275) q[0];
sx q[0];
rz(2.9344015) q[0];
rz(-1.0317135) q[1];
sx q[1];
rz(-1.1528015) q[1];
sx q[1];
rz(-1.656104) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88865796) q[0];
sx q[0];
rz(-2.2307848) q[0];
sx q[0];
rz(2.071991) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1319415) q[2];
sx q[2];
rz(-1.3232987) q[2];
sx q[2];
rz(-2.0405318) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.9511764) q[1];
sx q[1];
rz(-0.37280478) q[1];
sx q[1];
rz(-0.35721161) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6909825) q[3];
sx q[3];
rz(-0.73420364) q[3];
sx q[3];
rz(2.7199573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.1032054) q[2];
sx q[2];
rz(-2.3296671) q[2];
sx q[2];
rz(2.0337598) q[2];
rz(0.92249089) q[3];
sx q[3];
rz(-1.4100217) q[3];
sx q[3];
rz(0.24423519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
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
rz(1.51782) q[0];
sx q[0];
rz(-0.51384846) q[0];
sx q[0];
rz(-2.9129831) q[0];
rz(-0.27433431) q[1];
sx q[1];
rz(-2.3286596) q[1];
sx q[1];
rz(2.6913604) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58670973) q[0];
sx q[0];
rz(-2.3392117) q[0];
sx q[0];
rz(-0.35361617) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5374372) q[2];
sx q[2];
rz(-0.22936996) q[2];
sx q[2];
rz(-2.0663886) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6321559) q[1];
sx q[1];
rz(-0.56427279) q[1];
sx q[1];
rz(2.6585101) q[1];
rz(-pi) q[2];
rz(-2.8725876) q[3];
sx q[3];
rz(-1.9719567) q[3];
sx q[3];
rz(1.8393733) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.7776514) q[2];
sx q[2];
rz(-0.55299091) q[2];
sx q[2];
rz(0.027776329) q[2];
rz(-0.76644301) q[3];
sx q[3];
rz(-1.981363) q[3];
sx q[3];
rz(2.4465223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9195093) q[0];
sx q[0];
rz(-1.6351901) q[0];
sx q[0];
rz(-0.62244225) q[0];
rz(0.70603236) q[1];
sx q[1];
rz(-2.249554) q[1];
sx q[1];
rz(-0.58696729) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3764899) q[0];
sx q[0];
rz(-2.4060632) q[0];
sx q[0];
rz(2.2549948) q[0];
rz(-pi) q[1];
rz(2.6640011) q[2];
sx q[2];
rz(-2.6048663) q[2];
sx q[2];
rz(-1.1004741) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.7365954) q[1];
sx q[1];
rz(-2.6301503) q[1];
sx q[1];
rz(0.14528017) q[1];
rz(-pi) q[2];
rz(3.1312084) q[3];
sx q[3];
rz(-2.2828013) q[3];
sx q[3];
rz(2.4807348) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.065757699) q[2];
sx q[2];
rz(-1.4328052) q[2];
sx q[2];
rz(0.59448057) q[2];
rz(2.8822656) q[3];
sx q[3];
rz(-1.2649053) q[3];
sx q[3];
rz(-1.2172786) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.038430564) q[0];
sx q[0];
rz(-2.2362464) q[0];
sx q[0];
rz(-0.578798) q[0];
rz(1.3102866) q[1];
sx q[1];
rz(-1.2625887) q[1];
sx q[1];
rz(0.8126747) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90234251) q[0];
sx q[0];
rz(-0.50241155) q[0];
sx q[0];
rz(-0.30186304) q[0];
x q[1];
rz(1.9695639) q[2];
sx q[2];
rz(-0.59894717) q[2];
sx q[2];
rz(2.2676165) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.002269) q[1];
sx q[1];
rz(-1.2096757) q[1];
sx q[1];
rz(-2.4206672) q[1];
rz(-pi) q[2];
rz(-2.1401494) q[3];
sx q[3];
rz(-2.0096931) q[3];
sx q[3];
rz(0.90745181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.7909214) q[2];
sx q[2];
rz(-1.2988043) q[2];
sx q[2];
rz(0.92362967) q[2];
rz(-2.1212497) q[3];
sx q[3];
rz(-0.89710051) q[3];
sx q[3];
rz(-3.0237107) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(1.0608805) q[0];
sx q[0];
rz(-1.7055644) q[0];
sx q[0];
rz(-1.4870148) q[0];
rz(3.0912073) q[1];
sx q[1];
rz(-0.81704187) q[1];
sx q[1];
rz(-0.64687669) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9084307) q[0];
sx q[0];
rz(-1.1457503) q[0];
sx q[0];
rz(0.90954291) q[0];
x q[1];
rz(-2.0503309) q[2];
sx q[2];
rz(-2.1796436) q[2];
sx q[2];
rz(3.0225168) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.17740956) q[1];
sx q[1];
rz(-0.46176592) q[1];
sx q[1];
rz(-2.6564471) q[1];
rz(-2.991363) q[3];
sx q[3];
rz(-1.7818799) q[3];
sx q[3];
rz(0.37841132) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6369624) q[2];
sx q[2];
rz(-1.5550104) q[2];
sx q[2];
rz(-0.75616765) q[2];
rz(1.470083) q[3];
sx q[3];
rz(-0.78592891) q[3];
sx q[3];
rz(-2.3362931) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.037755448) q[0];
sx q[0];
rz(-1.2498195) q[0];
sx q[0];
rz(-1.1081498) q[0];
rz(-2.8773819) q[1];
sx q[1];
rz(-1.6725531) q[1];
sx q[1];
rz(-0.60233751) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7717404) q[0];
sx q[0];
rz(-2.8338972) q[0];
sx q[0];
rz(-1.0956531) q[0];
rz(-2.4265044) q[2];
sx q[2];
rz(-2.9093766) q[2];
sx q[2];
rz(0.84399904) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.332676) q[1];
sx q[1];
rz(-0.68854587) q[1];
sx q[1];
rz(1.6469514) q[1];
x q[2];
rz(-1.8014158) q[3];
sx q[3];
rz(-0.93276927) q[3];
sx q[3];
rz(-2.7122081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.63310528) q[2];
sx q[2];
rz(-2.5999531) q[2];
sx q[2];
rz(-0.87149054) q[2];
rz(-1.2095215) q[3];
sx q[3];
rz(-2.8612374) q[3];
sx q[3];
rz(-2.5476294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7021983) q[0];
sx q[0];
rz(-1.3855423) q[0];
sx q[0];
rz(2.6227797) q[0];
rz(-2.5595472) q[1];
sx q[1];
rz(-2.0379635) q[1];
sx q[1];
rz(1.4720974) q[1];
rz(-1.1956566) q[2];
sx q[2];
rz(-2.2715501) q[2];
sx q[2];
rz(3.0987433) q[2];
rz(1.7409848) q[3];
sx q[3];
rz(-2.2522829) q[3];
sx q[3];
rz(-1.9808148) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
