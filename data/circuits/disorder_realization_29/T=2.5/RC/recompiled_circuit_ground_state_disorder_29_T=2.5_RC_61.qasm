OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.86201) q[0];
sx q[0];
rz(-0.53505889) q[0];
sx q[0];
rz(0.23393272) q[0];
rz(0.84199953) q[1];
sx q[1];
rz(-1.7276126) q[1];
sx q[1];
rz(1.6230621) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46759638) q[0];
sx q[0];
rz(-2.5270871) q[0];
sx q[0];
rz(2.5835681) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1351003) q[2];
sx q[2];
rz(-0.47347906) q[2];
sx q[2];
rz(1.4896637) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.78511158) q[1];
sx q[1];
rz(-0.82725305) q[1];
sx q[1];
rz(0.77036304) q[1];
rz(-pi) q[2];
rz(-0.47187658) q[3];
sx q[3];
rz(-2.1778244) q[3];
sx q[3];
rz(-0.31479731) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.6003549) q[2];
sx q[2];
rz(-1.7011832) q[2];
sx q[2];
rz(2.3360628) q[2];
rz(2.7496998) q[3];
sx q[3];
rz(-1.3561748) q[3];
sx q[3];
rz(0.75683561) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5517752) q[0];
sx q[0];
rz(-0.68142319) q[0];
sx q[0];
rz(-2.7031194) q[0];
rz(-1.27502) q[1];
sx q[1];
rz(-1.6908815) q[1];
sx q[1];
rz(-0.10800392) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7399044) q[0];
sx q[0];
rz(-1.5956912) q[0];
sx q[0];
rz(1.4598926) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.503452) q[2];
sx q[2];
rz(-0.9305939) q[2];
sx q[2];
rz(1.3253044) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.7234791) q[1];
sx q[1];
rz(-1.8239841) q[1];
sx q[1];
rz(3.0550356) q[1];
rz(-1.0803079) q[3];
sx q[3];
rz(-1.401859) q[3];
sx q[3];
rz(-0.47034697) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.3652304) q[2];
sx q[2];
rz(-2.6680816) q[2];
sx q[2];
rz(-3.0253809) q[2];
rz(-2.727437) q[3];
sx q[3];
rz(-1.073758) q[3];
sx q[3];
rz(-2.3381086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47495833) q[0];
sx q[0];
rz(-1.3131498) q[0];
sx q[0];
rz(-0.83589244) q[0];
rz(-1.6235141) q[1];
sx q[1];
rz(-1.6267136) q[1];
sx q[1];
rz(-1.9035043) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9758751) q[0];
sx q[0];
rz(-1.200186) q[0];
sx q[0];
rz(-1.0281282) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.324685) q[2];
sx q[2];
rz(-1.2370584) q[2];
sx q[2];
rz(2.1527803) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.1819396) q[1];
sx q[1];
rz(-1.5132071) q[1];
sx q[1];
rz(-1.4532386) q[1];
rz(2.3457621) q[3];
sx q[3];
rz(-0.39271564) q[3];
sx q[3];
rz(-0.3317197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.8572924) q[2];
sx q[2];
rz(-2.9967873) q[2];
sx q[2];
rz(-2.4227552) q[2];
rz(2.5607064) q[3];
sx q[3];
rz(-1.418117) q[3];
sx q[3];
rz(-0.7287997) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-0.523664) q[0];
sx q[0];
rz(-1.5476462) q[0];
sx q[0];
rz(-1.1267598) q[0];
rz(-1.2976546) q[1];
sx q[1];
rz(-1.6512066) q[1];
sx q[1];
rz(1.0430956) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8617423) q[0];
sx q[0];
rz(-1.6011097) q[0];
sx q[0];
rz(1.8413196) q[0];
x q[1];
rz(1.9773433) q[2];
sx q[2];
rz(-2.5669328) q[2];
sx q[2];
rz(0.80917796) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.8759497) q[1];
sx q[1];
rz(-1.7247206) q[1];
sx q[1];
rz(1.2624082) q[1];
rz(-pi) q[2];
rz(1.6564441) q[3];
sx q[3];
rz(-2.2243735) q[3];
sx q[3];
rz(-1.9854922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.2780693) q[2];
sx q[2];
rz(-1.7570644) q[2];
sx q[2];
rz(2.1045904) q[2];
rz(2.1841124) q[3];
sx q[3];
rz(-1.0629531) q[3];
sx q[3];
rz(-2.3069042) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(1.0015513) q[0];
sx q[0];
rz(-0.20714864) q[0];
sx q[0];
rz(3.0539883) q[0];
rz(2.7032848) q[1];
sx q[1];
rz(-2.0594845) q[1];
sx q[1];
rz(1.8278488) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7676108) q[0];
sx q[0];
rz(-2.0519407) q[0];
sx q[0];
rz(2.8767881) q[0];
rz(1.6694267) q[2];
sx q[2];
rz(-1.6079788) q[2];
sx q[2];
rz(-1.1929026) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.524595) q[1];
sx q[1];
rz(-1.4819078) q[1];
sx q[1];
rz(-1.7434381) q[1];
rz(1.8270709) q[3];
sx q[3];
rz(-1.4087143) q[3];
sx q[3];
rz(-2.7173998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.4801415) q[2];
sx q[2];
rz(-1.9037312) q[2];
sx q[2];
rz(0.56337774) q[2];
rz(-1.2207458) q[3];
sx q[3];
rz(-1.7505587) q[3];
sx q[3];
rz(-2.5105072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0193943) q[0];
sx q[0];
rz(-0.65827426) q[0];
sx q[0];
rz(-0.011938183) q[0];
rz(-2.7519233) q[1];
sx q[1];
rz(-0.7020815) q[1];
sx q[1];
rz(0.24519244) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3936886) q[0];
sx q[0];
rz(-0.98391279) q[0];
sx q[0];
rz(-2.6176207) q[0];
x q[1];
rz(0.25845627) q[2];
sx q[2];
rz(-1.788013) q[2];
sx q[2];
rz(2.1803792) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.17055146) q[1];
sx q[1];
rz(-1.028113) q[1];
sx q[1];
rz(-1.0566864) q[1];
rz(-pi) q[2];
rz(-2.7390476) q[3];
sx q[3];
rz(-1.3005101) q[3];
sx q[3];
rz(-0.2337993) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.7616854) q[2];
sx q[2];
rz(-1.3096389) q[2];
sx q[2];
rz(-1.3690108) q[2];
rz(2.6109429) q[3];
sx q[3];
rz(-0.68415087) q[3];
sx q[3];
rz(1.0859038) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1928007) q[0];
sx q[0];
rz(-1.3230319) q[0];
sx q[0];
rz(0.2970933) q[0];
rz(1.6311215) q[1];
sx q[1];
rz(-0.62989569) q[1];
sx q[1];
rz(0.43103257) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98497903) q[0];
sx q[0];
rz(-2.7811858) q[0];
sx q[0];
rz(0.63215881) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9210429) q[2];
sx q[2];
rz(-0.23633453) q[2];
sx q[2];
rz(2.5823809) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3869074) q[1];
sx q[1];
rz(-1.0656989) q[1];
sx q[1];
rz(-2.7688857) q[1];
rz(-2.0484974) q[3];
sx q[3];
rz(-0.58148958) q[3];
sx q[3];
rz(-3.1149816) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.3057574) q[2];
sx q[2];
rz(-2.367986) q[2];
sx q[2];
rz(2.8744899) q[2];
rz(-3.109572) q[3];
sx q[3];
rz(-1.1705541) q[3];
sx q[3];
rz(3.035868) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0778462) q[0];
sx q[0];
rz(-1.8505322) q[0];
sx q[0];
rz(-0.63999501) q[0];
rz(0.85482875) q[1];
sx q[1];
rz(-1.6565485) q[1];
sx q[1];
rz(1.4124426) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0542568) q[0];
sx q[0];
rz(-2.1001513) q[0];
sx q[0];
rz(2.0482778) q[0];
rz(-1.475263) q[2];
sx q[2];
rz(-2.6474806) q[2];
sx q[2];
rz(1.6803368) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.0113653) q[1];
sx q[1];
rz(-1.9107358) q[1];
sx q[1];
rz(1.4814427) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9962003) q[3];
sx q[3];
rz(-1.914905) q[3];
sx q[3];
rz(-0.95208012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.5414446) q[2];
sx q[2];
rz(-1.4358127) q[2];
sx q[2];
rz(2.7195462) q[2];
rz(-2.6680434) q[3];
sx q[3];
rz(-0.62255064) q[3];
sx q[3];
rz(-2.9464909) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7710829) q[0];
sx q[0];
rz(-2.2042553) q[0];
sx q[0];
rz(1.7899845) q[0];
rz(-0.31556684) q[1];
sx q[1];
rz(-1.4548929) q[1];
sx q[1];
rz(1.1531166) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2446072) q[0];
sx q[0];
rz(-0.073052064) q[0];
sx q[0];
rz(0.79985072) q[0];
rz(-1.5222237) q[2];
sx q[2];
rz(-2.1056386) q[2];
sx q[2];
rz(-0.64316197) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.8421346) q[1];
sx q[1];
rz(-2.149035) q[1];
sx q[1];
rz(1.4186354) q[1];
rz(-pi) q[2];
rz(1.9824636) q[3];
sx q[3];
rz(-1.921706) q[3];
sx q[3];
rz(-1.0699492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.45712581) q[2];
sx q[2];
rz(-1.9292597) q[2];
sx q[2];
rz(2.7678164) q[2];
rz(-1.9155546) q[3];
sx q[3];
rz(-2.2235179) q[3];
sx q[3];
rz(-1.8019684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95947295) q[0];
sx q[0];
rz(-2.4612893) q[0];
sx q[0];
rz(1.8208338) q[0];
rz(-2.2044115) q[1];
sx q[1];
rz(-1.3359741) q[1];
sx q[1];
rz(2.6722867) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3166817) q[0];
sx q[0];
rz(-1.4535722) q[0];
sx q[0];
rz(3.0646695) q[0];
x q[1];
rz(0.98913828) q[2];
sx q[2];
rz(-1.7301705) q[2];
sx q[2];
rz(-1.6362658) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.7368245) q[1];
sx q[1];
rz(-1.4206845) q[1];
sx q[1];
rz(-1.4895358) q[1];
rz(-pi) q[2];
rz(2.8402249) q[3];
sx q[3];
rz(-1.3346938) q[3];
sx q[3];
rz(0.88168854) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4361973) q[2];
sx q[2];
rz(-0.86279482) q[2];
sx q[2];
rz(-1.9311284) q[2];
rz(0.56810275) q[3];
sx q[3];
rz(-1.3848687) q[3];
sx q[3];
rz(-0.97584045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9217011) q[0];
sx q[0];
rz(-2.2911063) q[0];
sx q[0];
rz(-1.6794857) q[0];
rz(0.89314356) q[1];
sx q[1];
rz(-1.8291263) q[1];
sx q[1];
rz(-1.6222454) q[1];
rz(2.91165) q[2];
sx q[2];
rz(-1.2886921) q[2];
sx q[2];
rz(-3.1293426) q[2];
rz(-1.0004956) q[3];
sx q[3];
rz(-1.199493) q[3];
sx q[3];
rz(-1.5293157) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
