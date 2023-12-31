OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.3611276) q[0];
sx q[0];
rz(-0.68929231) q[0];
sx q[0];
rz(2.8110992) q[0];
rz(-2.8117872) q[1];
sx q[1];
rz(-2.2916315) q[1];
sx q[1];
rz(2.4324774) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.015895122) q[0];
sx q[0];
rz(-2.1791434) q[0];
sx q[0];
rz(2.7817821) q[0];
rz(-pi) q[1];
rz(-1.8843295) q[2];
sx q[2];
rz(-1.6013813) q[2];
sx q[2];
rz(-2.7262296) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.658537) q[1];
sx q[1];
rz(-0.72209789) q[1];
sx q[1];
rz(-2.2621821) q[1];
rz(-pi) q[2];
rz(-2.764774) q[3];
sx q[3];
rz(-1.0816649) q[3];
sx q[3];
rz(2.7917142) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.4101397) q[2];
sx q[2];
rz(-2.7316766) q[2];
sx q[2];
rz(1.5343792) q[2];
rz(0.93506995) q[3];
sx q[3];
rz(-1.3053223) q[3];
sx q[3];
rz(-2.3527761) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4193029) q[0];
sx q[0];
rz(-3.0794444) q[0];
sx q[0];
rz(0.62227917) q[0];
rz(0.17624804) q[1];
sx q[1];
rz(-1.9269678) q[1];
sx q[1];
rz(-2.2252749) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30442552) q[0];
sx q[0];
rz(-2.2080748) q[0];
sx q[0];
rz(-0.62563719) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7158521) q[2];
sx q[2];
rz(-2.3999891) q[2];
sx q[2];
rz(-1.2765826) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.898145) q[1];
sx q[1];
rz(-1.2444278) q[1];
sx q[1];
rz(1.2136202) q[1];
rz(-pi) q[2];
x q[2];
rz(0.022577062) q[3];
sx q[3];
rz(-1.1728247) q[3];
sx q[3];
rz(0.28385362) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9006485) q[2];
sx q[2];
rz(-0.99664656) q[2];
sx q[2];
rz(2.3201578) q[2];
rz(3.1243096) q[3];
sx q[3];
rz(-2.343943) q[3];
sx q[3];
rz(-0.78330529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
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
rz(0.18773742) q[0];
sx q[0];
rz(-1.563235) q[0];
sx q[0];
rz(2.1333372) q[0];
rz(0.035765212) q[1];
sx q[1];
rz(-1.4859896) q[1];
sx q[1];
rz(-0.52454138) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5265822) q[0];
sx q[0];
rz(-0.16631642) q[0];
sx q[0];
rz(1.389857) q[0];
x q[1];
rz(2.4630765) q[2];
sx q[2];
rz(-1.3635474) q[2];
sx q[2];
rz(-3.0845272) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.14568612) q[1];
sx q[1];
rz(-1.6325163) q[1];
sx q[1];
rz(-1.593959) q[1];
x q[2];
rz(-1.1509622) q[3];
sx q[3];
rz(-1.435278) q[3];
sx q[3];
rz(-2.6672222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.77164578) q[2];
sx q[2];
rz(-1.5180072) q[2];
sx q[2];
rz(1.4952205) q[2];
rz(-1.8203991) q[3];
sx q[3];
rz(-1.4097872) q[3];
sx q[3];
rz(2.7041919) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8111073) q[0];
sx q[0];
rz(-2.2225668) q[0];
sx q[0];
rz(-0.088949732) q[0];
rz(0.51070172) q[1];
sx q[1];
rz(-0.82534868) q[1];
sx q[1];
rz(2.451992) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3084429) q[0];
sx q[0];
rz(-1.2763378) q[0];
sx q[0];
rz(-1.740728) q[0];
rz(-pi) q[1];
rz(1.7961411) q[2];
sx q[2];
rz(-1.0190522) q[2];
sx q[2];
rz(-1.6793041) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.62024414) q[1];
sx q[1];
rz(-0.84940956) q[1];
sx q[1];
rz(2.0116624) q[1];
rz(-pi) q[2];
rz(-1.425399) q[3];
sx q[3];
rz(-1.0905915) q[3];
sx q[3];
rz(-2.8499545) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.4758063) q[2];
sx q[2];
rz(-0.67725724) q[2];
sx q[2];
rz(2.518667) q[2];
rz(1.1359435) q[3];
sx q[3];
rz(-2.9562852) q[3];
sx q[3];
rz(0.46666551) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2899807) q[0];
sx q[0];
rz(-1.8988134) q[0];
sx q[0];
rz(0.583453) q[0];
rz(1.9955697) q[1];
sx q[1];
rz(-1.6505516) q[1];
sx q[1];
rz(-1.4978283) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72901112) q[0];
sx q[0];
rz(-0.81775613) q[0];
sx q[0];
rz(1.9304995) q[0];
x q[1];
rz(-0.79767144) q[2];
sx q[2];
rz(-1.6399709) q[2];
sx q[2];
rz(-2.5728512) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.0651107) q[1];
sx q[1];
rz(-1.5434192) q[1];
sx q[1];
rz(-2.4379424) q[1];
rz(-pi) q[2];
rz(-1.4086401) q[3];
sx q[3];
rz(-2.080653) q[3];
sx q[3];
rz(1.4354524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.65537611) q[2];
sx q[2];
rz(-1.5950173) q[2];
sx q[2];
rz(-1.5636469) q[2];
rz(-0.90562138) q[3];
sx q[3];
rz(-0.27035299) q[3];
sx q[3];
rz(-0.63846987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49801302) q[0];
sx q[0];
rz(-1.4587412) q[0];
sx q[0];
rz(-0.046982732) q[0];
rz(-2.9934096) q[1];
sx q[1];
rz(-2.2427185) q[1];
sx q[1];
rz(-1.4354338) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88046056) q[0];
sx q[0];
rz(-1.4229703) q[0];
sx q[0];
rz(2.3418505) q[0];
rz(-pi) q[1];
rz(-3.0827423) q[2];
sx q[2];
rz(-1.0219136) q[2];
sx q[2];
rz(2.5318052) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6378577) q[1];
sx q[1];
rz(-1.4459472) q[1];
sx q[1];
rz(-1.1605074) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.1141838) q[3];
sx q[3];
rz(-1.8084744) q[3];
sx q[3];
rz(-2.4159367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.0662213) q[2];
sx q[2];
rz(-2.1263289) q[2];
sx q[2];
rz(3.0701239) q[2];
rz(1.4525157) q[3];
sx q[3];
rz(-1.6966597) q[3];
sx q[3];
rz(0.92648363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-2.8554095) q[0];
sx q[0];
rz(-1.3278642) q[0];
sx q[0];
rz(2.563971) q[0];
rz(1.2795992) q[1];
sx q[1];
rz(-1.7956693) q[1];
sx q[1];
rz(1.0095899) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5492064) q[0];
sx q[0];
rz(-1.9027332) q[0];
sx q[0];
rz(0.50360002) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3426443) q[2];
sx q[2];
rz(-1.9462799) q[2];
sx q[2];
rz(2.0827039) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.91612591) q[1];
sx q[1];
rz(-1.9088609) q[1];
sx q[1];
rz(2.613693) q[1];
rz(-pi) q[2];
x q[2];
rz(0.61543492) q[3];
sx q[3];
rz(-0.86343599) q[3];
sx q[3];
rz(-2.5497041) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.0916831) q[2];
sx q[2];
rz(-2.78806) q[2];
sx q[2];
rz(1.1996777) q[2];
rz(2.6692634) q[3];
sx q[3];
rz(-1.4940558) q[3];
sx q[3];
rz(-1.002731) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49884477) q[0];
sx q[0];
rz(-1.6413178) q[0];
sx q[0];
rz(0.2391267) q[0];
rz(-0.7827951) q[1];
sx q[1];
rz(-2.6233964) q[1];
sx q[1];
rz(-0.4447287) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77138222) q[0];
sx q[0];
rz(-2.133773) q[0];
sx q[0];
rz(-0.51008205) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3808001) q[2];
sx q[2];
rz(-2.3173601) q[2];
sx q[2];
rz(0.21812083) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0900314) q[1];
sx q[1];
rz(-1.4434442) q[1];
sx q[1];
rz(0.8831555) q[1];
rz(-pi) q[2];
rz(2.8430311) q[3];
sx q[3];
rz(-2.7136554) q[3];
sx q[3];
rz(0.096979389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.42177054) q[2];
sx q[2];
rz(-2.3193216) q[2];
sx q[2];
rz(-2.3262809) q[2];
rz(2.1067965) q[3];
sx q[3];
rz(-1.5650322) q[3];
sx q[3];
rz(-1.8243272) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3141044) q[0];
sx q[0];
rz(-2.1056471) q[0];
sx q[0];
rz(-0.28717336) q[0];
rz(2.9526967) q[1];
sx q[1];
rz(-2.4177528) q[1];
sx q[1];
rz(-2.8093991) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9438292) q[0];
sx q[0];
rz(-1.8837067) q[0];
sx q[0];
rz(-2.3373332) q[0];
rz(-pi) q[1];
rz(0.046499649) q[2];
sx q[2];
rz(-2.2691155) q[2];
sx q[2];
rz(-0.38682129) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6576782) q[1];
sx q[1];
rz(-2.6473443) q[1];
sx q[1];
rz(0.56234042) q[1];
rz(-pi) q[2];
rz(-1.9751777) q[3];
sx q[3];
rz(-1.0746733) q[3];
sx q[3];
rz(2.8465084) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2087848) q[2];
sx q[2];
rz(-1.00495) q[2];
sx q[2];
rz(2.3020111) q[2];
rz(1.2906637) q[3];
sx q[3];
rz(-1.3495812) q[3];
sx q[3];
rz(2.4850142) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.055450913) q[0];
sx q[0];
rz(-2.3738326) q[0];
sx q[0];
rz(1.0593876) q[0];
rz(2.1620031) q[1];
sx q[1];
rz(-1.1971808) q[1];
sx q[1];
rz(0.25451452) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8043038) q[0];
sx q[0];
rz(-0.24233195) q[0];
sx q[0];
rz(-1.2605577) q[0];
rz(-0.43912402) q[2];
sx q[2];
rz(-0.52999485) q[2];
sx q[2];
rz(-2.8119171) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0490129) q[1];
sx q[1];
rz(-1.6515459) q[1];
sx q[1];
rz(-2.8156274) q[1];
rz(-pi) q[2];
rz(2.9322846) q[3];
sx q[3];
rz(-0.44692398) q[3];
sx q[3];
rz(-2.350654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.0467726) q[2];
sx q[2];
rz(-0.54868713) q[2];
sx q[2];
rz(1.0478896) q[2];
rz(-1.3607599) q[3];
sx q[3];
rz(-0.69793099) q[3];
sx q[3];
rz(2.5361983) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1142674) q[0];
sx q[0];
rz(-2.0226759) q[0];
sx q[0];
rz(-0.080060536) q[0];
rz(-0.36021532) q[1];
sx q[1];
rz(-1.4604912) q[1];
sx q[1];
rz(2.1866658) q[1];
rz(1.04503) q[2];
sx q[2];
rz(-1.7797995) q[2];
sx q[2];
rz(-1.9893653) q[2];
rz(2.1847235) q[3];
sx q[3];
rz(-1.7603684) q[3];
sx q[3];
rz(-2.6852222) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
