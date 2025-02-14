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
rz(-0.57617968) q[0];
sx q[0];
rz(-1.4644858) q[0];
sx q[0];
rz(-0.65281868) q[0];
rz(-0.44261143) q[1];
sx q[1];
rz(-2.0191329) q[1];
sx q[1];
rz(5.6607487) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63875893) q[0];
sx q[0];
rz(-1.2593102) q[0];
sx q[0];
rz(0.12249662) q[0];
x q[1];
rz(2.6735252) q[2];
sx q[2];
rz(-0.65535802) q[2];
sx q[2];
rz(1.5214024) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.5099401) q[1];
sx q[1];
rz(-1.3871683) q[1];
sx q[1];
rz(2.0765523) q[1];
rz(-pi) q[2];
rz(1.7951598) q[3];
sx q[3];
rz(-2.5208559) q[3];
sx q[3];
rz(0.94514314) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.5114674) q[2];
sx q[2];
rz(-1.7897391) q[2];
sx q[2];
rz(2.5959065) q[2];
rz(0.68583471) q[3];
sx q[3];
rz(-0.67982173) q[3];
sx q[3];
rz(2.1849476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0551374) q[0];
sx q[0];
rz(-2.4392023) q[0];
sx q[0];
rz(2.0461244) q[0];
rz(0.99758863) q[1];
sx q[1];
rz(-1.2350524) q[1];
sx q[1];
rz(1.9445317) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72916146) q[0];
sx q[0];
rz(-2.1502156) q[0];
sx q[0];
rz(0.30381911) q[0];
rz(-pi) q[1];
rz(0.62792553) q[2];
sx q[2];
rz(-2.2596779) q[2];
sx q[2];
rz(2.9297926) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.61651308) q[1];
sx q[1];
rz(-2.5933806) q[1];
sx q[1];
rz(1.0031149) q[1];
x q[2];
rz(-2.8124146) q[3];
sx q[3];
rz(-1.7478895) q[3];
sx q[3];
rz(-0.30325034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.41584388) q[2];
sx q[2];
rz(-1.3943322) q[2];
sx q[2];
rz(1.7572629) q[2];
rz(1.3437126) q[3];
sx q[3];
rz(-1.5680771) q[3];
sx q[3];
rz(1.869092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62190732) q[0];
sx q[0];
rz(-0.78751957) q[0];
sx q[0];
rz(0.4271048) q[0];
rz(-1.9097795) q[1];
sx q[1];
rz(-1.6424664) q[1];
sx q[1];
rz(2.5274091) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6525567) q[0];
sx q[0];
rz(-1.8985735) q[0];
sx q[0];
rz(-1.2851738) q[0];
rz(-pi) q[1];
rz(2.9236058) q[2];
sx q[2];
rz(-2.3828016) q[2];
sx q[2];
rz(1.9050777) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.79518049) q[1];
sx q[1];
rz(-1.507057) q[1];
sx q[1];
rz(-1.6883114) q[1];
rz(1.3282458) q[3];
sx q[3];
rz(-1.0772675) q[3];
sx q[3];
rz(-2.7427615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.68982879) q[2];
sx q[2];
rz(-0.58388766) q[2];
sx q[2];
rz(0.27099398) q[2];
rz(-0.48218918) q[3];
sx q[3];
rz(-0.8420344) q[3];
sx q[3];
rz(1.2838001) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86376205) q[0];
sx q[0];
rz(-1.8067124) q[0];
sx q[0];
rz(-0.41896391) q[0];
rz(0.22467443) q[1];
sx q[1];
rz(-1.2089665) q[1];
sx q[1];
rz(-0.077542543) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2231295) q[0];
sx q[0];
rz(-1.0498199) q[0];
sx q[0];
rz(-0.33921241) q[0];
rz(-pi) q[1];
rz(-2.3691142) q[2];
sx q[2];
rz(-1.315552) q[2];
sx q[2];
rz(-0.30568631) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.3540707) q[1];
sx q[1];
rz(-1.4284627) q[1];
sx q[1];
rz(0.44146367) q[1];
rz(-pi) q[2];
x q[2];
rz(0.61309149) q[3];
sx q[3];
rz(-0.35260751) q[3];
sx q[3];
rz(-0.72474721) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-3.0858687) q[2];
sx q[2];
rz(-0.33129498) q[2];
sx q[2];
rz(-1.9191939) q[2];
rz(-2.870765) q[3];
sx q[3];
rz(-1.9041678) q[3];
sx q[3];
rz(0.45473155) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
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
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5643352) q[0];
sx q[0];
rz(-2.1463647) q[0];
sx q[0];
rz(-1.7876392) q[0];
rz(2.1063781) q[1];
sx q[1];
rz(-1.2151006) q[1];
sx q[1];
rz(-1.530102) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.852762) q[0];
sx q[0];
rz(-1.0438598) q[0];
sx q[0];
rz(-2.2599392) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8132854) q[2];
sx q[2];
rz(-1.7264778) q[2];
sx q[2];
rz(-0.98634431) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.6616231) q[1];
sx q[1];
rz(-1.9723867) q[1];
sx q[1];
rz(1.5127403) q[1];
rz(-pi) q[2];
rz(2.7306385) q[3];
sx q[3];
rz(-1.5180713) q[3];
sx q[3];
rz(-0.26782521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.80663854) q[2];
sx q[2];
rz(-0.97794509) q[2];
sx q[2];
rz(0.063892603) q[2];
rz(2.4334) q[3];
sx q[3];
rz(-1.6876561) q[3];
sx q[3];
rz(1.3341058) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(2.719249) q[0];
sx q[0];
rz(-0.023119211) q[0];
sx q[0];
rz(-2.0656021) q[0];
rz(3.0883582) q[1];
sx q[1];
rz(-2.2884171) q[1];
sx q[1];
rz(0.3784953) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9747978) q[0];
sx q[0];
rz(-2.85436) q[0];
sx q[0];
rz(-0.47606456) q[0];
rz(1.4827864) q[2];
sx q[2];
rz(-1.529379) q[2];
sx q[2];
rz(2.8633022) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.4617745) q[1];
sx q[1];
rz(-1.312612) q[1];
sx q[1];
rz(-1.3854909) q[1];
rz(-pi) q[2];
rz(-0.58121292) q[3];
sx q[3];
rz(-0.85678116) q[3];
sx q[3];
rz(2.9670144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.9516912) q[2];
sx q[2];
rz(-1.4611763) q[2];
sx q[2];
rz(1.0339197) q[2];
rz(-2.2377491) q[3];
sx q[3];
rz(-2.240182) q[3];
sx q[3];
rz(3.097539) q[3];
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
rz(pi/2) q[0];
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
rz(-1.1164301) q[0];
sx q[0];
rz(-1.4465541) q[0];
sx q[0];
rz(0.59301162) q[0];
rz(3.016839) q[1];
sx q[1];
rz(-1.9826823) q[1];
sx q[1];
rz(0.4932901) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.54087) q[0];
sx q[0];
rz(-0.42779765) q[0];
sx q[0];
rz(-0.85036253) q[0];
rz(-3.0346625) q[2];
sx q[2];
rz(-2.9122346) q[2];
sx q[2];
rz(1.3170674) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.3354937) q[1];
sx q[1];
rz(-2.4579398) q[1];
sx q[1];
rz(-0.70913507) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7760963) q[3];
sx q[3];
rz(-2.8474244) q[3];
sx q[3];
rz(0.95137596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.7571681) q[2];
sx q[2];
rz(-1.1606471) q[2];
sx q[2];
rz(1.1473131) q[2];
rz(-0.82289639) q[3];
sx q[3];
rz(-1.9673248) q[3];
sx q[3];
rz(2.4790922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
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
rz(0.61854521) q[0];
sx q[0];
rz(-1.86919) q[0];
sx q[0];
rz(-2.7230895) q[0];
rz(3.0283527) q[1];
sx q[1];
rz(-1.6721882) q[1];
sx q[1];
rz(2.0638827) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54544577) q[0];
sx q[0];
rz(-1.8592863) q[0];
sx q[0];
rz(3.122807) q[0];
rz(0.32025614) q[2];
sx q[2];
rz(-0.82112193) q[2];
sx q[2];
rz(-0.94515726) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7188765) q[1];
sx q[1];
rz(-0.36564974) q[1];
sx q[1];
rz(1.7083733) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6304794) q[3];
sx q[3];
rz(-0.91509295) q[3];
sx q[3];
rz(-2.7969517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.8354127) q[2];
sx q[2];
rz(-0.70859185) q[2];
sx q[2];
rz(1.6861247) q[2];
rz(-2.9742187) q[3];
sx q[3];
rz(-2.6264329) q[3];
sx q[3];
rz(-1.0719871) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6587332) q[0];
sx q[0];
rz(-1.3128244) q[0];
sx q[0];
rz(2.7400548) q[0];
rz(0.43831476) q[1];
sx q[1];
rz(-2.3500748) q[1];
sx q[1];
rz(-1.6848791) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12109893) q[0];
sx q[0];
rz(-1.5803845) q[0];
sx q[0];
rz(1.6047948) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4859547) q[2];
sx q[2];
rz(-1.3961424) q[2];
sx q[2];
rz(-1.8675592) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.5324133) q[1];
sx q[1];
rz(-1.0858156) q[1];
sx q[1];
rz(-2.3663125) q[1];
rz(-pi) q[2];
rz(0.73897408) q[3];
sx q[3];
rz(-0.65118507) q[3];
sx q[3];
rz(-0.84924752) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.7598286) q[2];
sx q[2];
rz(-0.20918736) q[2];
sx q[2];
rz(-2.4299202) q[2];
rz(-0.034218637) q[3];
sx q[3];
rz(-2.4331369) q[3];
sx q[3];
rz(-1.9862407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
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
rz(-1.3592247) q[0];
sx q[0];
rz(-1.490626) q[0];
sx q[0];
rz(-0.79886287) q[0];
rz(2.0955739) q[1];
sx q[1];
rz(-1.1963528) q[1];
sx q[1];
rz(-0.23652133) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2940148) q[0];
sx q[0];
rz(-1.6315407) q[0];
sx q[0];
rz(3.1192256) q[0];
rz(-pi) q[1];
rz(-0.6467713) q[2];
sx q[2];
rz(-2.1097217) q[2];
sx q[2];
rz(-2.981271) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.0276841) q[1];
sx q[1];
rz(-2.5792173) q[1];
sx q[1];
rz(-0.67090583) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0048851) q[3];
sx q[3];
rz(-0.71562222) q[3];
sx q[3];
rz(2.2608044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.025042621) q[2];
sx q[2];
rz(-2.5849085) q[2];
sx q[2];
rz(-0.64209783) q[2];
rz(2.5721278) q[3];
sx q[3];
rz(-2.6537708) q[3];
sx q[3];
rz(-0.08629442) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4702598) q[0];
sx q[0];
rz(-1.3339806) q[0];
sx q[0];
rz(0.95680923) q[0];
rz(1.9500465) q[1];
sx q[1];
rz(-1.5622495) q[1];
sx q[1];
rz(-1.539485) q[1];
rz(-3.0189162) q[2];
sx q[2];
rz(-0.79570607) q[2];
sx q[2];
rz(-0.88655587) q[2];
rz(-2.8068918) q[3];
sx q[3];
rz(-1.1664433) q[3];
sx q[3];
rz(0.56952624) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
