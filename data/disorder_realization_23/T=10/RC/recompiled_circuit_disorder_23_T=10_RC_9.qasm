OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.8795348) q[0];
sx q[0];
rz(-1.4095925) q[0];
sx q[0];
rz(-1.4341266) q[0];
rz(0.6342451) q[1];
sx q[1];
rz(-2.5399962) q[1];
sx q[1];
rz(-0.4184202) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19791662) q[0];
sx q[0];
rz(-1.1063873) q[0];
sx q[0];
rz(-0.29505131) q[0];
rz(-pi) q[1];
rz(1.7069874) q[2];
sx q[2];
rz(-1.4601267) q[2];
sx q[2];
rz(1.9905123) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.13222209) q[1];
sx q[1];
rz(-1.8144061) q[1];
sx q[1];
rz(2.0396712) q[1];
rz(-pi) q[2];
rz(2.9207346) q[3];
sx q[3];
rz(-2.0005895) q[3];
sx q[3];
rz(-0.62648279) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.819954) q[2];
sx q[2];
rz(-1.3691838) q[2];
sx q[2];
rz(-2.3036172) q[2];
rz(2.6485802) q[3];
sx q[3];
rz(-2.8686782) q[3];
sx q[3];
rz(0.078991927) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3943966) q[0];
sx q[0];
rz(-2.4173739) q[0];
sx q[0];
rz(-1.863742) q[0];
rz(-2.9648119) q[1];
sx q[1];
rz(-1.8272094) q[1];
sx q[1];
rz(-2.7094254) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9958187) q[0];
sx q[0];
rz(-2.5403025) q[0];
sx q[0];
rz(-0.50253089) q[0];
rz(-pi) q[1];
rz(-1.1439267) q[2];
sx q[2];
rz(-1.2032713) q[2];
sx q[2];
rz(1.9232242) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.4914815) q[1];
sx q[1];
rz(-2.0108129) q[1];
sx q[1];
rz(-0.12932175) q[1];
rz(-0.024859419) q[3];
sx q[3];
rz(-1.0059352) q[3];
sx q[3];
rz(-0.75627518) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.54923487) q[2];
sx q[2];
rz(-1.2640415) q[2];
sx q[2];
rz(2.6548927) q[2];
rz(1.7633847) q[3];
sx q[3];
rz(-1.2599726) q[3];
sx q[3];
rz(0.53282213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.10087) q[0];
sx q[0];
rz(-0.7464872) q[0];
sx q[0];
rz(-0.41734636) q[0];
rz(1.4886645) q[1];
sx q[1];
rz(-2.5960943) q[1];
sx q[1];
rz(-2.6352077) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34844549) q[0];
sx q[0];
rz(-1.2198824) q[0];
sx q[0];
rz(-0.1883513) q[0];
rz(0.26489139) q[2];
sx q[2];
rz(-2.3893917) q[2];
sx q[2];
rz(-1.0571935) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.80458927) q[1];
sx q[1];
rz(-2.5924006) q[1];
sx q[1];
rz(-0.6188436) q[1];
rz(0.18249986) q[3];
sx q[3];
rz(-1.2068519) q[3];
sx q[3];
rz(-0.58992995) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.69904077) q[2];
sx q[2];
rz(-0.46135819) q[2];
sx q[2];
rz(-0.59147269) q[2];
rz(2.5555723) q[3];
sx q[3];
rz(-1.932671) q[3];
sx q[3];
rz(1.7104141) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1699003) q[0];
sx q[0];
rz(-2.5132892) q[0];
sx q[0];
rz(0.83918321) q[0];
rz(3.1160141) q[1];
sx q[1];
rz(-0.69568101) q[1];
sx q[1];
rz(1.5930088) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4913113) q[0];
sx q[0];
rz(-0.45931739) q[0];
sx q[0];
rz(-1.3643054) q[0];
rz(-pi) q[1];
rz(2.3392802) q[2];
sx q[2];
rz(-2.1244086) q[2];
sx q[2];
rz(2.4137036) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.81739391) q[1];
sx q[1];
rz(-0.84016582) q[1];
sx q[1];
rz(0.52449951) q[1];
rz(-2.7110093) q[3];
sx q[3];
rz(-1.3123371) q[3];
sx q[3];
rz(-2.5653198) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.30535355) q[2];
sx q[2];
rz(-0.88399115) q[2];
sx q[2];
rz(3.0419066) q[2];
rz(0.95885197) q[3];
sx q[3];
rz(-1.8226263) q[3];
sx q[3];
rz(1.3249741) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5754159) q[0];
sx q[0];
rz(-1.382099) q[0];
sx q[0];
rz(-0.25594041) q[0];
rz(-2.6804965) q[1];
sx q[1];
rz(-1.0436811) q[1];
sx q[1];
rz(-0.76006132) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5998659) q[0];
sx q[0];
rz(-1.6875629) q[0];
sx q[0];
rz(1.4670502) q[0];
rz(1.4270093) q[2];
sx q[2];
rz(-0.41848768) q[2];
sx q[2];
rz(0.33868044) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.0415503) q[1];
sx q[1];
rz(-1.6987545) q[1];
sx q[1];
rz(3.0138739) q[1];
rz(-pi) q[2];
x q[2];
rz(0.92442583) q[3];
sx q[3];
rz(-1.2741718) q[3];
sx q[3];
rz(0.38277205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.57050675) q[2];
sx q[2];
rz(-1.7692302) q[2];
sx q[2];
rz(2.467353) q[2];
rz(-0.21480602) q[3];
sx q[3];
rz(-0.45682296) q[3];
sx q[3];
rz(0.017344346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6102585) q[0];
sx q[0];
rz(-1.6711618) q[0];
sx q[0];
rz(-1.9624788) q[0];
rz(-0.20482652) q[1];
sx q[1];
rz(-2.3463459) q[1];
sx q[1];
rz(-2.0746322) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93766312) q[0];
sx q[0];
rz(-1.5563037) q[0];
sx q[0];
rz(-3.1209164) q[0];
x q[1];
rz(-0.42491575) q[2];
sx q[2];
rz(-1.1669461) q[2];
sx q[2];
rz(2.0770819) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.8300007) q[1];
sx q[1];
rz(-2.316906) q[1];
sx q[1];
rz(0.064706116) q[1];
rz(-pi) q[2];
rz(2.3451869) q[3];
sx q[3];
rz(-1.4561597) q[3];
sx q[3];
rz(-2.93626) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.71010464) q[2];
sx q[2];
rz(-1.8523214) q[2];
sx q[2];
rz(0.55994326) q[2];
rz(0.7263178) q[3];
sx q[3];
rz(-0.30877078) q[3];
sx q[3];
rz(-2.8360951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2840435) q[0];
sx q[0];
rz(-0.54656583) q[0];
sx q[0];
rz(-1.7204826) q[0];
rz(-0.20206085) q[1];
sx q[1];
rz(-1.7077363) q[1];
sx q[1];
rz(2.2834159) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9439529) q[0];
sx q[0];
rz(-1.4240992) q[0];
sx q[0];
rz(0.12978817) q[0];
rz(-pi) q[1];
x q[1];
rz(0.045072149) q[2];
sx q[2];
rz(-1.1606693) q[2];
sx q[2];
rz(2.8571667) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.2724175) q[1];
sx q[1];
rz(-0.032748001) q[1];
sx q[1];
rz(-0.48519022) q[1];
rz(1.2248366) q[3];
sx q[3];
rz(-2.1145027) q[3];
sx q[3];
rz(-2.791415) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.8537366) q[2];
sx q[2];
rz(-0.47984543) q[2];
sx q[2];
rz(1.8161592) q[2];
rz(0.89007968) q[3];
sx q[3];
rz(-1.1471014) q[3];
sx q[3];
rz(-2.1530698) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6778075) q[0];
sx q[0];
rz(-0.57254922) q[0];
sx q[0];
rz(0.37471399) q[0];
rz(0.97887865) q[1];
sx q[1];
rz(-0.68190494) q[1];
sx q[1];
rz(-1.3495061) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5193072) q[0];
sx q[0];
rz(-1.95682) q[0];
sx q[0];
rz(-2.4849154) q[0];
x q[1];
rz(-1.5937514) q[2];
sx q[2];
rz(-0.85283961) q[2];
sx q[2];
rz(0.1644451) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.3186036) q[1];
sx q[1];
rz(-2.9530596) q[1];
sx q[1];
rz(1.8382501) q[1];
x q[2];
rz(-2.7262444) q[3];
sx q[3];
rz(-0.58598622) q[3];
sx q[3];
rz(1.8893482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.65138856) q[2];
sx q[2];
rz(-2.3888402) q[2];
sx q[2];
rz(-2.6728969) q[2];
rz(1.9474585) q[3];
sx q[3];
rz(-1.8374551) q[3];
sx q[3];
rz(1.3635925) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63012183) q[0];
sx q[0];
rz(-0.90634316) q[0];
sx q[0];
rz(1.2619031) q[0];
rz(-2.966554) q[1];
sx q[1];
rz(-1.1418399) q[1];
sx q[1];
rz(-1.6040241) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3106829) q[0];
sx q[0];
rz(-1.3498107) q[0];
sx q[0];
rz(1.2587147) q[0];
rz(-pi) q[1];
rz(1.6749009) q[2];
sx q[2];
rz(-1.1541919) q[2];
sx q[2];
rz(0.57550752) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.21737145) q[1];
sx q[1];
rz(-1.6842168) q[1];
sx q[1];
rz(-2.4356615) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0314625) q[3];
sx q[3];
rz(-1.2445645) q[3];
sx q[3];
rz(2.5903451) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.039915446) q[2];
sx q[2];
rz(-2.1890409) q[2];
sx q[2];
rz(-2.3802479) q[2];
rz(0.90041655) q[3];
sx q[3];
rz(-2.5420928) q[3];
sx q[3];
rz(3.0925687) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3392357) q[0];
sx q[0];
rz(-2.673322) q[0];
sx q[0];
rz(2.9246869) q[0];
rz(0.63198173) q[1];
sx q[1];
rz(-1.4914373) q[1];
sx q[1];
rz(0.95473081) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.023708658) q[0];
sx q[0];
rz(-0.81239359) q[0];
sx q[0];
rz(0.448416) q[0];
x q[1];
rz(2.5108699) q[2];
sx q[2];
rz(-0.64881697) q[2];
sx q[2];
rz(-2.0231252) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6001544) q[1];
sx q[1];
rz(-0.46240515) q[1];
sx q[1];
rz(-2.9701783) q[1];
rz(0.93654376) q[3];
sx q[3];
rz(-1.9359971) q[3];
sx q[3];
rz(-1.4432424) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.1251936) q[2];
sx q[2];
rz(-1.9026326) q[2];
sx q[2];
rz(-2.1968502) q[2];
rz(0.38481209) q[3];
sx q[3];
rz(-2.0231569) q[3];
sx q[3];
rz(2.184536) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64086296) q[0];
sx q[0];
rz(-0.50518112) q[0];
sx q[0];
rz(1.5541979) q[0];
rz(0.8846994) q[1];
sx q[1];
rz(-0.90507602) q[1];
sx q[1];
rz(-0.25837635) q[1];
rz(-0.83512886) q[2];
sx q[2];
rz(-1.0799115) q[2];
sx q[2];
rz(-0.66839914) q[2];
rz(-1.1273884) q[3];
sx q[3];
rz(-1.3848806) q[3];
sx q[3];
rz(1.0564907) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
