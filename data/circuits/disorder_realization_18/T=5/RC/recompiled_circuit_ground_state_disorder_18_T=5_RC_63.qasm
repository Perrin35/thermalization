OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.9590149) q[0];
sx q[0];
rz(-2.9380517) q[0];
sx q[0];
rz(-1.8939053) q[0];
rz(-1.7893451) q[1];
sx q[1];
rz(-2.4440553) q[1];
sx q[1];
rz(-2.4434659) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20839918) q[0];
sx q[0];
rz(-2.0740447) q[0];
sx q[0];
rz(-2.0397951) q[0];
rz(1.5641937) q[2];
sx q[2];
rz(-1.1175795) q[2];
sx q[2];
rz(-2.1389769) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.2064849) q[1];
sx q[1];
rz(-0.84151959) q[1];
sx q[1];
rz(2.7323099) q[1];
x q[2];
rz(3.1269253) q[3];
sx q[3];
rz(-1.8258182) q[3];
sx q[3];
rz(1.6153112) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.84833604) q[2];
sx q[2];
rz(-1.9170599) q[2];
sx q[2];
rz(-2.2621034) q[2];
rz(-2.4073811) q[3];
sx q[3];
rz(-2.0101533) q[3];
sx q[3];
rz(2.3122299) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5010928) q[0];
sx q[0];
rz(-1.9556029) q[0];
sx q[0];
rz(-2.1511141) q[0];
rz(-0.99539202) q[1];
sx q[1];
rz(-1.6973015) q[1];
sx q[1];
rz(2.676414) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2110307) q[0];
sx q[0];
rz(-2.6829236) q[0];
sx q[0];
rz(-0.42210292) q[0];
x q[1];
rz(0.37613531) q[2];
sx q[2];
rz(-2.1636164) q[2];
sx q[2];
rz(1.7629634) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.2485741) q[1];
sx q[1];
rz(-2.0573528) q[1];
sx q[1];
rz(2.3793329) q[1];
rz(-2.8797014) q[3];
sx q[3];
rz(-1.3757816) q[3];
sx q[3];
rz(-0.075500114) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.5140932) q[2];
sx q[2];
rz(-1.9072615) q[2];
sx q[2];
rz(-2.9631183) q[2];
rz(-0.56525362) q[3];
sx q[3];
rz(-1.7491128) q[3];
sx q[3];
rz(2.5843411) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31767118) q[0];
sx q[0];
rz(-1.3452106) q[0];
sx q[0];
rz(0.70096651) q[0];
rz(2.7340381) q[1];
sx q[1];
rz(-1.9302255) q[1];
sx q[1];
rz(2.0661381) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4385745) q[0];
sx q[0];
rz(-2.10963) q[0];
sx q[0];
rz(-1.1807022) q[0];
rz(-pi) q[1];
rz(-0.22547743) q[2];
sx q[2];
rz(-2.2144711) q[2];
sx q[2];
rz(-0.89356092) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6523931) q[1];
sx q[1];
rz(-1.4402188) q[1];
sx q[1];
rz(1.9741535) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8355708) q[3];
sx q[3];
rz(-1.522101) q[3];
sx q[3];
rz(2.8935561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.0740697) q[2];
sx q[2];
rz(-1.8068376) q[2];
sx q[2];
rz(1.4074116) q[2];
rz(0.47809005) q[3];
sx q[3];
rz(-2.7985088) q[3];
sx q[3];
rz(0.34496719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57426977) q[0];
sx q[0];
rz(-1.747921) q[0];
sx q[0];
rz(-2.0950914) q[0];
rz(0.25431713) q[1];
sx q[1];
rz(-2.3260702) q[1];
sx q[1];
rz(-2.8291124) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82332084) q[0];
sx q[0];
rz(-1.7004564) q[0];
sx q[0];
rz(1.983485) q[0];
rz(-pi) q[1];
rz(2.4813969) q[2];
sx q[2];
rz(-1.9593628) q[2];
sx q[2];
rz(1.6903433) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.1487138) q[1];
sx q[1];
rz(-1.8853469) q[1];
sx q[1];
rz(-1.3158362) q[1];
x q[2];
rz(-1.0572079) q[3];
sx q[3];
rz(-1.2826589) q[3];
sx q[3];
rz(2.6954034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.3089402) q[2];
sx q[2];
rz(-0.71844429) q[2];
sx q[2];
rz(0.18860513) q[2];
rz(-0.23379937) q[3];
sx q[3];
rz(-2.2457687) q[3];
sx q[3];
rz(-2.5578267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69498649) q[0];
sx q[0];
rz(-1.8291031) q[0];
sx q[0];
rz(-3.1332698) q[0];
rz(-0.43909973) q[1];
sx q[1];
rz(-1.3689901) q[1];
sx q[1];
rz(1.2459374) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3942892) q[0];
sx q[0];
rz(-0.4957333) q[0];
sx q[0];
rz(0.41252211) q[0];
rz(-pi) q[1];
rz(-2.9925084) q[2];
sx q[2];
rz(-1.8051762) q[2];
sx q[2];
rz(-0.54985505) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4250278) q[1];
sx q[1];
rz(-1.887117) q[1];
sx q[1];
rz(1.1891868) q[1];
x q[2];
rz(0.9160568) q[3];
sx q[3];
rz(-1.3406357) q[3];
sx q[3];
rz(-3.0482231) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.7370854) q[2];
sx q[2];
rz(-0.046701996) q[2];
sx q[2];
rz(-1.5076293) q[2];
rz(2.5493933) q[3];
sx q[3];
rz(-1.3785572) q[3];
sx q[3];
rz(0.73970214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1208039) q[0];
sx q[0];
rz(-1.4494267) q[0];
sx q[0];
rz(0.26594308) q[0];
rz(-1.2159011) q[1];
sx q[1];
rz(-0.78548702) q[1];
sx q[1];
rz(1.9507834) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62248373) q[0];
sx q[0];
rz(-1.568123) q[0];
sx q[0];
rz(0.83774211) q[0];
rz(-1.6295678) q[2];
sx q[2];
rz(-2.588039) q[2];
sx q[2];
rz(0.89828459) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.2289575) q[1];
sx q[1];
rz(-1.5785145) q[1];
sx q[1];
rz(0.75386366) q[1];
x q[2];
rz(2.9407752) q[3];
sx q[3];
rz(-2.8614223) q[3];
sx q[3];
rz(-2.7602989) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.33209458) q[2];
sx q[2];
rz(-1.3846493) q[2];
sx q[2];
rz(-0.14796999) q[2];
rz(-2.68908) q[3];
sx q[3];
rz(-1.0871525) q[3];
sx q[3];
rz(0.40542671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1389403) q[0];
sx q[0];
rz(-1.8754706) q[0];
sx q[0];
rz(1.7953605) q[0];
rz(3.1345308) q[1];
sx q[1];
rz(-2.4738753) q[1];
sx q[1];
rz(-0.5074358) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1088828) q[0];
sx q[0];
rz(-2.0883955) q[0];
sx q[0];
rz(0.93982307) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9865737) q[2];
sx q[2];
rz(-1.9117725) q[2];
sx q[2];
rz(0.91779681) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.9964357) q[1];
sx q[1];
rz(-1.5801589) q[1];
sx q[1];
rz(2.7476176) q[1];
x q[2];
rz(0.0094311992) q[3];
sx q[3];
rz(-1.8430897) q[3];
sx q[3];
rz(1.1444397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.70600447) q[2];
sx q[2];
rz(-2.0861349) q[2];
sx q[2];
rz(-2.9662507) q[2];
rz(0.38241479) q[3];
sx q[3];
rz(-1.9491111) q[3];
sx q[3];
rz(2.6800938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14313702) q[0];
sx q[0];
rz(-1.0112421) q[0];
sx q[0];
rz(1.3125516) q[0];
rz(3.0328499) q[1];
sx q[1];
rz(-2.5332632) q[1];
sx q[1];
rz(-0.67799062) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16081339) q[0];
sx q[0];
rz(-0.49404544) q[0];
sx q[0];
rz(1.52397) q[0];
rz(-pi) q[1];
rz(-0.74900519) q[2];
sx q[2];
rz(-1.555948) q[2];
sx q[2];
rz(2.2211563) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.99876632) q[1];
sx q[1];
rz(-1.8590392) q[1];
sx q[1];
rz(-1.0267797) q[1];
rz(-pi) q[2];
rz(0.59321065) q[3];
sx q[3];
rz(-0.81045356) q[3];
sx q[3];
rz(0.56014251) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.8678191) q[2];
sx q[2];
rz(-1.5742233) q[2];
sx q[2];
rz(-1.9847974) q[2];
rz(2.6371238) q[3];
sx q[3];
rz(-1.0320832) q[3];
sx q[3];
rz(-2.5696136) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2061763) q[0];
sx q[0];
rz(-1.6125866) q[0];
sx q[0];
rz(-0.61844283) q[0];
rz(-0.84856021) q[1];
sx q[1];
rz(-2.6237374) q[1];
sx q[1];
rz(0.79577622) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.764641) q[0];
sx q[0];
rz(-1.5466154) q[0];
sx q[0];
rz(-2.0410246) q[0];
x q[1];
rz(-0.82473849) q[2];
sx q[2];
rz(-1.5360263) q[2];
sx q[2];
rz(0.29037133) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.8164715) q[1];
sx q[1];
rz(-0.42342788) q[1];
sx q[1];
rz(-2.2864443) q[1];
x q[2];
rz(1.2451781) q[3];
sx q[3];
rz(-2.7005356) q[3];
sx q[3];
rz(-2.8056322) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.15647469) q[2];
sx q[2];
rz(-1.6343583) q[2];
sx q[2];
rz(1.9462684) q[2];
rz(-0.34504238) q[3];
sx q[3];
rz(-2.2796567) q[3];
sx q[3];
rz(-0.8663469) q[3];
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
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7095551) q[0];
sx q[0];
rz(-2.9032752) q[0];
sx q[0];
rz(-3.0639783) q[0];
rz(-1.9084825) q[1];
sx q[1];
rz(-1.0401007) q[1];
sx q[1];
rz(2.3495823) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38538489) q[0];
sx q[0];
rz(-1.0288219) q[0];
sx q[0];
rz(2.3535924) q[0];
rz(-pi) q[1];
rz(-0.77113232) q[2];
sx q[2];
rz(-0.22093102) q[2];
sx q[2];
rz(0.67459092) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.7866871) q[1];
sx q[1];
rz(-2.6514158) q[1];
sx q[1];
rz(-1.2398943) q[1];
rz(-pi) q[2];
x q[2];
rz(0.11802063) q[3];
sx q[3];
rz(-2.3343918) q[3];
sx q[3];
rz(-1.0420052) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.21976694) q[2];
sx q[2];
rz(-2.2302088) q[2];
sx q[2];
rz(2.0683973) q[2];
rz(-0.24751599) q[3];
sx q[3];
rz(-1.211834) q[3];
sx q[3];
rz(0.016935067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76219227) q[0];
sx q[0];
rz(-1.8397377) q[0];
sx q[0];
rz(2.4158438) q[0];
rz(2.2629867) q[1];
sx q[1];
rz(-0.85522063) q[1];
sx q[1];
rz(-1.9488889) q[1];
rz(1.0868308) q[2];
sx q[2];
rz(-1.7499583) q[2];
sx q[2];
rz(-1.3365895) q[2];
rz(-2.2976919) q[3];
sx q[3];
rz(-1.6132728) q[3];
sx q[3];
rz(-2.1584395) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
