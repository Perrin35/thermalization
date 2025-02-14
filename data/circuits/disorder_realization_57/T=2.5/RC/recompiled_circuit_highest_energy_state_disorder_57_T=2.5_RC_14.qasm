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
rz(1.3996742) q[0];
sx q[0];
rz(1.4631441) q[0];
sx q[0];
rz(14.027496) q[0];
rz(-0.48179102) q[1];
sx q[1];
rz(-0.81463373) q[1];
sx q[1];
rz(0.93946594) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50892539) q[0];
sx q[0];
rz(-1.5611751) q[0];
sx q[0];
rz(-1.6144572) q[0];
rz(-pi) q[1];
x q[1];
rz(0.58161503) q[2];
sx q[2];
rz(-2.5293859) q[2];
sx q[2];
rz(-0.78919995) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.11096) q[1];
sx q[1];
rz(-0.041555066) q[1];
sx q[1];
rz(1.3135733) q[1];
rz(-pi) q[2];
rz(-1.788673) q[3];
sx q[3];
rz(-1.9697726) q[3];
sx q[3];
rz(0.50673317) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.7751094) q[2];
sx q[2];
rz(-1.4668239) q[2];
sx q[2];
rz(0.62466204) q[2];
rz(2.1810253) q[3];
sx q[3];
rz(-1.4550236) q[3];
sx q[3];
rz(-2.2666986) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7333154) q[0];
sx q[0];
rz(-0.69077078) q[0];
sx q[0];
rz(0.1917924) q[0];
rz(-2.5674112) q[1];
sx q[1];
rz(-1.5230813) q[1];
sx q[1];
rz(-2.7343553) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8438661) q[0];
sx q[0];
rz(-2.9653984) q[0];
sx q[0];
rz(-1.5017444) q[0];
rz(-pi) q[1];
rz(-3.0767976) q[2];
sx q[2];
rz(-1.7470834) q[2];
sx q[2];
rz(0.48688146) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.15528725) q[1];
sx q[1];
rz(-2.1006507) q[1];
sx q[1];
rz(1.3572567) q[1];
rz(-0.089669946) q[3];
sx q[3];
rz(-1.0700135) q[3];
sx q[3];
rz(1.1550241) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.9521956) q[2];
sx q[2];
rz(-1.1789097) q[2];
sx q[2];
rz(0.65323812) q[2];
rz(2.2630528) q[3];
sx q[3];
rz(-2.0432751) q[3];
sx q[3];
rz(-3.0285335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1427796) q[0];
sx q[0];
rz(-1.7072562) q[0];
sx q[0];
rz(-2.0913731) q[0];
rz(0.6160008) q[1];
sx q[1];
rz(-1.7237677) q[1];
sx q[1];
rz(2.2284257) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28024689) q[0];
sx q[0];
rz(-1.6771835) q[0];
sx q[0];
rz(1.7411355) q[0];
rz(-pi) q[1];
rz(3.0229125) q[2];
sx q[2];
rz(-1.8556229) q[2];
sx q[2];
rz(-0.26155805) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.84308594) q[1];
sx q[1];
rz(-1.32846) q[1];
sx q[1];
rz(2.0451115) q[1];
x q[2];
rz(-2.764934) q[3];
sx q[3];
rz(-0.91729255) q[3];
sx q[3];
rz(-2.3953826) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.22127613) q[2];
sx q[2];
rz(-0.90711275) q[2];
sx q[2];
rz(0.059180666) q[2];
rz(-2.3024998) q[3];
sx q[3];
rz(-1.9234761) q[3];
sx q[3];
rz(2.2742417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3686104) q[0];
sx q[0];
rz(-2.1324069) q[0];
sx q[0];
rz(2.7918145) q[0];
rz(-1.3444258) q[1];
sx q[1];
rz(-0.58660048) q[1];
sx q[1];
rz(0.1127359) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7736334) q[0];
sx q[0];
rz(-1.9122352) q[0];
sx q[0];
rz(-0.48078791) q[0];
rz(-pi) q[1];
rz(0.91103986) q[2];
sx q[2];
rz(-1.8000812) q[2];
sx q[2];
rz(-2.2509172) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.73600125) q[1];
sx q[1];
rz(-1.8718157) q[1];
sx q[1];
rz(-0.89717866) q[1];
rz(-0.76072089) q[3];
sx q[3];
rz(-0.84555999) q[3];
sx q[3];
rz(-1.4694422) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3939487) q[2];
sx q[2];
rz(-1.5564206) q[2];
sx q[2];
rz(0.057223884) q[2];
rz(-1.2563489) q[3];
sx q[3];
rz(-2.5689954) q[3];
sx q[3];
rz(-0.73067874) q[3];
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
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8771862) q[0];
sx q[0];
rz(-1.974396) q[0];
sx q[0];
rz(2.4315244) q[0];
rz(2.8801628) q[1];
sx q[1];
rz(-1.5855887) q[1];
sx q[1];
rz(2.1905621) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20999058) q[0];
sx q[0];
rz(-1.8977802) q[0];
sx q[0];
rz(0.47173421) q[0];
x q[1];
rz(-2.0567953) q[2];
sx q[2];
rz(-2.1446781) q[2];
sx q[2];
rz(-1.4950372) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.1360395) q[1];
sx q[1];
rz(-2.4085718) q[1];
sx q[1];
rz(-0.3049149) q[1];
rz(-pi) q[2];
x q[2];
rz(0.099154648) q[3];
sx q[3];
rz(-2.2014849) q[3];
sx q[3];
rz(-2.658228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.6422673) q[2];
sx q[2];
rz(-2.817481) q[2];
sx q[2];
rz(0.30936852) q[2];
rz(-1.0950836) q[3];
sx q[3];
rz(-2.0131922) q[3];
sx q[3];
rz(-0.40422082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9139022) q[0];
sx q[0];
rz(-0.18276437) q[0];
sx q[0];
rz(-0.90232724) q[0];
rz(-3.0790216) q[1];
sx q[1];
rz(-1.8726655) q[1];
sx q[1];
rz(1.8960457) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.093687261) q[0];
sx q[0];
rz(-1.7249247) q[0];
sx q[0];
rz(2.4612294) q[0];
x q[1];
rz(-0.056137265) q[2];
sx q[2];
rz(-1.0318163) q[2];
sx q[2];
rz(-2.7201729) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.9334792) q[1];
sx q[1];
rz(-1.1565349) q[1];
sx q[1];
rz(0.53908344) q[1];
x q[2];
rz(2.4687103) q[3];
sx q[3];
rz(-1.1816634) q[3];
sx q[3];
rz(-2.0418039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7237599) q[2];
sx q[2];
rz(-1.87169) q[2];
sx q[2];
rz(-1.343824) q[2];
rz(1.8580681) q[3];
sx q[3];
rz(-0.44461861) q[3];
sx q[3];
rz(2.4166687) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.035324) q[0];
sx q[0];
rz(-3.1043053) q[0];
sx q[0];
rz(-0.37621793) q[0];
rz(0.47209921) q[1];
sx q[1];
rz(-2.4596228) q[1];
sx q[1];
rz(-0.45904407) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7721586) q[0];
sx q[0];
rz(-0.87352814) q[0];
sx q[0];
rz(-2.5791772) q[0];
rz(-pi) q[1];
x q[1];
rz(0.44957513) q[2];
sx q[2];
rz(-1.4130985) q[2];
sx q[2];
rz(1.1166935) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.3891148) q[1];
sx q[1];
rz(-1.6387617) q[1];
sx q[1];
rz(-0.1631921) q[1];
x q[2];
rz(-1.5581896) q[3];
sx q[3];
rz(-0.9731826) q[3];
sx q[3];
rz(-1.2659092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(3.0692733) q[2];
sx q[2];
rz(-0.71162144) q[2];
sx q[2];
rz(0.35815987) q[2];
rz(2.2530796) q[3];
sx q[3];
rz(-1.0212967) q[3];
sx q[3];
rz(-2.2108938) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.096980378) q[0];
sx q[0];
rz(-2.8428069) q[0];
sx q[0];
rz(0.97691798) q[0];
rz(2.3878429) q[1];
sx q[1];
rz(-1.6981373) q[1];
sx q[1];
rz(-1.5026198) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96732765) q[0];
sx q[0];
rz(-2.2696724) q[0];
sx q[0];
rz(1.4682659) q[0];
rz(0.66960488) q[2];
sx q[2];
rz(-1.9960224) q[2];
sx q[2];
rz(-2.2086672) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.0148919) q[1];
sx q[1];
rz(-2.3897244) q[1];
sx q[1];
rz(0.45104646) q[1];
rz(2.9108293) q[3];
sx q[3];
rz(-0.9597336) q[3];
sx q[3];
rz(0.83702528) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.94540191) q[2];
sx q[2];
rz(-1.6795029) q[2];
sx q[2];
rz(-2.2958344) q[2];
rz(2.3614007) q[3];
sx q[3];
rz(-1.3581685) q[3];
sx q[3];
rz(0.63854027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82275259) q[0];
sx q[0];
rz(-1.0718811) q[0];
sx q[0];
rz(2.7986797) q[0];
rz(-1.0036537) q[1];
sx q[1];
rz(-1.9099078) q[1];
sx q[1];
rz(-0.71663219) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6824237) q[0];
sx q[0];
rz(-2.3083516) q[0];
sx q[0];
rz(2.6618945) q[0];
x q[1];
rz(-2.6478397) q[2];
sx q[2];
rz(-1.2896364) q[2];
sx q[2];
rz(-1.8914831) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.20880675) q[1];
sx q[1];
rz(-2.1541671) q[1];
sx q[1];
rz(1.8933184) q[1];
x q[2];
rz(1.9699205) q[3];
sx q[3];
rz(-2.2033327) q[3];
sx q[3];
rz(0.61033953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.0922682) q[2];
sx q[2];
rz(-1.6246139) q[2];
sx q[2];
rz(2.1410904) q[2];
rz(2.3246824) q[3];
sx q[3];
rz(-1.2714081) q[3];
sx q[3];
rz(-2.7774155) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7526806) q[0];
sx q[0];
rz(-2.2081544) q[0];
sx q[0];
rz(-2.3626784) q[0];
rz(-2.3498416) q[1];
sx q[1];
rz(-1.9154895) q[1];
sx q[1];
rz(2.7955999) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8824005) q[0];
sx q[0];
rz(-2.0171875) q[0];
sx q[0];
rz(1.4784228) q[0];
rz(-1.6093639) q[2];
sx q[2];
rz(-0.30353217) q[2];
sx q[2];
rz(-2.8736049) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.78948632) q[1];
sx q[1];
rz(-2.5946684) q[1];
sx q[1];
rz(1.5032306) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7818809) q[3];
sx q[3];
rz(-0.87950686) q[3];
sx q[3];
rz(0.018010294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.9927696) q[2];
sx q[2];
rz(-1.537348) q[2];
sx q[2];
rz(2.0796622) q[2];
rz(-2.3289833) q[3];
sx q[3];
rz(-1.9989697) q[3];
sx q[3];
rz(0.44212166) q[3];
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
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84520311) q[0];
sx q[0];
rz(-2.3136105) q[0];
sx q[0];
rz(-1.4985341) q[0];
rz(-2.8614112) q[1];
sx q[1];
rz(-1.9277086) q[1];
sx q[1];
rz(2.3452506) q[1];
rz(2.7325999) q[2];
sx q[2];
rz(-1.6074642) q[2];
sx q[2];
rz(1.4516914) q[2];
rz(0.82874684) q[3];
sx q[3];
rz(-1.8846877) q[3];
sx q[3];
rz(1.7170513) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
