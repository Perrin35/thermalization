OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.5947333) q[0];
sx q[0];
rz(-1.5164627) q[0];
sx q[0];
rz(0.2642785) q[0];
rz(2.1677986) q[1];
sx q[1];
rz(-1.9314613) q[1];
sx q[1];
rz(2.4063453) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7289294) q[0];
sx q[0];
rz(-2.465415) q[0];
sx q[0];
rz(0.23763188) q[0];
rz(0.67970694) q[2];
sx q[2];
rz(-1.9787111) q[2];
sx q[2];
rz(-2.9349875) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.6890251) q[1];
sx q[1];
rz(-1.2847632) q[1];
sx q[1];
rz(0.16098117) q[1];
rz(-0.58044503) q[3];
sx q[3];
rz(-1.5942897) q[3];
sx q[3];
rz(-1.98222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.66951093) q[2];
sx q[2];
rz(-1.3005723) q[2];
sx q[2];
rz(1.1038587) q[2];
rz(-1.8707229) q[3];
sx q[3];
rz(-1.2277675) q[3];
sx q[3];
rz(2.8675573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(2.1141777) q[0];
sx q[0];
rz(-0.52863055) q[0];
sx q[0];
rz(0.43637481) q[0];
rz(2.6787058) q[1];
sx q[1];
rz(-1.0375689) q[1];
sx q[1];
rz(-2.8754821) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3811831) q[0];
sx q[0];
rz(-1.1927483) q[0];
sx q[0];
rz(-0.80947431) q[0];
rz(1.3219464) q[2];
sx q[2];
rz(-2.5517533) q[2];
sx q[2];
rz(0.8283386) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.7539566) q[1];
sx q[1];
rz(-2.1093183) q[1];
sx q[1];
rz(0.83420475) q[1];
rz(-pi) q[2];
rz(2.9588685) q[3];
sx q[3];
rz(-2.1696739) q[3];
sx q[3];
rz(3.0298508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.6767072) q[2];
sx q[2];
rz(-1.2767982) q[2];
sx q[2];
rz(-2.6300988) q[2];
rz(0.809508) q[3];
sx q[3];
rz(-1.6102689) q[3];
sx q[3];
rz(0.28765837) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70616102) q[0];
sx q[0];
rz(-1.5937188) q[0];
sx q[0];
rz(0.92873746) q[0];
rz(-1.4061032) q[1];
sx q[1];
rz(-0.7000674) q[1];
sx q[1];
rz(1.3471289) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12362326) q[0];
sx q[0];
rz(-0.53640134) q[0];
sx q[0];
rz(-0.53039741) q[0];
x q[1];
rz(2.9834619) q[2];
sx q[2];
rz(-1.2889382) q[2];
sx q[2];
rz(-2.6513197) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.0128847) q[1];
sx q[1];
rz(-1.8599659) q[1];
sx q[1];
rz(0.94990001) q[1];
x q[2];
rz(1.7473162) q[3];
sx q[3];
rz(-0.71478292) q[3];
sx q[3];
rz(0.96283462) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.9946263) q[2];
sx q[2];
rz(-1.7549843) q[2];
sx q[2];
rz(-1.5220801) q[2];
rz(-2.8772723) q[3];
sx q[3];
rz(-2.1051354) q[3];
sx q[3];
rz(-0.49595293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3494444) q[0];
sx q[0];
rz(-1.9932207) q[0];
sx q[0];
rz(2.175892) q[0];
rz(2.4194338) q[1];
sx q[1];
rz(-1.637371) q[1];
sx q[1];
rz(2.5818363) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9352919) q[0];
sx q[0];
rz(-1.1456523) q[0];
sx q[0];
rz(-1.6598808) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.70978769) q[2];
sx q[2];
rz(-0.94883942) q[2];
sx q[2];
rz(-0.70993916) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.2985845) q[1];
sx q[1];
rz(-2.1857939) q[1];
sx q[1];
rz(-3.0857012) q[1];
rz(-0.23493725) q[3];
sx q[3];
rz(-2.3384691) q[3];
sx q[3];
rz(-2.0149751) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.7136148) q[2];
sx q[2];
rz(-1.5087936) q[2];
sx q[2];
rz(-1.7170061) q[2];
rz(-2.8811841) q[3];
sx q[3];
rz(-1.7491165) q[3];
sx q[3];
rz(-0.4310472) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.150862) q[0];
sx q[0];
rz(-1.6435511) q[0];
sx q[0];
rz(0.36636233) q[0];
rz(-1.5953966) q[1];
sx q[1];
rz(-2.5876744) q[1];
sx q[1];
rz(0.34367925) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2491964) q[0];
sx q[0];
rz(-2.2049892) q[0];
sx q[0];
rz(0.87974324) q[0];
rz(-pi) q[1];
rz(0.91243773) q[2];
sx q[2];
rz(-2.3235333) q[2];
sx q[2];
rz(2.22544) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.57203612) q[1];
sx q[1];
rz(-1.0735638) q[1];
sx q[1];
rz(-3.0576502) q[1];
x q[2];
rz(-3.120317) q[3];
sx q[3];
rz(-0.92704232) q[3];
sx q[3];
rz(1.6568041) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.7498103) q[2];
sx q[2];
rz(-2.2822773) q[2];
sx q[2];
rz(-3.0878477) q[2];
rz(1.404445) q[3];
sx q[3];
rz(-2.6440547) q[3];
sx q[3];
rz(-0.29156175) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6546201) q[0];
sx q[0];
rz(-0.64990652) q[0];
sx q[0];
rz(0.75575954) q[0];
rz(-0.02515633) q[1];
sx q[1];
rz(-2.2143366) q[1];
sx q[1];
rz(-2.8818534) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4070968) q[0];
sx q[0];
rz(-1.8844814) q[0];
sx q[0];
rz(-1.19019) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2345215) q[2];
sx q[2];
rz(-1.5521126) q[2];
sx q[2];
rz(-2.3980354) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.427634) q[1];
sx q[1];
rz(-0.40363388) q[1];
sx q[1];
rz(2.1778818) q[1];
rz(0.017976956) q[3];
sx q[3];
rz(-2.020105) q[3];
sx q[3];
rz(-2.9177005) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.48866895) q[2];
sx q[2];
rz(-0.80604625) q[2];
sx q[2];
rz(-3.0409813) q[2];
rz(0.18209022) q[3];
sx q[3];
rz(-0.87825769) q[3];
sx q[3];
rz(1.7939059) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
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
rz(-1.1473734) q[0];
sx q[0];
rz(-1.7213151) q[0];
sx q[0];
rz(2.8314262) q[0];
rz(-2.639333) q[1];
sx q[1];
rz(-0.52400932) q[1];
sx q[1];
rz(2.535634) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3340671) q[0];
sx q[0];
rz(-1.7958607) q[0];
sx q[0];
rz(-2.5146757) q[0];
x q[1];
rz(0.29166834) q[2];
sx q[2];
rz(-1.3783611) q[2];
sx q[2];
rz(0.33484909) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1389321) q[1];
sx q[1];
rz(-1.9117038) q[1];
sx q[1];
rz(-2.1041811) q[1];
rz(-pi) q[2];
rz(-0.53020729) q[3];
sx q[3];
rz(-1.0876417) q[3];
sx q[3];
rz(2.9999441) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.24017748) q[2];
sx q[2];
rz(-1.9577953) q[2];
sx q[2];
rz(-0.85285464) q[2];
rz(1.7715706) q[3];
sx q[3];
rz(-1.4586689) q[3];
sx q[3];
rz(0.20496932) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2296427) q[0];
sx q[0];
rz(-2.5456173) q[0];
sx q[0];
rz(-1.6802616) q[0];
rz(-1.7386859) q[1];
sx q[1];
rz(-2.1673514) q[1];
sx q[1];
rz(-3.0775552) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94992276) q[0];
sx q[0];
rz(-0.64282286) q[0];
sx q[0];
rz(-1.2547917) q[0];
x q[1];
rz(-2.089558) q[2];
sx q[2];
rz(-0.56781893) q[2];
sx q[2];
rz(1.0713112) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.0361955) q[1];
sx q[1];
rz(-0.35696778) q[1];
sx q[1];
rz(3.0371975) q[1];
rz(-pi) q[2];
rz(-1.6927034) q[3];
sx q[3];
rz(-0.31414437) q[3];
sx q[3];
rz(-1.7243408) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-3.0503851) q[2];
sx q[2];
rz(-0.63082266) q[2];
sx q[2];
rz(-1.5861661) q[2];
rz(-2.2533916) q[3];
sx q[3];
rz(-1.9017838) q[3];
sx q[3];
rz(2.2122038) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14426194) q[0];
sx q[0];
rz(-2.0781131) q[0];
sx q[0];
rz(-1.4319179) q[0];
rz(-0.56888467) q[1];
sx q[1];
rz(-2.6061997) q[1];
sx q[1];
rz(-1.127839) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9147946) q[0];
sx q[0];
rz(-0.1495805) q[0];
sx q[0];
rz(0.10766715) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.608689) q[2];
sx q[2];
rz(-2.1295296) q[2];
sx q[2];
rz(1.2677873) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9319218) q[1];
sx q[1];
rz(-1.5823963) q[1];
sx q[1];
rz(-1.8124707) q[1];
rz(-pi) q[2];
rz(0.3663775) q[3];
sx q[3];
rz(-1.0421703) q[3];
sx q[3];
rz(-0.60929326) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.613712) q[2];
sx q[2];
rz(-0.85313672) q[2];
sx q[2];
rz(2.9679427) q[2];
rz(-2.8052143) q[3];
sx q[3];
rz(-1.9189546) q[3];
sx q[3];
rz(0.39150795) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7062475) q[0];
sx q[0];
rz(-2.5543537) q[0];
sx q[0];
rz(1.6760814) q[0];
rz(-0.82410518) q[1];
sx q[1];
rz(-1.5287639) q[1];
sx q[1];
rz(0.5724268) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51047072) q[0];
sx q[0];
rz(-2.0513751) q[0];
sx q[0];
rz(-0.62453385) q[0];
rz(-1.053439) q[2];
sx q[2];
rz(-2.6247019) q[2];
sx q[2];
rz(0.042339485) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(3.010667) q[1];
sx q[1];
rz(-1.3548684) q[1];
sx q[1];
rz(2.8972577) q[1];
x q[2];
rz(-2.4329348) q[3];
sx q[3];
rz(-1.7950247) q[3];
sx q[3];
rz(-1.538016) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.3616025) q[2];
sx q[2];
rz(-2.6670691) q[2];
sx q[2];
rz(2.6043716) q[2];
rz(1.0572664) q[3];
sx q[3];
rz(-0.89151645) q[3];
sx q[3];
rz(2.4479772) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.854241) q[0];
sx q[0];
rz(-1.1653405) q[0];
sx q[0];
rz(-1.5821138) q[0];
rz(-2.6782425) q[1];
sx q[1];
rz(-2.2644823) q[1];
sx q[1];
rz(1.5092441) q[1];
rz(2.2255185) q[2];
sx q[2];
rz(-1.4929885) q[2];
sx q[2];
rz(2.3103466) q[2];
rz(0.061999576) q[3];
sx q[3];
rz(-2.0580895) q[3];
sx q[3];
rz(2.5930391) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
