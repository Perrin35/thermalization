OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.20004162) q[0];
sx q[0];
rz(3.3257421) q[0];
sx q[0];
rz(10.781547) q[0];
rz(1.4461203) q[1];
sx q[1];
rz(-1.5674011) q[1];
sx q[1];
rz(-2.8653436) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7295099) q[0];
sx q[0];
rz(-2.9449521) q[0];
sx q[0];
rz(1.6565408) q[0];
x q[1];
rz(-3.0577781) q[2];
sx q[2];
rz(-3.028125) q[2];
sx q[2];
rz(-2.6423315) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.4893814) q[1];
sx q[1];
rz(-1.6141574) q[1];
sx q[1];
rz(1.5549589) q[1];
rz(-pi) q[2];
rz(-1.6365693) q[3];
sx q[3];
rz(-1.2322174) q[3];
sx q[3];
rz(-2.8985646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.9487298) q[2];
sx q[2];
rz(-0.36878815) q[2];
sx q[2];
rz(-0.31664872) q[2];
rz(-0.51490319) q[3];
sx q[3];
rz(-3.1334936) q[3];
sx q[3];
rz(0.81075877) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0136451) q[0];
sx q[0];
rz(-0.34334308) q[0];
sx q[0];
rz(-3.1105594) q[0];
rz(1.4399928) q[1];
sx q[1];
rz(-2.3411935) q[1];
sx q[1];
rz(1.7354785) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7714376) q[0];
sx q[0];
rz(-1.9906128) q[0];
sx q[0];
rz(0.19924723) q[0];
rz(-pi) q[1];
x q[1];
rz(0.038766301) q[2];
sx q[2];
rz(-1.3666461) q[2];
sx q[2];
rz(0.37224712) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.5819491) q[1];
sx q[1];
rz(-0.20155263) q[1];
sx q[1];
rz(-2.6879716) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3651016) q[3];
sx q[3];
rz(-1.6776976) q[3];
sx q[3];
rz(2.2872146) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.7180353) q[2];
sx q[2];
rz(-0.2090629) q[2];
sx q[2];
rz(-2.5162856) q[2];
rz(-1.5705671) q[3];
sx q[3];
rz(-0.26206854) q[3];
sx q[3];
rz(3.1077793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5728979) q[0];
sx q[0];
rz(-2.608572) q[0];
sx q[0];
rz(-2.0030588) q[0];
rz(2.2472527) q[1];
sx q[1];
rz(-3.1411451) q[1];
sx q[1];
rz(-2.0785887) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.244449) q[0];
sx q[0];
rz(-0.70310601) q[0];
sx q[0];
rz(-1.1887111) q[0];
rz(-pi) q[1];
rz(2.393542) q[2];
sx q[2];
rz(-1.7123634) q[2];
sx q[2];
rz(-0.54439616) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6394807) q[1];
sx q[1];
rz(-1.6422137) q[1];
sx q[1];
rz(1.7014496) q[1];
x q[2];
rz(1.7785709) q[3];
sx q[3];
rz(-1.7881743) q[3];
sx q[3];
rz(0.3849963) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.2669617) q[2];
sx q[2];
rz(-2.3571372) q[2];
sx q[2];
rz(1.1355012) q[2];
rz(1.1442432) q[3];
sx q[3];
rz(-0.091322986) q[3];
sx q[3];
rz(1.7448366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3369694) q[0];
sx q[0];
rz(-0.40358821) q[0];
sx q[0];
rz(-1.1176156) q[0];
rz(-3.0506253) q[1];
sx q[1];
rz(-0.03799835) q[1];
sx q[1];
rz(1.5100347) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6081597) q[0];
sx q[0];
rz(-1.3638487) q[0];
sx q[0];
rz(0.98790913) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.48520605) q[2];
sx q[2];
rz(-0.23113568) q[2];
sx q[2];
rz(0.71454465) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.9659057) q[1];
sx q[1];
rz(-1.693778) q[1];
sx q[1];
rz(-1.8781652) q[1];
x q[2];
rz(-0.7078739) q[3];
sx q[3];
rz(-2.3972569) q[3];
sx q[3];
rz(-2.3553948) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.1893922) q[2];
sx q[2];
rz(-0.11710937) q[2];
sx q[2];
rz(1.2516775) q[2];
rz(0.39024726) q[3];
sx q[3];
rz(-0.016642112) q[3];
sx q[3];
rz(2.9515631) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2802884) q[0];
sx q[0];
rz(-2.5087293) q[0];
sx q[0];
rz(1.2907775) q[0];
rz(3.1351008) q[1];
sx q[1];
rz(-0.26632729) q[1];
sx q[1];
rz(-1.3351306) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4075942) q[0];
sx q[0];
rz(-1.4346135) q[0];
sx q[0];
rz(-0.34466095) q[0];
x q[1];
rz(-3.0090121) q[2];
sx q[2];
rz(-2.6040316) q[2];
sx q[2];
rz(-0.81576306) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8346429) q[1];
sx q[1];
rz(-1.4473302) q[1];
sx q[1];
rz(-2.929652) q[1];
rz(-pi) q[2];
rz(1.9919768) q[3];
sx q[3];
rz(-1.4011746) q[3];
sx q[3];
rz(-2.6785525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.73162762) q[2];
sx q[2];
rz(-0.16606398) q[2];
sx q[2];
rz(-0.49327332) q[2];
rz(2.0491681) q[3];
sx q[3];
rz(-0.0041996669) q[3];
sx q[3];
rz(-2.5815729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4527556) q[0];
sx q[0];
rz(-1.1863363) q[0];
sx q[0];
rz(1.3415927) q[0];
rz(1.5635368) q[1];
sx q[1];
rz(-0.23461239) q[1];
sx q[1];
rz(-2.2685307) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4463947) q[0];
sx q[0];
rz(-1.4322011) q[0];
sx q[0];
rz(2.6202484) q[0];
rz(-pi) q[1];
rz(-2.701917) q[2];
sx q[2];
rz(-1.623325) q[2];
sx q[2];
rz(0.66212468) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(3.0256589) q[1];
sx q[1];
rz(-1.4486989) q[1];
sx q[1];
rz(3.0662838) q[1];
rz(1.3033684) q[3];
sx q[3];
rz(-1.7828327) q[3];
sx q[3];
rz(2.7515281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.0688429) q[2];
sx q[2];
rz(-1.170697) q[2];
sx q[2];
rz(1.076131) q[2];
rz(-0.38232803) q[3];
sx q[3];
rz(-3.101442) q[3];
sx q[3];
rz(0.7676777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
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
rz(-0.024993984) q[0];
sx q[0];
rz(-3.0285663) q[0];
sx q[0];
rz(0.22079994) q[0];
rz(2.2258672) q[1];
sx q[1];
rz(-2.974739) q[1];
sx q[1];
rz(1.6545666) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2325848) q[0];
sx q[0];
rz(-0.56031432) q[0];
sx q[0];
rz(2.4057897) q[0];
x q[1];
rz(1.9369097) q[2];
sx q[2];
rz(-1.7559474) q[2];
sx q[2];
rz(-2.2529035) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.88363632) q[1];
sx q[1];
rz(-0.90054926) q[1];
sx q[1];
rz(0.61780091) q[1];
x q[2];
rz(1.4288032) q[3];
sx q[3];
rz(-1.4513223) q[3];
sx q[3];
rz(0.17585299) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.23358341) q[2];
sx q[2];
rz(-3.0685232) q[2];
sx q[2];
rz(1.118534) q[2];
rz(1.7782878) q[3];
sx q[3];
rz(-0.076863591) q[3];
sx q[3];
rz(-1.793687) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8175279) q[0];
sx q[0];
rz(-2.9370566) q[0];
sx q[0];
rz(1.2943785) q[0];
rz(1.4140465) q[1];
sx q[1];
rz(-0.050948016) q[1];
sx q[1];
rz(2.9185413) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4572499) q[0];
sx q[0];
rz(-0.97829223) q[0];
sx q[0];
rz(-2.1564583) q[0];
x q[1];
rz(1.7086008) q[2];
sx q[2];
rz(-1.3827822) q[2];
sx q[2];
rz(1.6013946) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.9893045) q[1];
sx q[1];
rz(-2.3803601) q[1];
sx q[1];
rz(-0.95639317) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.11745061) q[3];
sx q[3];
rz(-0.87908213) q[3];
sx q[3];
rz(2.1827616) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.602542) q[2];
sx q[2];
rz(-1.9571303) q[2];
sx q[2];
rz(-2.7089673) q[2];
rz(2.1116347) q[3];
sx q[3];
rz(-3.0410788) q[3];
sx q[3];
rz(-2.4242134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3887727) q[0];
sx q[0];
rz(-2.6322375) q[0];
sx q[0];
rz(2.9260522) q[0];
rz(-2.0095677) q[1];
sx q[1];
rz(-2.3339381) q[1];
sx q[1];
rz(-1.4707123) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0316353) q[0];
sx q[0];
rz(-0.67207974) q[0];
sx q[0];
rz(2.0874778) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1342513) q[2];
sx q[2];
rz(-1.2255395) q[2];
sx q[2];
rz(-0.79011142) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.017478213) q[1];
sx q[1];
rz(-0.062659293) q[1];
sx q[1];
rz(-2.6397328) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3683794) q[3];
sx q[3];
rz(-0.63211956) q[3];
sx q[3];
rz(-2.8912197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1620764) q[2];
sx q[2];
rz(-0.039204892) q[2];
sx q[2];
rz(2.6230679) q[2];
rz(2.3833852) q[3];
sx q[3];
rz(-0.77986217) q[3];
sx q[3];
rz(1.6583748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-2.6766847) q[0];
sx q[0];
rz(-1.9187036) q[0];
sx q[0];
rz(-1.884888) q[0];
rz(-2.8769809) q[1];
sx q[1];
rz(-3.1393317) q[1];
sx q[1];
rz(1.1983926) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4751202) q[0];
sx q[0];
rz(-1.571093) q[0];
sx q[0];
rz(-0.0065684321) q[0];
rz(-3.0727238) q[2];
sx q[2];
rz(-1.9039717) q[2];
sx q[2];
rz(-2.888713) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.2923515) q[1];
sx q[1];
rz(-1.3363943) q[1];
sx q[1];
rz(-1.5567907) q[1];
rz(-pi) q[2];
rz(1.7400473) q[3];
sx q[3];
rz(-1.5334606) q[3];
sx q[3];
rz(-2.7558655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.6780702) q[2];
sx q[2];
rz(-0.0040201298) q[2];
sx q[2];
rz(-3.1275952) q[2];
rz(2.8409581) q[3];
sx q[3];
rz(-0.0036573452) q[3];
sx q[3];
rz(1.5699056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32776253) q[0];
sx q[0];
rz(-1.7131085) q[0];
sx q[0];
rz(-1.499927) q[0];
rz(-0.035540237) q[1];
sx q[1];
rz(-2.7062369) q[1];
sx q[1];
rz(-2.8027262) q[1];
rz(-0.28066228) q[2];
sx q[2];
rz(-0.71692313) q[2];
sx q[2];
rz(3.0871152) q[2];
rz(0.20094674) q[3];
sx q[3];
rz(-2.3280716) q[3];
sx q[3];
rz(0.18542326) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
