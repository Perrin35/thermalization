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
rz(-2.3081245) q[0];
sx q[0];
rz(3.9552116) q[0];
sx q[0];
rz(7.8501346) q[0];
rz(-1.1235224) q[1];
sx q[1];
rz(-0.82668537) q[1];
sx q[1];
rz(-2.6348662) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6228765) q[0];
sx q[0];
rz(-0.55326792) q[0];
sx q[0];
rz(3.1068463) q[0];
rz(-pi) q[1];
rz(2.0000524) q[2];
sx q[2];
rz(-2.5876293) q[2];
sx q[2];
rz(-2.9316154) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.0953137) q[1];
sx q[1];
rz(-1.6254402) q[1];
sx q[1];
rz(-2.8846963) q[1];
rz(-pi) q[2];
rz(1.5208901) q[3];
sx q[3];
rz(-2.4766233) q[3];
sx q[3];
rz(-1.7902399) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.43164051) q[2];
sx q[2];
rz(-1.9441354) q[2];
sx q[2];
rz(2.6006827) q[2];
rz(-1.8139808) q[3];
sx q[3];
rz(-1.6874467) q[3];
sx q[3];
rz(-0.66837627) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3427908) q[0];
sx q[0];
rz(-1.3725766) q[0];
sx q[0];
rz(1.0569245) q[0];
rz(2.7150555) q[1];
sx q[1];
rz(-1.5142454) q[1];
sx q[1];
rz(0.0043409745) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26165923) q[0];
sx q[0];
rz(-0.78908217) q[0];
sx q[0];
rz(2.1270149) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0356194) q[2];
sx q[2];
rz(-0.7053203) q[2];
sx q[2];
rz(-0.83503228) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.8651826) q[1];
sx q[1];
rz(-1.4917871) q[1];
sx q[1];
rz(-3.1311796) q[1];
x q[2];
rz(3.0375157) q[3];
sx q[3];
rz(-0.53734578) q[3];
sx q[3];
rz(1.4817103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.1987552) q[2];
sx q[2];
rz(-1.402907) q[2];
sx q[2];
rz(2.7665561) q[2];
rz(-0.047920553) q[3];
sx q[3];
rz(-1.5791357) q[3];
sx q[3];
rz(-0.16938773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38309836) q[0];
sx q[0];
rz(-0.57621813) q[0];
sx q[0];
rz(-1.5118442) q[0];
rz(-1.0668628) q[1];
sx q[1];
rz(-1.8424415) q[1];
sx q[1];
rz(0.80924928) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.209037) q[0];
sx q[0];
rz(-1.5309628) q[0];
sx q[0];
rz(1.6074384) q[0];
rz(-pi) q[1];
x q[1];
rz(0.079945076) q[2];
sx q[2];
rz(-0.71942753) q[2];
sx q[2];
rz(1.5647174) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.0195649) q[1];
sx q[1];
rz(-2.3246336) q[1];
sx q[1];
rz(-2.8754068) q[1];
rz(-pi) q[2];
rz(3.1339614) q[3];
sx q[3];
rz(-0.3618218) q[3];
sx q[3];
rz(-2.1280888) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.2201436) q[2];
sx q[2];
rz(-0.67110944) q[2];
sx q[2];
rz(1.6371833) q[2];
rz(-0.97357059) q[3];
sx q[3];
rz(-0.74876553) q[3];
sx q[3];
rz(-1.0934632) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1188537) q[0];
sx q[0];
rz(-1.5596507) q[0];
sx q[0];
rz(-2.7795025) q[0];
rz(1.3478966) q[1];
sx q[1];
rz(-1.3840414) q[1];
sx q[1];
rz(-2.0787584) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84447538) q[0];
sx q[0];
rz(-1.3985635) q[0];
sx q[0];
rz(-2.9976241) q[0];
rz(-pi) q[1];
rz(1.8824018) q[2];
sx q[2];
rz(-2.5172293) q[2];
sx q[2];
rz(-1.7758689) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.98403) q[1];
sx q[1];
rz(-1.2783935) q[1];
sx q[1];
rz(-1.2127905) q[1];
x q[2];
rz(-1.7881259) q[3];
sx q[3];
rz(-1.8193075) q[3];
sx q[3];
rz(-2.0686764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.7245543) q[2];
sx q[2];
rz(-0.78196708) q[2];
sx q[2];
rz(-3.0843201) q[2];
rz(0.96632424) q[3];
sx q[3];
rz(-0.94760197) q[3];
sx q[3];
rz(1.2053325) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9178977) q[0];
sx q[0];
rz(-0.88858336) q[0];
sx q[0];
rz(0.59820557) q[0];
rz(-1.573645) q[1];
sx q[1];
rz(-1.6697829) q[1];
sx q[1];
rz(0.15984687) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.94598818) q[0];
sx q[0];
rz(-2.3538082) q[0];
sx q[0];
rz(-0.24290084) q[0];
x q[1];
rz(-0.27879981) q[2];
sx q[2];
rz(-0.21950102) q[2];
sx q[2];
rz(2.1967099) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.7084078) q[1];
sx q[1];
rz(-1.1918648) q[1];
sx q[1];
rz(-2.4772732) q[1];
rz(-pi) q[2];
rz(-1.615404) q[3];
sx q[3];
rz(-1.6470223) q[3];
sx q[3];
rz(1.3520698) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.9827031) q[2];
sx q[2];
rz(-0.44124678) q[2];
sx q[2];
rz(-2.5510447) q[2];
rz(0.96747893) q[3];
sx q[3];
rz(-1.5891821) q[3];
sx q[3];
rz(-1.1837186) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9443053) q[0];
sx q[0];
rz(-0.66216457) q[0];
sx q[0];
rz(-2.1024607) q[0];
rz(0.93152535) q[1];
sx q[1];
rz(-0.68382278) q[1];
sx q[1];
rz(2.0328124) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0174579) q[0];
sx q[0];
rz(-2.0724247) q[0];
sx q[0];
rz(2.5162016) q[0];
rz(2.336727) q[2];
sx q[2];
rz(-2.1535843) q[2];
sx q[2];
rz(2.1104167) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0045864) q[1];
sx q[1];
rz(-0.69939628) q[1];
sx q[1];
rz(-0.57199332) q[1];
rz(-pi) q[2];
x q[2];
rz(2.118947) q[3];
sx q[3];
rz(-1.7860869) q[3];
sx q[3];
rz(1.5529274) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.80591566) q[2];
sx q[2];
rz(-0.47199619) q[2];
sx q[2];
rz(0.3042039) q[2];
rz(2.975587) q[3];
sx q[3];
rz(-1.7793572) q[3];
sx q[3];
rz(-1.5267052) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-2.6840927) q[0];
sx q[0];
rz(-0.2245716) q[0];
sx q[0];
rz(-1.4373047) q[0];
rz(2.8999088) q[1];
sx q[1];
rz(-1.4983404) q[1];
sx q[1];
rz(-0.93745747) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8305377) q[0];
sx q[0];
rz(-1.8094581) q[0];
sx q[0];
rz(-2.412459) q[0];
rz(-pi) q[1];
x q[1];
rz(1.676399) q[2];
sx q[2];
rz(-2.8594976) q[2];
sx q[2];
rz(-2.8959236) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.7852206) q[1];
sx q[1];
rz(-0.8360148) q[1];
sx q[1];
rz(-2.4874175) q[1];
x q[2];
rz(0.99284692) q[3];
sx q[3];
rz(-1.9521023) q[3];
sx q[3];
rz(0.15814885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.9044056) q[2];
sx q[2];
rz(-1.097647) q[2];
sx q[2];
rz(0.22906765) q[2];
rz(1.1711052) q[3];
sx q[3];
rz(-1.3874715) q[3];
sx q[3];
rz(2.4598725) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0197297) q[0];
sx q[0];
rz(-0.64058146) q[0];
sx q[0];
rz(-1.4135452) q[0];
rz(-1.9238663) q[1];
sx q[1];
rz(-1.4378589) q[1];
sx q[1];
rz(0.53122836) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1517881) q[0];
sx q[0];
rz(-1.2453016) q[0];
sx q[0];
rz(-0.20533104) q[0];
rz(0.81905535) q[2];
sx q[2];
rz(-1.2651288) q[2];
sx q[2];
rz(-1.7421987) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.8063691) q[1];
sx q[1];
rz(-2.1861736) q[1];
sx q[1];
rz(-1.022831) q[1];
x q[2];
rz(1.2770428) q[3];
sx q[3];
rz(-1.862502) q[3];
sx q[3];
rz(-0.13495378) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.5765877) q[2];
sx q[2];
rz(-1.4385303) q[2];
sx q[2];
rz(-0.79565221) q[2];
rz(0.82120454) q[3];
sx q[3];
rz(-0.20902769) q[3];
sx q[3];
rz(0.30278277) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
sx q[3];
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
rz(-2.8177556) q[0];
sx q[0];
rz(-0.19446401) q[0];
sx q[0];
rz(1.4870462) q[0];
rz(-1.6200292) q[1];
sx q[1];
rz(-1.227102) q[1];
sx q[1];
rz(-1.0618658) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.056039) q[0];
sx q[0];
rz(-2.0889707) q[0];
sx q[0];
rz(-2.0033625) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7369607) q[2];
sx q[2];
rz(-2.8434357) q[2];
sx q[2];
rz(2.1782137) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7567021) q[1];
sx q[1];
rz(-1.9918009) q[1];
sx q[1];
rz(-2.9744991) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9412334) q[3];
sx q[3];
rz(-1.9780156) q[3];
sx q[3];
rz(-0.20663156) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.1294641) q[2];
sx q[2];
rz(-1.1295372) q[2];
sx q[2];
rz(-0.10936347) q[2];
rz(-0.28956413) q[3];
sx q[3];
rz(-1.3907631) q[3];
sx q[3];
rz(0.65350986) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43082032) q[0];
sx q[0];
rz(-1.9662974) q[0];
sx q[0];
rz(1.9868504) q[0];
rz(0.42354241) q[1];
sx q[1];
rz(-1.7083028) q[1];
sx q[1];
rz(2.4636041) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7586939) q[0];
sx q[0];
rz(-0.80307559) q[0];
sx q[0];
rz(-2.4944958) q[0];
rz(0.26099796) q[2];
sx q[2];
rz(-1.3512865) q[2];
sx q[2];
rz(-1.0936119) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.0748234) q[1];
sx q[1];
rz(-0.69157234) q[1];
sx q[1];
rz(2.2645386) q[1];
rz(0.20991169) q[3];
sx q[3];
rz(-0.60214199) q[3];
sx q[3];
rz(-1.9937552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.28422156) q[2];
sx q[2];
rz(-2.325433) q[2];
sx q[2];
rz(0.59874272) q[2];
rz(1.6019999) q[3];
sx q[3];
rz(-1.9339823) q[3];
sx q[3];
rz(-1.0369176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67958647) q[0];
sx q[0];
rz(-0.34272598) q[0];
sx q[0];
rz(1.7525679) q[0];
rz(-0.012808375) q[1];
sx q[1];
rz(-0.3777596) q[1];
sx q[1];
rz(-1.6526745) q[1];
rz(0.92973443) q[2];
sx q[2];
rz(-1.4818746) q[2];
sx q[2];
rz(-1.6043639) q[2];
rz(1.5690702) q[3];
sx q[3];
rz(-0.55894096) q[3];
sx q[3];
rz(-0.10645549) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
