OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.3333617) q[0];
sx q[0];
rz(-0.97167492) q[0];
sx q[0];
rz(-1.4810286) q[0];
rz(-0.12227585) q[1];
sx q[1];
rz(3.0552157) q[1];
sx q[1];
rz(9.4019158) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.011933) q[0];
sx q[0];
rz(-1.7190022) q[0];
sx q[0];
rz(0.11797842) q[0];
rz(-pi) q[1];
rz(-0.51486751) q[2];
sx q[2];
rz(-1.1866743) q[2];
sx q[2];
rz(1.1615953) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.80427985) q[1];
sx q[1];
rz(-2.7645281) q[1];
sx q[1];
rz(1.9878597) q[1];
x q[2];
rz(-0.25306563) q[3];
sx q[3];
rz(-2.1742637) q[3];
sx q[3];
rz(0.30482182) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.7227398) q[2];
sx q[2];
rz(-2.3143694) q[2];
sx q[2];
rz(-2.375405) q[2];
rz(-2.9700759) q[3];
sx q[3];
rz(-2.4092716) q[3];
sx q[3];
rz(0.58656251) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3294753) q[0];
sx q[0];
rz(-2.7024039) q[0];
sx q[0];
rz(-3.0549333) q[0];
rz(-0.86241972) q[1];
sx q[1];
rz(-2.5984867) q[1];
sx q[1];
rz(-0.08509732) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21762411) q[0];
sx q[0];
rz(-2.1930165) q[0];
sx q[0];
rz(-1.5505962) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1727824) q[2];
sx q[2];
rz(-1.3170871) q[2];
sx q[2];
rz(-1.3041376) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.3955235) q[1];
sx q[1];
rz(-1.6568526) q[1];
sx q[1];
rz(1.6834016) q[1];
rz(2.8220903) q[3];
sx q[3];
rz(-1.1211683) q[3];
sx q[3];
rz(3.086123) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.26248419) q[2];
sx q[2];
rz(-1.5905249) q[2];
sx q[2];
rz(1.4651728) q[2];
rz(1.1761459) q[3];
sx q[3];
rz(-0.90548235) q[3];
sx q[3];
rz(2.8951077) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6572606) q[0];
sx q[0];
rz(-0.13467877) q[0];
sx q[0];
rz(2.235967) q[0];
rz(2.8088645) q[1];
sx q[1];
rz(-0.85488027) q[1];
sx q[1];
rz(-1.7571626) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8850088) q[0];
sx q[0];
rz(-0.48261595) q[0];
sx q[0];
rz(0.17871876) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.92817523) q[2];
sx q[2];
rz(-0.60545063) q[2];
sx q[2];
rz(1.6124992) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.82921925) q[1];
sx q[1];
rz(-0.81554283) q[1];
sx q[1];
rz(1.5401328) q[1];
rz(-pi) q[2];
rz(1.0043886) q[3];
sx q[3];
rz(-1.1096138) q[3];
sx q[3];
rz(-0.084463488) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.8973792) q[2];
sx q[2];
rz(-0.34003568) q[2];
sx q[2];
rz(1.3282233) q[2];
rz(-1.7437079) q[3];
sx q[3];
rz(-0.54026794) q[3];
sx q[3];
rz(3.0298997) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0760647) q[0];
sx q[0];
rz(-2.9632443) q[0];
sx q[0];
rz(-0.67101014) q[0];
rz(1.3440075) q[1];
sx q[1];
rz(-1.1616511) q[1];
sx q[1];
rz(2.9342594) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3509388) q[0];
sx q[0];
rz(-1.6121943) q[0];
sx q[0];
rz(-0.008965094) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.030427293) q[2];
sx q[2];
rz(-1.7446339) q[2];
sx q[2];
rz(2.0174842) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.3934717) q[1];
sx q[1];
rz(-1.3337743) q[1];
sx q[1];
rz(2.9841828) q[1];
rz(-pi) q[2];
x q[2];
rz(3.1158833) q[3];
sx q[3];
rz(-2.6298012) q[3];
sx q[3];
rz(0.6230841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.8308782) q[2];
sx q[2];
rz(-1.0881311) q[2];
sx q[2];
rz(-2.3332398) q[2];
rz(0.80777848) q[3];
sx q[3];
rz(-0.72179663) q[3];
sx q[3];
rz(0.19367735) q[3];
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
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5508674) q[0];
sx q[0];
rz(-0.71722513) q[0];
sx q[0];
rz(2.2461058) q[0];
rz(0.26516178) q[1];
sx q[1];
rz(-2.3064955) q[1];
sx q[1];
rz(-2.6079544) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95414872) q[0];
sx q[0];
rz(-1.8549518) q[0];
sx q[0];
rz(0.28676333) q[0];
rz(-pi) q[1];
rz(0.76782121) q[2];
sx q[2];
rz(-1.4066833) q[2];
sx q[2];
rz(-2.5196911) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.67709778) q[1];
sx q[1];
rz(-0.2650731) q[1];
sx q[1];
rz(-2.5527843) q[1];
rz(-pi) q[2];
rz(1.2920462) q[3];
sx q[3];
rz(-2.4933443) q[3];
sx q[3];
rz(-1.3384782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.9020033) q[2];
sx q[2];
rz(-1.5782372) q[2];
sx q[2];
rz(2.5732102) q[2];
rz(1.9249632) q[3];
sx q[3];
rz(-0.47074461) q[3];
sx q[3];
rz(2.7980428) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
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
rz(2.6333273) q[0];
sx q[0];
rz(-0.87843043) q[0];
sx q[0];
rz(-1.9653962) q[0];
rz(1.5559224) q[1];
sx q[1];
rz(-0.95247477) q[1];
sx q[1];
rz(-2.1369381) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83197901) q[0];
sx q[0];
rz(-1.2011315) q[0];
sx q[0];
rz(1.2020822) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4899646) q[2];
sx q[2];
rz(-1.1751375) q[2];
sx q[2];
rz(-2.253502) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(3.1130044) q[1];
sx q[1];
rz(-0.60739809) q[1];
sx q[1];
rz(-2.7007156) q[1];
rz(-pi) q[2];
rz(-1.9846356) q[3];
sx q[3];
rz(-1.5817747) q[3];
sx q[3];
rz(-1.4423086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.4040318) q[2];
sx q[2];
rz(-2.1514251) q[2];
sx q[2];
rz(-0.44815865) q[2];
rz(-2.8435977) q[3];
sx q[3];
rz(-0.69152504) q[3];
sx q[3];
rz(-0.34860778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(0.77666831) q[0];
sx q[0];
rz(-1.8402599) q[0];
sx q[0];
rz(-0.65814322) q[0];
rz(1.2843885) q[1];
sx q[1];
rz(-2.6874266) q[1];
sx q[1];
rz(-2.4694494) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2015398) q[0];
sx q[0];
rz(-1.7367898) q[0];
sx q[0];
rz(0.0075017651) q[0];
x q[1];
rz(-1.1291814) q[2];
sx q[2];
rz(-1.1814983) q[2];
sx q[2];
rz(-2.0814975) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0969095) q[1];
sx q[1];
rz(-1.9367366) q[1];
sx q[1];
rz(-2.6198322) q[1];
rz(-pi) q[2];
rz(1.4268731) q[3];
sx q[3];
rz(-0.87464273) q[3];
sx q[3];
rz(-1.5211943) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5333574) q[2];
sx q[2];
rz(-2.50864) q[2];
sx q[2];
rz(-2.7427924) q[2];
rz(0.82459015) q[3];
sx q[3];
rz(-1.7493068) q[3];
sx q[3];
rz(0.59857541) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44052112) q[0];
sx q[0];
rz(-2.7109787) q[0];
sx q[0];
rz(-2.9833802) q[0];
rz(2.054706) q[1];
sx q[1];
rz(-2.465076) q[1];
sx q[1];
rz(-2.3349082) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54803941) q[0];
sx q[0];
rz(-2.002945) q[0];
sx q[0];
rz(-0.12950626) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9358747) q[2];
sx q[2];
rz(-1.6763902) q[2];
sx q[2];
rz(-0.92668698) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.403277) q[1];
sx q[1];
rz(-2.3936831) q[1];
sx q[1];
rz(2.7110389) q[1];
rz(-pi) q[2];
rz(1.673647) q[3];
sx q[3];
rz(-0.30234435) q[3];
sx q[3];
rz(-0.24149382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.57758254) q[2];
sx q[2];
rz(-1.6489886) q[2];
sx q[2];
rz(0.99964833) q[2];
rz(3.0380761) q[3];
sx q[3];
rz(-3.0045356) q[3];
sx q[3];
rz(-1.1243533) q[3];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0209811) q[0];
sx q[0];
rz(-2.3636901) q[0];
sx q[0];
rz(2.4840684) q[0];
rz(-2.855037) q[1];
sx q[1];
rz(-2.2189238) q[1];
sx q[1];
rz(0.34067571) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4083456) q[0];
sx q[0];
rz(-0.78662965) q[0];
sx q[0];
rz(2.0287081) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9129487) q[2];
sx q[2];
rz(-0.59235448) q[2];
sx q[2];
rz(2.5387788) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.68742467) q[1];
sx q[1];
rz(-1.7523645) q[1];
sx q[1];
rz(-1.0086235) q[1];
x q[2];
rz(0.70298985) q[3];
sx q[3];
rz(-0.25745108) q[3];
sx q[3];
rz(1.5233745) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.75491536) q[2];
sx q[2];
rz(-1.9415386) q[2];
sx q[2];
rz(2.6436451) q[2];
rz(-0.30161101) q[3];
sx q[3];
rz(-2.7409654) q[3];
sx q[3];
rz(2.6224459) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72934812) q[0];
sx q[0];
rz(-0.059818581) q[0];
sx q[0];
rz(-0.27858946) q[0];
rz(2.5623698) q[1];
sx q[1];
rz(-0.93943739) q[1];
sx q[1];
rz(0.07671193) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4519276) q[0];
sx q[0];
rz(-1.2566914) q[0];
sx q[0];
rz(-3.0118224) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9701091) q[2];
sx q[2];
rz(-2.137261) q[2];
sx q[2];
rz(0.4085853) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.0753239) q[1];
sx q[1];
rz(-1.9362209) q[1];
sx q[1];
rz(-1.5931904) q[1];
rz(0.46080188) q[3];
sx q[3];
rz(-0.19690234) q[3];
sx q[3];
rz(2.868696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.3672093) q[2];
sx q[2];
rz(-0.70720208) q[2];
sx q[2];
rz(-2.6860766) q[2];
rz(0.43141836) q[3];
sx q[3];
rz(-0.12532561) q[3];
sx q[3];
rz(3.07807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1800304) q[0];
sx q[0];
rz(-1.4008235) q[0];
sx q[0];
rz(0.93749198) q[0];
rz(-2.5554399) q[1];
sx q[1];
rz(-1.502232) q[1];
sx q[1];
rz(-1.4486817) q[1];
rz(-1.0473245) q[2];
sx q[2];
rz(-2.6101255) q[2];
sx q[2];
rz(2.5892467) q[2];
rz(2.7145731) q[3];
sx q[3];
rz(-1.2357161) q[3];
sx q[3];
rz(-1.2824035) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
