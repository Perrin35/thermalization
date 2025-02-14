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
rz(0.91520619) q[0];
sx q[0];
rz(2.6829166) q[0];
sx q[0];
rz(12.248326) q[0];
rz(7.6492352) q[1];
sx q[1];
rz(0.69861689) q[1];
sx q[1];
rz(3.0787992) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.013115766) q[0];
sx q[0];
rz(-0.72762353) q[0];
sx q[0];
rz(-0.81650556) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7032825) q[2];
sx q[2];
rz(-0.93899124) q[2];
sx q[2];
rz(0.28011986) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5972728) q[1];
sx q[1];
rz(-1.8493127) q[1];
sx q[1];
rz(-2.0803948) q[1];
rz(-pi) q[2];
rz(-2.1413689) q[3];
sx q[3];
rz(-2.4834342) q[3];
sx q[3];
rz(0.47129813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.723145) q[2];
sx q[2];
rz(-1.3157088) q[2];
sx q[2];
rz(-2.1212228) q[2];
rz(-2.4127035) q[3];
sx q[3];
rz(-0.43292361) q[3];
sx q[3];
rz(-1.5322022) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6272524) q[0];
sx q[0];
rz(-1.1662551) q[0];
sx q[0];
rz(3.1033206) q[0];
rz(0.12380869) q[1];
sx q[1];
rz(-2.6688711) q[1];
sx q[1];
rz(1.5708057) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1403303) q[0];
sx q[0];
rz(-1.5902983) q[0];
sx q[0];
rz(-3.131585) q[0];
rz(-0.10199396) q[2];
sx q[2];
rz(-0.77384678) q[2];
sx q[2];
rz(0.71600658) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.99421895) q[1];
sx q[1];
rz(-3.1404853) q[1];
sx q[1];
rz(1.963041) q[1];
rz(-1.6917211) q[3];
sx q[3];
rz(-2.0837415) q[3];
sx q[3];
rz(-2.4642022) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.6191285) q[2];
sx q[2];
rz(-1.3294486) q[2];
sx q[2];
rz(2.5627356) q[2];
rz(2.0959496) q[3];
sx q[3];
rz(-2.3561616) q[3];
sx q[3];
rz(-2.3845909) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4724562) q[0];
sx q[0];
rz(-2.0643015) q[0];
sx q[0];
rz(-2.6876167) q[0];
rz(-0.16481915) q[1];
sx q[1];
rz(-2.3655128) q[1];
sx q[1];
rz(-1.4705315) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1810581) q[0];
sx q[0];
rz(-1.1832602) q[0];
sx q[0];
rz(-0.18523943) q[0];
x q[1];
rz(-2.1039338) q[2];
sx q[2];
rz(-0.59269822) q[2];
sx q[2];
rz(1.3190086) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.78782234) q[1];
sx q[1];
rz(-2.2332237) q[1];
sx q[1];
rz(-1.1819928) q[1];
x q[2];
rz(-2.4824247) q[3];
sx q[3];
rz(-1.9999095) q[3];
sx q[3];
rz(1.8076713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.37956023) q[2];
sx q[2];
rz(-1.6635868) q[2];
sx q[2];
rz(1.3564159) q[2];
rz(0.49946579) q[3];
sx q[3];
rz(-1.5288589) q[3];
sx q[3];
rz(-2.0904026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22719638) q[0];
sx q[0];
rz(-0.90407404) q[0];
sx q[0];
rz(-1.0585349) q[0];
rz(1.1610441) q[1];
sx q[1];
rz(-0.5608905) q[1];
sx q[1];
rz(-0.11722359) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40848428) q[0];
sx q[0];
rz(-0.54382864) q[0];
sx q[0];
rz(-0.075790719) q[0];
rz(-pi) q[1];
rz(0.19660825) q[2];
sx q[2];
rz(-0.93485445) q[2];
sx q[2];
rz(1.3351118) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.83671727) q[1];
sx q[1];
rz(-2.3099897) q[1];
sx q[1];
rz(-1.9518714) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6113308) q[3];
sx q[3];
rz(-0.70928364) q[3];
sx q[3];
rz(-1.6146631) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.78493541) q[2];
sx q[2];
rz(-1.4484118) q[2];
sx q[2];
rz(1.7735205) q[2];
rz(1.8562227) q[3];
sx q[3];
rz(-2.0338438) q[3];
sx q[3];
rz(0.72803298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0990937) q[0];
sx q[0];
rz(-0.92785257) q[0];
sx q[0];
rz(1.4121144) q[0];
rz(2.1389351) q[1];
sx q[1];
rz(-2.899677) q[1];
sx q[1];
rz(-1.5589421) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5413645) q[0];
sx q[0];
rz(-1.1523655) q[0];
sx q[0];
rz(2.5264242) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4232403) q[2];
sx q[2];
rz(-1.6334772) q[2];
sx q[2];
rz(2.2557392) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.69935) q[1];
sx q[1];
rz(-2.1210056) q[1];
sx q[1];
rz(3.1139598) q[1];
rz(-0.094314625) q[3];
sx q[3];
rz(-1.2504761) q[3];
sx q[3];
rz(-1.8220779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.5530508) q[2];
sx q[2];
rz(-1.43247) q[2];
sx q[2];
rz(-2.7480965) q[2];
rz(2.1816175) q[3];
sx q[3];
rz(-2.4656651) q[3];
sx q[3];
rz(-2.1737607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2734964) q[0];
sx q[0];
rz(-3.0143026) q[0];
sx q[0];
rz(-2.9582276) q[0];
rz(-2.6496437) q[1];
sx q[1];
rz(-2.349647) q[1];
sx q[1];
rz(1.9221745) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3122092) q[0];
sx q[0];
rz(-1.4903194) q[0];
sx q[0];
rz(-0.089437251) q[0];
rz(-0.030410847) q[2];
sx q[2];
rz(-1.6492427) q[2];
sx q[2];
rz(0.65785656) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.26039429) q[1];
sx q[1];
rz(-1.9250084) q[1];
sx q[1];
rz(-0.43448914) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0991092) q[3];
sx q[3];
rz(-1.6460643) q[3];
sx q[3];
rz(-1.4510297) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.1313974) q[2];
sx q[2];
rz(-1.7051899) q[2];
sx q[2];
rz(-3.0592226) q[2];
rz(2.2149337) q[3];
sx q[3];
rz(-2.4121273) q[3];
sx q[3];
rz(-1.5026708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6019186) q[0];
sx q[0];
rz(-0.32783666) q[0];
sx q[0];
rz(-0.71640054) q[0];
rz(-2.7377103) q[1];
sx q[1];
rz(-1.5733746) q[1];
sx q[1];
rz(-1.7049888) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88750171) q[0];
sx q[0];
rz(-1.9483101) q[0];
sx q[0];
rz(0.15477263) q[0];
rz(2.613343) q[2];
sx q[2];
rz(-1.6852323) q[2];
sx q[2];
rz(-0.50811646) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.0246504) q[1];
sx q[1];
rz(-0.41100178) q[1];
sx q[1];
rz(0.25115369) q[1];
rz(-pi) q[2];
rz(0.8248803) q[3];
sx q[3];
rz(-2.2268104) q[3];
sx q[3];
rz(0.069930596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.7776103) q[2];
sx q[2];
rz(-2.1668375) q[2];
sx q[2];
rz(-0.065356143) q[2];
rz(0.86723793) q[3];
sx q[3];
rz(-1.6403653) q[3];
sx q[3];
rz(1.8170099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48677483) q[0];
sx q[0];
rz(-1.4210533) q[0];
sx q[0];
rz(-2.4105893) q[0];
rz(2.3664318) q[1];
sx q[1];
rz(-2.8072) q[1];
sx q[1];
rz(2.7269272) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0867942) q[0];
sx q[0];
rz(-2.2761184) q[0];
sx q[0];
rz(2.2948625) q[0];
rz(2.7887949) q[2];
sx q[2];
rz(-0.55856201) q[2];
sx q[2];
rz(-2.0143904) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.2036671) q[1];
sx q[1];
rz(-1.8216743) q[1];
sx q[1];
rz(2.832042) q[1];
x q[2];
rz(1.7348792) q[3];
sx q[3];
rz(-1.7030431) q[3];
sx q[3];
rz(-1.2175206) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9097462) q[2];
sx q[2];
rz(-0.48603386) q[2];
sx q[2];
rz(-3.0492142) q[2];
rz(-2.3653024) q[3];
sx q[3];
rz(-1.4624566) q[3];
sx q[3];
rz(0.20364729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48731503) q[0];
sx q[0];
rz(-1.646811) q[0];
sx q[0];
rz(3.1238632) q[0];
rz(-2.0954633) q[1];
sx q[1];
rz(-1.4132063) q[1];
sx q[1];
rz(2.0808992) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47155505) q[0];
sx q[0];
rz(-0.070138358) q[0];
sx q[0];
rz(0.20151968) q[0];
rz(1.1938042) q[2];
sx q[2];
rz(-0.62219884) q[2];
sx q[2];
rz(0.54905984) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.6010161) q[1];
sx q[1];
rz(-0.81588826) q[1];
sx q[1];
rz(0.23758446) q[1];
x q[2];
rz(-0.10194998) q[3];
sx q[3];
rz(-0.52604874) q[3];
sx q[3];
rz(-0.40204429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.9866508) q[2];
sx q[2];
rz(-0.81601802) q[2];
sx q[2];
rz(0.6443392) q[2];
rz(-2.2996357) q[3];
sx q[3];
rz(-1.569845) q[3];
sx q[3];
rz(0.90698609) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3802721) q[0];
sx q[0];
rz(-0.43571061) q[0];
sx q[0];
rz(-2.6897588) q[0];
rz(1.0093581) q[1];
sx q[1];
rz(-1.6296856) q[1];
sx q[1];
rz(-1.5712646) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6752351) q[0];
sx q[0];
rz(-1.5769099) q[0];
sx q[0];
rz(0.17019043) q[0];
rz(-pi) q[1];
rz(-1.3469996) q[2];
sx q[2];
rz(-0.48529709) q[2];
sx q[2];
rz(1.3027625) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.81262368) q[1];
sx q[1];
rz(-1.4419893) q[1];
sx q[1];
rz(1.9616496) q[1];
x q[2];
rz(-1.5630003) q[3];
sx q[3];
rz(-0.8484133) q[3];
sx q[3];
rz(-2.9786199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.62341225) q[2];
sx q[2];
rz(-1.0717012) q[2];
sx q[2];
rz(-1.5706459) q[2];
rz(3.0359641) q[3];
sx q[3];
rz(-0.94625866) q[3];
sx q[3];
rz(-2.6251729) q[3];
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
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8138206) q[0];
sx q[0];
rz(-2.1229424) q[0];
sx q[0];
rz(0.13737296) q[0];
rz(2.9227921) q[1];
sx q[1];
rz(-1.7411502) q[1];
sx q[1];
rz(2.2314744) q[1];
rz(-1.1458746) q[2];
sx q[2];
rz(-2.6617202) q[2];
sx q[2];
rz(-0.94750994) q[2];
rz(-1.2050592) q[3];
sx q[3];
rz(-0.94612056) q[3];
sx q[3];
rz(-2.003086) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
