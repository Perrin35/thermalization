OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.5443213) q[0];
sx q[0];
rz(-2.7379524) q[0];
sx q[0];
rz(0.37024745) q[0];
rz(-2.9397842) q[1];
sx q[1];
rz(-1.8887853) q[1];
sx q[1];
rz(-1.7226146) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1365294) q[0];
sx q[0];
rz(-2.8423474) q[0];
sx q[0];
rz(-2.5314999) q[0];
rz(2.6909157) q[2];
sx q[2];
rz(-2.5833231) q[2];
sx q[2];
rz(-2.6435341) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.287392) q[1];
sx q[1];
rz(-0.55254793) q[1];
sx q[1];
rz(-1.6709177) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2526413) q[3];
sx q[3];
rz(-0.83359026) q[3];
sx q[3];
rz(-2.6300501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.7449164) q[2];
sx q[2];
rz(-1.6884721) q[2];
sx q[2];
rz(-2.7837226) q[2];
rz(-0.19168028) q[3];
sx q[3];
rz(-0.43281698) q[3];
sx q[3];
rz(2.5884957) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1382004) q[0];
sx q[0];
rz(-1.1509742) q[0];
sx q[0];
rz(-0.068280846) q[0];
rz(1.0066907) q[1];
sx q[1];
rz(-0.12193646) q[1];
sx q[1];
rz(2.4904747) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2383645) q[0];
sx q[0];
rz(-0.67175409) q[0];
sx q[0];
rz(2.8147459) q[0];
rz(-pi) q[1];
rz(-1.5468555) q[2];
sx q[2];
rz(-2.0964453) q[2];
sx q[2];
rz(-2.0999694) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.5755641) q[1];
sx q[1];
rz(-2.063942) q[1];
sx q[1];
rz(-0.50599392) q[1];
x q[2];
rz(-0.39204709) q[3];
sx q[3];
rz(-0.52467504) q[3];
sx q[3];
rz(-0.80313659) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6887001) q[2];
sx q[2];
rz(-2.9850027) q[2];
sx q[2];
rz(-2.2382656) q[2];
rz(-0.7545169) q[3];
sx q[3];
rz(-1.4947596) q[3];
sx q[3];
rz(2.8675458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84905255) q[0];
sx q[0];
rz(-0.64340574) q[0];
sx q[0];
rz(-2.6480411) q[0];
rz(1.4986528) q[1];
sx q[1];
rz(-0.40619266) q[1];
sx q[1];
rz(2.1123871) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4714067) q[0];
sx q[0];
rz(-2.0051415) q[0];
sx q[0];
rz(1.9451408) q[0];
rz(2.2649293) q[2];
sx q[2];
rz(-0.52774094) q[2];
sx q[2];
rz(-1.9805679) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.88941105) q[1];
sx q[1];
rz(-1.2671789) q[1];
sx q[1];
rz(1.3510515) q[1];
x q[2];
rz(3.0164099) q[3];
sx q[3];
rz(-1.0297965) q[3];
sx q[3];
rz(-2.560905) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.9323953) q[2];
sx q[2];
rz(-2.5019427) q[2];
sx q[2];
rz(-1.8445245) q[2];
rz(2.5629937) q[3];
sx q[3];
rz(-1.9208627) q[3];
sx q[3];
rz(-0.67563081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9349174) q[0];
sx q[0];
rz(-0.6466372) q[0];
sx q[0];
rz(-2.7329965) q[0];
rz(-1.714255) q[1];
sx q[1];
rz(-1.4988377) q[1];
sx q[1];
rz(2.2600007) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0094385) q[0];
sx q[0];
rz(-0.75320019) q[0];
sx q[0];
rz(0.25235812) q[0];
rz(-pi) q[1];
rz(2.4741715) q[2];
sx q[2];
rz(-2.1758658) q[2];
sx q[2];
rz(0.24250008) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.91729212) q[1];
sx q[1];
rz(-0.91887337) q[1];
sx q[1];
rz(-0.93445458) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.41739695) q[3];
sx q[3];
rz(-1.5098803) q[3];
sx q[3];
rz(0.062373769) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.0980229) q[2];
sx q[2];
rz(-1.6514401) q[2];
sx q[2];
rz(2.0289452) q[2];
rz(-2.5860795) q[3];
sx q[3];
rz(-1.7881309) q[3];
sx q[3];
rz(1.6120733) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9920138) q[0];
sx q[0];
rz(-1.1554138) q[0];
sx q[0];
rz(-1.7768815) q[0];
rz(-2.8240906) q[1];
sx q[1];
rz(-0.96173871) q[1];
sx q[1];
rz(-3.024335) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8506354) q[0];
sx q[0];
rz(-0.98941776) q[0];
sx q[0];
rz(1.4145538) q[0];
rz(-pi) q[1];
x q[1];
rz(0.60301493) q[2];
sx q[2];
rz(-1.7479959) q[2];
sx q[2];
rz(3.035383) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.458141) q[1];
sx q[1];
rz(-1.7675753) q[1];
sx q[1];
rz(-1.7276006) q[1];
x q[2];
rz(-1.6860784) q[3];
sx q[3];
rz(-0.59024631) q[3];
sx q[3];
rz(0.78199996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.34510288) q[2];
sx q[2];
rz(-1.3619276) q[2];
sx q[2];
rz(-1.7618746) q[2];
rz(-1.1896677) q[3];
sx q[3];
rz(-2.9822571) q[3];
sx q[3];
rz(-3.0670847) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18774408) q[0];
sx q[0];
rz(-1.501361) q[0];
sx q[0];
rz(-0.70621079) q[0];
rz(2.0300991) q[1];
sx q[1];
rz(-0.62875426) q[1];
sx q[1];
rz(0.10791735) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88718016) q[0];
sx q[0];
rz(-2.2590056) q[0];
sx q[0];
rz(1.956091) q[0];
x q[1];
rz(0.040060476) q[2];
sx q[2];
rz(-1.5821579) q[2];
sx q[2];
rz(-0.60919112) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.2006827) q[1];
sx q[1];
rz(-0.74146491) q[1];
sx q[1];
rz(-0.42592589) q[1];
x q[2];
rz(2.1351486) q[3];
sx q[3];
rz(-0.77277771) q[3];
sx q[3];
rz(0.3533065) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.32249054) q[2];
sx q[2];
rz(-2.1671447) q[2];
sx q[2];
rz(-1.1090013) q[2];
rz(-1.2935982) q[3];
sx q[3];
rz(-1.785991) q[3];
sx q[3];
rz(0.10722815) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7727707) q[0];
sx q[0];
rz(-1.4819551) q[0];
sx q[0];
rz(-0.014904508) q[0];
rz(-2.7203454) q[1];
sx q[1];
rz(-2.0887471) q[1];
sx q[1];
rz(0.79963911) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9124303) q[0];
sx q[0];
rz(-1.8254571) q[0];
sx q[0];
rz(-2.4486662) q[0];
x q[1];
rz(-2.8268379) q[2];
sx q[2];
rz(-1.8416648) q[2];
sx q[2];
rz(-2.1115007) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.0886503) q[1];
sx q[1];
rz(-2.9859516) q[1];
sx q[1];
rz(2.774653) q[1];
rz(-3.0040967) q[3];
sx q[3];
rz(-2.1318448) q[3];
sx q[3];
rz(1.4631127) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.556276) q[2];
sx q[2];
rz(-0.63627807) q[2];
sx q[2];
rz(-0.9220534) q[2];
rz(1.8317892) q[3];
sx q[3];
rz(-1.9366879) q[3];
sx q[3];
rz(0.95782763) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(-0.17469445) q[0];
sx q[0];
rz(-2.0860724) q[0];
sx q[0];
rz(-2.877537) q[0];
rz(-1.7407725) q[1];
sx q[1];
rz(-1.6903279) q[1];
sx q[1];
rz(0.31731269) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2469651) q[0];
sx q[0];
rz(-1.6435701) q[0];
sx q[0];
rz(-1.9586246) q[0];
rz(-pi) q[1];
rz(0.1605026) q[2];
sx q[2];
rz(-1.2264894) q[2];
sx q[2];
rz(0.80489327) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6416157) q[1];
sx q[1];
rz(-0.14650211) q[1];
sx q[1];
rz(-2.4661857) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4044912) q[3];
sx q[3];
rz(-2.2095223) q[3];
sx q[3];
rz(2.2685662) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.7060966) q[2];
sx q[2];
rz(-2.2075682) q[2];
sx q[2];
rz(1.6395817) q[2];
rz(-2.8912985) q[3];
sx q[3];
rz(-1.7264265) q[3];
sx q[3];
rz(-1.1423473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89649993) q[0];
sx q[0];
rz(-2.8332906) q[0];
sx q[0];
rz(1.4244351) q[0];
rz(0.4793438) q[1];
sx q[1];
rz(-1.4762676) q[1];
sx q[1];
rz(0.11553484) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.696366) q[0];
sx q[0];
rz(-1.0123024) q[0];
sx q[0];
rz(0.88748705) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4268482) q[2];
sx q[2];
rz(-1.7822767) q[2];
sx q[2];
rz(-1.5874869) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.1852473) q[1];
sx q[1];
rz(-2.6223409) q[1];
sx q[1];
rz(0.96038702) q[1];
x q[2];
rz(-1.450591) q[3];
sx q[3];
rz(-2.2634441) q[3];
sx q[3];
rz(-1.8409178) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.66403786) q[2];
sx q[2];
rz(-1.1096191) q[2];
sx q[2];
rz(1.7158562) q[2];
rz(-1.4005631) q[3];
sx q[3];
rz(-1.0481342) q[3];
sx q[3];
rz(-2.8588296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7857159) q[0];
sx q[0];
rz(-2.4079005) q[0];
sx q[0];
rz(0.26915959) q[0];
rz(-2.0458938) q[1];
sx q[1];
rz(-0.91047374) q[1];
sx q[1];
rz(1.7620618) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3740765) q[0];
sx q[0];
rz(-2.9449468) q[0];
sx q[0];
rz(1.0023414) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.079897957) q[2];
sx q[2];
rz(-1.1546087) q[2];
sx q[2];
rz(1.4518567) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.97314944) q[1];
sx q[1];
rz(-1.7440376) q[1];
sx q[1];
rz(0.030862191) q[1];
rz(1.2606603) q[3];
sx q[3];
rz(-1.2536067) q[3];
sx q[3];
rz(-0.96084259) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.049008869) q[2];
sx q[2];
rz(-2.9276586) q[2];
sx q[2];
rz(1.5157549) q[2];
rz(-1.1670636) q[3];
sx q[3];
rz(-1.5853106) q[3];
sx q[3];
rz(-1.0313755) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2172858) q[0];
sx q[0];
rz(-1.3047682) q[0];
sx q[0];
rz(-2.7346942) q[0];
rz(2.6869607) q[1];
sx q[1];
rz(-1.1063207) q[1];
sx q[1];
rz(2.8938821) q[1];
rz(1.8998014) q[2];
sx q[2];
rz(-2.6458089) q[2];
sx q[2];
rz(-2.1363346) q[2];
rz(1.7975939) q[3];
sx q[3];
rz(-2.8759225) q[3];
sx q[3];
rz(0.21951036) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
