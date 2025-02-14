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
rz(-1.6725809) q[0];
sx q[0];
rz(-2.2218158) q[0];
sx q[0];
rz(-0.65499175) q[0];
rz(-1.4870149) q[1];
sx q[1];
rz(-1.4227285) q[1];
sx q[1];
rz(1.3230327) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.060299035) q[0];
sx q[0];
rz(-0.64612389) q[0];
sx q[0];
rz(0.52559121) q[0];
rz(0.31794117) q[2];
sx q[2];
rz(-1.7934718) q[2];
sx q[2];
rz(-0.26938619) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.2046584) q[1];
sx q[1];
rz(-1.2950674) q[1];
sx q[1];
rz(3.0380121) q[1];
rz(0.72812702) q[3];
sx q[3];
rz(-1.1525197) q[3];
sx q[3];
rz(0.74467105) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.677864) q[2];
sx q[2];
rz(-0.68631309) q[2];
sx q[2];
rz(0.65832552) q[2];
rz(2.5074734) q[3];
sx q[3];
rz(-1.442797) q[3];
sx q[3];
rz(1.770795) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.051801) q[0];
sx q[0];
rz(-2.8103516) q[0];
sx q[0];
rz(0.077089699) q[0];
rz(-2.5568621) q[1];
sx q[1];
rz(-1.3405776) q[1];
sx q[1];
rz(-1.1999406) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3886007) q[0];
sx q[0];
rz(-1.7677099) q[0];
sx q[0];
rz(-2.0373175) q[0];
x q[1];
rz(-1.4064155) q[2];
sx q[2];
rz(-1.5417678) q[2];
sx q[2];
rz(2.8433702) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.6226095) q[1];
sx q[1];
rz(-0.72777343) q[1];
sx q[1];
rz(0.98811291) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6317815) q[3];
sx q[3];
rz(-1.012523) q[3];
sx q[3];
rz(1.3577485) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.8889019) q[2];
sx q[2];
rz(-0.98068792) q[2];
sx q[2];
rz(0.67548951) q[2];
rz(0.0096970079) q[3];
sx q[3];
rz(-2.8986425) q[3];
sx q[3];
rz(-0.01827904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6460687) q[0];
sx q[0];
rz(-1.6697474) q[0];
sx q[0];
rz(-0.75463265) q[0];
rz(-3.1309639) q[1];
sx q[1];
rz(-1.3898712) q[1];
sx q[1];
rz(-1.0108112) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0389688) q[0];
sx q[0];
rz(-1.083923) q[0];
sx q[0];
rz(-2.7953933) q[0];
x q[1];
rz(1.6057683) q[2];
sx q[2];
rz(-2.4596697) q[2];
sx q[2];
rz(0.65820314) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.6200808) q[1];
sx q[1];
rz(-0.82640582) q[1];
sx q[1];
rz(-1.6537731) q[1];
rz(-pi) q[2];
rz(0.76356865) q[3];
sx q[3];
rz(-2.0393622) q[3];
sx q[3];
rz(2.6811622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.2406771) q[2];
sx q[2];
rz(-2.7132576) q[2];
sx q[2];
rz(0.50673103) q[2];
rz(-1.2725376) q[3];
sx q[3];
rz(-1.7839909) q[3];
sx q[3];
rz(-0.27819628) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5946567) q[0];
sx q[0];
rz(-1.2242067) q[0];
sx q[0];
rz(2.0261672) q[0];
rz(-1.2660654) q[1];
sx q[1];
rz(-1.5645809) q[1];
sx q[1];
rz(2.26684) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47622555) q[0];
sx q[0];
rz(-1.3194808) q[0];
sx q[0];
rz(3.0500924) q[0];
x q[1];
rz(-0.91712894) q[2];
sx q[2];
rz(-1.5857595) q[2];
sx q[2];
rz(-1.5749168) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.529244) q[1];
sx q[1];
rz(-1.1643641) q[1];
sx q[1];
rz(-1.2807349) q[1];
x q[2];
rz(2.277209) q[3];
sx q[3];
rz(-1.4339851) q[3];
sx q[3];
rz(-0.19696008) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.11883417) q[2];
sx q[2];
rz(-0.47670445) q[2];
sx q[2];
rz(1.9258707) q[2];
rz(-1.7440965) q[3];
sx q[3];
rz(-1.6179061) q[3];
sx q[3];
rz(1.7106748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27595156) q[0];
sx q[0];
rz(-3.0396099) q[0];
sx q[0];
rz(-1.4707461) q[0];
rz(1.3784846) q[1];
sx q[1];
rz(-2.3536317) q[1];
sx q[1];
rz(-1.7313622) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84274693) q[0];
sx q[0];
rz(-2.0618467) q[0];
sx q[0];
rz(0.23570717) q[0];
rz(-pi) q[1];
x q[1];
rz(1.330014) q[2];
sx q[2];
rz(-2.2301939) q[2];
sx q[2];
rz(1.8383775) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.8157364) q[1];
sx q[1];
rz(-1.5724025) q[1];
sx q[1];
rz(-2.8291507) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.89680421) q[3];
sx q[3];
rz(-2.6072579) q[3];
sx q[3];
rz(0.33269445) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.1269647) q[2];
sx q[2];
rz(-2.2829599) q[2];
sx q[2];
rz(-1.8178168) q[2];
rz(1.3592367) q[3];
sx q[3];
rz(-1.2991354) q[3];
sx q[3];
rz(-2.7544379) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1112261) q[0];
sx q[0];
rz(-1.68196) q[0];
sx q[0];
rz(-3.1193745) q[0];
rz(-2.4413595) q[1];
sx q[1];
rz(-1.2513688) q[1];
sx q[1];
rz(0.95692316) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8533405) q[0];
sx q[0];
rz(-1.7159675) q[0];
sx q[0];
rz(2.4397544) q[0];
x q[1];
rz(1.7526723) q[2];
sx q[2];
rz(-1.3564912) q[2];
sx q[2];
rz(-0.3265115) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.14084223) q[1];
sx q[1];
rz(-0.4943119) q[1];
sx q[1];
rz(-2.1782081) q[1];
rz(-2.2150008) q[3];
sx q[3];
rz(-0.53172382) q[3];
sx q[3];
rz(0.34429541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.9243098) q[2];
sx q[2];
rz(-2.0192912) q[2];
sx q[2];
rz(1.699532) q[2];
rz(-1.1066655) q[3];
sx q[3];
rz(-0.60551867) q[3];
sx q[3];
rz(-1.1520011) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42665136) q[0];
sx q[0];
rz(-1.8639257) q[0];
sx q[0];
rz(-2.6746124) q[0];
rz(-1.3106208) q[1];
sx q[1];
rz(-0.95342031) q[1];
sx q[1];
rz(2.7511645) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.098179558) q[0];
sx q[0];
rz(-1.5382086) q[0];
sx q[0];
rz(0.14049732) q[0];
x q[1];
rz(-1.9127185) q[2];
sx q[2];
rz(-0.31154267) q[2];
sx q[2];
rz(2.3780413) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.5891287) q[1];
sx q[1];
rz(-0.93163449) q[1];
sx q[1];
rz(0.57560779) q[1];
rz(-2.1834063) q[3];
sx q[3];
rz(-1.4496007) q[3];
sx q[3];
rz(-2.6874128) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.9499669) q[2];
sx q[2];
rz(-1.2167296) q[2];
sx q[2];
rz(2.1245655) q[2];
rz(0.93938604) q[3];
sx q[3];
rz(-2.2685969) q[3];
sx q[3];
rz(2.2123607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4262714) q[0];
sx q[0];
rz(-2.4763698) q[0];
sx q[0];
rz(1.2128879) q[0];
rz(2.5384278) q[1];
sx q[1];
rz(-1.1319356) q[1];
sx q[1];
rz(-2.718198) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0608262) q[0];
sx q[0];
rz(-2.2696324) q[0];
sx q[0];
rz(-0.10159512) q[0];
rz(-0.076831623) q[2];
sx q[2];
rz(-1.4203826) q[2];
sx q[2];
rz(-0.39932775) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.5880206) q[1];
sx q[1];
rz(-1.3979853) q[1];
sx q[1];
rz(-0.2134448) q[1];
rz(-3.0984108) q[3];
sx q[3];
rz(-0.16999741) q[3];
sx q[3];
rz(1.377533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6508871) q[2];
sx q[2];
rz(-0.31279534) q[2];
sx q[2];
rz(-2.1631964) q[2];
rz(-2.4466416) q[3];
sx q[3];
rz(-1.2000822) q[3];
sx q[3];
rz(1.8470701) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-3.0438743) q[0];
sx q[0];
rz(-2.9053423) q[0];
sx q[0];
rz(-0.043721113) q[0];
rz(-1.1737191) q[1];
sx q[1];
rz(-0.81467384) q[1];
sx q[1];
rz(2.3497605) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39003644) q[0];
sx q[0];
rz(-1.1564213) q[0];
sx q[0];
rz(0.66585559) q[0];
rz(2.1570683) q[2];
sx q[2];
rz(-2.7862264) q[2];
sx q[2];
rz(-0.1600858) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.14230141) q[1];
sx q[1];
rz(-0.55146927) q[1];
sx q[1];
rz(-2.2427008) q[1];
x q[2];
rz(-0.47747647) q[3];
sx q[3];
rz(-1.2315742) q[3];
sx q[3];
rz(-0.71996688) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.23058471) q[2];
sx q[2];
rz(-0.30969301) q[2];
sx q[2];
rz(1.2370375) q[2];
rz(0.48219314) q[3];
sx q[3];
rz(-2.2174406) q[3];
sx q[3];
rz(-2.6528416) q[3];
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
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1173387) q[0];
sx q[0];
rz(-0.92702213) q[0];
sx q[0];
rz(-1.8035969) q[0];
rz(0.85211873) q[1];
sx q[1];
rz(-1.9133277) q[1];
sx q[1];
rz(1.1504014) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43088461) q[0];
sx q[0];
rz(-1.0477433) q[0];
sx q[0];
rz(0.56351785) q[0];
rz(-pi) q[1];
x q[1];
rz(3.030376) q[2];
sx q[2];
rz(-1.2555946) q[2];
sx q[2];
rz(-2.1325796) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.51440367) q[1];
sx q[1];
rz(-2.7233988) q[1];
sx q[1];
rz(-1.2209284) q[1];
rz(-2.7546553) q[3];
sx q[3];
rz(-2.134857) q[3];
sx q[3];
rz(1.8534941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.368025) q[2];
sx q[2];
rz(-1.0057534) q[2];
sx q[2];
rz(0.31039882) q[2];
rz(-0.6271022) q[3];
sx q[3];
rz(-2.1576364) q[3];
sx q[3];
rz(1.9703126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2720168) q[0];
sx q[0];
rz(-1.5675114) q[0];
sx q[0];
rz(-1.5935224) q[0];
rz(-1.2177474) q[1];
sx q[1];
rz(-1.4532614) q[1];
sx q[1];
rz(-0.90167602) q[1];
rz(0.2097585) q[2];
sx q[2];
rz(-2.0611613) q[2];
sx q[2];
rz(0.47917889) q[2];
rz(1.6321833) q[3];
sx q[3];
rz(-1.8922378) q[3];
sx q[3];
rz(0.031382244) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
