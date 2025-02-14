OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.6731113) q[0];
sx q[0];
rz(-3.0335479) q[0];
sx q[0];
rz(-0.98969069) q[0];
rz(-1.5762848) q[1];
sx q[1];
rz(-1.58374) q[1];
sx q[1];
rz(-1.6852112) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7262048) q[0];
sx q[0];
rz(-1.9378462) q[0];
sx q[0];
rz(-2.996654) q[0];
rz(1.5338495) q[2];
sx q[2];
rz(-1.4375028) q[2];
sx q[2];
rz(-2.5306866) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.8238358) q[1];
sx q[1];
rz(-1.8537304) q[1];
sx q[1];
rz(2.1499277) q[1];
x q[2];
rz(1.5037886) q[3];
sx q[3];
rz(-2.5634804) q[3];
sx q[3];
rz(1.029795) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.1504537) q[2];
sx q[2];
rz(-3.1296215) q[2];
sx q[2];
rz(-1.0452622) q[2];
rz(-2.1446877) q[3];
sx q[3];
rz(-0.0054587047) q[3];
sx q[3];
rz(-1.3158984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.588722) q[0];
sx q[0];
rz(-1.2405688) q[0];
sx q[0];
rz(1.3528104) q[0];
rz(-3.1006587) q[1];
sx q[1];
rz(-1.9238238) q[1];
sx q[1];
rz(1.5997684) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7517325) q[0];
sx q[0];
rz(-0.19987389) q[0];
sx q[0];
rz(0.78011192) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.1253042) q[2];
sx q[2];
rz(-1.5495346) q[2];
sx q[2];
rz(2.4580815) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.92506986) q[1];
sx q[1];
rz(-2.1669186) q[1];
sx q[1];
rz(0.03742569) q[1];
rz(-pi) q[2];
rz(-2.1435817) q[3];
sx q[3];
rz(-1.9011455) q[3];
sx q[3];
rz(2.9444061) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.6359977) q[2];
sx q[2];
rz(-3.1090241) q[2];
sx q[2];
rz(0.51844281) q[2];
rz(2.017766) q[3];
sx q[3];
rz(-0.80945194) q[3];
sx q[3];
rz(-0.43784416) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4951303) q[0];
sx q[0];
rz(-3.0914682) q[0];
sx q[0];
rz(-1.5396402) q[0];
rz(-2.4165972) q[1];
sx q[1];
rz(-0.031818964) q[1];
sx q[1];
rz(2.4620788) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2452263) q[0];
sx q[0];
rz(-0.99649444) q[0];
sx q[0];
rz(-2.0942874) q[0];
x q[1];
rz(1.5660921) q[2];
sx q[2];
rz(-1.3563547) q[2];
sx q[2];
rz(-3.0115779) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.0371548) q[1];
sx q[1];
rz(-1.3111115) q[1];
sx q[1];
rz(1.2472769) q[1];
rz(-pi) q[2];
x q[2];
rz(0.98285003) q[3];
sx q[3];
rz(-1.790739) q[3];
sx q[3];
rz(1.4877121) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.2346377) q[2];
sx q[2];
rz(-0.24181557) q[2];
sx q[2];
rz(0.50092906) q[2];
rz(2.6856954) q[3];
sx q[3];
rz(-0.026345043) q[3];
sx q[3];
rz(2.946089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2723715) q[0];
sx q[0];
rz(-0.048373241) q[0];
sx q[0];
rz(-0.79917556) q[0];
rz(2.8436106) q[1];
sx q[1];
rz(-0.25845343) q[1];
sx q[1];
rz(2.2395649) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.395754) q[0];
sx q[0];
rz(-2.3180048) q[0];
sx q[0];
rz(1.3831148) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.600014) q[2];
sx q[2];
rz(-1.4409833) q[2];
sx q[2];
rz(2.1101348) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.23556701) q[1];
sx q[1];
rz(-0.14793554) q[1];
sx q[1];
rz(-1.6355913) q[1];
rz(3.0864363) q[3];
sx q[3];
rz(-1.6481019) q[3];
sx q[3];
rz(-2.5173924) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.67348376) q[2];
sx q[2];
rz(-3.1099042) q[2];
sx q[2];
rz(1.4704977) q[2];
rz(2.9521613) q[3];
sx q[3];
rz(-0.048308689) q[3];
sx q[3];
rz(-0.76907492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22537941) q[0];
sx q[0];
rz(-0.12752859) q[0];
sx q[0];
rz(2.7498229) q[0];
rz(-2.0562992) q[1];
sx q[1];
rz(-3.1317874) q[1];
sx q[1];
rz(2.772803) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.987326) q[0];
sx q[0];
rz(-1.0079103) q[0];
sx q[0];
rz(1.994654) q[0];
rz(-pi) q[1];
rz(-2.872345) q[2];
sx q[2];
rz(-1.1675861) q[2];
sx q[2];
rz(2.0079812) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.60978598) q[1];
sx q[1];
rz(-1.577425) q[1];
sx q[1];
rz(3.1408491) q[1];
rz(1.3372793) q[3];
sx q[3];
rz(-0.87461014) q[3];
sx q[3];
rz(2.1748469) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.5607249) q[2];
sx q[2];
rz(-3.0256425) q[2];
sx q[2];
rz(1.3690534) q[2];
rz(-3.0154058) q[3];
sx q[3];
rz(-2.7044665) q[3];
sx q[3];
rz(-2.0807467) q[3];
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
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5996025) q[0];
sx q[0];
rz(-2.3805711) q[0];
sx q[0];
rz(-1.5916995) q[0];
rz(-2.1688993) q[1];
sx q[1];
rz(-0.21316554) q[1];
sx q[1];
rz(-2.1125643) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21287928) q[0];
sx q[0];
rz(-2.166011) q[0];
sx q[0];
rz(1.9363602) q[0];
rz(-2.6112399) q[2];
sx q[2];
rz(-2.8674922) q[2];
sx q[2];
rz(2.0567187) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.2756535) q[1];
sx q[1];
rz(-1.5703619) q[1];
sx q[1];
rz(-1.666953) q[1];
x q[2];
rz(-1.738284) q[3];
sx q[3];
rz(-1.5351908) q[3];
sx q[3];
rz(1.6990341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.782393) q[2];
sx q[2];
rz(-1.7155557) q[2];
sx q[2];
rz(-2.7101809) q[2];
rz(2.1443478) q[3];
sx q[3];
rz(-3.1044208) q[3];
sx q[3];
rz(0.55716151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7998841) q[0];
sx q[0];
rz(-0.48483098) q[0];
sx q[0];
rz(1.0579911) q[0];
rz(-0.82408389) q[1];
sx q[1];
rz(-3.1415756) q[1];
sx q[1];
rz(-2.3242059) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0213172) q[0];
sx q[0];
rz(-2.6861354) q[0];
sx q[0];
rz(1.0767471) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5777052) q[2];
sx q[2];
rz(-1.5525609) q[2];
sx q[2];
rz(1.697714) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.73002023) q[1];
sx q[1];
rz(-1.3294833) q[1];
sx q[1];
rz(3.1012721) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.75202077) q[3];
sx q[3];
rz(-1.9746426) q[3];
sx q[3];
rz(2.2026317) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.99007964) q[2];
sx q[2];
rz(-1.238287) q[2];
sx q[2];
rz(1.5129169) q[2];
rz(-0.40142909) q[3];
sx q[3];
rz(-0.028086834) q[3];
sx q[3];
rz(2.1872971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3681188) q[0];
sx q[0];
rz(-0.064726949) q[0];
sx q[0];
rz(-1.7777959) q[0];
rz(3.0996481) q[1];
sx q[1];
rz(-3.0105803) q[1];
sx q[1];
rz(1.0036453) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.100228) q[0];
sx q[0];
rz(-1.4778839) q[0];
sx q[0];
rz(0.63444699) q[0];
rz(-1.5940158) q[2];
sx q[2];
rz(-2.8041511) q[2];
sx q[2];
rz(-2.0793123) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(3.0692301) q[1];
sx q[1];
rz(-1.3389999) q[1];
sx q[1];
rz(-2.9332471) q[1];
rz(2.4859927) q[3];
sx q[3];
rz(-1.3274945) q[3];
sx q[3];
rz(-0.52025627) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.7046788) q[2];
sx q[2];
rz(-3.0824326) q[2];
sx q[2];
rz(1.3979647) q[2];
rz(0.49020234) q[3];
sx q[3];
rz(-0.043488113) q[3];
sx q[3];
rz(-1.1787666) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.040084664) q[0];
sx q[0];
rz(-3.0136216) q[0];
sx q[0];
rz(-2.9301933) q[0];
rz(1.118411) q[1];
sx q[1];
rz(-3.1299997) q[1];
sx q[1];
rz(1.5107907) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3113041) q[0];
sx q[0];
rz(-1.2823866) q[0];
sx q[0];
rz(2.8187241) q[0];
x q[1];
rz(2.0076934) q[2];
sx q[2];
rz(-2.1497576) q[2];
sx q[2];
rz(-1.7086264) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.37247) q[1];
sx q[1];
rz(-1.5536397) q[1];
sx q[1];
rz(-2.9869798) q[1];
x q[2];
rz(-2.9752067) q[3];
sx q[3];
rz(-1.659044) q[3];
sx q[3];
rz(-1.4180753) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.3596892) q[2];
sx q[2];
rz(-0.017904559) q[2];
sx q[2];
rz(-0.27931279) q[2];
rz(0.8638047) q[3];
sx q[3];
rz(-3.1371208) q[3];
sx q[3];
rz(-2.4590676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.913468) q[0];
sx q[0];
rz(-2.0371912) q[0];
sx q[0];
rz(-1.8260691) q[0];
rz(-0.77519351) q[1];
sx q[1];
rz(-2.9029791) q[1];
sx q[1];
rz(1.388789) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42487803) q[0];
sx q[0];
rz(-1.8535063) q[0];
sx q[0];
rz(-1.4155875) q[0];
x q[1];
rz(1.3451683) q[2];
sx q[2];
rz(-1.1755845) q[2];
sx q[2];
rz(1.164142) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.5701616) q[1];
sx q[1];
rz(-1.5716911) q[1];
sx q[1];
rz(-1.570163) q[1];
rz(-pi) q[2];
rz(-0.78755112) q[3];
sx q[3];
rz(-2.5803704) q[3];
sx q[3];
rz(0.40555996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.76643252) q[2];
sx q[2];
rz(-3.1258686) q[2];
sx q[2];
rz(-1.6895705) q[2];
rz(0.25894138) q[3];
sx q[3];
rz(-2.9539234) q[3];
sx q[3];
rz(1.1567206) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5158841) q[0];
sx q[0];
rz(-2.4242171) q[0];
sx q[0];
rz(-1.7051359) q[0];
rz(1.4868078) q[1];
sx q[1];
rz(-0.27324067) q[1];
sx q[1];
rz(-2.9437093) q[1];
rz(3.1080053) q[2];
sx q[2];
rz(-2.941249) q[2];
sx q[2];
rz(-2.8810101) q[2];
rz(1.6428357) q[3];
sx q[3];
rz(-0.34769736) q[3];
sx q[3];
rz(-3.0394239) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
