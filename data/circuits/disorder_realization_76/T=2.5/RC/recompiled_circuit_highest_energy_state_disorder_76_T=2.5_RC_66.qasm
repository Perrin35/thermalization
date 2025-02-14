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
rz(-2.2657356) q[0];
sx q[0];
rz(-2.1729204) q[0];
sx q[0];
rz(0.92585603) q[0];
rz(0.61697382) q[1];
sx q[1];
rz(8.7595688) q[1];
sx q[1];
rz(8.1002357) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61875859) q[0];
sx q[0];
rz(-1.5069739) q[0];
sx q[0];
rz(-0.60019779) q[0];
rz(-pi) q[1];
rz(0.1466897) q[2];
sx q[2];
rz(-0.58977276) q[2];
sx q[2];
rz(-1.9751994) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.263773) q[1];
sx q[1];
rz(-0.78697453) q[1];
sx q[1];
rz(-2.469648) q[1];
rz(-pi) q[2];
rz(0.76530568) q[3];
sx q[3];
rz(-2.0767143) q[3];
sx q[3];
rz(1.4368865) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.2609743) q[2];
sx q[2];
rz(-1.7558492) q[2];
sx q[2];
rz(-2.3495038) q[2];
rz(0.875862) q[3];
sx q[3];
rz(-3.0328817) q[3];
sx q[3];
rz(0.66332269) q[3];
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
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51139128) q[0];
sx q[0];
rz(-1.317861) q[0];
sx q[0];
rz(2.200101) q[0];
rz(2.1425653) q[1];
sx q[1];
rz(-2.2143054) q[1];
sx q[1];
rz(-0.082854465) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1076528) q[0];
sx q[0];
rz(-1.6710738) q[0];
sx q[0];
rz(-2.4058008) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8993334) q[2];
sx q[2];
rz(-1.9177116) q[2];
sx q[2];
rz(-1.0936979) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.283047) q[1];
sx q[1];
rz(-2.0783922) q[1];
sx q[1];
rz(1.0495484) q[1];
x q[2];
rz(-0.96534713) q[3];
sx q[3];
rz(-1.0167398) q[3];
sx q[3];
rz(-0.29747552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.89722172) q[2];
sx q[2];
rz(-1.7873849) q[2];
sx q[2];
rz(1.2574035) q[2];
rz(-0.8463549) q[3];
sx q[3];
rz(-0.1592764) q[3];
sx q[3];
rz(-2.4634821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1988679) q[0];
sx q[0];
rz(-1.6833479) q[0];
sx q[0];
rz(-2.9123836) q[0];
rz(-2.9275059) q[1];
sx q[1];
rz(-0.87132088) q[1];
sx q[1];
rz(2.7600938) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5435181) q[0];
sx q[0];
rz(-1.9207354) q[0];
sx q[0];
rz(0.66611163) q[0];
rz(-0.83865954) q[2];
sx q[2];
rz(-1.7146972) q[2];
sx q[2];
rz(0.14140192) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.4134564) q[1];
sx q[1];
rz(-1.6237139) q[1];
sx q[1];
rz(0.028832988) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.34708628) q[3];
sx q[3];
rz(-1.5659837) q[3];
sx q[3];
rz(-2.6322685) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.7233589) q[2];
sx q[2];
rz(-2.8853719) q[2];
sx q[2];
rz(-2.5733433) q[2];
rz(1.4189643) q[3];
sx q[3];
rz(-1.4293554) q[3];
sx q[3];
rz(-2.0645781) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(0.028246183) q[0];
sx q[0];
rz(-2.009511) q[0];
sx q[0];
rz(-2.8908253) q[0];
rz(-0.084550683) q[1];
sx q[1];
rz(-1.0944347) q[1];
sx q[1];
rz(-1.2006522) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4290743) q[0];
sx q[0];
rz(-1.2116644) q[0];
sx q[0];
rz(0.85690367) q[0];
rz(-pi) q[1];
rz(-2.8941706) q[2];
sx q[2];
rz(-1.1411925) q[2];
sx q[2];
rz(-1.0216433) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.1657287) q[1];
sx q[1];
rz(-1.4714676) q[1];
sx q[1];
rz(0.71290675) q[1];
x q[2];
rz(0.26504604) q[3];
sx q[3];
rz(-2.0816021) q[3];
sx q[3];
rz(3.0731415) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.4667929) q[2];
sx q[2];
rz(-1.1226706) q[2];
sx q[2];
rz(-1.8505081) q[2];
rz(-0.42168266) q[3];
sx q[3];
rz(-0.45160523) q[3];
sx q[3];
rz(1.5765367) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52294937) q[0];
sx q[0];
rz(-0.72294253) q[0];
sx q[0];
rz(1.6408386) q[0];
rz(-0.067642637) q[1];
sx q[1];
rz(-1.5875971) q[1];
sx q[1];
rz(-1.1281475) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2826506) q[0];
sx q[0];
rz(-0.19154405) q[0];
sx q[0];
rz(2.9454977) q[0];
rz(1.0447776) q[2];
sx q[2];
rz(-1.8785718) q[2];
sx q[2];
rz(0.36164944) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.6349244) q[1];
sx q[1];
rz(-1.6841869) q[1];
sx q[1];
rz(0.054111295) q[1];
rz(1.2913843) q[3];
sx q[3];
rz(-1.3456723) q[3];
sx q[3];
rz(1.0042508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.03380123) q[2];
sx q[2];
rz(-1.1148323) q[2];
sx q[2];
rz(2.4647253) q[2];
rz(-0.90773165) q[3];
sx q[3];
rz(-1.9151442) q[3];
sx q[3];
rz(2.7270253) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9229729) q[0];
sx q[0];
rz(-2.3899879) q[0];
sx q[0];
rz(0.76939097) q[0];
rz(-1.9006624) q[1];
sx q[1];
rz(-0.85378328) q[1];
sx q[1];
rz(-2.9262537) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.365959) q[0];
sx q[0];
rz(-1.3833519) q[0];
sx q[0];
rz(0.98926095) q[0];
rz(-3.1245232) q[2];
sx q[2];
rz(-1.2400377) q[2];
sx q[2];
rz(1.1137373) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.0508381) q[1];
sx q[1];
rz(-1.8148481) q[1];
sx q[1];
rz(0.65907101) q[1];
rz(-1.9086544) q[3];
sx q[3];
rz(-2.2400899) q[3];
sx q[3];
rz(0.026642628) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.26663366) q[2];
sx q[2];
rz(-1.6338438) q[2];
sx q[2];
rz(2.5249262) q[2];
rz(2.941361) q[3];
sx q[3];
rz(-2.4112406) q[3];
sx q[3];
rz(-2.0088137) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.013833372) q[0];
sx q[0];
rz(-2.5795689) q[0];
sx q[0];
rz(0.01509893) q[0];
rz(-0.12241441) q[1];
sx q[1];
rz(-1.8465123) q[1];
sx q[1];
rz(2.1133568) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1249753) q[0];
sx q[0];
rz(-0.65799558) q[0];
sx q[0];
rz(-0.090975002) q[0];
x q[1];
rz(-1.8625268) q[2];
sx q[2];
rz(-2.5030067) q[2];
sx q[2];
rz(-1.2995468) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.15785698) q[1];
sx q[1];
rz(-2.7089556) q[1];
sx q[1];
rz(3.0745201) q[1];
rz(-pi) q[2];
rz(-2.8260303) q[3];
sx q[3];
rz(-1.3129819) q[3];
sx q[3];
rz(-0.97168789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.60468173) q[2];
sx q[2];
rz(-1.9313507) q[2];
sx q[2];
rz(0.49986419) q[2];
rz(-0.31575051) q[3];
sx q[3];
rz(-0.67086589) q[3];
sx q[3];
rz(-2.1226814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71378088) q[0];
sx q[0];
rz(-1.2232057) q[0];
sx q[0];
rz(-0.94386238) q[0];
rz(1.7367412) q[1];
sx q[1];
rz(-1.8753139) q[1];
sx q[1];
rz(2.2199383) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93225828) q[0];
sx q[0];
rz(-2.8186322) q[0];
sx q[0];
rz(-2.9684116) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8692851) q[2];
sx q[2];
rz(-1.3278074) q[2];
sx q[2];
rz(1.1557494) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7870175) q[1];
sx q[1];
rz(-1.0598618) q[1];
sx q[1];
rz(-2.4035011) q[1];
rz(-pi) q[2];
rz(0.0042762386) q[3];
sx q[3];
rz(-1.0579526) q[3];
sx q[3];
rz(-2.9073334) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.9965808) q[2];
sx q[2];
rz(-1.8700446) q[2];
sx q[2];
rz(-1.5550295) q[2];
rz(-1.0541213) q[3];
sx q[3];
rz(-1.8127706) q[3];
sx q[3];
rz(1.4012977) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15040511) q[0];
sx q[0];
rz(-2.0441971) q[0];
sx q[0];
rz(2.4990668) q[0];
rz(-1.4379028) q[1];
sx q[1];
rz(-0.75690401) q[1];
sx q[1];
rz(-2.9978851) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3966345) q[0];
sx q[0];
rz(-1.2681343) q[0];
sx q[0];
rz(-2.55654) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2177198) q[2];
sx q[2];
rz(-2.0028186) q[2];
sx q[2];
rz(-0.96521711) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.70348008) q[1];
sx q[1];
rz(-2.2128155) q[1];
sx q[1];
rz(1.9025406) q[1];
rz(-0.71227422) q[3];
sx q[3];
rz(-2.5872719) q[3];
sx q[3];
rz(-1.2620776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.6045427) q[2];
sx q[2];
rz(-1.9436676) q[2];
sx q[2];
rz(0.26270467) q[2];
rz(-0.25775868) q[3];
sx q[3];
rz(-2.7596605) q[3];
sx q[3];
rz(-0.8530544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2951374) q[0];
sx q[0];
rz(-1.4895804) q[0];
sx q[0];
rz(-1.9211796) q[0];
rz(2.659761) q[1];
sx q[1];
rz(-2.316663) q[1];
sx q[1];
rz(-0.89734546) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37384826) q[0];
sx q[0];
rz(-1.4486827) q[0];
sx q[0];
rz(-0.51885817) q[0];
x q[1];
rz(0.32081066) q[2];
sx q[2];
rz(-1.6751373) q[2];
sx q[2];
rz(-3.0719245) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.5612302) q[1];
sx q[1];
rz(-1.3127316) q[1];
sx q[1];
rz(-3.0901872) q[1];
rz(-pi) q[2];
rz(-1.7190211) q[3];
sx q[3];
rz(-2.5276042) q[3];
sx q[3];
rz(-0.84166354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.4974978) q[2];
sx q[2];
rz(-0.92659014) q[2];
sx q[2];
rz(-0.11478718) q[2];
rz(0.74350205) q[3];
sx q[3];
rz(-1.8758352) q[3];
sx q[3];
rz(2.6928597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7262309) q[0];
sx q[0];
rz(-2.3804433) q[0];
sx q[0];
rz(-0.95809715) q[0];
rz(1.505898) q[1];
sx q[1];
rz(-2.0151357) q[1];
sx q[1];
rz(-1.9912079) q[1];
rz(0.74581292) q[2];
sx q[2];
rz(-1.0619166) q[2];
sx q[2];
rz(-0.3084736) q[2];
rz(-2.791849) q[3];
sx q[3];
rz(-2.6833224) q[3];
sx q[3];
rz(0.53862017) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
